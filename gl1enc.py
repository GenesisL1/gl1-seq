#!/usr/bin/env python3
"""
gl1enc.py — Minimal-dependency encoder/decoder for:
  • 2-bit nucleotides (A,C,G,T; optional U->T canonicalization)
  • 4-bit nucleotides (SAM/BAM NT16; optional U-preserving table)
  • 6-bit proteins (20 aa + * X B Z J U O -)

Features:
  • FASTA / FASTQ / RAW input (auto-detect)
  • Multi-record container (.gl1): one .gl1 file can store many sequences (FASTA entries / FASTQ reads)
  • Chunking for large sequences (chromosome-scale)
  • Decode single file OR decode many files by pointing to a folder (recursive optional)
  • Optional filename prefixing when combining outputs

Dependencies: Python stdlib only.

-------------------------------------------------------
GL1F v1 container (simple, self-contained, fast IO)
-------------------------------------------------------
FILE HEADER (20 bytes):
  <4s H H I Q
    magic=b'GL1F'
    version=1
    reserved=0
    record_count (u32)
    reserved_q=0

RECORD HEADER (32 bytes) + name bytes:
  <4s H B B I Q I H H I
    magic=b'GL1E'
    version=1
    kind: 1=nucleotide, 2=protein
    codec: 1=NUC2, 2=NUC4_BAM, 3=NUC4_U, 16=PROT6
    flags (u32)
    length_symbols (u64)
    chunk_symbols (u32)   # 0 => unchunked (single chunk per record)
    name_len (u16)
    reserved (u16)
    chunk_count (u32)
  name bytes (UTF-8)

CHUNKS:
  for each chunk:
    <I I  (chunk_symbols, chunk_bytes)
    chunk_payload bytes

Flags (u32):
  bit0: HAD_U
  bit1: HAD_LOWERCASE
  bit2: FORCED_MODE
  bit3: INPUT_FASTQ
  bit4: NUC_U2T_APPLIED (informational)
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import struct
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Tuple

# ----------------------------
# Constants / IDs
# ----------------------------

MAGIC_FILE = b"GL1F"
MAGIC_REC  = b"GL1E"
FILE_VERSION = 1
REC_VERSION = 1

KIND_NUC = 1
KIND_PROT = 2

CODEC_NUC2 = 1
CODEC_NUC4_BAM = 2
CODEC_NUC4_U = 3
CODEC_PROT6 = 16

FLAG_HAD_U = 1 << 0
FLAG_HAD_LOWERCASE = 1 << 1
FLAG_FORCED_MODE = 1 << 2
FLAG_INPUT_FASTQ = 1 << 3
FLAG_NUC_U2T_APPLIED = 1 << 4

FILE_HDR_FMT = "<4sHHIQ"     # magic, ver, reserved, record_count, reserved_q
REC_HDR_FMT  = "<4sHBBIQIHHI"  # magic, ver, kind, codec, flags, length, chunk_symbols, name_len, reserved, chunk_count
CHUNK_HDR_FMT = "<II"

FILE_HDR_SIZE = struct.calcsize(FILE_HDR_FMT)
REC_HDR_SIZE  = struct.calcsize(REC_HDR_FMT)
CHUNK_HDR_SIZE = struct.calcsize(CHUNK_HDR_FMT)

# ----------------------------
# Alphabets (for inference)
# ----------------------------

# Nucleotide IUPAC + '=' (BAM), plus U
NUC_SET = set("ACGTURYSWKMBDHVN=")
# Protein letters + ambiguity + stop + gap
PROT_SET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-")
# Letters that strongly imply protein (not used in nucleotide IUPAC set)
PROT_STRONG = set("EFILPQOJZX*")

# ----------------------------
# Code tables
# ----------------------------

# NUC2: A=0, C=1, G=2, T=3 (U treated as T in NUC2 encoding)
_NUC2_CODE: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}

# SAM/BAM NT16: "=ACMGRSVTWYHKDBN" -> 0..15
_NUC4_BAM_ORDER = "=ACMGRSVTWYHKDBN"
_NUC4_BAM_CODE: Dict[str, int] = {ch: i for i, ch in enumerate(_NUC4_BAM_ORDER)}

# NUC4_U: preserve U by mapping U -> 0 (uses '=' slot); keep other BAM codes.
_NUC4_U_CODE: Dict[str, int] = dict(_NUC4_BAM_CODE)
_NUC4_U_CODE["U"] = 0
_NUC4_U_CODE["="] = 0  # tolerate '='

# PROT6: 0..19 are 20 canonical amino acids in fixed order
_AA20 = "ACDEFGHIKLMNPQRSTVWY"
_PROT6_CODE: Dict[str, int] = {aa: i for i, aa in enumerate(_AA20)}
_PROT6_CODE.update({
    "*": 20,
    "X": 21,
    "B": 22,  # D or N
    "Z": 23,  # E or Q
    "J": 24,  # I or L
    "U": 25,  # selenocysteine
    "O": 26,  # pyrrolysine
    "-": 27,  # gap
})

# Reverse maps
_NUC2_REV = b"ACGT"
_NUC4_BAM_REV = bytes(_NUC4_BAM_ORDER, "ascii")
_NUC4_U_REV = bytearray(_NUC4_BAM_REV)
_NUC4_U_REV[0] = ord("U")
_NUC4_U_REV = bytes(_NUC4_U_REV)

_PROT6_REV = bytearray(b"X" * 64)
for ch, code in _PROT6_CODE.items():
    _PROT6_REV[code] = ord(ch)
_PROT6_REV = bytes(_PROT6_REV)

def _make_translate_table(mapping: Dict[str, int], default: int = 255) -> bytes:
    tbl = bytearray([default] * 256)
    for ch, code in mapping.items():
        tbl[ord(ch)] = code
        lo = ch.lower()
        if lo != ch:
            tbl[ord(lo)] = code
    return bytes(tbl)

NUC2_TRANS = _make_translate_table(_NUC2_CODE)
# NUC4_BAM cannot represent U; we map U->T unless strict mode rejects it
NUC4_BAM_TRANS = _make_translate_table({**_NUC4_BAM_CODE, "U": _NUC4_BAM_CODE["T"]})
NUC4_U_TRANS   = _make_translate_table(_NUC4_U_CODE)

PROT6_TRANS = _make_translate_table(_PROT6_CODE)

# Fast decode LUTs
_NUC2_BYTE_TO_4 = [None] * 256
for b in range(256):
    _NUC2_BYTE_TO_4[b] = bytes((
        _NUC2_REV[(b >> 6) & 3],
        _NUC2_REV[(b >> 4) & 3],
        _NUC2_REV[(b >> 2) & 3],
        _NUC2_REV[b & 3],
    ))

_NUC4_BYTE_TO_2_BAM = [None] * 256
_NUC4_BYTE_TO_2_U   = [None] * 256
for b in range(256):
    hi = (b >> 4) & 0xF
    lo = b & 0xF
    _NUC4_BYTE_TO_2_BAM[b] = bytes((_NUC4_BAM_REV[hi], _NUC4_BAM_REV[lo]))
    _NUC4_BYTE_TO_2_U[b]   = bytes((_NUC4_U_REV[hi], _NUC4_U_REV[lo]))

# ----------------------------
# Plans
# ----------------------------

@dataclass
class RecordPlan:
    name: str
    length: int
    kind: int
    codec: int
    flags: int
    chunk_symbols: int
    chunk_count: int

# ----------------------------
# Input formats
# ----------------------------

class InputFormat:
    FASTA = "fasta"
    FASTQ = "fastq"
    RAW = "raw"
    AUTO = "auto"

def _open_text(path: str) -> io.TextIOBase:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")

def detect_format(path: str, forced: str) -> str:
    if forced != InputFormat.AUTO:
        return forced
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                return InputFormat.FASTA
            if s.startswith("@"):
                return InputFormat.FASTQ
            return InputFormat.RAW
    return InputFormat.RAW

def _clean_seq_line(line: str) -> str:
    return "".join(ch for ch in line.strip() if not ch.isspace())

def stream_events(handle: io.TextIOBase, fmt: str) -> Iterator[Tuple[str, Optional[str]]]:
    """
    Yields ('record_start', name), ('seq', chunk), ('record_end', None)
    """
    if fmt == InputFormat.FASTA:
        name = None
        for line in handle:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield ("record_end", None)
                header = line[1:].strip()
                nm = header.split()[0] if header else "seq"
                name = nm
                yield ("record_start", name)
            else:
                s = _clean_seq_line(line)
                if s:
                    yield ("seq", s)
        if name is not None:
            yield ("record_end", None)
        return

    if fmt == InputFormat.FASTQ:
        while True:
            header = handle.readline()
            if not header:
                break
            header = header.rstrip("\n\r")
            if not header:
                continue
            if not header.startswith("@"):
                # degrade to raw
                yield ("record_start", "seq1")
                yield ("seq", _clean_seq_line(header))
                for line in handle:
                    s = _clean_seq_line(line)
                    if s:
                        yield ("seq", s)
                yield ("record_end", None)
                return

            name = header[1:].strip().split()[0] if header[1:].strip() else "seq"
            yield ("record_start", name)

            # sequence until '+'
            seq_len = 0
            while True:
                line = handle.readline()
                if not line:
                    break
                if line.startswith("+"):
                    break
                s = _clean_seq_line(line)
                if s:
                    seq_len += len(s)
                    yield ("seq", s)

            # skip quality lines (length seq_len)
            remaining = seq_len
            while remaining > 0:
                q = handle.readline()
                if not q:
                    break
                qline = q.rstrip("\n\r")
                remaining -= len(qline)

            yield ("record_end", None)
        return

    # RAW
    yield ("record_start", "seq1")
    for line in handle:
        s = _clean_seq_line(line)
        if s:
            yield ("seq", s)
    yield ("record_end", None)

# ----------------------------
# Inference / codec selection
# ----------------------------

def infer_kind(letters_upper: set) -> int:
    letters = set(letters_upper)
    letters.discard("-")
    if letters & PROT_STRONG:
        return KIND_PROT
    if letters <= NUC_SET:
        return KIND_NUC
    if letters <= PROT_SET:
        return KIND_PROT
    invalid = sorted(ch for ch in letters if ch not in NUC_SET and ch not in PROT_SET)
    raise ValueError(f"Unrecognized symbols in sequence: {invalid[:30]}{'...' if len(invalid) > 30 else ''}")

def choose_codec(
    letters_upper: set,
    forced_mode: str,
    nuc4_table: str,
    strict: bool,
    had_lower: bool,
    had_u: bool,
    unknown_to_x: bool,
) -> Tuple[int, int, int]:
    """
    Returns (kind, codec, flags).
    forced_mode: auto|2bit|4bit|6bit
    """
    flags = 0
    if had_lower:
        flags |= FLAG_HAD_LOWERCASE
    if had_u:
        flags |= FLAG_HAD_U
    if forced_mode != "auto":
        flags |= FLAG_FORCED_MODE

    # Forced modes override auto inference of kind/codec.
    if forced_mode == "2bit":
        kind = KIND_NUC
        codec = CODEC_NUC2
        allowed = set("ACGTU")
        if not letters_upper <= allowed:
            bad = sorted(letters_upper - allowed)
            raise ValueError(f"--encode 2bit requested but found non-ACGTU symbols: {bad[:30]}")
        if had_u:
            if strict:
                raise ValueError("--encode 2bit in strict mode cannot represent U; use --encode 4bit or disable --strict")
            flags |= FLAG_NUC_U2T_APPLIED
        return kind, codec, flags

    if forced_mode == "4bit":
        kind = KIND_NUC
        if nuc4_table == "auto":
            nuc4_table = "u" if had_u else "bam"
        codec = CODEC_NUC4_BAM if nuc4_table == "bam" else CODEC_NUC4_U
        # Validate nucleotide-ish letters
        if not letters_upper <= NUC_SET:
            bad = sorted(letters_upper - NUC_SET)
            raise ValueError(f"--encode 4bit requested but found non-nucleotide symbols: {bad[:30]}")
        if codec == CODEC_NUC4_BAM and strict and had_u:
            raise ValueError("--encode 4bit --nuc4-table bam in strict mode cannot represent U; use --nuc4-table u")
        return kind, codec, flags

    if forced_mode == "6bit":
        kind = KIND_PROT
        codec = CODEC_PROT6
        # If unknown_to_x, allow A-Z plus '*' '-' (others error). Otherwise require PROT_SET.
        if unknown_to_x:
            for ch in letters_upper:
                if ch == "*" or ch == "-":
                    continue
                if not ("A" <= ch <= "Z"):
                    raise ValueError(f"--encode 6bit: invalid protein character {repr(ch)} (non A-Z/*/-)")
        else:
            if not letters_upper <= PROT_SET:
                bad = sorted(letters_upper - PROT_SET)
                raise ValueError(f"--encode 6bit requested but found non-protein symbols: {bad[:30]}")
        return kind, codec, flags

    # AUTO:
    kind = infer_kind(letters_upper)

    if kind == KIND_PROT:
        return KIND_PROT, CODEC_PROT6, flags

    # nucleotide auto codec
    core = set("ACGT")
    core_u = set("ACGTU")

    if letters_upper <= core:
        return KIND_NUC, CODEC_NUC2, flags

    if letters_upper <= core_u and had_u:
        # preserve U by default
        if nuc4_table in ("auto", "u"):
            return KIND_NUC, CODEC_NUC4_U, flags
        return KIND_NUC, CODEC_NUC4_BAM, flags

    # ambiguity => 4-bit
    if nuc4_table == "auto":
        return KIND_NUC, (CODEC_NUC4_U if had_u else CODEC_NUC4_BAM), flags
    return KIND_NUC, (CODEC_NUC4_BAM if nuc4_table == "bam" else CODEC_NUC4_U), flags

# ----------------------------
# Encoding primitives
# ----------------------------

def enc_nuc2_bytes(seq_bytes_upper: bytes) -> bytes:
    codes = seq_bytes_upper.translate(NUC2_TRANS)
    if 255 in codes:
        bad = sorted(set(chr(ch) for ch in seq_bytes_upper if NUC2_TRANS[ch] == 255))
        raise ValueError(f"NUC2 cannot encode symbols: {bad[:30]}")
    n = len(codes)
    out = bytearray((n + 3) // 4)
    mv = codes
    i = 0
    j = 0
    n4 = n - (n % 4)
    while i < n4:
        out[j] = (mv[i] << 6) | (mv[i+1] << 4) | (mv[i+2] << 2) | mv[i+3]
        i += 4
        j += 1
    rem = n - n4
    if rem:
        b0 = mv[n4] << 6
        b1 = (mv[n4+1] << 4) if rem > 1 else 0
        b2 = (mv[n4+2] << 2) if rem > 2 else 0
        out[j] = b0 | b1 | b2
    return bytes(out)

def enc_nuc4_bytes(seq_bytes_upper: bytes, table: str, strict: bool) -> bytes:
    trans = NUC4_BAM_TRANS if table == "bam" else NUC4_U_TRANS
    codes = seq_bytes_upper.translate(trans)
    if 255 in codes:
        bad = sorted(set(chr(ch) for ch in seq_bytes_upper if trans[ch] == 255))
        raise ValueError(f"NUC4 cannot encode symbols: {bad[:30]}")
    if table == "bam" and strict and (b"U" in seq_bytes_upper):
        raise ValueError("NUC4 bam table cannot represent U in strict mode; use --nuc4-table u")
    n = len(codes)
    out = bytearray((n + 1) // 2)
    mv = codes
    i = 0
    j = 0
    n2 = n - (n % 2)
    while i < n2:
        out[j] = (mv[i] << 4) | mv[i+1]
        i += 2
        j += 1
    if n % 2:
        out[j] = (mv[n-1] << 4)
    return bytes(out)

def enc_prot6_bytes(seq_bytes_upper: bytes, unknown_to_x: bool) -> bytes:
    codes = seq_bytes_upper.translate(PROT6_TRANS)
    if 255 in codes:
        if unknown_to_x:
            arr = bytearray(codes)
            xcode = _PROT6_CODE["X"]
            for i, ch in enumerate(seq_bytes_upper):
                if PROT6_TRANS[ch] == 255:
                    arr[i] = xcode
            codes = bytes(arr)
        else:
            bad = sorted(set(chr(ch) for ch in seq_bytes_upper if PROT6_TRANS[ch] == 255))
            raise ValueError(f"PROT6 cannot encode symbols: {bad[:30]}")
    out = bytearray()
    acc = 0
    acc_bits = 0
    for c in codes:
        acc = (acc << 6) | (c & 0x3F)
        acc_bits += 6
        while acc_bits >= 8:
            acc_bits -= 8
            out.append((acc >> acc_bits) & 0xFF)
            acc &= (1 << acc_bits) - 1 if acc_bits else 0
    if acc_bits:
        out.append((acc << (8 - acc_bits)) & 0xFF)
    return bytes(out)

# ----------------------------
# Decoding primitives
# ----------------------------

def dec_nuc2_chunk(chunk_bytes: bytes, symbols: int) -> bytes:
    out = bytearray()
    for b in chunk_bytes:
        out.extend(_NUC2_BYTE_TO_4[b])
    return bytes(out[:symbols])

def dec_nuc4_chunk(chunk_bytes: bytes, symbols: int, table: str) -> bytes:
    out = bytearray()
    lut = _NUC4_BYTE_TO_2_BAM if table == "bam" else _NUC4_BYTE_TO_2_U
    for b in chunk_bytes:
        out.extend(lut[b])
    return bytes(out[:symbols])

def dec_prot6_chunk(chunk_bytes: bytes, symbols: int) -> bytes:
    out = bytearray()
    acc = 0
    acc_bits = 0
    produced = 0
    for b in chunk_bytes:
        acc = (acc << 8) | b
        acc_bits += 8
        while acc_bits >= 6 and produced < symbols:
            acc_bits -= 6
            code = (acc >> acc_bits) & 0x3F
            out.append(_PROT6_REV[code])
            produced += 1
            acc &= (1 << acc_bits) - 1 if acc_bits else 0
        if produced >= symbols:
            break
    return bytes(out)

# ----------------------------
# Scan pass
# ----------------------------

def scan_records(
    path: str,
    fmt: str,
    forced_mode: str,
    nuc4_table: str,
    strict: bool,
    unknown_to_x: bool,
) -> Tuple[List[RecordPlan], str]:
    detected = detect_format(path, fmt)
    plans: List[RecordPlan] = []

    current_name: Optional[str] = None
    letters_upper: set = set()
    length = 0
    had_lower = False
    had_u = False

    with _open_text(path) as f:
        for ev, val in stream_events(f, detected):
            if ev == "record_start":
                current_name = val or "seq"
                letters_upper = set()
                length = 0
                had_lower = False
                had_u = False

            elif ev == "seq":
                s = val or ""
                if not s:
                    continue
                u = s.upper()
                if u != s:
                    had_lower = True
                if "U" in u:
                    had_u = True
                letters_upper |= set(u)
                length += len(u)

            elif ev == "record_end":
                if current_name is None:
                    continue
                kind, codec, flags = choose_codec(
                    letters_upper,
                    forced_mode=forced_mode,
                    nuc4_table=nuc4_table,
                    strict=strict,
                    had_lower=had_lower,
                    had_u=had_u,
                    unknown_to_x=unknown_to_x,
                )
                plans.append(RecordPlan(
                    name=current_name,
                    length=length,
                    kind=kind,
                    codec=codec,
                    flags=flags,
                    chunk_symbols=0,
                    chunk_count=1,
                ))
                current_name = None

    return plans, detected

# ----------------------------
# Encode: write GL1F
# ----------------------------

def encode_file(
    in_path: str,
    out_path: str,
    fmt: str,
    forced_mode: str,
    chunk_symbols: int,
    nuc4_table: str,
    strict: bool,
    unknown_to_x: bool,
):
    plans, detected_fmt = scan_records(in_path, fmt, forced_mode, nuc4_table, strict, unknown_to_x)

    for p in plans:
        if chunk_symbols and chunk_symbols > 0:
            p.chunk_symbols = int(chunk_symbols)
            p.chunk_count = (p.length + p.chunk_symbols - 1) // p.chunk_symbols
        else:
            p.chunk_symbols = 0
            p.chunk_count = 1

    with open(out_path, "wb") as out:
        out.write(struct.pack(FILE_HDR_FMT, MAGIC_FILE, FILE_VERSION, 0, len(plans), 0))

        with _open_text(in_path) as f:
            events = stream_events(f, detected_fmt)
            plan_i = -1
            active: Optional[RecordPlan] = None
            buf = bytearray()
            is_fastq = (detected_fmt == InputFormat.FASTQ)

            for ev, val in events:
                if ev == "record_start":
                    plan_i += 1
                    if plan_i >= len(plans):
                        raise RuntimeError("Input changed between scan and encode (more records found).")
                    active = plans[plan_i]
                    buf.clear()

                    name_bytes = active.name.encode("utf-8")
                    flags = active.flags | (FLAG_INPUT_FASTQ if is_fastq else 0)

                    out.write(struct.pack(
                        REC_HDR_FMT,
                        MAGIC_REC,
                        REC_VERSION,
                        active.kind,
                        active.codec,
                        flags,
                        active.length,
                        active.chunk_symbols,
                        len(name_bytes),
                        0,
                        active.chunk_count,
                    ))
                    out.write(name_bytes)

                elif ev == "seq":
                    if active is None:
                        continue
                    s = val or ""
                    if not s:
                        continue
                    b = s.upper().encode("ascii", "strict")
                    buf.extend(b)

                    if active.chunk_symbols and active.chunk_symbols > 0:
                        cs = active.chunk_symbols
                        while len(buf) >= cs:
                            part = bytes(buf[:cs])
                            del buf[:cs]
                            _write_chunk(out, active, part, nuc4_table, strict, unknown_to_x)

                elif ev == "record_end":
                    if active is None:
                        continue

                    if active.chunk_symbols and active.chunk_symbols > 0:
                        if buf:
                            _write_chunk(out, active, bytes(buf), nuc4_table, strict, unknown_to_x)
                            buf.clear()
                    else:
                        _write_chunk(out, active, bytes(buf), nuc4_table, strict, unknown_to_x)
                        buf.clear()
                    active = None

def _write_chunk(
    out: io.BufferedWriter,
    plan: RecordPlan,
    seq_bytes_upper: bytes,
    nuc4_table: str,
    strict: bool,
    unknown_to_x: bool,
):
    sym = len(seq_bytes_upper)
    if sym == 0:
        enc = b""
    else:
        if plan.kind == KIND_NUC:
            if plan.codec == CODEC_NUC2:
                enc = enc_nuc2_bytes(seq_bytes_upper)
            elif plan.codec == CODEC_NUC4_BAM:
                enc = enc_nuc4_bytes(seq_bytes_upper, table="bam", strict=strict)
            elif plan.codec == CODEC_NUC4_U:
                enc = enc_nuc4_bytes(seq_bytes_upper, table="u", strict=strict)
            else:
                raise ValueError(f"Unknown nucleotide codec {plan.codec}")
        else:
            if plan.codec != CODEC_PROT6:
                raise ValueError(f"Unknown protein codec {plan.codec}")
            enc = enc_prot6_bytes(seq_bytes_upper, unknown_to_x=unknown_to_x)

    out.write(struct.pack(CHUNK_HDR_FMT, sym, len(enc)))
    if enc:
        out.write(enc)

# ----------------------------
# Decode (single file) — write to handle
# ----------------------------

def decode_one_gl1(
    in_path: str,
    out_handle: io.TextIOBase,
    out_fmt: str,
    wrap: int,
    only_name: Optional[str],
    name_prefix: str = "",
):
    with open(in_path, "rb") as f:
        magic, ver, _, rec_count, _ = struct.unpack(FILE_HDR_FMT, f.read(FILE_HDR_SIZE))
        if magic != MAGIC_FILE or ver != FILE_VERSION:
            raise ValueError(f"Not a GL1F v{FILE_VERSION} file: {in_path}")

        for _ in range(rec_count):
            hdr = f.read(REC_HDR_SIZE)
            if len(hdr) != REC_HDR_SIZE:
                raise EOFError("Truncated record header")

            (rmagic, rver, kind, codec, flags, length, chunk_symbols, name_len, _, chunk_count) = struct.unpack(REC_HDR_FMT, hdr)
            if rmagic != MAGIC_REC or rver != REC_VERSION:
                raise ValueError("Bad GL1E record header")

            name = f.read(name_len).decode("utf-8", "replace")
            want = (only_name is None) or (name == only_name)

            if want and out_fmt == "fasta":
                out_handle.write(f">{name_prefix}{name}\n")

            line_buf = bytearray()
            nuc4_table = "bam" if codec == CODEC_NUC4_BAM else "u"

            for _ci in range(chunk_count):
                ch = f.read(CHUNK_HDR_SIZE)
                if len(ch) != CHUNK_HDR_SIZE:
                    raise EOFError("Truncated chunk header")
                sym, blen = struct.unpack(CHUNK_HDR_FMT, ch)
                payload = f.read(blen)
                if len(payload) != blen:
                    raise EOFError("Truncated chunk payload")

                if not want:
                    continue

                if kind == KIND_NUC:
                    if codec == CODEC_NUC2:
                        decoded = dec_nuc2_chunk(payload, sym)
                    elif codec in (CODEC_NUC4_BAM, CODEC_NUC4_U):
                        decoded = dec_nuc4_chunk(payload, sym, table=nuc4_table)
                    else:
                        raise ValueError(f"Unknown nucleotide codec {codec}")
                else:
                    if codec != CODEC_PROT6:
                        raise ValueError(f"Unknown protein codec {codec}")
                    decoded = dec_prot6_chunk(payload, sym)

                if out_fmt == "raw":
                    out_handle.write(decoded.decode("ascii"))
                    continue

                if wrap <= 0:
                    out_handle.write(decoded.decode("ascii"))
                    continue

                line_buf.extend(decoded)
                while len(line_buf) >= wrap:
                    out_handle.write(line_buf[:wrap].decode("ascii"))
                    out_handle.write("\n")
                    del line_buf[:wrap]

            if want and out_fmt == "fasta":
                if wrap > 0:
                    if line_buf:
                        out_handle.write(line_buf.decode("ascii"))
                        out_handle.write("\n")
                else:
                    out_handle.write("\n")

# ----------------------------
# Multi-input expansion
# ----------------------------

def expand_gl1_inputs(inputs: List[str], recursive: bool) -> List[str]:
    files: List[str] = []
    for p in inputs:
        if os.path.isdir(p):
            if recursive:
                for root, _, fnames in os.walk(p):
                    for fn in fnames:
                        if fn.endswith(".gl1"):
                            files.append(os.path.join(root, fn))
            else:
                for fn in os.listdir(p):
                    if fn.endswith(".gl1"):
                        files.append(os.path.join(p, fn))
        else:
            files.append(p)
    # de-dup and sort
    files = sorted(set(files))
    return files

# ----------------------------
# Decode dispatcher
# ----------------------------

def decode_many(
    in_paths: List[str],
    out_path: str,
    out_fmt: str,
    wrap: int,
    only_name: Optional[str],
    recursive: bool,
    outdir: Optional[str],
    prefix_file: bool,
):
    files = expand_gl1_inputs(in_paths, recursive=recursive)
    if not files:
        raise SystemExit("No input files found.")

    ext = ".fa" if out_fmt == "fasta" else ".txt"

    # Per-file outputs if outdir is provided OR out_path is an existing directory
    per_file = False
    if outdir is not None:
        per_file = True
        target_dir = outdir
    elif out_path != "-" and os.path.isdir(out_path):
        per_file = True
        target_dir = out_path
    else:
        target_dir = None

    if per_file:
        if out_path == "-":
            raise SystemExit("When decoding multiple files to per-file outputs, --out cannot be stdout ('-').")
        os.makedirs(target_dir, exist_ok=True)
        for fpath in files:
            stem = os.path.splitext(os.path.basename(fpath))[0]
            out_file = os.path.join(target_dir, stem + ext)
            with open(out_file, "wt", encoding="utf-8", newline="\n") as oh:
                prefix = f"{stem}|" if prefix_file else ""
                decode_one_gl1(fpath, oh, out_fmt, wrap, only_name, name_prefix=prefix)
        return

    # Combined output into one file or stdout
    if out_path == "-":
        oh = sys.stdout
        close = False
    else:
        oh = open(out_path, "wt", encoding="utf-8", newline="\n")
        close = True

    try:
        for fpath in files:
            stem = os.path.splitext(os.path.basename(fpath))[0]
            prefix = f"{stem}|" if prefix_file else ""
            decode_one_gl1(fpath, oh, out_fmt, wrap, only_name, name_prefix=prefix)
    finally:
        if close:
            oh.close()

# ----------------------------
# Info (supports multiple files/dirs)
# ----------------------------

def info_one(in_path: str):
    with open(in_path, "rb") as f:
        magic, ver, _, rec_count, _ = struct.unpack(FILE_HDR_FMT, f.read(FILE_HDR_SIZE))
        if magic != MAGIC_FILE:
            raise ValueError(f"Not a GL1F file: {in_path}")
        print(f"{in_path}: GL1F v{ver} records={rec_count}")
        for _ in range(rec_count):
            hdr = f.read(REC_HDR_SIZE)
            (rmagic, rver, kind, codec, flags, length, chunk_symbols, name_len, _, chunk_count) = struct.unpack(REC_HDR_FMT, hdr)
            name = f.read(name_len).decode("utf-8", "replace")

            for _ci in range(chunk_count):
                sym, blen = struct.unpack(CHUNK_HDR_FMT, f.read(CHUNK_HDR_SIZE))
                f.seek(blen, io.SEEK_CUR)

            kind_s = "NUC" if kind == KIND_NUC else "PROT"
            codec_s = {
                CODEC_NUC2: "NUC2",
                CODEC_NUC4_BAM: "NUC4_BAM",
                CODEC_NUC4_U: "NUC4_U",
                CODEC_PROT6: "PROT6",
            }.get(codec, str(codec))
            chunk_s = str(chunk_symbols) if chunk_symbols else "none"
            print(f"  - {name}: kind={kind_s} codec={codec_s} length={length} chunk={chunk_s} chunks={chunk_count} flags=0x{flags:08x}")

def info_many(in_paths: List[str], recursive: bool):
    files = expand_gl1_inputs(in_paths, recursive=recursive)
    if not files:
        raise SystemExit("No input files found.")
    for fpath in files:
        info_one(fpath)

# ----------------------------
# CLI
# ----------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="gl1enc", description="Encode/decode GL1 (2b nuc, 4b nuc, 6b prot) with chunking.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_enc = sub.add_parser("encode", help="Encode FASTA/FASTQ/raw into .gl1 (GL1F container).")
    ap_enc.add_argument("-i", "--in", dest="inp", required=True, help="Input file path (FASTA/FASTQ/raw). .gz supported. Use '-' for stdin.")
    ap_enc.add_argument("-o", "--out", dest="out", required=True, help="Output .gl1 file path.")
    ap_enc.add_argument("--format", choices=[InputFormat.AUTO, InputFormat.FASTA, InputFormat.FASTQ, InputFormat.RAW], default=InputFormat.AUTO)
    ap_enc.add_argument("--encode", dest="mode", choices=["auto", "2bit", "4bit", "6bit"], default="auto",
                        help="Force encoding mode (auto chooses per-record).")
    ap_enc.add_argument("--chunk", type=int, default=0, help="Chunk size in symbols (bases/residues). 0 = single chunk per record.")
    ap_enc.add_argument("--nuc4-table", choices=["auto", "bam", "u"], default="auto",
                        help="For 4-bit nucleotides: 'bam' uses SAM/BAM NT16 (U->T), 'u' preserves U as code 0.")
    ap_enc.add_argument("--strict", action="store_true",
                        help="Strict validation: reject symbols not representable in chosen codec (e.g., U in bam table).")
    ap_enc.add_argument("--unknown-to-x", action="store_true",
                        help="Protein: map unknown letters to X instead of error (also affects --encode 6bit validation).")

    ap_dec = sub.add_parser("decode", help="Decode .gl1 back to FASTA or raw. Supports files and folders.")
    ap_dec.add_argument("-i", "--in", dest="inp", nargs="+", required=True,
                        help="Input .gl1 file(s) and/or folder(s). If a folder is given, all *.gl1 inside are decoded.")
    ap_dec.add_argument("-o", "--out", dest="out", default="-",
                        help="Output path. '-' = stdout. If multiple inputs: provide a file path to combine, OR a directory to write per-file outputs.")
    ap_dec.add_argument("--outdir", default=None,
                        help="Directory for per-file outputs (overrides --out). One output per input .gl1 file.")
    ap_dec.add_argument("--recursive", action="store_true", help="When an input is a directory, scan subfolders for *.gl1.")
    ap_dec.add_argument("--to", choices=["fasta", "raw"], default="fasta")
    ap_dec.add_argument("--wrap", type=int, default=60, help="FASTA line wrap length (0 = no wrap).")
    ap_dec.add_argument("--record", default=None, help="Decode only the record with this exact name (per input file).")
    ap_dec.add_argument("--prefix-file", action="store_true",
                        help="When combining outputs from multiple input files, prefix FASTA record names with '<file>|' to avoid collisions.")

    ap_info = sub.add_parser("info", help="Print summary of .gl1 file(s). Supports folders.")
    ap_info.add_argument("-i", "--in", dest="inp", nargs="+", required=True, help="Input .gl1 file(s) and/or folder(s).")
    ap_info.add_argument("--recursive", action="store_true", help="When an input is a directory, scan subfolders for *.gl1.")

    args = ap.parse_args(argv)

    if args.cmd == "encode":
        if args.inp == "-":
            with tempfile.NamedTemporaryFile("wb", delete=False) as tf:
                tmp_path = tf.name
                tf.write(sys.stdin.buffer.read())
            try:
                encode_file(tmp_path, args.out, args.format, args.mode, args.chunk, args.nuc4_table, args.strict, args.unknown_to_x)
            finally:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
        else:
            encode_file(args.inp, args.out, args.format, args.mode, args.chunk, args.nuc4_table, args.strict, args.unknown_to_x)
        return 0

    if args.cmd == "decode":
        decode_many(
            in_paths=args.inp,
            out_path=args.out,
            out_fmt=args.to,
            wrap=args.wrap,
            only_name=args.record,
            recursive=args.recursive,
            outdir=args.outdir,
            prefix_file=args.prefix_file,
        )
        return 0

    if args.cmd == "info":
        info_many(args.inp, recursive=args.recursive)
        return 0

    return 2

if __name__ == "__main__":
    raise SystemExit(main())
