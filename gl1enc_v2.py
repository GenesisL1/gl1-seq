"""
gl1enc.py — GL1 Sequence Codec Engine (reference implementation)
==============================================================

This module defines the *canonical*, deterministic way GL1 represents biological
sequences as bytes for hashing, storage, and on-chain/off-chain parity.

Design goals
------------
- Deterministic and Solidity-portable rules (no locale dependencies).
- Canonical normalization (uppercase + ASCII whitespace stripped) so identical inputs
  always hash to identical bytes.
- Multiple compact encodings (2/4/6-bit) plus ASCII fallback:
    (0) DNA2: 2-bit-per-base packing for pure A/C/G/T sequences.
    (2) DNA4: 4-bit (nibble) encoding for IUPAC nucleotide ambiguity codes (+ optional gap),
              using bitmask semantics (A=1,C=2,G=4,T=8, ambiguities are ORs).
    (3) SIXBIT: 6-bit packed encoding for a fixed 64-character alphabet (ALPH64),
                useful for protein strings and other uppercase symbolic sequences.
    (1) ASCII: raw uppercase ASCII bytes fallback for everything else.

Why keep ASCII even with SIXBIT?
--------------------------------
SIXBIT only supports a fixed 64-char alphabet. ASCII is the "future proof" fallback
that can represent any normalized ASCII sequence (including symbols not in ALPH64).

Versioning
----------
- __version__ is semantic versioning for the Python module.
- The encoding *format* itself is versioned via GL1ENC_FORMAT.
  If you ever change bit packing, normalization, canonical selection rules,
  or the fixed alphabets/mappings, bump GL1ENC_FORMAT (breaking change).

Public API
----------
- SeqEncoding (IntEnum)
- EncodedSequence (dataclass)
- normalize_seq
- is_pure_acgt / is_dna4_eligible / is_sixbit_eligible
- choose_canonical_encoding
- dna2_byte_len / dna4_byte_len / sixbit_byte_len
- pack_* / unpack_* and validate_* helpers
- encode_seq / decode_seq
- canonical_encode / canonical_decode
- reverse_complement
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Dict, Tuple

__all__ = [
    "GL1ENC_FORMAT",
    "__version__",
    "SeqEncoding",
    "EncodedSequence",
    "normalize_seq",
    "is_pure_acgt",
    "is_dna4_eligible",
    "is_sixbit_eligible",
    "choose_canonical_encoding",
    "dna2_byte_len",
    "dna4_byte_len",
    "sixbit_byte_len",
    "pack_dna2",
    "unpack_dna2",
    "validate_dna2_bytes",
    "pack_dna4",
    "unpack_dna4",
    "validate_dna4_bytes",
    "pack_sixbit",
    "unpack_sixbit",
    "validate_sixbit_bytes",
    "encode_seq",
    "decode_seq",
    "canonical_encode",
    "canonical_decode",
    "reverse_complement",
]

__version__ = "2.0.0"

# IMPORTANT: Treat this as the wire-format identifier.
GL1ENC_FORMAT = "GL1ENCv2"


class SeqEncoding(IntEnum):
    """
    Canonical encoding identifiers (uint8 on-chain).

    IMPORTANT: IDs 0 and 1 are kept compatible with GL1ENCv1:
      0 = DNA2
      1 = ASCII

    New in GL1ENCv2:
      2 = DNA4 (IUPAC nibble / bitmask)
      3 = SIXBIT (packed 6-bit, ALPH64)
    """

    DNA2 = 0
    ASCII = 1
    DNA4 = 2
    SIXBIT = 3

    @classmethod
    def from_id(cls, enc_id: int) -> "SeqEncoding":
        try:
            return cls(enc_id)
        except ValueError as e:
            raise ValueError(f"Unknown SeqEncoding id: {enc_id}") from e


def normalize_seq(seq: str) -> str:
    """
    Canonical normalization applied before any encoding.

    Rules (must match Solidity implementation):
    - strip() leading/trailing whitespace
    - uppercase
    - remove ASCII spaces and \\n \\r \\t anywhere in the string

    Notes:
    - We deliberately do NOT remove all Unicode whitespace; inputs should be ASCII.
    - If non-ASCII characters remain, ASCII encoding will fail.
    """
    return (
        seq.strip()
        .upper()
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .replace("\t", "")
    )


# ---------------------------------------------------------------------
# DNA2 (2-bit) encoding
# ---------------------------------------------------------------------

# DNA2 mapping: A=00, C=01, G=10, T=11
_DNA2_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}
_DNA2_RMAP = {v: k for k, v in _DNA2_MAP.items()}


def is_pure_acgt(seq: str, *, assume_normalized: bool = False) -> bool:
    """
    True iff the (normalized) sequence consists only of A, C, G, T characters.

    This is the eligibility condition for DNA2 encoding.
    """
    s = seq if assume_normalized else normalize_seq(seq)
    for ch in s:
        if ch not in "ACGT":
            return False
    return True


def dna2_byte_len(length: int) -> int:
    """Bytes required to store `length` bases in DNA2 (4 bases/byte)."""
    if length < 0:
        raise ValueError("length must be >= 0")
    return (length + 3) // 4


def pack_dna2(seq: str) -> bytes:
    """
    Pack a pure A/C/G/T sequence into 2-bit DNA2 bytes.

    Bit layout within each byte:
      base0 -> bits 7..6 (MSB)
      base1 -> bits 5..4
      base2 -> bits 3..2
      base3 -> bits 1..0 (LSB)

    If the final chunk has fewer than 4 bases, remaining 2-bit slots are padded with 0 (A).
    """
    s = normalize_seq(seq)
    if not is_pure_acgt(s, assume_normalized=True):
        raise ValueError("DNA2 packing requires only A,C,G,T after normalization")

    out = bytearray()
    for i in range(0, len(s), 4):
        b = 0
        for j in range(4):
            if i + j < len(s):
                code = _DNA2_MAP[s[i + j]]
            else:
                code = 0  # padding (A)
            shift = (3 - j) * 2
            b |= (code & 0b11) << shift
        out.append(b)
    return bytes(out)


def validate_dna2_bytes(data: bytes, length: int) -> None:
    """
    Validate canonical DNA2 bytes for a given base length.

    Checks:
    - byte length matches dna2_byte_len(length)
    - unused (padded) 2-bit slots in the final byte are zero
    """
    if length < 0:
        raise ValueError("length must be >= 0")
    exp = dna2_byte_len(length)
    if len(data) != exp:
        raise ValueError(f"DNA2 byte length mismatch: got {len(data)} expected {exp}")

    rem = length % 4
    if rem == 0 or exp == 0:
        return

    last = data[-1]
    unused_slots = 4 - rem
    # unused bits are the least-significant 2*unused_slots bits
    mask = (1 << (2 * unused_slots)) - 1
    if (last & mask) != 0:
        raise ValueError("Non-canonical DNA2 padding bits: expected 0s in unused slots")


def unpack_dna2(data: bytes, length: int) -> str:
    """
    Decode DNA2 bytes into an A/C/G/T string, using `length` to ignore padding.
    """
    validate_dna2_bytes(data, length)

    out_chars = []
    n = 0
    for byte in data:
        for j in range(4):
            if n >= length:
                break
            shift = (3 - j) * 2
            code = (byte >> shift) & 0b11
            out_chars.append(_DNA2_RMAP[code])
            n += 1
        if n >= length:
            break
    return "".join(out_chars)


# ---------------------------------------------------------------------
# DNA4 (4-bit nibble / IUPAC) encoding
# ---------------------------------------------------------------------

# DNA4 is defined as a 4-bit *bitmask* of {A,C,G,T}:
# A=0001, C=0010, G=0100, T=1000
# Ambiguity codes are OR-combinations. We also assign GAP '-' = 0000.
#
# This yields a very useful property:
#   base_matches(code, base) <=> (mask(base) & mask(code)) != 0
#
# Mapping uses the standard IUPAC letters for nucleotides.
_DNA4_CHAR_TO_NIBBLE: Dict[str, int] = {
    "-": 0x0,  # gap / unknown/empty slot
    "A": 0x1,
    "C": 0x2,
    "G": 0x4,
    "T": 0x8,
    "M": 0x3,  # A or C
    "R": 0x5,  # A or G
    "S": 0x6,  # C or G
    "V": 0x7,  # A or C or G
    "W": 0x9,  # A or T
    "Y": 0xA,  # C or T
    "H": 0xB,  # A or C or T
    "K": 0xC,  # G or T
    "D": 0xD,  # A or G or T
    "B": 0xE,  # C or G or T
    "N": 0xF,  # A or C or G or T
}
_DNA4_NIBBLE_TO_CHAR: Dict[int, str] = {v: k for k, v in _DNA4_CHAR_TO_NIBBLE.items()}


def is_dna4_eligible(seq: str, *, assume_normalized: bool = False) -> bool:
    """
    True iff every character is encodable in DNA4.

    DNA4 supports: A,C,G,T and IUPAC ambiguity codes {MRWSYKVHDBN} plus '-' (gap).
    It intentionally does NOT include 'U' (RNA). If 'U' appears, SIXBIT/ASCII is used.
    """
    s = seq if assume_normalized else normalize_seq(seq)
    for ch in s:
        if ch not in _DNA4_CHAR_TO_NIBBLE:
            return False
    return True


def dna4_byte_len(length: int) -> int:
    """Bytes required to store `length` symbols in DNA4 (2 symbols/byte)."""
    if length < 0:
        raise ValueError("length must be >= 0")
    return (length + 1) // 2


def pack_dna4(seq: str) -> bytes:
    """
    Pack a DNA4-eligible sequence into 4-bit nibble bytes.

    Layout:
      symbol0 -> high nibble (bits 7..4)
      symbol1 -> low nibble  (bits 3..0)

    If length is odd, final low nibble is padded with 0 (GAP '-').
    """
    s = normalize_seq(seq)
    if not is_dna4_eligible(s, assume_normalized=True):
        raise ValueError("DNA4 packing requires only DNA4 alphabet symbols after normalization")

    out = bytearray()
    i = 0
    while i < len(s):
        hi = _DNA4_CHAR_TO_NIBBLE[s[i]] & 0xF
        lo = 0
        if i + 1 < len(s):
            lo = _DNA4_CHAR_TO_NIBBLE[s[i + 1]] & 0xF
        # if odd length, lo stays 0 as padding
        out.append((hi << 4) | lo)
        i += 2
    return bytes(out)


def validate_dna4_bytes(data: bytes, length: int) -> None:
    """
    Validate canonical DNA4 bytes for a given symbol length.

    Checks:
    - byte length matches dna4_byte_len(length)
    - if length is odd, the padded (low) nibble of the last byte must be 0
    - all nibbles used by `length` symbols map to defined DNA4 characters
    """
    if length < 0:
        raise ValueError("length must be >= 0")
    exp = dna4_byte_len(length)
    if len(data) != exp:
        raise ValueError(f"DNA4 byte length mismatch: got {len(data)} expected {exp}")

    if exp == 0:
        return

    # validate used nibbles
    for idx, byte in enumerate(data):
        hi = (byte >> 4) & 0xF
        lo = byte & 0xF
        sym_index = idx * 2
        if sym_index < length:
            if hi not in _DNA4_NIBBLE_TO_CHAR:
                raise ValueError(f"Invalid DNA4 high nibble 0x{hi:X} at byte index {idx}")
        if sym_index + 1 < length:
            if lo not in _DNA4_NIBBLE_TO_CHAR:
                raise ValueError(f"Invalid DNA4 low nibble 0x{lo:X} at byte index {idx}")
        else:
            # padded nibble must be zero
            if length % 2 == 1 and lo != 0:
                raise ValueError("Non-canonical DNA4 padding nibble: expected 0x0 in last low nibble")


def unpack_dna4(data: bytes, length: int) -> str:
    """
    Decode DNA4 bytes into a string, using `length` to ignore padding.
    """
    validate_dna4_bytes(data, length)

    out_chars = []
    n = 0
    for byte in data:
        if n >= length:
            break
        hi = (byte >> 4) & 0xF
        out_chars.append(_DNA4_NIBBLE_TO_CHAR[hi])
        n += 1
        if n >= length:
            break
        lo = byte & 0xF
        out_chars.append(_DNA4_NIBBLE_TO_CHAR[lo])
        n += 1
    return "".join(out_chars)


# ---------------------------------------------------------------------
# SIXBIT (6-bit packed) encoding
# ---------------------------------------------------------------------

# Fixed 64-character alphabet for SIXBIT.
# Must be exactly 64 chars; index is the 6-bit code.
SIXBIT_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-*._:;,|/\\+=()[]{}<>#$%&?!@^"
assert len(SIXBIT_ALPHABET) == 64, "SIXBIT_ALPHABET must be exactly 64 characters"

_SIXBIT_CHAR_TO_CODE: Dict[str, int] = {ch: i for i, ch in enumerate(SIXBIT_ALPHABET)}
_SIXBIT_CODE_TO_CHAR: Dict[int, str] = {i: ch for i, ch in enumerate(SIXBIT_ALPHABET)}


def is_sixbit_eligible(seq: str, *, assume_normalized: bool = False) -> bool:
    """
    True iff every character is encodable in SIXBIT (ALPH64).

    SIXBIT supports:
      - uppercase A-Z
      - digits 0-9
      - and the symbol set: -*._:;,|/\\+=()[]{}<>#$%&?!@^

    (Whitespace is removed by normalization.)
    """
    s = seq if assume_normalized else normalize_seq(seq)
    for ch in s:
        if ch not in _SIXBIT_CHAR_TO_CODE:
            return False
    return True


def sixbit_byte_len(length: int) -> int:
    """Bytes required to store `length` SIXBIT symbols (ceil(length*6/8))."""
    if length < 0:
        raise ValueError("length must be >= 0")
    return (length * 6 + 7) // 8


def pack_sixbit(seq: str) -> bytes:
    """
    Pack a SIXBIT-eligible sequence into bytes.

    Packing is MSB-first bitstream concatenation of 6-bit codes, padded with 0 bits at the end.
    """
    s = normalize_seq(seq)
    if not is_sixbit_eligible(s, assume_normalized=True):
        raise ValueError("SIXBIT packing requires only ALPH64 symbols after normalization")

    out = bytearray()
    buf = 0
    bits = 0
    for ch in s:
        code = _SIXBIT_CHAR_TO_CODE[ch] & 0x3F
        buf = (buf << 6) | code
        bits += 6
        while bits >= 8:
            bits -= 8
            out.append((buf >> bits) & 0xFF)
            # keep only remaining bits in buf to avoid huge ints
            buf &= (1 << bits) - 1 if bits > 0 else 0

    if bits > 0:
        # pad remaining bits with zeros on the right (LSB)
        out.append((buf << (8 - bits)) & 0xFF)

    return bytes(out)


def validate_sixbit_bytes(data: bytes, length: int) -> None:
    """
    Validate canonical SIXBIT bytes for a given symbol length.

    Checks:
    - byte length matches sixbit_byte_len(length)
    - unused padding bits in the final byte (if any) must be 0
    - all extracted 6-bit codes are within 0..63 (always true) (alphabet validity checked on decode)
    """
    if length < 0:
        raise ValueError("length must be >= 0")

    exp = sixbit_byte_len(length)
    if len(data) != exp:
        raise ValueError(f"SIXBIT byte length mismatch: got {len(data)} expected {exp}")

    if exp == 0:
        return

    total_bits = length * 6
    rem = total_bits % 8
    if rem != 0:
        pad_bits = 8 - rem
        last = data[-1]
        mask = (1 << pad_bits) - 1  # low pad_bits must be zero
        if (last & mask) != 0:
            raise ValueError("Non-canonical SIXBIT padding bits: expected 0s in unused LSBs")


def unpack_sixbit(data: bytes, length: int) -> str:
    """
    Decode SIXBIT bytes into a string, using `length` to ignore padding.
    """
    validate_sixbit_bytes(data, length)

    out_chars = []
    buf = 0
    bits = 0
    i = 0
    for byte in data:
        buf = (buf << 8) | byte
        bits += 8
        while bits >= 6 and i < length:
            bits -= 6
            code = (buf >> bits) & 0x3F
            out_chars.append(_SIXBIT_CODE_TO_CHAR[code])
            i += 1
        # keep only remaining bits
        buf &= (1 << bits) - 1 if bits > 0 else 0
        if i >= length:
            break

    if i != length:
        raise ValueError("SIXBIT decode underflow: not enough data for requested length")

    return "".join(out_chars)


# ---------------------------------------------------------------------
# Unified encode/decode and canonical selection
# ---------------------------------------------------------------------

def choose_canonical_encoding(seq: str) -> SeqEncoding:
    """
    Canonical encoding choice (GL1ENCv2).

    Rule order (deterministic):
      1) If normalized sequence is NON-EMPTY and pure A/C/G/T -> DNA2
      2) Else if all symbols are DNA4 alphabet -> DNA4
      3) Else if all symbols are SIXBIT ALPH64 -> SIXBIT
      4) Else -> ASCII

    The empty string is encoded as ASCII (stable convention).
    """
    s = normalize_seq(seq)
    if not s:
        return SeqEncoding.ASCII
    if is_pure_acgt(s, assume_normalized=True):
        return SeqEncoding.DNA2
    if is_dna4_eligible(s, assume_normalized=True):
        return SeqEncoding.DNA4
    if is_sixbit_eligible(s, assume_normalized=True):
        return SeqEncoding.SIXBIT
    return SeqEncoding.ASCII


def encode_seq(seq: str, encoding: SeqEncoding) -> Tuple[bytes, int]:
    """
    Encode a sequence into (bytes, length_in_symbols).

    The returned length is ALWAYS the number of characters in the normalized sequence.
    """
    s = normalize_seq(seq)

    if encoding == SeqEncoding.DNA2:
        data = pack_dna2(s)
        return data, len(s)

    if encoding == SeqEncoding.DNA4:
        data = pack_dna4(s)
        return data, len(s)

    if encoding == SeqEncoding.SIXBIT:
        data = pack_sixbit(s)
        return data, len(s)

    if encoding == SeqEncoding.ASCII:
        try:
            data = s.encode("ascii")
        except UnicodeEncodeError as e:
            raise ValueError("ASCII encoding requires ASCII-only characters after normalization") from e
        return data, len(s)

    raise ValueError(f"Unknown encoding: {encoding}")


def decode_seq(data: bytes, length: int, encoding: SeqEncoding) -> str:
    """
    Decode bytes back into the normalized sequence string.

    Canonical strictness:
    - DNA2: validate_dna2_bytes(data, length)
    - DNA4: validate_dna4_bytes(data, length)
    - SIXBIT: validate_sixbit_bytes(data, length)
    - ASCII: require len(data) == length and bytes are valid ASCII
    """
    if length < 0:
        raise ValueError("length must be >= 0")

    if encoding == SeqEncoding.DNA2:
        return unpack_dna2(data, length)

    if encoding == SeqEncoding.DNA4:
        return unpack_dna4(data, length)

    if encoding == SeqEncoding.SIXBIT:
        return unpack_sixbit(data, length)

    if encoding == SeqEncoding.ASCII:
        if len(data) != length:
            raise ValueError(f"ASCII byte length mismatch: got {len(data)} expected {length}")
        try:
            return data.decode("ascii")
        except UnicodeDecodeError as e:
            raise ValueError("ASCII decoding failed; bytes are not valid ASCII") from e

    raise ValueError(f"Unknown encoding: {encoding}")


@dataclass(frozen=True)
class EncodedSequence:
    """
    Self-contained encoded sequence.

    - encoding: SeqEncoding (uint8 on-chain)
    - length: number of symbols AFTER normalization
    - data: encoded bytes (may include canonical padding)
    """

    encoding: SeqEncoding
    length: int
    data: bytes

    def validate(self) -> None:
        """Validate canonical structure for this encoding."""
        if self.length < 0:
            raise ValueError("length must be >= 0")
        if self.encoding == SeqEncoding.DNA2:
            validate_dna2_bytes(self.data, self.length)
        elif self.encoding == SeqEncoding.DNA4:
            validate_dna4_bytes(self.data, self.length)
        elif self.encoding == SeqEncoding.SIXBIT:
            validate_sixbit_bytes(self.data, self.length)
        elif self.encoding == SeqEncoding.ASCII:
            if len(self.data) != self.length:
                raise ValueError(f"ASCII byte length mismatch: got {len(self.data)} expected {self.length}")
        else:
            raise ValueError(f"Unknown encoding: {self.encoding}")

    def decode(self) -> str:
        """Decode back into a normalized string."""
        return decode_seq(self.data, self.length, self.encoding)


def canonical_encode(seq: str) -> EncodedSequence:
    """Canonically encode a sequence using choose_canonical_encoding()."""
    enc = choose_canonical_encoding(seq)
    data, length = encode_seq(seq, enc)
    out = EncodedSequence(enc, length, data)
    out.validate()
    return out


def canonical_decode(encoded: EncodedSequence) -> str:
    """Decode an EncodedSequence back into its normalized string."""
    encoded.validate()
    return encoded.decode()


# ---------------------------------------------------------------------
# Reverse complement (strand helper)
# ---------------------------------------------------------------------

# IUPAC complement table (uppercase). Unrecognized characters remain unchanged.
_IUPAC_COMPLEMENT = {
    "-": "-",  # gap stays gap
    "*": "*",
    ".": ".",
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",  # if RNA appears
    "R": "Y",  # A/G <-> C/T
    "Y": "R",
    "S": "S",  # G/C
    "W": "W",  # A/T
    "K": "M",  # G/T <-> A/C
    "M": "K",
    "B": "V",  # C/G/T <-> A/C/G
    "D": "H",  # A/G/T <-> A/C/T
    "H": "D",
    "V": "B",
    "N": "N",
}
_COMPLEMENT_TRANS = str.maketrans(_IUPAC_COMPLEMENT)


def reverse_complement(seq: str) -> str:
    """
    Reverse-complement the sequence after normalization.

    - Normalizes input first (GL1 rules).
    - Applies IUPAC complements where defined.
    - Leaves unknown characters unchanged (but still reverses).
    """
    s = normalize_seq(seq)
    return s.translate(_COMPLEMENT_TRANS)[::-1]


# ---------------------------------------------------------------------
# Self-test / vectors
# ---------------------------------------------------------------------

def _selftest() -> None:
    # DNA2 pack/unpack vectors
    assert pack_dna2("ACGT") == bytes([0x1B])
    assert unpack_dna2(bytes([0x1B]), 4) == "ACGT"
    assert pack_dna2("T") == bytes([0xC0])
    assert unpack_dna2(bytes([0xC0]), 1) == "T"
    assert pack_dna2("ACGTAC") == bytes([0x1B, 0x10])
    assert unpack_dna2(bytes([0x1B, 0x10]), 6) == "ACGTAC"

    # DNA4 vectors
    assert pack_dna4("A") == bytes([0x10])  # A then pad(0)
    assert unpack_dna4(bytes([0x10]), 1) == "A"
    assert pack_dna4("ACGT") == bytes([0x12, 0x48])
    assert unpack_dna4(bytes([0x12, 0x48]), 4) == "ACGT"
    assert pack_dna4("N-") == bytes([0xF0])
    assert unpack_dna4(bytes([0xF0]), 2) == "N-"

    # SIXBIT vectors
    assert pack_sixbit("ABCD") == bytes([0x00, 0x10, 0x83])
    assert unpack_sixbit(bytes([0x00, 0x10, 0x83]), 4) == "ABCD"

    # canonical selection
    assert choose_canonical_encoding("ACGT") == SeqEncoding.DNA2
    assert choose_canonical_encoding("ACGN") == SeqEncoding.DNA4
    assert choose_canonical_encoding("MKWVTFISLL") == SeqEncoding.SIXBIT
    assert choose_canonical_encoding("") == SeqEncoding.ASCII

    # reverse complement
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAGT") == "ACTT"

    print("gl1enc selftest: OK")


if __name__ == "__main__":
    _selftest()
