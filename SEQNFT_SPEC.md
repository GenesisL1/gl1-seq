# SeqNFT + GL1 Unified Pipeline Specification (v2, GL1ENCv2)

This document specifies a **single pipeline** that runs with identical semantics in:

- **Off-chain mode** (local files, object storage, IPFS, cloud)
- **On-chain mode** (smart contracts storing bytes or anchoring content-addressed pointers)

It is built on top of the **GL1 format family** and **GL1ENCv2** codec engine:

- `.gl1`  — **GL1F v1**: chunked payload container
- `.gl1x` — **GL1X v2**: random access + metadata (mask, exceptions, hashes, Merkle)
- `.gl1q` — **GL1Q v2**: FASTQ qualities with TOC (optional)
- `.gl1i` — **GL1I v1**: minimizer index with buckets (optional)

The key additional concept is **SeqNFT objects + commits**, which make GL1 **history-aware**
(git-like) and **NFT-native**.

---

## 0. What changed from v1

v2 updates the pipeline to match **`gl1enc_v2.py`** (GL1ENCv2):

- Codec IDs are aligned with `SeqEncoding`:
  - `0 = DNA2` (2-bit A/C/G/T)
  - `1 = ASCII` (byte-per-symbol fallback)
  - `2 = DNA4` (4-bit IUPAC nucleotides + `-` gap)
  - `3 = SIXBIT` (6-bit, 64-symbol alphabet)

- **NUC2 “exceptions” remain supported**, but they are now defined as:
  - payload uses **DNA2** over a base sequence where non-ACGT positions are replaced (e.g. `A`)
  - exceptions store the *true* symbol using **DNA4 nibble codes** (0..15)

- The old `nuc4_table=bam/u` distinction is removed in v2 (GL1ENCv2’s DNA4 mapping includes `-` as code `0` and does **not** encode `U` separately).
  - If your data contains `U` (RNA), you may:
    - use ASCII codec, or
    - pre-normalize `U→T` before minting (this changes semantics/hashes and must be explicit in your workflow).

---

## 1. Goals

1) **Identical semantics on-chain and off-chain**  
   Same hashes, same roots, same record/chunk layout, same region extraction output.

2) **Future-proof “git-like” history for sequences**  
   Commits form a DAG; new commits reuse unchanged chunks; only modified chunks become new objects/blobs.

3) **NFT-native hierarchy**  
   - Chromosome/Record = parent NFT  
   - Chunk = child NFT  
   - Assembly/Commit = root NFTs (versioned dataset)

4) **Two deployment modes**
   - **Fully on-chain**: chunk bytes are stored directly in contracts (or in a data contract).
   - **Anchor/IPFS**: contracts store only `(cid, sha256, byte_length)`; clients fetch bytes and verify.

5) Ready for:
   - scalable region queries (“view”)
   - scalable search (“minimizer index”)
   - deterministic derivations (variants → new commits)
   - alignment-ready expansion (future objects referencing blocks + indexes)

---

## 2. GL1ENCv2 (codec engine)

### 2.1 Encodings

GL1ENCv2 defines four canonical encodings:

| Encoding | ID | Bits/symbol | Intended use |
|---|---:|---:|---|
| DNA2 | 0 | 2 | `A,C,G,T` only |
| ASCII | 1 | 8 | fallback, preserves any ASCII symbol |
| DNA4 | 2 | 4 | IUPAC nucleotides + `-` gap |
| SIXBIT | 3 | 6 | protein alphabet + extra symbols |

### 2.2 DNA4 mapping (nibble → char)

DNA4 is defined by a 4-bit mask (IUPAC bitmask) with `-` gap in code 0:

```
0  '-'   (gap)
1  'A'
2  'C'
3  'M'   (A|C)
4  'G'
5  'R'   (A|G)
6  'S'   (C|G)
7  'V'   (A|C|G)
8  'T'
9  'W'   (A|T)
10 'Y'   (C|T)
11 'H'   (A|C|T)
12 'K'   (G|T)
13 'D'   (A|G|T)
14 'B'   (C|G|T)
15 'N'   (A|C|G|T)
```

This is exactly the mapping in `gl1enc_v2.py`.

---

## 3. GL1F v1 — `.gl1` main container

GL1F stores a sequence dataset as:

- file header
- N records
- each record contains:
  - record header + name
  - chunk list: each chunk is `CHUNK_HDR + payload`

### 3.1 Binary layout

All integer fields are **little-endian**.

#### File header (`FILE_HDR_FMT = "<4sHHII"`)

| Field | Type | Meaning |
|---|---|---|
| magic | 4s | `b"GL1F"` |
| version | u16 | `1` |
| flags | u16 | reserved (0) |
| record_count | u32 | number of records |
| reserved | u32 | reserved (0) |

#### Record header (`REC_HDR_FMT = "<4sHBBHIIHHI"`)

| Field | Type | Meaning |
|---|---|---|
| magic | 4s | `b"GL1R"` |
| version | u16 | `1` |
| kind | u8 | `1=nucleotide`, `2=protein` |
| codec | u8 | `0 DNA2`, `1 ASCII`, `2 DNA4`, `3 SIXBIT` |
| flags | u16 | see §3.3 |
| length | u32 | total symbols in record |
| chunk_symbols | u32 | symbols per chunk (last chunk may be shorter) |
| name_len | u16 | bytes of UTF-8 name |
| reserved | u16 | reserved |
| chunk_count | u32 | number of chunks |
| name bytes | name_len | UTF-8 |

#### Chunk header (`CHUNK_HDR_FMT = "<II"`)

| Field | Type | Meaning |
|---|---|---|
| sym_count | u32 | number of symbols in this chunk |
| payload_len | u32 | bytes in payload |
| payload | bytes | codec-dependent |

### 3.2 Payload semantics

The chunk payload is encoded by `codec`:

- `codec=DNA2`: GL1ENCv2 `pack_dna2` over an **ACGT-only base sequence**.
  - If the record uses exceptions (see §4.3), non-ACGT symbols are restored by the exception list in GL1X / SeqNFT chunk objects.
- `codec=DNA4`: GL1ENCv2 `pack_dna4` over IUPAC+`-`.
- `codec=SIXBIT`: GL1ENCv2 `pack_sixbit` over SIXBIT alphabet.
- `codec=ASCII`: raw uppercase ASCII bytes (`len(payload)==sym_count`).

### 3.3 Record flags

Flags are a u16 bitfield:

- `0x0001 FLAG_HAS_MASK` — record has mask intervals (lowercase runs)
- `0x0002 FLAG_HAS_EXCP` — record uses exceptions (meaningful with `codec=DNA2`)
- `0x0004 FLAG_HAS_QUAL` — qualities exist in `.gl1q`

---

## 4. GL1X v2 — `.gl1x` random access + metadata

GL1X adds:

- record offsets for random access into `.gl1`
- per-chunk offsets
- mask intervals (absolute coordinates)
- exception points (absolute coordinates)
- per-chunk `sha256(payload)` hashes
- Merkle root over chunk hashes

### 4.1 Header (`X_HDR_FMT = "<4sHII"`)

| Field | Type | Meaning |
|---|---|---|
| magic | 4s | `b"GL1X"` |
| version | u16 | `2` |
| record_count | u32 | number of records |
| reserved | u32 | reserved |

### 4.2 Record blocks

For each record:

| Field | Type |
|---|---|
| name_len | u16 |
| name | bytes |
| record_off | u64 |
| kind | u8 |
| codec | u8 |
| flags | u16 |
| length | u32 |
| chunk_symbols | u32 |
| chunk_count | u32 |
| mask_count | u32 |
| excp_count | u32 |
| chunk_offs | chunk_count × u64 |
| mask_intervals | mask_count × (u32 start, u32 len) |
| exceptions | excp_count × (u32 pos, u8 code) |
| chunk_hashes | chunk_count × 32 bytes |
| merkle_root | 32 bytes |

### 4.3 Exceptions (DNA2 + sparse DNA4 codes)

When a nucleotide record is stored as `codec=DNA2` but contains non-ACGT symbols:

- payload stores DNA2 of a base sequence where the non-ACGT symbols are replaced (typically `A`)
- `exceptions` restores the true symbol at those positions:
  - `code` is a **DNA4 nibble** (0..15) from §2.2
  - decode applies exceptions after unpacking DNA2

This keeps “mostly ACGT” genomes near 2-bit size without losing ambiguity.

---

## 5. GL1Q v2 — `.gl1q` qualities (optional)

GL1Q stores FASTQ quality strings per record, chunked identically to `.gl1`.

### 5.1 Header (`Q_HDR_FMT = "<4sHHQII"`)

| Field | Type | Meaning |
|---|---|---|
| magic | 4s | `b"GL1Q"` |
| version | u16 | `2` |
| reserved | u16 | 0 |
| toc_off | u64 | file offset of TOC block |
| record_count | u32 | records with qualities |
| reserved2 | u32 | 0 |

### 5.2 Record block

At offset `toc[name]`:

| Field | Type |
|---|---|
| name_len | u16 |
| name | bytes |
| meta | `<III>` = total_len, chunk_symbols, chunk_count |
| chunks | chunk_count × (Q_CHUNK_HDR + zbytes) |

Chunk header (`Q_CHUNK_HDR_FMT="<II"`):

| Field | Type |
|---|---|
| sym_count | u32 |
| zlen | u32 |
| zbytes | zlen bytes |

The uncompressed bytes are raw FASTQ qualities for that chunk.

### 5.3 TOC

At `toc_off`:

Header: `QTOC_HDR_FMT="<4sHI"` with magic `b"QTOC"`, version `1`, record_count.

Then entries:

- name_len (u16), name bytes, offset (u64)

---

## 6. GL1I v1 — `.gl1i` minimizer index (optional)

GL1I stores **minimizer postings** for fast seed hits.

This reference pipeline builds GL1I for **nucleotide records only** by default:

- only kmers containing A/C/G/T are indexed (ambiguous kmers are skipped)

### 6.1 Header (`I_HDR_FMT="<4sHBBHIII"`)

| Field | Type | Meaning |
|---|---|---|
| magic | 4s | `b"GL1I"` |
| version | u16 | `1` |
| k | u8 | k-mer length |
| w | u8 | minimizer window length (in kmers) |
| flags | u16 | bits: canonical / nuc / prot |
| bucket_count | u32 | number of buckets |
| record_count | u32 | number of records in table |
| reserved | u32 | 0 |

Flags:

- `I_FLAG_CANON (0x1)` canonical (reverse-complement min) for nucleotide keys
- `I_FLAG_NUC (0x2)` nucleotide keys (2-bit)
- `I_FLAG_PROT (0x4)` protein keys (6-bit) — reserved for future

### 6.2 Records table

`record_count` entries:

- name_len (u16), name bytes, length (u32)

### 6.3 Bucket directory

`bucket_count` entries (`I_BUCKET_DIR_FMT="<QII"`):

- offset (u64) — where bucket block begins
- key_count (u32)
- reserved (u32)

### 6.4 Bucket block

At `offset`:

- key_count (u32)
- repeated `key_count` times:
  - key (u64)
  - posting_count (u32)
  - postings: posting_count × (rec_id u32, pos u32)

Bucket is chosen as `bucket = key % bucket_count`.

---

## 7. SeqNFT object model (content-addressed)

All objects are immutable and addressed by `sha256` of a **canonical encoding** of the object core.

### 7.1 Blob

Raw bytes addressed by `sha256(bytes)`.

Used for:

- GL1 chunk blobs (`CHUNK_HDR + payload`)
- compressed quality chunk bytes (`Q_CHUNK_HDR + zbytes`)
- `.gl1i` bytes

### 7.2 Chunk object (`type="chunk"`)

Represents a single chunk + chunk-local metadata.

Core fields:

- `chunk_blob`: blob id for exact GL1 bytes `CHUNK_HDR + payload`
- `payload_sha256`: sha256(payload) (leaf for seq Merkle root)
- `mask`: list of `(start,len)` **within chunk**
- `excp`: list of `(pos, code)` **within chunk** (DNA4 nibble codes)
- `qual_chunk` (optional): blob id for `Q_CHUNK_HDR + zbytes`

### 7.3 Chromosome/Record object (`type="chrom"`)

Record manifest = **parent NFT**.

Core fields:

- `name`, `kind`, `codec`, `flags`, `length`, `chunk_symbols`
- `chunks`: ordered list of chunk object ids
- integrity roots:
  - `seq_payload_root`: Merkle over `sha256(payload)` leaves (GL1X-compatible)
  - `mask_root`, `excp_root`, `qual_root`: Merkle over per-chunk meta leaves (pipeline-level)

### 7.4 Assembly object (`type="assembly"`)

Dataset manifest.

Core fields:

- `assembly`: name/version string
- `chromosomes`: list of chromosome object ids
- `assembly_root`: Merkle over chromosome ids (pipeline-level)

### 7.5 Commit object (`type="commit"`)

Git-like commit node; can be minted as a “Version NFT”.

Core fields:

- `tree`: assembly id
- `parents`: list of parent commit ids
- `author`, `message`, `time`
- `tags`: structured metadata (pipeline params, input filename, etc.)

---

## 8. Canonical hashing

### 8.1 Canonical JSON

Canonical object encoding (v2) is JSON with:

- sorted keys
- separators `(",", ":")`
- UTF-8

Then:

```
object_id = sha256( b"seqnft:" + type + b"\0" + canonical_json(core) )
```

### 8.2 Merkle

- leaf hashes are 32-byte sha256 digests
- parent = sha256(left || right)
- if odd count, duplicate last (“dup-last”)

---

## 9. NFT mapping (on-chain)

### 9.1 Token IDs

Token IDs SHOULD be the 256-bit value of the object id (stored as `bytes32` or `uint256`).

### 9.2 Hierarchy

- ChromosomeNFT token id == chromosome object id
- ChunkNFT token id == chunk object id
- ChromosomeNFT “children” are `chunks[]`

### 9.3 Deployment modes

- **onchain**: ChunkNFT stores bytes directly (or stores blob id that resolves onchain).
- **ipfs/anchor**: ChunkNFT stores `(cid, sha256, byte_length)`.

In both cases, `sha256` is authoritative.

---

## 10. Evolution via commits

Because record manifests are lists of chunk ids:

- mutations affecting only a region rewrite only overlapping chunks
- all other chunk ids are reused
- a new chromosome id and assembly id result (hash changes), but storage is incremental

Each evolutionary step is a commit; history is replayable.

---

## 11. Summary

This spec makes GL1 into a **history-aware, NFT-native sequence system**:

- GL1ENCv2 remains the payload codec layer.
- GL1F/GL1X/GL1Q/GL1I provide file compatibility and integrity.
- SeqNFT objects add Git-like commits, on-chain deployment, and explicit integrity roots.
