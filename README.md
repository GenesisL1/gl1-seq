# GL1ENCv2 and SeqNFT pipeline — canonical, on-chain ready genomics

This repository pairs two production-oriented references:

- **`gl1enc_v2.py`** — the canonical GL1 sequence codec (2/4/6‑bit + ASCII) with deterministic normalization, padding rules, and validation.
- **`seqnft_pipeline.py`** — a dual-mode SeqNFT simulation pipeline that uses GL1 encodings to build content-addressed, Merkleized payloads that stay identical on-chain and off-chain.

Both modules emphasize scientific reproducibility (explicit byte layouts, fixed alphabets, deterministic hashing) and on-chain practicality (uint8 codec IDs, canonical padding, and length-aware payloads).

---

## Repository components

- **`gl1enc_v2.py`** — reference implementation of GL1ENCv2 with DNA2/DNA4/SIXBIT/ASCII encoders, strict validators, and reverse-complement utilities.
- **`seqnft_pipeline.py`** — pipeline that constructs GL1F/GL1X/GL1Q/GL1I containers, hashes objects with deterministic JSON canonicalization, and supports both “fully on-chain bytes” and “anchor/IPFS” deployments.
- **`GL1ENCv2_SPEC.md`** and **`SEQNFT_SPEC.md`** — narrative specifications complementing the code.

---

## GL1ENCv2 (gl1enc_v2.py) in detail

GL1ENCv2 is the normalization and encoding engine used by all GL1 payloads. It supplies both the Python API and the wire-format identifiers required for Solidity parity.

### Canonical normalization

- Strip leading/trailing whitespace, uppercase, and remove ASCII spaces, tabs, and newlines to obtain a deterministic `normalized_seq`.
- Normalization is deliberately minimal so the same rules can be replicated in Solidity without locale surprises.

### Encodings and IDs

| Encoding | ID (uint8) | Bits / symbol | Intended payloads |
|---|---:|---:|---|
| **DNA2** | 0 | 2 | Strict A/C/G/T sequences with 4 bases per byte (00/01/10/11 mapping). |
| **ASCII** | 1 | 8 | Future-proof fallback for any normalized ASCII string. |
| **DNA4** | 2 | 4 | IUPAC nucleotide bitmasks + gap, with canonical padding in the low nibble. |
| **SIXBIT** | 3 | 6 | Fixed 64-character ALPH64 alphabet for proteins and symbolic strings. |

Canonical selection prefers the most compact eligible encoding: non-empty A/C/G/T → DNA2; otherwise DNA4 → SIXBIT → ASCII; empty strings default to ASCII.

### Deterministic packing and validation

- **DNA2** packs four bases per byte; unused slots pad with `00`. Validators enforce expected length and zeroed padding.
- **DNA4** packs two symbols per byte with a required `0x0` pad nibble for odd lengths.
- **SIXBIT** streams 6-bit codes MSB-first and requires trailing zero bits; alphabet indices are fixed by `ALPH64` order.
- **ASCII** simply requires byte length to equal normalized length and ASCII validity.

All encoders return `(bytes, length_in_symbols)`; the `EncodedSequence` dataclass stores `encoding`, `length`, and `data`, and can self-validate or decode. The helper `reverse_complement` operates on normalized input and supports the full IUPAC complement set.

### Minimal usage example

```python
from gl1enc_v2 import canonical_encode, canonical_decode

enc = canonical_encode("acgtN-")      # chooses DNA4, pads low nibble if needed
seq = canonical_decode(enc)            # -> "ACGTN-"
```

This workflow is safe for on-chain use because the codec IDs, padding rules, and byte lengths are deterministic and validated before decoding.

---

## SeqNFT dual-mode pipeline (seqnft_pipeline.py)

The pipeline emulates the on-chain SeqNFT system locally while preserving byte-for-byte parity with Solidity-friendly layouts. It binds GL1ENCv2 to higher-level objects and container files.

### Data model and containers

- **GL1F** (main) files chunk encoded sequences into records with explicit codec IDs (`DNA2/ASCII/DNA4/SIXBIT`) and symbol lengths.
- **GL1X** sidecars (v2) add random access tables, chunk hashes, and Merkle roots for verification without re-reading payloads.
- **GL1Q** (v2) stores quality strings with its own table of contents; **GL1I** (v1) stores minimizer indices with per-bucket directories and k-mer flags.
- Record flags capture lowercase masks, exception masks for DNA2 payloads, and attached qualities.

### Content addressing and determinism

- Uses canonical JSON (sorted keys, compact separators) plus SHA-256 to derive IDs for assemblies, records, and chunks; the same logic can be swapped to keccak256 for Solidity without changing canonical bytes.
- Merkle roots duplicate the last leaf when needed, matching common on-chain tree semantics.
- All payloads depend on `gl1enc_v2` for canonical encoding; length fields and padding checks prevent ambiguous encodings.

### Deployment modes

- **Fully on-chain**: store GL1F/GL1X bytes directly in the “chain simulator,” retaining complete payloads and hashes.
- **Anchor/IPFS**: store only CID-like pointers with SHA-256 anchors, keeping GL1 hashes verifiable while moving bulk payloads off-chain.

These modes let the same chunked objects travel between L1 contracts, L2 blobs, IPFS, or conventional storage without breaking equivalence.

---

## Why GL1 for on-chain genomics?

- **Canonical byte commitments**: deterministic normalization, explicit padding rules, and fixed alphabets remove alternate encodings that could break Merkle proofs or smart-contract verification.
- **Bit-efficient encodings**: 2/4/6-bit packing keeps calldata, blobs, and state writes small while retaining exact symbolic sequences.
- **Length-coupled metadata**: every payload carries symbol counts alongside bytes, ensuring safe slicing, chunking, and reverse-complement operations on-chain.
- **Composable containers**: GL1F/GL1X/GL1Q/GL1I bundle payloads, quality scores, and indices while exposing hashes suitable for commitments, proofs, or bridge attestations.
- **Cross-environment parity**: the same code paths run locally and in Solidity (ID compatibility and padding checks), reducing the surface for consensus bugs.

---

## Comparison with common genomics formats (on-chain focus)

| Format / container | Scope | Determinism & canonicalization | Compression / packing | Mutation & chunking | On-chain readiness | GL1 advantage |
|---|---|---|---|---|---|---|
| **GL1 (GL1ENCv2 + GL1F/GL1X/GL1Q/GL1I)** | Encoded sequences, qualities, indices with Merkle roots | Full normalization + fixed alphabets; validators reject non-canonical padding | 2/4/6-bit codecs with explicit padding; optional gzip for sidecars | Native chunk headers + per-chunk hashes; reverse-complement helper in core codec | Codec IDs are uint8; lengths stored explicitly; SHA-256/keccak-ready; IPFS or on-chain bytes interchangeable | Purpose-built for deterministic proofs, low calldata, and safe smart-contract decoding |
| FASTA | Plain sequences | No canonical whitespace rules; mixed casing common | ASCII only; whitespace variable | No standard chunking | Needs external hashing; ambiguous whitespace harms commitments | GL1 enforces whitespace stripping and compact encodings |
| FASTQ | Sequences + qualities | Header and line wrapping not standardized | ASCII bases + ASCII qualities; size scales linearly | Record boundaries implicit, not chunked | Requires preprocessing to normalize; large calldata | GL1Q pairs qualities with deterministic TOC and chunk sizes |
| BAM/CRAM | Aligned reads | Header order and auxiliary tags can vary; CRAM codec versions shift | Binary with compression; CRAM depends on reference | Record blocks compress contextually; edits require recompression | Hard to validate deterministically on-chain; heavy dependencies | GL1F chunking + fixed-bit codecs are simpler to verify and refit on-chain |
| VCF/BCF | Variant calls | Header order and INFO ordering vary; normalization is workflow-specific | Text (VCF) or binary (BCF) with variable INFO/FORMAT | Record streaming but no Merkleized chunk hashes | Requires bespoke schemas on-chain; INFO parsing costly | GL1 encodes only canonical sequences + masks, leaving variant semantics to higher layers |
| GFA | Assembly graphs | Depends on line ordering and path notation | Textual; no fixed packing | Graph updates change many lines | Not byte-stable for commitments | GL1 separates immutable sequence payloads from graph semantics |
| Generic IPFS blobs | Arbitrary content | Determinism depends on producer | Depends on content | No semantic chunking | On-chain contracts cannot interpret payload safely | GL1 defines interpretable bytes plus chunk-level hashes, making IPFS anchors verifiable |

---

## Versioning and Solidity alignment

- `GL1ENC_FORMAT = "GL1ENCv2"` identifies the wire format; the Python module is versioned with `__version__ = "2.0.0"`. Any change to normalization, alphabets, packing, or canonical selection must bump the format string.
- Encoding IDs are stable (`0=DNA2`, `1=ASCII`, `2=DNA4`, `3=SIXBIT`) and are intended to be stored as uint8 on-chain.
- Solidity ports should avoid normalizing user-provided strings; expect already-normalized bytes plus a length for DNA2/DNA4/SIXBIT, and enforce zeroed padding in final bytes/nibbles.

---

## Key reference vectors

- DNA2: `ACGT` → `0x1B` (len 4); `T` → `0xC0` (len 1); `ACGTAC` → `0x1B10` (len 6).
- DNA4: `A` → `0x10` (len 1); `ACGT` → `0x1248` (len 4); `N-` → `0xF0` (len 2).
- SIXBIT: `ABCD` → `0x00 0x10 0x83` (len 4) using ALPH64 indices.
- Canonical choice: `ACGT` → DNA2; `ACGN` → DNA4; `MKWVTFISLL` → SIXBIT; empty string → ASCII.
- Reverse complement: `AAGT` → `ACTT` after normalization.

These vectors must remain consistent across Python, Solidity, and any storage layer that uses GL1 commitments.
