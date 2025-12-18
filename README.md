# gl1enc — ultra-light 2-bit / 4-bit / 6-bit sequence encoder/decoder

`gl1enc.py` is a **single-file**, **zero-dependency** Python tool for encoding and decoding biological sequences into a compact **`.gl1` container**.

It supports:

- **NUC2 (2-bit nucleotides):** `A C G T` (optionally accepts `U` and stores it as `T` unless strict)
- **NUC4 (4-bit nucleotides, SAM/BAM NT16):** `= A C M G R S V T W Y H K D B N` (IUPAC ambiguity codes)
- **NUC4_U (4-bit nucleotides with explicit U):** preserves RNA `U` by using code `0` for `U`
- **PROT6 (6-bit proteins):** 20 amino acids + `* X B Z J U O -`

This is intended as a practical **binary payload layer** that is fast to encode/decode and easy to integrate into indexing, storage, chunk/NFT pipelines, and Web2/Web3 distribution.

---

## Table of contents

- [Benefits](#benefits)
  - [Compact storage](#compact-storage)
  - [Fast IO + large sequence support](#fast-io--large-sequence-support)
  - [Interoperability-friendly](#interoperability-friendly)
  - [Practical workflows](#practical-workflows)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input formats](#input-formats)
  - [`--format auto` (default)](#--format-auto-default)
  - [FASTA](#fasta)
  - [FASTQ](#fastq)
  - [RAW](#raw)
  - [Gzip input](#gzip-input)
- [Output: what is a `.gl1` file?](#output-what-is-a-gl1-file)
- [Encoding mode selection (AUTO)](#encoding-mode-selection-auto)
  - [Protein vs nucleotide inference](#protein-vs-nucleotide-inference)
  - [Nucleotide codec in AUTO](#nucleotide-codec-in-auto)
- [Command reference](#command-reference)
  - [`encode`](#1-encode)
  - [`decode`](#2-decode)
  - [`info`](#3-info)
- [Examples (realistic)](#examples-realistic)
- [Performance tips](#performance-tips)
- [Limitations (current)](#limitations-current)
- [File naming](#file-naming)
- [License (Public Domain)](#license-public-domain)

---

## Benefits

### Compact storage

- **2-bit nucleotides:** ~4× smaller than ASCII sequence text  
- **4-bit nucleotides:** ~2× smaller than ASCII while supporting ambiguity  
- **6-bit proteins:** ~25% smaller than ASCII proteins  

### Fast IO + large sequence support

- Chunking keeps memory bounded even for whole chromosomes.
- Decoding is LUT-based for nucleotides (very fast).

### Interoperability-friendly

- 4-bit nucleotide mode mirrors the **SAM/BAM NT16 symbol ordering** (`=ACMGRSVTWYHKDBN`), so it’s easy to map to/from HTS tooling internals.
- Output can be decoded back to normal FASTA anytime.

### Practical workflows

- One `.gl1` file can store **multiple sequences** (multi-FASTA / FASTQ reads).
- The tool can **decode a single file**, **many files**, or **an entire folder** (recursive optional).

---

## Installation

No installation needed. Just copy the script.

```bash
chmod +x gl1enc.py
./gl1enc.py --help
```

**Requirements:** Python 3.8+ (stdlib only)

---

## Quick start

### Encode a large FASTA (recommended chunking)

```bash
python3 gl1enc.py encode -i chr1.fa -o chr1.gl1 --chunk 1000000
python3 gl1enc.py info -i chr1.gl1
```

### Decode back to FASTA

```bash
python3 gl1enc.py decode -i chr1.gl1 -o chr1_out.fa --to fasta --wrap 60
```

### Decode a whole folder (one FASTA per `.gl1`)

```bash
python3 gl1enc.py decode -i ./gl1_folder -o ./decoded --to fasta --wrap 60
```

### Decode a folder recursively

```bash
python3 gl1enc.py decode -i ./gl1_folder -o ./decoded --recursive
```

### Combine many `.gl1` files into one FASTA (avoid name collisions)

```bash
python3 gl1enc.py decode -i ./gl1_folder --recursive -o all.fa --prefix-file
```

---

## Input formats

### `--format auto` (default)

Auto-detects by the first meaningful line:

- `>` → **FASTA**
- `@` → **FASTQ**
- otherwise → **RAW**

### FASTA

Multi-record FASTA becomes multiple records inside one `.gl1` container.

Example:

```text
>seq1
ACGTACGT...
>seq2
ACGTNRY...
```

➡️ `out.gl1` contains 2 records: `seq1`, `seq2`

### FASTQ

Each read becomes a record; qualities are skipped (sequence only).

### RAW

Whole file becomes a single record named `seq1` (whitespace removed).

### Gzip input

If input ends with `.gz`, it’s read transparently.

---

## Output: what is a `.gl1` file?

A `.gl1` is a container that can hold one or many records (sequences). Each record has:

- name
- kind: nucleotide or protein
- codec: NUC2 / NUC4_BAM / NUC4_U / PROT6
- length
- chunk size and chunk count
- chunk payload bytes

Chunking is internal. Decoding automatically concatenates chunks back into the original sequence.

---

## Encoding mode selection (AUTO)

### Protein vs nucleotide inference

Heuristic:

- If the record contains characters that strongly imply protein (e.g. `E,F,I,L,P,Q,O,J,Z,X,*`) → protein
- Otherwise, if it fits nucleotide + IUPAC alphabet → nucleotide

If you have edge cases, force the mode with `--encode`.

### Nucleotide codec in AUTO

- Only `A,C,G,T` → **NUC2**
- If ambiguity is present (`N,R,Y,S,W,K,M,B,D,H,V,=`) → **NUC4**
- If `U` is present:
  - AUTO prefers **NUC4_U** (preserve `U`) unless you force `--nuc4-table bam`

---

## Command reference

### 1) `encode`

Encode FASTA/FASTQ/RAW into `.gl1`.

```bash
python3 gl1enc.py encode -i INPUT -o OUTPUT.gl1 [options]
```

#### Flags

- `-i, --in PATH`  
  Input file path (FASTA/FASTQ/RAW). `.gz` supported.  
  Use `-` for stdin (spooled to a temp file because encoding is two-pass in auto mode).

- `-o, --out PATH`  
  Output `.gl1` file.

- `--format {auto,fasta,fastq,raw}`  
  Force input parser.

- `--encode {auto,2bit,4bit,6bit}`  
  Force encoding mode:
  - `2bit` → NUC2 (only `A,C,G,T,U`)
  - `4bit` → NUC4 (IUPAC; uses `--nuc4-table`)
  - `6bit` → PROT6
  - `auto` → per-record decision

- `--chunk N`  
  Chunk size in symbols (bases/residues).  
  `0` = unchunked (one chunk per record; can use lots of RAM for big sequences)

  Recommended:
  - whole chromosome: `--chunk 1000000` to `--chunk 4000000`
  - protein sets: `--chunk 200000`

- `--nuc4-table {auto,bam,u}`  
  Controls 4-bit nucleotide mapping:
  - `bam`: SAM/BAM NT16 (`=ACMGRSVTWYHKDBN`). U is mapped to T unless `--strict`.
  - `u`: preserve U using code 0.
  - `auto`: choose `u` if U appears, else `bam`.

- `--strict`  
  Strict validation (fail on anything not representable). Examples:
  - `--encode 2bit --strict` fails if U appears
  - `--encode 4bit --nuc4-table bam --strict` fails if U appears

- `--unknown-to-x`  
  Protein mode: map unknown letters to X instead of error.  
  Also relaxes `--encode 6bit` validation to accept any A–Z (mapped to X if unknown).

---

### 2) `decode`

Decode `.gl1` to FASTA or raw. Supports single file, multiple files, and folders.

```bash
python3 gl1enc.py decode -i INPUT [INPUT2 ...] -o OUT [options]
```

#### Inputs

`-i/--in` accepts:

- one or more `.gl1` files
- directories containing `.gl1` files

#### Outputs

You have two output patterns:

**A) Combine all decoded output into one file (or stdout)**

```bash
python3 gl1enc.py decode -i a.gl1 b.gl1 -o combined.fa --prefix-file
```

**B) One output per input `.gl1`**

```bash
python3 gl1enc.py decode -i ./folder -o ./decoded
# or explicitly:
python3 gl1enc.py decode -i ./folder --outdir ./decoded
```

#### Flags

- `-i, --in ...`  
  List of files/folders.

- `-o, --out PATH`  
  `-` = stdout (default).  
  If you pass a directory path, it writes one output per input file.

- `--outdir DIR`  
  Explicit directory for per-file outputs (overrides `--out`).

- `--recursive`  
  If inputs include directories, scan subfolders too.

- `--to {fasta,raw}`  
  Output format:
  - `fasta` = headers + sequence
  - `raw` = sequence only (no headers)

- `--wrap N`  
  FASTA line wrap length (default 60).  
  `0` = no wrap.

- `--record NAME`  
  Decode only the record with this exact name (applies per input `.gl1` file).

- `--prefix-file`  
  When combining outputs from multiple input files, prefix record names with `<file>|` to avoid collisions.

---

### 3) `info`

Summarize `.gl1` file(s). Supports folders.

```bash
python3 gl1enc.py info -i file.gl1
python3 gl1enc.py info -i ./folder --recursive
```

Shows per record:

- kind (NUC/PROT)
- codec (NUC2 / NUC4_BAM / NUC4_U / PROT6)
- length
- chunk size + chunk count
- flags

---

## Examples (realistic)

### Multi-FASTA: two sequences, different codecs inside `one.gl1`

If `seq1` is pure A/C/G/T and `seq2` has ambiguity codes:

```bash
python3 gl1enc.py encode -i two_seqs.fa -o two_seqs.gl1 --chunk 500
python3 gl1enc.py info -i two_seqs.gl1
```

Expected:

- `seq1` → NUC2, chunks = 10 (for 5000 bases, chunk 500)
- `seq2` → NUC4_*, chunks = 10

### Force 2-bit (only A/C/G/T/U allowed)

```bash
python3 gl1enc.py encode -i clean.fa -o clean.gl1 --encode 2bit --chunk 1000000
```

### Force 4-bit BAM table (U→T unless strict)

```bash
python3 gl1enc.py encode -i asm.fa -o asm.gl1 --encode 4bit --nuc4-table bam
```

### Preserve U in RNA explicitly

```bash
python3 gl1enc.py encode -i rna.fa -o rna.gl1 --encode 4bit --nuc4-table u
```

### Protein one

```bash
python3 gl1enc.py encode -i proteins.faa -o proteins.gl1 --encode 6bit --chunk 200000
```

---

## Performance tips

- For chromosome-scale data: always use `--chunk` (e.g., 1,000,000).
- Avoid `--chunk 0` for large sequences (it buffers the entire record).
- AUTO mode is two-pass (scan then encode). If you need single-pass streaming, you must force `--encode` and we can add a streaming-only path.

---

## Limitations (current)

- Lowercase masking is not preserved (presence is recorded in flags only).
- FASTQ qualities are not stored (sequence only).
- Auto protein/nucleotide detection is heuristic; for edge cases use `--encode`.

---

## File naming

Filenames do not need to encode codec info because `info` is the source of truth.

If you enforce a single codec per file (forced mode), optional human hints are fine:

- `chr1__nuc2.gl1`
- `asm__nuc4bam.gl1`
- `transcripts__nuc4u.gl1`
- `uniprot__prot6.gl1`

---

## License (Public Domain)

This project is dedicated to the public domain by Bioinformatics Company, using **The Unlicense**.

See `LICENSE` for the full text.
