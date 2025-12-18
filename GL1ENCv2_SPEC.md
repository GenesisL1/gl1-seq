# GL1ENC — GL1 Sequence Codec Engine Specification (GL1ENCv2)

This document specifies **GL1ENCv2**, the canonical encoding used by GL1 to convert
sequence strings into bytes for:

- hashing (Blob IDs, Merkle leaves, commit roots),
- storage (on‑chain bytes or IPFS payload bytes),
- strict on‑chain/off‑chain parity (Solidity must match these rules exactly).

This spec corresponds to **`gl1enc.py` module version 2.0.0**.

---

## 1) Scope

GL1ENC defines how a sequence string (DNA/RNA/protein/symbolic sequences) is:

1. **normalized** into a canonical uppercase, whitespace-free string, and
2. **encoded** into bytes using one of these encodings:

| Encoding | ID (uint8) | Bits / symbol | Typical use |
|---|---:|---:|---|
| DNA2 | 0 | 2 | ultra‑compact A/C/G/T genomes |
| DNA4 | 2 | 4 | IUPAC nucleotide ambiguity codes (fast bitmask semantics) |
| SIXBIT | 3 | 6 | proteins and other symbolic strings within a fixed 64‑char alphabet |
| ASCII | 1 | 8 | future‑proof fallback for any normalized ASCII |

> GL1ENC does **not** define hashing (keccak256) itself. It only defines the exact byte representation.

---

## 2) Versioning and identifiers

- Encoding wire-format identifier: `GL1ENC_FORMAT = "GL1ENCv2"`
- Python module semantic version: `__version__ = "2.0.0"`

Any breaking change to **normalization**, **canonical selection order**, **alphabets**, or **bit packing** requires bumping `GL1ENC_FORMAT`.

---

## 3) Normalization

Before any encoding, the input sequence is normalized.

### 3.1 Algorithm

Given input `seq`:

1. `seq = seq.strip()`
2. `seq = seq.upper()`
3. Remove the following characters **anywhere** in the string:
   - space `' '`
   - newline `'\n'`
   - carriage return `'\r'`
   - tab `'\t'`

The result is `normalized_seq`.

### 3.2 Notes

- Normalization is intentionally simple and Solidity-portable.
- Non-ASCII characters are **not** removed; they will cause ASCII encoding to fail.
- **Symbol length** used elsewhere in GL1 is `L = len(normalized_seq)`.

---

## 4) Canonical encoding selection

Given normalized sequence `s`:

1. If `s` is **non-empty** and consists only of `A`, `C`, `G`, `T`: choose **DNA2** (ID 0)
2. Else if every char is within the **DNA4 alphabet**: choose **DNA4** (ID 2)
3. Else if every char is within the **SIXBIT alphabet (ALPH64)**: choose **SIXBIT** (ID 3)
4. Else choose **ASCII** (ID 1)

The empty string `""` is encoded as **ASCII**.

This order is deterministic and always prefers the most compact eligible encoding.

---

## 5) DNA2 encoding (ID 0, 2-bit)

### 5.1 Eligibility

DNA2 may be used only if the normalized sequence is:
- non-empty, and
- contains only: `A C G T`

### 5.2 Base-to-bit mapping

| Base | Bits | Value |
|---|---|---:|
| A | `00` | 0 |
| C | `01` | 1 |
| G | `10` | 2 |
| T | `11` | 3 |

### 5.3 Packing order (4 bases per byte)

Within each output byte:

- base0 → bits 7..6 (MSB)
- base1 → bits 5..4
- base2 → bits 3..2
- base3 → bits 1..0 (LSB)

If the final group has fewer than 4 bases, remaining 2-bit slots are padded with `00` (A).

### 5.4 Byte length

For base length `L`:

```
byte_len = ceil(L / 4) = (L + 3) // 4
```

### 5.5 Canonical validity (MUST)

A DNA2 payload is canonical iff:
1. `len(data_bytes) == byte_len`
2. Unused padded 2-bit slots in the final byte are all `00`.

---

## 6) DNA4 encoding (ID 2, 4-bit nibble / IUPAC bitmask)

DNA4 is designed for **efficient ambiguity-aware matching**.

### 6.1 Semantics: bitmask over {A,C,G,T}

Each symbol encodes a 4-bit mask:

- `A = 0001 (0x1)`
- `C = 0010 (0x2)`
- `G = 0100 (0x4)`
- `T = 1000 (0x8)`

Ambiguity codes are bitwise ORs of the above.
A special `GAP` symbol `'-'` is encoded as `0000 (0x0)`.

Matching property:
- A concrete base `b` matches ambiguity symbol `x` iff `(mask(b) & mask(x)) != 0`.

### 6.2 Alphabet and nibble values

| Char | Meaning | Nibble (hex) | Bits |
|---|---|---:|---|
| `-` | GAP | 0x0 | 0000 |
| `A` | A | 0x1 | 0001 |
| `C` | C | 0x2 | 0010 |
| `G` | G | 0x4 | 0100 |
| `T` | T | 0x8 | 1000 |
| `M` | A or C | 0x3 | 0011 |
| `R` | A or G | 0x5 | 0101 |
| `S` | C or G | 0x6 | 0110 |
| `V` | A or C or G | 0x7 | 0111 |
| `W` | A or T | 0x9 | 1001 |
| `Y` | C or T | 0xA | 1010 |
| `H` | A or C or T | 0xB | 1011 |
| `K` | G or T | 0xC | 1100 |
| `D` | A or G or T | 0xD | 1101 |
| `B` | C or G or T | 0xE | 1110 |
| `N` | A or C or G or T | 0xF | 1111 |

Notes:
- DNA4 intentionally does **not** include `U` (RNA). If `U` appears, SIXBIT or ASCII is used.

### 6.3 Packing order (2 symbols per byte)

For symbols `(s0, s1)`:

- `s0` nibble is stored in the **high nibble** (bits 7..4)
- `s1` nibble is stored in the **low nibble** (bits 3..0)

If `L` is odd, the final low nibble is padded with `0x0` (GAP).

### 6.4 Byte length

For symbol length `L`:

```
byte_len = ceil(L / 2) = (L + 1) // 2
```

### 6.5 Canonical validity (MUST)

A DNA4 payload is canonical iff:
1. `len(data_bytes) == byte_len`
2. If `L` is odd, the last byte’s low nibble is `0x0`
3. Every used nibble is one of the values listed above.

---

## 7) SIXBIT encoding (ID 3, 6-bit packed ALPH64)

SIXBIT is for symbolic sequences that fit into a fixed 64-character alphabet. It is useful for protein sequences and many other uppercase symbolic strings.

### 7.1 ALPH64 definition (MUST)

The alphabet is an ordered string of exactly 64 characters:

```
ALPH64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-*._:;,|/\\+=()[]{}<>#$%&?!@^"
```

The 6-bit code of a character is its **index** in this string (0..63).

### 7.2 Eligibility

SIXBIT may be used iff every character of the normalized sequence is in ALPH64.

### 7.3 Packing order (bitstream)

SIXBIT concatenates 6-bit codes MSB-first into a bitstream, then outputs bytes.

- Codes are appended in sequence order.
- The final byte is padded with **0 bits** (on the right / LSBs) if necessary.

### 7.4 Byte length

For symbol length `L`:

```
byte_len = ceil((L * 6) / 8) = (L*6 + 7) // 8
```

### 7.5 Canonical validity (MUST)

A SIXBIT payload is canonical iff:
1. `len(data_bytes) == byte_len`
2. Unused padding bits in the final byte (if any) are all zero.

---

## 8) ASCII encoding (ID 1, 8-bit)

### 8.1 Rules

The normalized sequence is encoded as ASCII:

- `data_bytes = normalized_seq.encode("ascii")`
- `length = len(normalized_seq)`

### 8.2 Canonical validity (MUST)

An ASCII payload is canonical iff:
1. `len(data_bytes) == length`
2. `data_bytes` are valid ASCII bytes

ASCII is the “future‑proof” fallback when a symbol is not in DNA2/DNA4/SIXBIT alphabets.

---

## 9) Reverse complement helper

GL1ENC defines a reverse-complement helper for strand handling.

### 9.1 Normalization-first

Reverse complement is defined over the **normalized sequence**.

### 9.2 Complement mapping (IUPAC)

Mapping (uppercase):

- `A ↔ T`
- `C ↔ G`
- `R ↔ Y`
- `K ↔ M`
- `B ↔ V`
- `D ↔ H`
- `S ↔ S`, `W ↔ W`, `N ↔ N`
- `U → A` (if present; note `U` is not DNA4-eligible but is supported by the helper)
- `- → -`, `* → *`, `. → .`

Unknown characters are left unchanged.

### 9.3 Operation

```
reverse_complement(seq) = reverse( complement( normalize(seq) ) )
```

---

## 10) Test vectors (MUST)

These vectors must match in Python and Solidity implementations.

### 10.1 DNA2 vectors

1. `seq = "ACGT"`
   - normalized: `"ACGT"`
   - bytes: `0x1b`
   - length: `4`

2. `seq = "T"`
   - bytes: `0xc0`
   - length: `1`

3. `seq = "ACGTAC"`
   - bytes: `0x1b10`
   - length: `6`

### 10.2 DNA4 vectors

1. `seq = "A"`
   - bytes: `0x10` (A in high nibble, pad 0 in low nibble)
   - length: `1`

2. `seq = "ACGT"`
   - A=0x1, C=0x2 → byte0 `0x12`
   - G=0x4, T=0x8 → byte1 `0x48`
   - bytes: `0x1248`
   - length: `4`

3. `seq = "N-"`
   - N=0xF, '-'=0x0 → byte `0xF0`
   - bytes: `0xF0`
   - length: `2`

### 10.3 SIXBIT vectors

Using ALPH64 index codes (A=0,B=1,C=2,D=3):

- `seq = "ABCD"`
- codes: `0,1,2,3`
- packed bytes: `0x001083`
- length: `4`

### 10.4 Canonical selection vectors

- `"ACGT"` → DNA2
- `"ACGN"` → DNA4
- `"MKWVTFISLL"` → SIXBIT
- `""` → ASCII

### 10.5 Reverse complement vector

- `"AAGT"` → `"ACTT"`

---

## 11) Solidity porting notes

- If you accept string input on-chain, normalization is expensive.
  Most designs will accept *already-normalized bytes* (or store bytes).
- Always store `length` alongside encoded bytes for DNA2/DNA4/SIXBIT.
- Enforce canonical padding bits/nibbles to prevent alternate encodings of the same logical string.
- Keep encoding IDs stable forever. If you ever need a new mapping/alphabet, introduce a new format version.
