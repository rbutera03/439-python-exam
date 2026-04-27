# DSP439 Python Exam - K-mer Analyzer
This repository contains a command-line DNA k-mer analyzer and a full pytest test suite.
The analyzer reads DNA sequences from a text file, validates them, counts k-mers, tracks
which nucleotide appears immediately after each k-mer, and writes aggregated results to an output file.

## Repository Structure
- `kmer_analyzer.py`: Main implementation and CLI entry point.
- `test_kmer_analyzer.py`: Unit tests for validation, counting logic, output formatting, and end-to-end behavior.
- `README.md`: Project documentation and usage instructions.

## How It Works
The workflow in `kmer_analyzer.py` is:

1. Parse CLI arguments: input file path, `k` value, and output file path.
2. Read each sequence line from the input file.
3. Normalize each sequence to uppercase.
4. Validate each sequence:
   - Sequence must be a string.
   - `k` must be a positive integer (not boolean).
   - Sequence length must be at least `k`.
   - Sequence can contain only `A`, `C`, `G`, `T`.
5. Skip invalid sequences with a warning.
6. For valid sequences, compute k-mer counts with following-character context.
7. Aggregate counts across all valid lines.
8. Write sorted results to the output file.

## Input Format
The input file should contain one DNA sequence per line, for example:

```text
ATGCGT
ATGATG
ATXG
atgc
```

- Invalid lines (like `ATXG`) are skipped.
- Lowercase letters are accepted because sequences are normalized to uppercase.

## Output Format
The output file contains one line per k-mer:

```text
<kmer> <total_count> <next_char_1>:<count> <next_char_2>:<count> ...
```

Example:

```text
AT 3 C:1 G:2
TG 2 A:1 C:1
```

Notes:
- K-mers are sorted alphabetically.
- Next characters for each k-mer are sorted alphabetically.
- Only k-mers that have a following character are counted.

## Run the Analyzer
From the repository root:

```bash
python kmer_analyzer.py <input_file> <k> <output_file>
```

Example:

```bash
python kmer_analyzer.py sequences.txt 2 output.txt
```

## Run the Tests
Install pytest (if needed), then run:

```bash
pytest -q
```

This executes tests for:
- Sequence validation rules
- K-mer update and counting behavior
- Output formatting and sorting
- Main CLI aggregation across multiple sequences

## AI Usage
AI assistance was used mainly for creating and expanding the test suite in
`test_kmer_analyzer.py`, since writing many similar edge-case and output-format
checks was the most repetitive and tedious part of this project.

I also used AI as a review pass for bug fixes I made in `kmer_analyzer.py` to
help confirm I had not missed anything crucial, especially around logic and
edge-case handling.
