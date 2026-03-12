# IggNition

**Ultra-fast nucleotide-level Aho numbering for antibody variable domains.**

IggNition replaces ANARCI/HMMER with a purpose-built Rust aligner against pre-numbered human germline V/D/J genes. It assigns [Aho scheme](https://doi.org/10.1006/jmbi.2001.4662) positions to every nucleotide in an antibody variable domain sequence, producing a fixed-length coordinate frame suitable for repertoire analysis and AbLLM training.

| | ANARCI (Python/HMMER) | **IgNition** |
|---|---|---|
| Throughput (single core) | ~4 seq/s | **>50,000 seq/s** |
| Language | Python + HMMER | Rust |
| Output | amino-acid level | **nucleotide level** |
| Parallelism | external (shell) | built-in (Rayon) |

## Installation

```bash
pip install iggnition
```

Pre-compiled wheels are available for Linux (x86\_64, aarch64), macOS (Apple Silicon + Intel), and Windows.

## Quick Start

### Single sequence

```python
import iggnition

df = iggnition.run(
    nt_seq="CAGGTGCAGCTGGTGCAGTCTGGAGCT...",
    aa_seq="QVQLVQSGAE...",
)
print(df)
# shape: (447, 7)  ← 149 Aho positions × 3 nucleotides for a heavy chain
# ┌─────────────┬───────┬─────────────┬──────────────┬───────────────┬────────────┬───────────┐
# │ sequence_id ┆ chain ┆ nt_position ┆ aho_position ┆ codon_position┆ nucleotide ┆ amino_acid│
# │ ---         ┆ ---   ┆ ---         ┆ ---          ┆ ---           ┆ ---        ┆ ---       │
# │ u32         ┆ str   ┆ u32         ┆ u32          ┆ u32           ┆ str        ┆ str       │
```

### AIRR-format DataFrame (single chain)

```python
import polars as pl
import iggnition

df = pl.read_csv("airr_table.tsv", separator="\t")

results, errors = iggnition.run(
    df,
    nt_col="sequence",
    aa_col="sequence_aa",
    locus_col="locus",   # e.g. "IGH", "IGK", "IGL"
)
```

### Paired heavy + light (PairPlex-style)

```python
results, errors = iggnition.run(
    df,
    paired=True,
    nt_col_heavy="sequence:0",
    aa_col_heavy="sequence_aa:0",
    nt_col_light="sequence:1",
    aa_col_light="sequence_aa:1",
)
```

### File paths

```python
# FASTA → DataFrame
results, errors = iggnition.run("input.fasta")

# Parquet → Parquet
iggnition.run("input.parquet", output="numbered.parquet")

# TSV → TSV
iggnition.run("input.tsv", output="numbered.tsv")
```

### Output formats

```python
# Per-codon (one row per Aho position)
df = iggnition.run(nt_seq=nt, aa_seq=aa, per_codon=True)
# columns: sequence_id, chain, aho_position, codon, amino_acid

# Wide format (one row per sequence, positional columns)
df = iggnition.run(nt_seq=nt, aa_seq=aa, wide=True)
# columns: sequence_id, H_nt_1, H_nt_2, ..., H_nt_447
```

## CLI

```bash
# FASTA → TSV (stdout)
iggnition run input.fasta

# FASTA → TSV (file)
iggnition run input.fasta output.tsv

# Parquet → Parquet
iggnition run input.parquet output.parquet

# TSV (AIRR) → TSV, per-codon
iggnition run input.tsv output.tsv --per-codon

# Wide format, 8 threads
iggnition run input.fasta output.tsv --wide --threads 8
```

### CLI options

| Option | Default | Description |
|--------|---------|-------------|
| `--per-codon` | off | One row per codon instead of per nucleotide |
| `--wide` | off | Pivot to wide format |
| `--nt-col` | `sequence` | NT column name (TSV/Parquet) |
| `--aa-col` | `sequence_aa` | AA column name |
| `--nt-col-heavy` | `sequence:0` | Heavy chain NT column |
| `--aa-col-heavy` | `sequence_aa:0` | Heavy chain AA column |
| `--nt-col-light` | `sequence:1` | Light chain NT column |
| `--aa-col-light` | `sequence_aa:1` | Light chain AA column |
| `--no-aa` | off | Auto-translate NT (fallback, emits warning) |
| `--threads` | all cores | Rayon worker threads |
| `--chunk-size` | 10000 | Sequences per processing chunk |

## Output Schema

### Per-nucleotide (default)

| Column | Type | Description |
|--------|------|-------------|
| `sequence_id` | u32 | Row index from input |
| `chain` | str | `"H"`, `"K"`, or `"L"` |
| `nt_position` | u32 | Absolute nucleotide position (1-based) |
| `aho_position` | u32 | Aho amino acid position (1-based) |
| `codon_position` | u32 | Position within codon (1, 2, 3) |
| `nucleotide` | str | `A`/`T`/`G`/`C` or `-` for gaps |
| `amino_acid` | str | Single-letter AA or `-` for gaps |

### Aho position ranges

| Chain | Max Aho position | Max NT positions |
|-------|-----------------|-----------------|
| H | 149 | 447 |
| K | 148 | 444 |
| L | 148 | 444 |

## How It Works

1. All human V/D/J germline genes are pre-numbered with Aho positions and **embedded in the binary at compile time** — no external database files.
2. For each query: find the closest germline via Needleman-Wunsch alignment (amino acid level), then transfer Aho positions from the germline to the query.
3. Map each occupied Aho position to its codon from the nucleotide sequence: position `N` → nucleotides `(N-1)*3+1`, `(N-1)*3+2`, `(N-1)*3+3`. Unoccupied positions become gaps (`-`).
4. **Rayon** parallelises across sequences in the batch; the Python GIL is never held during alignment.

## Building from Source

```bash
git clone https://github.com/bnemoz/ignition
cd ignition
pip install maturin
maturin develop --release
```

## License

MIT. If you use IgNition in published research, please cite accordingly.

## References

- Honegger, A. & Plückthun, A. (2001). Yet another numbering scheme for immunoglobulin variable domains. *J Mol Biol*, 309(3), 657–670.
- Dunbar, J. & Deane, C.M. (2016). ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics*, 32(2), 298–300.
