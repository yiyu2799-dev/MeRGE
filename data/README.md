# Data Preparation

This guide describes the input format and preprocessing tools for **MeRGE (MEthylation Reconstruction with Graph Edge-features)**. Run the commands from the repository root. Paths beginning with `/path/to/` are placeholders for local data directories; the multiline examples use Bash syntax.

## Datasets and availability

The loaders support the following datasets through a shared array format. Dataset names are labels passed with `--dataset_name`; they do not select preprocessing rules or training parameters automatically.

| Dataset | Source | Input availability and supported protocol |
| --- | --- | --- |
| Hemato | [GSE87197](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87197) | MambaCpG-format inputs; chromosome holdout. |
| Neuron-Mouse | [GSE97179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97179) | MambaCpG-format inputs; chromosome holdout. |
| Neuron-Homo | [GSE97179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97179) | MambaCpG-format inputs; chromosome holdout. |
| LG_ACCPU | [HumanCellEpigenomeAtlas metadata](https://huggingface.co/datasets/zhoujt1994/HumanCellEpigenomeAtlas_metadata) | The documented experiment uses a chromosome 18 subset and a within-chromosome positional split. MeRGE-ready arrays and raw-data preprocessing for this subset are not distributed in this repository. |

The existing external processed-data source for Hemato, Neuron-Mouse, and Neuron-Homo is the [MambaCpG Zenodo release](https://doi.org/10.5281/zenodo.15571992). This is a third-party release, not a MeRGE data DOI. Obtain the matching DNA sequence, methylation, and position files from that source and check their format below before conversion. This repository does not contain datasets, pretrained checkpoints, or a download URL for a complete MeRGE processed-data package. Additional datasets can use the same format, but no additional dataset-specific preprocessing pipeline or reproduction configuration is supplied.

## Source arrays and conversion

The converter accepts three NPZ archives with matching chromosome keys:

| Archive | Array stored under each chromosome key |
| --- | --- |
| `<Dataset>_X.npz` | A one-dimensional encoded reference DNA sequence. Tokens `0`, `1`, `2`, and `3` retain the source encoding; all other values become `4` (unknown/padding). |
| `<Dataset>_y.npz` | A two-dimensional **CpG site × cell** methylation matrix with `0` = unmethylated, `1` = methylated, and `-1` = missing. |
| `<Dataset>_pos.npz` | A one-dimensional integer coordinate array in the same site order as the methylation matrix. |

All chromosomes must use the same cell-column order and cell count. The converter processes chromosome keys shared by all three archives; check that this intersection includes every intended chromosome. The script preserves the source nucleotide encoding and does not infer a nucleotide-to-token mapping. Reference sequences and CpG coordinates must use the same genome assembly and coordinate convention.

For Hemato:

```bash
python scripts/convert_mambacpg_to_merge.py \
  --input_dir /path/to/data_source/Hemato \
  --x_file Hemato_X.npz \
  --y_file Hemato_y.npz \
  --pos_file Hemato_pos.npz \
  --dataset_name Hemato \
  --output_dir /path/to/data/Hemato \
  --window 201 \
  --pos_base auto \
  --val_chrom 5 \
  --test_chrom 10
```

Use the corresponding filenames and dataset label for Neuron-Mouse or Neuron-Homo. The DNA window must have a positive odd width and agree with the model's `--dna_window` setting. Windows are centered on the supplied CpG coordinate; positions outside the reference sequence are padded with token `4`.

`--pos_base auto` uses a heuristic based on sequence-pair frequencies in up to three chromosomes. Inspect its reported choice against the source data. If the coordinate convention is known, use `--pos_base 0` or `--pos_base 1` explicitly. The converter preserves original coordinates in `pos.npy` and records the convention in `metadata.json`; it subtracts `pos_base` when indexing the reference sequence.

Use a new output directory for each conversion: the converter writes its output arrays in place. Its `--val_chrom` and `--test_chrom` options record split metadata; they do not create separate train/validation/test files. Actual training splits are selected by the training arguments described below.

## Processed-data format

Each processed directory must provide aligned arrays:

| File | Shape and contents |
| --- | --- |
| `meth_matrix.npy` | `(n_sites, n_cells)`, integer methylation values in `{-1, 0, 1}`. |
| `dna_windows_centerC.npy` | `(n_sites, window)`, integer DNA tokens in `{0, 1, 2, 3, 4}`. The default window is `201`. |
| `pos.npy` | `(n_sites,)`, genomic coordinates using the convention recorded by `pos_base`. |
| `chrom.npy` | `(n_sites,)`, chromosome labels. The converter writes fixed-width Unicode labels without a `chr` prefix. Use a non-object NumPy dtype for memory-mapped loading. |
| `metadata.json` | Dataset provenance and array metadata, including `n_cells`, `window`, `pos_base`, and `chromosome_lengths`. |

Keep site rows aligned across all four arrays and preserve cell-column identities between training, evaluation, and imputation. Cell count alone cannot verify cell identity or order.

The default `reference_length` position normalization requires reference chromosome lengths. The converter writes `chromosome_lengths` from `len(X[chrom])`, so `X` must contain the matching reference sequence rather than only a list of CpG bases. It leaves `reference_assembly` unset; retain the actual assembly information with the dataset. For manually prepared arrays, supply the corresponding lengths in `metadata.json`. Training also accepts `--reference_lengths_file` pointing to a JSON mapping from chromosome labels to lengths; relative paths are resolved within `--processed_dir`. The preflight checker does not accept that override and checks lengths in `metadata.json`.

If reference lengths are missing, reconvert compatible source archives into a new directory with the supplied converter, or provide verified reference metadata. Do not change position normalization solely to bypass missing metadata when reproducing an existing run.

## Splitting protocols

### Chromosome holdout

The default protocol assigns chromosome 5 to validation, chromosome 10 to testing, and the remaining chromosomes to training. It is used in the documented Hemato, Neuron-Mouse, and Neuron-Homo workflow. Validation and test chromosomes must both exist and differ, and at least one other chromosome is needed for training. Use the same settings in every training stage:

```text
--split_mode chromosome_holdout --val_chrom 5 --test_chrom 10
```

These are arguments to `train.py`, not a standalone shell command. Training reads the command-line settings rather than applying `metadata.json["split_rule"]` automatically.

### Within-chromosome positional split

For the documented LG_ACCPU chromosome 18 subset, use:

```text
--split_mode within_chromosome --split_chrom 18 --split_fractions 0.8 0.1 0.1
```

Pass these arguments to `train.py` in every stage with `--dataset_name LG_ACCPU` and the corresponding processed-data path. The processed arrays must contain exactly chromosome 18; this mode does not filter a multi-chromosome dataset automatically. The loader sorts sites by genomic position, rejects duplicate positions, and partitions them into contiguous training, validation, and test regions by **CpG-site count**, not by physical chromosome length. All fractions must be positive, sum to one, and yield nonempty regions. Segment overlap context and local neighbors remain within their assigned split.

Evaluation and imputation reconstruct the split from the checkpoint configuration. Existing data directory names and historical checkpoint labels need not be renamed to use the standardized dataset name in new runs.

## Data validation

After converting a chromosome-holdout dataset:

```bash
python scripts/check_processed.py \
  --processed_dir /path/to/data/Hemato \
  --expected_window 201 \
  --val_chrom 5 \
  --test_chrom 10 \
  --require_reference_lengths
```

For already prepared LG_ACCPU chromosome 18 arrays:

```bash
python scripts/check_processed.py \
  --processed_dir /path/to/data/LG_ACCPU \
  --expected_window 201 \
  --split_mode within_chromosome \
  --split_chrom 18 \
  --split_fractions 0.8 0.1 0.1 \
  --require_reference_lengths
```

Add `--expected_cells` with the verified cell count for the specific dataset when available. The checker validates array dimensions, split requirements, and reference-length metadata, and samples rows for methylation and DNA token checks. It is a preflight check, not an exhaustive validation of every value, cell identity, genome assembly, or sequence alignment.
