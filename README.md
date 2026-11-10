# MeRGE

**MeRGE (MEthylation Reconstruction with Graph Edge-features)** is a graph neural network framework for imputing missing methylation states in sparse single-cell DNA methylation data. It combines a cell–CpG graph with DNA sequence, genomic position, cell identity, and multi-scale neighboring CpG information to predict methylation probabilities for unobserved cell–CpG pairs.

## Overview

The implementation contains four functional blocks:

1. **Input Block:** encodes CpG-centered DNA windows with a convolutional network, genomic coordinates with Fourier features, and cells with learnable identity embeddings. Observed binary methylation states provide graph edge features.
2. **GRAPE-based Global Graph Block:** performs edge-enhanced message passing on a bipartite cell–CpG graph to obtain global cell and CpG representations within each genomic segment. GRAPE refers to the underlying graph learning method.
3. **Multi-scale Neighbor Block:** aggregates same-cell neighboring CpGs at multiple scales using target-conditioned attention, methylation states, observation masks, and relative positions and distances.
4. **Output Block:** combines the target cell and CpG representations with local context in an MLP decoder. Sigmoid-transformed logits give methylation probabilities.

Training uses masked observed entries as prediction targets. The staged workflow first trains the global branch, then the local branch and expanded decoder, and finally fine-tunes the complete model. End-to-end training and ablation options are also exposed through `python train.py --help`.

## Repository Structure

```text
MeRGE/
├── README.md
├── LICENSE
├── requirements.txt
├── environment.yml
├── config.py                         # Model, data, and training defaults
├── train.py                          # Staged training entry point
├── evaluate.py                       # Evaluation of masked observed entries
├── impute.py                         # Prediction of originally missing entries
├── data/
│   ├── README.md                     # Input format and data preparation
│   ├── __init__.py
│   ├── dataset.py                    # Chromosome-holdout segments
│   ├── splitters.py                  # Genomic partition helpers
│   └── within_chromosome_dataset.py  # Contiguous positional partitions
├── models/
│   ├── __init__.py
│   ├── dna_encoder.py
│   ├── global_graph.py
│   ├── local_neighbor.py
│   └── model.py                      # Complete model
├── training/
│   ├── __init__.py
│   ├── checkpoints.py
│   ├── metrics.py
│   └── trainer.py
├── scripts/
│   ├── convert_mambacpg_to_merge.py   # Data conversion entry point
│   ├── convert_mambacpg_to_grapecpg.py # Legacy-compatible conversion entry point
│   └── check_processed.py            # Input validation
└── licenses/
    └── GRAPE_LICENSE                 # Upstream license notice
```

Processed datasets, pretrained checkpoints, and a test suite are not included in this checkout. Run the project directly from the repository root; it does not provide `setup.py` or `pyproject.toml` packaging.

## Installation

Clone the repository and run all subsequent commands from its root:

```bash
git clone https://github.com/yiyu2799-dev/MeRGE.git
cd MeRGE
conda env create -f environment.yml
conda activate merge
```

The environment specification pins Python 3.10.18, NumPy 1.26.4, SciPy 1.11.4, scikit-learn 1.3.2, PyTorch 2.0.1 with CUDA 11.8, and PyTorch Geometric 2.3.1. GPU execution requires a compatible NVIDIA driver. The entry points select CUDA when available and otherwise use the CPU.

Alternatively, in a Python 3.10 environment:

```bash
python -m pip install -r requirements.txt
```

`requirements.txt` additionally lists pytest; no tests are currently distributed in the repository. The examples below use Bash line continuations. On Windows PowerShell, put each command on one line or use PowerShell continuation syntax.

## Data Preparation

See [data/README.md](data/README.md) for data sources, conversion parameters, coordinate conventions, and validation commands. The supported input format consists of:

| File | Contents |
| --- | --- |
| `meth_matrix.npy` | CpG sites × cells; `0` unmethylated, `1` methylated, `-1` missing |
| `dna_windows_centerC.npy` | Integer-encoded, CpG-centered DNA windows aligned with the CpG rows |
| `pos.npy` | Genomic coordinate for each CpG row |
| `chrom.npy` | Chromosome identifier for each CpG row |
| `metadata.json` | Dataset metadata, including coordinate base and reference chromosome lengths |

Keep row alignment across arrays and cell-column order unchanged between training, evaluation, and imputation. Default position normalization requires reference chromosome lengths, supplied through metadata or the training option `--reference_lengths_file`.

The documented datasets and protocols are:

| Dataset | Preparation route | Partition protocol |
| --- | --- | --- |
| Hemato | Convert compatible MambaCpG NPZ files | Chromosome 5 validation; chromosome 10 test |
| Neuron-Mouse | Convert compatible MambaCpG NPZ files | Chromosome 5 validation; chromosome 10 test |
| Neuron-Homo | Convert compatible MambaCpG NPZ files | Chromosome 5 validation; chromosome 10 test |
| LG_ACCPU | Supply arrays in the format above; raw preprocessing and processed data are not included | Chromosome 18 subset; contiguous 80%/10%/10% split by CpG count |

These names identify datasets, not bundled example data. The loader accepts other datasets satisfying the same schema, but no additional dataset-specific preparation pipeline is provided. The repository does not currently supply a MeRGE processed-data release URL.

For example, convert and validate Hemato before training:

```bash
python scripts/convert_mambacpg_to_merge.py \
  --input_dir /path/to/data_source/Hemato \
  --x_file Hemato_X.npz \
  --y_file Hemato_y.npz \
  --pos_file Hemato_pos.npz \
  --dataset_name Hemato \
  --output_dir /path/to/data/Hemato \
  --window 201 --pos_base auto --val_chrom 5 --test_chrom 10

python scripts/check_processed.py \
  --processed_dir /path/to/data/Hemato \
  --expected_window 201 --val_chrom 5 --test_chrom 10 \
  --require_reference_lengths
```

All `/path/to/...` values are placeholders for local data directories. Confirm the input coordinate convention before conversion; `--pos_base auto` infers it from the supplied sequence and positions.

## Training

Defaults are defined in [config.py](config.py). The following commands use Hemato and the existing chromosome-holdout defaults. Each stage writes to `runs/<dataset_name>/<log_dir>/`.

### Stage 1: Global training

Train the global graph branch and its decoder with the local branch disabled:

```bash
python train.py \
  --dataset_name Hemato --processed_dir /path/to/data/Hemato \
  --segment_size 1024 --stage stage1 --log_dir stage1
```

### Stage 2: Local training

Load the best Stage-1 checkpoint. By default, the global branch is frozen, and the local branch and expanded decoder are trained. The default `d3` initialization preserves the Stage-1 decoder's global input columns.

```bash
python train.py \
  --dataset_name Hemato --processed_dir /path/to/data/Hemato \
  --segment_size 1024 --stage stage2 --log_dir stage2 \
  --load_checkpoint runs/Hemato/stage1/model_best_val_auroc.pt
```

### Stage 3: Joint fine-tuning

Load the best Stage-2 checkpoint and fine-tune all branches:

```bash
python train.py \
  --dataset_name Hemato --processed_dir /path/to/data/Hemato \
  --segment_size 1024 --stage stage3 --log_dir stage3 \
  --load_checkpoint runs/Hemato/stage2/model_best_val_auroc.pt
```

For Neuron-Mouse and Neuron-Homo, replace the dataset name and paths consistently. The previously documented segment sizes are 512 and 256, respectively; the parser does not select them automatically. Keep the same segment size, split settings, and relevant model options across stages.

For LG_ACCPU, use `--dataset_name LG_ACCPU` and its processed directory, and append `--split_mode within_chromosome --split_chrom 18 --split_fractions 0.8 0.1 0.1` to **each** training stage. Select the segment size from the original experiment configuration; no LG_ACCPU-specific segment-size preset is supplied. Adapt checkpoint paths to `runs/LG_ACCPU/<log_dir>/`.

Training records `args.json`, `cmd_input.txt`, `environment.json`, and `result_summary.json`. The best finite validation AUROC determines `model_best_val_auroc.pt`; early stopping separately monitors validation loss. Stage 3 and end-to-end training also evaluate the test split and write `final_test_metrics.json`.

## Evaluation

Evaluate masked observed entries using a checkpoint produced by this training code:

```bash
python evaluate.py \
  --checkpoint runs/Hemato/stage3/model_best_val_auroc.pt \
  --processed_dir /path/to/data/Hemato --split test \
  --output runs/Hemato/stage3/evaluation_test.json
```

Use `--split val` for validation. The model and split configuration are restored from the checkpoint; `--processed_dir` only overrides the data location. Reported fields include loss, accuracy, precision, recall, specificity, F1, MCC, AUROC, and AUPRC (computed as average precision). AUROC and AUPRC are undefined when the evaluated labels contain only one class. The script requires an embedded checkpoint `config`; a bare state dictionary is insufficient for this entry point.

### Imputation of missing entries

Predict entries originally encoded as `-1`:

```bash
python impute.py \
  --checkpoint runs/Hemato/stage3/model_best_val_auroc.pt \
  --processed_dir /path/to/data/Hemato --split all \
  --output_dir runs/Hemato/stage3/imputed
```

Imputation uses observed methylation entries as context and writes one NPZ shard per segment with missing targets, plus `manifest.json`. Each shard contains `global_site_index`, `cell_index`, `probability`, and `predicted_state` (threshold 0.5). These indices refer to the original input arrays. This is distinct from held-out-label evaluation and does not produce a dense imputed matrix. Use a separate output directory for each checkpoint and scope.

## Reproducibility

- Preserve the exact processed arrays, reference genome and coordinate convention, chromosome lengths, CpG row order, and cell-column order. Cell identity embeddings require the same cells in the same order; matching cell counts alone is insufficient.
- Record the Git revision, complete commands, dataset provenance, and generated run metadata. The distributed defaults and examples do not replace the original per-experiment configuration or supply missing datasets and checkpoints.
- The default seed is 1. Training masks are resampled each epoch by default; validation and test use fixed seeds with offsets of 100000 and 200000. GPU operations are not configured for strict deterministic execution, so bitwise reproducibility is not guaranteed.
- Chromosome holdout is the default. For a within-chromosome run, the chromosome and all three fractions must be explicit. Positional partitions, overlap context, and local neighborhoods remain within their respective splits.
- Keep preprocessing, masking, segmentation, architecture, and split options consistent across stages. Use new run directories to preserve previous results. Checkpoint selection uses validation AUROC, while early stopping uses validation loss.
- `MeRGEModel` is the public model name. The legacy `GrapeCpGModel` class and conversion-script path remain available for compatibility; this naming update does not change model tensors, state-dictionary keys, checkpoint format, or training defaults. Existing conda environments keep their original names; newly created environments use `merge`.

## Acknowledgements

The global graph block builds on [GRAPE](https://github.com/maxiaoba/GRAPE). Data conversion supports the processed NPZ format used by [MambaCpG](https://github.com/Lee-qiangzee/MambaCpG). These names refer to separate existing methods.

## License

MeRGE is distributed under the [MIT License](LICENSE). The upstream GRAPE notice is retained in [licenses/GRAPE_LICENSE](licenses/GRAPE_LICENSE).
