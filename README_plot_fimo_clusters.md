# plot_fimo_clusters.py

Plot clustered FIMO motif hits relative to a transcription start site (TSS).

This script takes FIMO output files, groups motif hits into families, identifies spatial clusters of motif density around a reference TSS, and generates both track-style plots and heatmaps summarizing motif enrichment patterns.

It is designed for downstream visualization and interpretation of motif-discovery results.

---

## Overview

This script:
- Reads FIMO TSV output files
- Assigns motifs to transcription factor families
- Detects clusters of motif hits relative to a TSS
- Outputs summary tables and publication-ready plots

The workflow is intended to be reproducible and flexible for promoter or enhancer motif analyses.

---

## Requirements

- Python ≥ 3.8

### Python dependencies

Typical dependencies include:
```bash
pip install pandas numpy matplotlib seaborn scipy
```

(The exact imports are defined in the script.)

---

## Input

### Required
- One or more **FIMO TSV output files**

FIMO files must contain standard columns such as:
- motif_id
- sequence_name
- start
- stop
- strand
- score
- p-value
- q-value

### Optional metadata
- Motif family definitions (embedded in the script via `FAMILY_MAP`)
- Custom TSS position (via command-line arguments)

---

## Usage

Run the script from the command line:

```bash
python plot_fimo_clusters.py \
    --fimo results/fimo.tsv \
    --outprefix gata4_motifs \
    --tss 0
```

Multiple FIMO files can be provided if supported by your analysis setup.

---

## Key Arguments

Common arguments include:

- `--fimo`  
  Path to one or more FIMO TSV files

- `--outprefix`  
  Prefix used for all output files

- `--tss`  
  Reference transcription start site position (default: 0)

- `--window`  
  Window size around the TSS to consider for clustering

- `--bins`  
  Number of bins used for heatmap visualization

- `--title`  
  Custom plot title

(See the script’s `argparse` section for the full and authoritative list.)

---

## Output

The script generates the following files:

- `<outprefix>_windows.csv`  
  Sliding window counts used for cluster detection

- `<outprefix>_clusters.csv`  
  Summary of detected motif clusters

- `<outprefix>_tracks.png` / `.pdf`  
  Motif hit tracks plotted relative to the TSS

- `<outprefix>_heatmap.png` / `.pdf`  
  Heatmap of motif density across bins

- `<outprefix>_cluster<N>_hits.tsv`  
  Individual TSV files containing all motif hits for each detected cluster

---

## Motif Families

Motifs are grouped into transcription factor families using predefined mappings in the script:

- `FAMILY_MAP`
- `FAMILY_ORDER`

You may edit these dictionaries to:
- Add new TF families
- Reassign motifs
- Control plotting order

---

## Customization

You may want to modify the script to:
- Change clustering thresholds
- Adjust window or bin sizes
- Customize color schemes
- Restrict analysis to specific motif families
- Integrate with upstream motif-selection pipelines

All parameters are centralized and documented in the code for ease of modification.

---

## Example Workflow

1. Run FIMO with a MEME motif file:
   ```bash
   fimo motifs.meme sequences.fasta
   ```

2. Plot clustered motif hits:
   ```bash
   python plot_fimo_clusters.py --fimo fimo.tsv --outprefix analysis1
   ```

3. Use the generated plots and tables for interpretation or representation.

---


## Author

Nicolas Noel

For questions or contributions, please open an issue or pull request.
