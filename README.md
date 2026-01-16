# The Structure of ApoB100 From Human Low-Density Lipoprotein

Code and data files to reproduce the custom image analysis and crosslink analysis from the paper.

## Citation

> [Paper citation - add full reference here]

## Repository Structure

```
├── data/                           # Crosslink data and cryo-EM class averages
│   ├── all_common_xlinks_small_and_large.csv   # Common crosslinks between datasets
│   ├── all_xlinks_large.csv                     # All crosslinks (large particles)
│   ├── all_xlinks_small.csv                     # All crosslinks (small particles)
│   ├── unique_xlinks_large.csv                  # Unique crosslinks (large)
│   ├── unique_xlinks_small.csv                  # Unique crosslinks (small)
│   └── LDL_200Classes.mrc.zip                   # Cryo-EM 2D class averages
├── python/                         # Python analysis scripts
│   ├── xlink_analysis_script.py                 # Main analysis script
│   └── xlink_analysis_functions.py              # Supporting functions
├── matlab/                         # MatLab analysis scripts
│   └── LDL_diameterANDeccentricity.m            # Morphological analysis
├── requirements.txt                # Python dependencies
└── LICENSE
```

## Data Files

The CSV files contain crosslinking data organized by residue pairs with:
- C-alpha distances for both apoB100 models
- Spectral counts
- Sequence distances
- Domain associations

The file `all_common_xlinks_small_and_large.csv` contains crosslinks found in common between the two independent datasets, with spectral counts averaged.

## Python Analysis

### Requirements

- Python 3.9+
- pandas
- matplotlib

### Installation

```bash
pip install -r requirements.txt
```

### Usage

Run from the repository root directory:

```bash
python python/xlink_analysis_script.py
```

This will reproduce the crosslinking analysis plots and tables from the paper:
- C-alpha Distance vs Spectral Count scatter plot
- Sequence Distance vs C-alpha Distance scatter plot
- Cumulative distribution of crosslink distances
- Domain statistics and summaries

## MatLab Analysis

### Requirements

- MatLab R2022a or later
- Image Processing Toolbox
- [ReadMRC](https://www.mathworks.com/matlabcentral/fileexchange/27021-imagic-mrc-dm-and-star-file-i-o) function from MATLAB File Exchange

### Usage

1. Download and unzip `data/LDL_200Classes.mrc.zip`
2. Download the ReadMRC script from the link above
3. Run `matlab/LDL_diameterANDeccentricity.m`

This performs morphological analysis of cryo-EM 2D class averages of LDL particles.

## System Information

Code was developed and tested on:
- MacBook Air M2, 24GB RAM
- macOS Sonoma 14.6.1
- Python 3.9.13
- MatLab R2022a with Image Processing Toolbox

Typical runtime: A few seconds for both Python and MatLab scripts.

## License

MIT License - see [LICENSE](LICENSE) file.
