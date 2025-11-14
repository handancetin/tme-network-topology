# tme-network-topology
This repository contains the complete computational pipeline and analysis code for the manuscript **"Multi-scale Network Topology Analysis Reveals Metabolic Reprogramming and Therapeutic Vulnerabilities in the Tumor Microenvironment"**.

```
tme-network-topology/
├── matlab/                   # MATLAB code for metabolic modeling
│   ├── functions/            # Custom MATLAB functions
│   ├── output/               # MATLAB analysis outputs
│   │   └── simulations/      # Flux sampling and knockout results
│   └── scripts/              # Main MATLAB analysis scripts
├── python/                   # Python notebooks for analysis and visualization
│   ├── output/               # Python analysis outputs
│   ├── fig1_gem_sampling.ipynb           # Context-specific models and sampling analysis
│   ├── fig2_knockout_survival.ipynb      # Gene knockout and survival analysis
│   ├── fig3_centrality_definitions.ipynb # Centrality measure definitions
│   ├── fig4_network_topology.ipynb       # Topological analysis
│   ├── fig5_role_transitions.ipynb       # Metabolite role transition analysis
│   └── fig6_graph_geometry.ipynb         # Multifractal and curvature analysis
└── r/                        # R code for scRNA-seq processing
│   └── output/               # R analysis outputs
├── figures/                  # Generated figures
│   ├── pdf/                  # PDF versions of all figures
│   └── png/                  # PNG versions of all figures
```

All main figures and supplementary figures can be regenerated using the scripts in Jupyter Notebooks.

## Citation
If you use this code or data, please cite:

```
[Citation will be added upon publication]
```
