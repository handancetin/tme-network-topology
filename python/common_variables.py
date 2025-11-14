
# Basic Python Libraries
import os
import re
import math
import glob
import ast
import warnings
warnings.filterwarnings('ignore', message='n_jobs value.*overridden to 1 by setting random_state')
warnings.filterwarnings('ignore', message='Input histogram consists of integer.*')
warnings.filterwarnings('ignore', message='.*c.* argument looks like a single numeric RGB or RGBA sequence.*')
                       

# Progress Tracking
from tqdm import tqdm

# Data Manipulation and Analysis
import numpy as np
import pandas as pd
import itertools
from itertools import combinations
from collections import OrderedDict, Counter, defaultdict

# Scientific Computing and Statistics
import scipy.stats as stats
from scipy import io
from scipy.sparse import csr_matrix
from scipy.spatial.distance import euclidean, jaccard, squareform
from scipy.stats import zscore, wasserstein_distance
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list, ward, set_link_color_palette
from statsmodels.stats.multitest import fdrcorrection

# Machine Learning
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score 
from sklearn.feature_selection import VarianceThreshold

# Visualization
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch, Rectangle, Polygon, Wedge
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from adjustText import adjust_text

# Network Analysis
import networkx as nx
import networkit as nk
from GraphRicciCurvature.OllivierRicci import OllivierRicci
from GraphRicciCurvature.FormanRicci import FormanRicci

# Dimensionality Reduction
import umap

# ===== CONFIGURATION CONSTANTS =====

markers_all = ['Malignant', 'FGFR2+', 'CD73+', 'ICAM1+', 'ICAM1-', 'FAP+', 'THBS1+', 'VCAN+', 'MARCO+']
markers_F   = ['FGFR2+', 'CD73+', 'ICAM1+', 'ICAM1-', 'FAP+']
markers_M   = ['THBS1+', 'VCAN+', 'MARCO+']

# Color definitions - use CAPS for constants
marker_colors = {
    "Malignant": "#8d8d8d",  
    "Epithelial": "#8d8d8d",  
    "FAP+": "#FD8D3D", 
    "ICAM1+": "#DE77AE",  
    "ICAM1-": "#FFC2E7",  
    "FGFR2+": "#C51C7D", 
    "CD73+": "#b6202b", 
    #"MFAP5+": "#FECC5C",  # myofibroblasts
    #"DES+": "#893776",  # myofibroblasts
    "THBS1+": "#90DEF9",  
    "VCAN+": "#00499b",  
    "MARCO+": "#007ef7"
}

marker_classes = {
    "Malignant":"Epithelial",   
    "FGFR2+":   "Fibroblast", 
    "ICAM1+":   "Fibroblast", 
    "FAP+":     "Fibroblast", 
    "ICAM1-":   "Fibroblast", 
    "CD73+":    "Fibroblast", 
    #"MFAP5+":   "Myofibroblast", 
    #"DES+":     "Myofibroblast", 
    "THBS1+":   "Macrophage",  
    "VCAN+":    "Macrophage",  
    "MARCO+":   "Macrophage"
}

# Subsystem conversion dictionary
subsystem_conversion_dict = {
    # Maintain these larger subsystems as they are
    "transport reactions": "transport reactions",
    "exchange/demand reactions": "exchange reactions", 
    "glycolysis / gluconeogenesis": "glycolysis and gluconeogenesis",
    "pentose phosphate pathway": "pentose phosphate pathway",
    "pyruvate metabolism": "pyruvate metabolism",
    "oxidative phosphorylation": "oxidative phosphorylation",
    "artificial reactions": "artificial reactions",
    "pool reactions": "artificial reactions",
    "sphingolipid metabolism": "sphingolipid metabolism",
    "glycerolipid metabolism": "glycerolipid metabolism", 

    # Beta Oxidation merging
    "beta oxidation of even-chain fatty acids (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of even-chain fatty acids (peroxisomal)": "fatty acid metabolism",
    "beta oxidation of odd-chain fatty acids (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of odd-chain fatty acids (peroxisomal)": "fatty acid metabolism",
    "beta oxidation of poly-unsaturated fatty acids (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of poly-unsaturated fatty acids (n-6) (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of poly-unsaturated fatty acids (n-9) (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of di-unsaturated fatty acids (n-6) (peroxisomal)": "fatty acid metabolism",
    "beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of branched-chain fatty acids (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of unsaturated fatty acids (n-7) (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of unsaturated fatty acids (n-7) (peroxisomal)": "fatty acid metabolism",
    "beta oxidation of unsaturated fatty acids (n-9) (mitochondrial)": "fatty acid metabolism",
    "beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)": "fatty acid metabolism",
    "beta oxidation of phytanic acid (peroxisomal)": "fatty acid metabolism",
    
    # Fatty Acid Metabolism 
    "fatty acid biosynthesis (even-chain)": "fatty acid metabolism",
    "fatty acid biosynthesis (odd-chain)": "fatty acid metabolism",
    "fatty acid biosynthesis (unsaturated)": "fatty acid metabolism",
    "fatty acid oxidation": "fatty acid metabolism",
    "omega-3 fatty acid metabolism": "fatty acid metabolism",
    "omega-6 fatty acid metabolism": "fatty acid metabolism",
    "linoleate metabolism": "fatty acid metabolism",
    "fatty acid activation (cytosolic)": "fatty acid metabolism", 
    "fatty acid activation (endoplasmic reticular)": "fatty acid metabolism",
    "fatty acid elongation (even-chain)": "fatty acid metabolism",
    "fatty acid elongation (odd-chain)": "fatty acid metabolism",
    "fatty acid desaturation (even-chain)": "fatty acid metabolism",
    "fatty acid desaturation (odd-chain)": "fatty acid metabolism",
    "fatty acid degradation": "fatty acid metabolism",
    
    # Carnitine Shuttle merging
    "carnitine shuttle (endoplasmic reticular)": "carnitine shuttle metabolism",
    "carnitine shuttle (mitochondrial)": "carnitine shuttle metabolism",
    "carnitine shuttle (cytosolic)": "carnitine shuttle metabolism",
    "carnitine shuttle (peroxisomal)": "carnitine shuttle metabolism",
    
    # Vitamin Metabolism merging
    "vitamin a metabolism": "vitamin metabolism",
    "vitamin b2 metabolism": "vitamin metabolism",
    "vitamin b6 metabolism": "vitamin metabolism",
    "vitamin b12 metabolism": "vitamin metabolism",
    "vitamin c metabolism": "vitamin metabolism",
    "vitamin d metabolism": "vitamin metabolism",
    "vitamin e metabolism": "vitamin metabolism",
    "retinol metabolism": "vitamin metabolism",
    "riboflavin metabolism": "vitamin metabolism",
    
    # Glycosphingolipid Metabolism merging
    "glycosphingolipid biosynthesis-lacto and neolacto series": "glycosphingolipid metabolism",
    "glycosphingolipid biosynthesis-ganglio series": "glycosphingolipid metabolism",
    "glycosphingolipid biosynthesis-globo series": "glycosphingolipid metabolism",
    
    # Cholesterol Metabolism merging
    "formation and hydrolysis of cholesterol esters": "cholesterol metabolism",
    "cholesterol biosynthesis 1 (bloch pathway)": "cholesterol biosynthesis",
    "cholesterol biosynthesis 2": "cholesterol biosynthesis",
    "cholesterol biosynthesis 3 (kandustch-russell pathway)": "cholesterol biosynthesis",
    
    # Amino Acid Metabolism merging
    "cysteine and methionine metabolism": "amino acid metabolism",
    "glycine, serine and threonine metabolism": "amino acid metabolism",
    "alanine, aspartate and glutamate metabolism": "amino acid metabolism",
    "valine, leucine, and isoleucine metabolism": "amino acid metabolism",
    "arginine and proline metabolism": "amino acid metabolism",
    "methionine metabolism": "amino acid metabolism",
    "lysine metabolism": "amino acid metabolism",
    "histidine metabolism": "amino acid metabolism",
    "tyrosine metabolism": "amino acid metabolism",
    "tryptophan metabolism": "amino acid metabolism",
    "phenylalanine, tyrosine and tryptophan biosynthesis": "amino acid metabolism",
    "beta-alanine metabolism": "amino acid metabolism",
    "phenylalanine metabolism": "amino acid metabolism",
    "metabolism of other amino acids": "amino acid metabolism",
    
    # Lipid Metabolism merging
    "eicosanoid metabolism": "lipid metabolism",
    "steroid metabolism": "lipid metabolism",
    "acylglycerides metabolism": "lipid metabolism",
    "prostaglandin biosynthesis": "lipid metabolism",
    "acyl-coa hydrolysis": "lipid metabolism",
    "other lipid metabolism": "lipid metabolism",
    "lipoic acid metabolism": "lipid metabolism",
    "triacylglycerol synthesis": "lipid metabolism",
    "inositol phosphate metabolism": "lipid metabolism",
    "phosphatidylinositol phosphate metabolism": "lipid metabolism",
    "phosphatidylinositol": "lipid metabolism",
    "ether lipid metabolism": "lipid metabolism",
    
    # Bile Acid Metabolism merging
    "bile acid biosynthesis": "bile acid metabolism",
    "bile acid recycling": "bile acid metabolism",
    
    # Carbohydrate Metabolism merging
    "starch and sucrose metabolism": "carbohydrate metabolism",
    "fructose and mannose metabolism": "carbohydrate metabolism",
    "galactose metabolism": "carbohydrate metabolism",
    "pentose and glucuronate interconversions": "carbohydrate metabolism",
    "ascorbate and aldarate metabolism": "carbohydrate metabolism",
    "glyoxylate and dicarboxylate metabolism": "carbohydrate metabolism",
    "propanoate and nicotinamide metabolism": "carbohydrate metabolism",
    "c5-branched dibasic acid metabolism": "carbohydrate metabolism",
    
    # Nucleotide and DNA/RNA merging
    "nucleotide metabolism": "nucleotide metabolism",
    "purine metabolism": "nucleotide metabolism",
    "pyrimidine metabolism": "nucleotide metabolism",
    "nicotinate metabolism": "nucleotide metabolism",
    "amino sugar and nucleotide sugar metabolism": "nucleotide metabolism",
    "aminoacyl-trna biosynthesis": "nucleotide metabolism",
    
    # TCA cycle related
    "tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism": "tricarboxylic acid cycle",
    
    # Hormone merging
    "androgen metabolism": "hormone metabolism",
    "estrogen metabolism": "hormone metabolism",
    "insect hormone biosynthesis": "hormone metabolism",
    
    # Protein and Peptide merging
    "peptide metabolism": "protein and peptide metabolism",
    "protein metabolism": "protein and peptide metabolism",
    "protein degradation": "protein and peptide metabolism",
    "protein assembly": "protein and peptide metabolism",
    "protein modification": "protein and peptide metabolism",
    
    # Sulfur Metabolism
    "keratin sulfide degradation": "sulfur metabolism",
    "keratan sulfate biosynthesis": "sulfur metabolism",
    "heparan sulfate degradation": "sulfur metabolism",
    "chondroitin / heparan sulfate biosynthesis": "sulfur metabolism",
    "chondroitin sulfate degradation": "sulfur metabolism",
    "taurine metabolism": "sulfur metabolism",
    "glutathione metabolism": "glutathione metabolism",
    "sulfur metabolism": "sulfur metabolism",
    
    # Xenobiotic Metabolism
    "xenobiotics metabolism": "xenobiotic metabolism",
    "drug metabolism": "xenobiotic metabolism",
    
    # Glycan Metabolism
    "o-glycan metabolism": "glycan metabolism",
    "n-glycan metabolism": "glycan metabolism",
    "glycosylphosphatidylinositol (gpi)-anchor biosynthesis": "glycan metabolism",
    "blood group biosynthesis": "glycan metabolism",
    
    # Coenzyme Metabolism
    "pantothenate and coa metabolism": "coenzyme metabolism",
    "porphyrin metabolism": "coenzyme metabolism",
    "biotin metabolism": "coenzyme metabolism",
    "folate metabolism": "coenzyme metabolism",
    "thiamine metabolism": "coenzyme metabolism",
    "ubiquinone and other terpenoid-quinone biosynthesis": "coenzyme metabolism",
    "ubiquinone synthesis": "coenzyme metabolism",
    "heme degradation": "coenzyme metabolism",

    # Arachidonic acid metabolism
    "leukotriene metabolism": "leukotriene and arachidonic acid metabolism",
    "arachidonic acid metabolism": "leukotriene and arachidonic acid metabolism",
    "\tarachidonic acid metabolism": "leukotriene and arachidonic acid metabolism",
 
    # Other not-so specific pathways
    "terpenoid backbone biosynthesis": "other reactions",
    "serotonin and melatonin biosynthesis": "other reactions",
    "inositate metabolism": "other reactions",
    "glucuronide biosynthesis": "other reactions",
    "glucocorticoid biosynthesis": "other reactions",
    "butanoate metabolism": "other reactions",
    "octane oxidation": "other reactions",
    "urea cycle": "other reactions", 
    "dietary fiber binding": "other reactions",
    "hippurate metabolism": "other reactions",
    "toluene degradation": "other reactions",
    "alkaloids biosynthesis": "other reactions",
    "miscellaneous": "other reactions",
    "isolated": "other reactions",
}

# Ordering for pathways (from differential results)
def getPathwayImportanceOrdering():
        return {
            'fatty acid metabolism': 4, 'carnitine shuttle metabolism': 4, 'leukotriene and arachidonic acid metabolism': 4,
            'amino acid metabolism': 4, 'nucleotide metabolism': 4, 'ros detoxification': 4, 'lipid metabolism': 4,
            'bile acid metabolism': 4, 'glutathione metabolism': 4, 'glycolysis and gluconeogenesis': 3,
            'coenzyme metabolism': 3, 'tricarboxylic acid cycle': 3, 'oxidative phosphorylation': 3,
            'hormone metabolism': 3, 'vitamin metabolism': 3, 'pentose phosphate pathway': 3,
            'carbohydrate metabolism': 2, 'sphingolipid metabolism': 2, 'glycerolipid metabolism': 2,
            'glycerophospholipid metabolism': 2, 'cholesterol metabolism': 2, 'cholesterol biosynthesis': 2,
            'fatty acid biosynthesis': 2, 'pyruvate metabolism': 2, 'propanoate metabolism': 2,
            'nicotinate and nicotinamide metabolism': 2, 'biopterin metabolism': 2, 'sulfur metabolism': 2,
            'protein and peptide metabolism': 2, 'glycan metabolism': 2, 'glycosphingolipid metabolism': 2,
            'keratan sulfate degradation': 1, 'xenobiotic metabolism': 1, 'other reactions': 1,
            'exchange reactions': 0, 'transport reactions': 0, 'artificial reactions': 0
}

# ===== UTILITY FUNCTIONS =====

def set_figure_style():
    sns.set_theme(style="ticks", rc={
        'axes.facecolor': 'white', 
        'grid.color': '0.97',
        'font.family': 'Arial',
        'font.size': 16,              
        'axes.labelsize': 18,       
        'axes.titlesize': 18,   
        'xtick.labelsize': 14,      
        'ytick.labelsize': 14, 
        'text.color': '#000000',
        'axes.labelcolor': '#000000',
        'xtick.color': '#000000',
        'ytick.color': '#000000', 
        'axes.edgecolor': '#000000',
        'axes.grid': True, 
        'lines.color': '#000000',     
        'grid.linestyle': '-'
    })

def save_figure(fig, filename, output_dir='../figures/', transparent=True):
    from pathlib import Path
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    png_dir = output_path / 'png'
    svg_dir = output_path / 'svg'
    pdf_dir = output_path / 'pdf'
    png_dir.mkdir(exist_ok=True)
    svg_dir.mkdir(exist_ok=True)
    pdf_dir.mkdir(exist_ok=True)
    
    # Save figures
    fig.savefig(png_dir / f"{filename}.png", dpi=300, bbox_inches="tight", format='png', transparent=transparent)
    fig.savefig(svg_dir / f"{filename}.svg", dpi=300, bbox_inches="tight", format='svg', transparent=transparent)
    fig.savefig(pdf_dir / f"{filename}.pdf", dpi=300, bbox_inches="tight", format='pdf', transparent=transparent)

    print(f"Figure saved: {filename}")

def get_base_paths():
    username = os.getenv("USER") or os.getenv("USERNAME")
    
    if username == "handan":  # iMac
        base_dir = "/Users/handan/Documents/CSBL/project-files/tme-network-topology"
    elif username == "hcetin":  # MacBook
        base_dir = "/Users/hcetin/Documents/CSBL/project-files/tme-network-topology"
    else:
        raise ValueError("Please specify base_dir path manually.")
    
    return {
        'base_dir': base_dir,
        'dir_binary_rxns': f"{base_dir}/matlab/output/binaryRxnsTable.csv",
        'dir_binary_mets': f"{base_dir}/matlab/output/binaryMetsTable.csv",
        'dir_binary_genes': f"{base_dir}/matlab/output/binaryGenesTable.csv",
        'dir_umap_data': f"{base_dir}/r/output/umap_pos_data.csv",
        'dir_samples': f"{base_dir}/matlab/output/samples",
        'dir_nodesedges': f"{base_dir}/matlab/output/graphs",
        'dir_simulations': f"{base_dir}/matlab/output/simulations",
        'dir_survival_data': f"{base_dir}/r/output_data"
    }

def apply_subsystem_conversion(df):
    df['subSystems'] = df['subSystems_og'].map(lambda x: subsystem_conversion_dict.get(x, x))
    return df

def load_model_tables():
    paths = get_base_paths()
    
    # Load datasets
    binary_rxns_df = pd.read_csv(paths['dir_binary_rxns'], header=0)
    binary_rxns_df['ReactionNames'] = binary_rxns_df['ReactionNames'].fillna(binary_rxns_df['Reactions'])
    binary_rxns_df['subSystems_og'] = binary_rxns_df['subSystems'].str.lower()
    
    binary_mets_df = pd.read_csv(paths['dir_binary_mets'], header=0)
    binary_mets_df['MetaboliteNames'] = binary_mets_df['MetaboliteNames'].fillna(binary_mets_df['Metabolites'])
    
    binary_genes_df = pd.read_csv(paths['dir_binary_genes'], header=0)
    binary_genes_df['GeneShortNames'] = binary_genes_df['GeneShortNames'].fillna(binary_genes_df['Genes'])
    
    # Apply subsystem conversion
    binary_rxns_df = apply_subsystem_conversion(binary_rxns_df)
    
    print("Data loading complete!")
    
    return {
        'binary_rxns': binary_rxns_df,
        'binary_mets': binary_mets_df,
        'binary_genes': binary_genes_df,
        'paths': paths
    }

class GradientLineHandler:
    """Handler to create gradient lines in matplotlib legends"""
    def __init__(self, colors, n_segments=10):
        self.colors = colors
        self.n_segments = n_segments

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        
        # Create line segments
        x = np.linspace(x0, x0 + width, self.n_segments)
        y = np.ones(self.n_segments) * (y0 + height/2)
        
        # Create line collection with gradient
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        # Create gradient colors
        cmap = LinearSegmentedColormap.from_list("custom", self.colors)
        colors = [cmap(i/self.n_segments) for i in range(self.n_segments-1)]
        
        lc = LineCollection(segments, colors=colors, linewidth=12)
        handlebox.add_artist(lc)
        return lc


# ===== INITIALIZATION =====
# Apply styling automatically when imported
set_figure_style()

# Load data automatically when imported
DATA = load_model_tables()

print("Data loaded into DATA dictionary.")