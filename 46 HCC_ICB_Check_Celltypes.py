
import os
import pandas as pd
from alpineer import io_utils
from ark.utils.plot_utils import cohort_cluster_plot, color_segmentation_by_stat
import natsort as ns
import ark.settings as settings
import colorcet as cc
import colorcet
from matplotlib import colormaps
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

path_to_data = r'...\Celltype total metadata.csv'
cells = pd.read_csv(path_to_data)
data = cells

custom_cmap = pd.DataFrame({
    'cell_meta_cluster': ["Heptocyte","Endothelial","CD8T","Macrophage","Bile Duct","immune_others",
"CD4T","Neutrophil","Fibroblast","B","Lymphatic","Unassign"],
    'color': ["#caf0f6","#f29f05","#007f4e","#6c584c","#72b043","#757575",
               "#8a508f","#057dcd","#dfe2fe","#ed008c","#e6c994","#4C4C4D"]
})

unique_fovs = data['fov'].unique()

cell_cluster.cell_meta_cluster.unique()
cell_cluster = data
cell_cluster = data[['fov','label','celltype']]
cell_cluster.rename(columns={'celltype': 'cell_meta_cluster'}, inplace=True)
   
image_dir = r'...\Split_channels'
seg_dir = r'...\deepcell_output'
save_dir =r'...\CellType Plot'
    
if not os.path.exists(save_dir):
    os.makedirs(save_dir) 

cohort_cluster_plot(
    fovs=unique_fovs,
    seg_dir=seg_dir,
    save_dir= save_dir,
    cell_data= cell_cluster,
    erode=True,
    fov_col=settings.FOV_ID,
    label_col=settings.CELL_LABEL,
    cluster_col=settings.CELL_TYPE,
    seg_suffix="_whole_cell.tiff",
    cmap=custom_cmap,
    display_fig=True,
    fig_file_type="png",
    figsize=(10, 10),
    dpi=300,
)






