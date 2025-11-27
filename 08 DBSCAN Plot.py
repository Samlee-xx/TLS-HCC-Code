
import os
import pandas as pd
from alpineer import io_utils
from ark.utils.plot_utils import cohort_cluster_plot, color_segmentation_by_stat
import natsort as ns
import ark.settings as settings
import colorcet as cc
import colorcet
from matplotlib import colormaps


for slide_name in slide:
    one_slide = data[data['slide'] == slide_name]
    unique_fovs = one_slide['fov'].unique()

    replace_dict = {0: 'cluster1', -1: 'others', 1: 'cluster2', 2: 'cluster3', 
                    3: 'cluster4', 4: 'cluster5', 5: 'cluster6', 
                    6: 'cluster7', 7: 'cluster8', 8: 'cluster9', 9: 'cluster10',
                    10: 'cluster11'}
    one_slide['cluster_label'] = one_slide['cluster_label'].replace(replace_dict)

    cell_cluster = one_slide
    cell_cluster = one_slide[['fov','label','cluster_label']]
    cell_cluster.rename(columns={'cluster_label': 'cell_meta_cluster'}, inplace=True)
    
    image_dir = r'...\02 Qupath Output\{}'.format(slide_name)
    seg_dir = r'...\03 Segmentation\{}\deepcell_output'.format(slide_name)
    save_dir = r'...\01 DBSCAN Plot 80\{}'.format(slide_name)
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir) 

    fovs = ns.natsorted(io_utils.list_folders(image_dir))
    cohort_cluster_plot(
        fovs=fovs,
        seg_dir=seg_dir,
        save_dir= save_dir,
        cell_data= cell_cluster,
        erode=True,
        fov_col=settings.FOV_ID,
        label_col=settings.CELL_LABEL,
        cluster_col=settings.CELL_TYPE,
        seg_suffix="_whole_cell.tiff",
        cmap="cet_glasbey_bw_minc_20",
        display_fig=True,
        fig_file_type="png",
        figsize=(10, 10),
        dpi=300,
    )
    
