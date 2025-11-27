

# -*- coding: utf-8 -*-
"""
Sam lee

This script is used to conduct segmentations of all fovs
2024/01/03

"""
import os
import warnings
import matplotlib.pyplot as plt
import numpy as np
import skimage.io as io
import xarray as xr
from alpineer import io_utils
from ark.segmentation import marker_quantification, segmentation_utils
from ark.utils import (deepcell_service_utils, example_dataset,
                       plot_utils)

base_dir = r""
file_path = ''
result_path = file_path

cell_table_dir = os.path.join(result_path, "cell_table")
deepcell_input_dir = os.path.join(result_path, "deepcell_input")
deepcell_output_dir = os.path.join(result_path, "deepcell_output")
deepcell_visualization_dir = os.path.join(result_path, "deepcell_visualization")

for directory in [cell_table_dir, deepcell_input_dir, deepcell_output_dir, deepcell_visualization_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)
        
# input data file path
tiff_dir = base_dir
fovs = io_utils.list_folders(tiff_dir)

nucs = ['DAPI']
mems = ['CD45','CD20','CD3e','CD8','CD4','CD14','CD66b','CD11c']

# generate and save deepcell input tiffs
# set img_sub_folder param to None if the image files in tiff_dir are not in a separate sub folder 
deepcell_service_utils.generate_deepcell_input(
    deepcell_input_dir,
    tiff_dir,
    nucs,
    mems,
    fovs,
    img_sub_folder=None
)

# Mesmer was trained on data acquired at 20X resolution. 
#If your image data was acquired at a different resolution,
# you will get the best performance by rescaling. 
#The rescale factor will increase or decrease the image resolution by the value you provide.
# For example, if you data was acquired at 10X, use a `rescale_factor` of 2.
# If your data was acquired at 60X resolution, use a `rescale_factor` of 0.33.

rescale_factor = 1.0
deepcell_service_utils.create_deepcell_output(deepcell_input_dir, 
                                              deepcell_output_dir, 
                                              fovs=fovs, 
                                              scale=rescale_factor)

# save the overlaid segmentation labels for each fov (these will not display, but will save in viz_dir)
segmentation_utils.save_segmentation_labels(
    segmentation_dir=deepcell_output_dir,
    data_dir=deepcell_input_dir,
    output_dir=deepcell_visualization_dir,
    fovs=io_utils.remove_file_extensions(fovs),
    channels=['nuclear_channel', 'membrane_channel']
)

# set to True to add nuclear cell properties to the expression matrix
nuclear_counts = False

# set to True to bypass expensive cell property calculations
# only cell label, size, and centroid will be extracted if True
fast_extraction = False

# now extract the segmented imaging data to create normalized and transformed expression matrices
# note that if you're loading your own dataset, please make sure all the imaging data is in the same folder
# with each fov given its own folder and all fovs having the same channels
#The bactch size means run rov numbers simutaneously
cell_table_size_normalized, cell_table_arcsinh_transformed = \
    marker_quantification.generate_cell_table(segmentation_dir=deepcell_output_dir,
                                              tiff_dir=tiff_dir,
                                              img_sub_folder=None,
                                              fovs=fovs,
                                              batch_size=5,
                                              nuclear_counts=nuclear_counts,
                                              fast_extraction=fast_extraction)

compression = None
cell_table_size_normalized.to_csv(os.path.join(cell_table_dir, 'cell_table_size_normalized.csv'),
                                  compression=compression, index=False)
cell_table_arcsinh_transformed.to_csv(os.path.join(cell_table_dir, 'cell_table_arcsinh_transformed.csv'),
                                      compression=compression, index=False)






