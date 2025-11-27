

data = pd.read_csv(r'...\TLS Positive fov metadata.csv')
data = data[['label','centroid.0','centroid.1','fov','celltype']]
data = data[~data['celltype'].isin(['Heptocyte'])]

unique_fovs = data['fov'].unique()
data['cluster_label'] = -2

for fov in unique_fovs:
    test_data = data[data['fov'] == fov]
    X = test_data[['centroid.0', 'centroid.1']].values
    dbscan = DBSCAN(eps=35, min_samples=50)
    labels = dbscan.fit_predict(X)
    data.loc[data['fov'] == fov, 'cluster_label'] = labels
    
replace_dict = {0: 'cluster1', -1: 'others', 1: 'cluster2', 2: 'cluster3', 
                    3: 'cluster4', 4: 'cluster5', 5: 'cluster6', 
                    6: 'cluster7', 7: 'cluster8', 8: 'cluster9', 9: 'cluster10',
                    10: 'cluster11',11:'cluster12',12:'cluster13',13:'cluster14',
                    14:'cluster15'}
data['cluster_label'] = data['cluster_label'].replace(replace_dict)
cell_cluster = data
cell_cluster = data[['fov','label','cluster_label']]
cell_cluster.rename(columns={'cluster_label': 'cell_meta_cluster'}, inplace=True)
    
image_dir = r'...\Split_channels'
seg_dir = r'...\Segmentation_Out\deepcell_output'
save_dir =r'...\Analysis\DBSCAN Plot\eps 35'
    
if not os.path.exists(save_dir):
    os.makedirs(save_dir) 

data.to_csv(r'...\Analysis\TLS FOV DBSCAN.csv')
custom_cmap = pd.DataFrame({
    'cell_meta_cluster': ['cluster1', 'cluster3', 'cluster4', 'cluster2', 'cluster5',
                          'cluster6', 'cluster7', 'cluster8', 'cluster9', 'cluster10', 
                          'cluster11', 'cluster12','others','cluster13','cluster14','cluster15'],
    'color':  ["#caf0f6","#f29f05","#007f4e","#6c584c","#72b043","#757575",
               "#8a508f","#057dcd","#dfe2fe","#ed008c","#e6c994","#4C4C4D",'#4C4C4D',
               '#4C4C4D','#4C4C4D','#4C4C4D']
})
fovs =  data['fov'].unique()

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
    cmap=custom_cmap,
    display_fig=True,
    fig_file_type="png",
    figsize=(10, 10),
    dpi=300,
)


