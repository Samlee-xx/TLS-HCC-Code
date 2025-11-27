

# -*- coding: utf-8 -*-
"""
20240617
Lishuai

@author: DELL
"""

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import networkx as nx
import matplotlib.cm as cm
data = pd.read_csv(r'E:\01-TLS\02-CODEX\11-TA TLS Analysis\TA_TLS_Metadata.csv')
tls_id = data['tls_id'].unique()

#123942_L_11_TLS0
#176552_L_2_TLS3
#136759_N_6_TLS11
#176552_T_2_TLS0
final_boundary = pd.DataFrame(columns=['Level', 'CellID'])
for tls_name in tqdm(tls_id):
    test = data[data['tls_id']== tls_name]
    test_coord = test[['centroid.0','centroid.1']].values
    celltypes = test[['celltype']].values
    cellid = test[['cellid']].values
    tri = Delaunay(test_coord)
    triangles = tri.simplices

    max_distance = 40
    for triangle in tri.simplices:
        for i in range(3):
            for j in range(i+1,3):
                if np.linalg.norm(test_coord[triangle[i]] - test_coord[triangle[j]]) < max_distance:
                    plt.plot([test_coord[triangle[i],0],test_coord[triangle[j],0]],
                             [test_coord[triangle[i],1],test_coord[triangle[j],1]], linewidth = 0.1)

    colors = plt.cm.tab10(np.arange(len(test['celltype'].unique())))

    for i, cell_type in enumerate(test['celltype'].unique()):
        plt.scatter(test[test['celltype']==cell_type]['centroid.0'],
                    test[test['celltype']==cell_type]['centroid.1'],
                    s=0.5, label = cell_type, color = colors[i])
    
    G = nx.Graph()
    # 遍历三角形并在距离小于max_distance时添加边
    for triangle in tri.simplices:
        for i in range(3):
            for j in range(i+1, 3):
                distance = np.linalg.norm(test_coord[triangle[i]] - test_coord[triangle[j]])
                if distance < max_distance:
                    G.add_edge(triangle[i], triangle[j])
    
    radius = 45 # 根据需要调整半径大小
    neighbor_count_threshold =15 # 根据需要调整邻居数量阈值
    # 计算每个点周围一定半径内的邻居数量
    neighbor_counts = []
    for i in G.nodes():
        neighbors_within_radius = [j for j in G.nodes() if i != j and np.linalg.norm(test_coord[i] - test_coord[j]) <= radius]
        neighbor_counts.append(len(neighbors_within_radius))

    # 根据半径内邻居数量筛选节点
    filtered_nodes_by_neighbor_count = [node for i, node in enumerate(G.nodes()) if neighbor_counts[i] < neighbor_count_threshold]

    filtered_nodes_by_neighbor_count = list(set(filtered_nodes_by_neighbor_count))
    boundary_nodes = filtered_nodes_by_neighbor_count
    final_boundary_nodes = [boundary_nodes]  # 存储不同级别的边界节点

    # 定义扩展的次数
    expansion_iterations = 100

    for _ in range(expansion_iterations):
        new_boundary_nodes = []
        for node in boundary_nodes:
            neighbors = list(G.neighbors(node))
            new_boundary_nodes.extend(neighbors)
            new_boundary_nodes = list(set(new_boundary_nodes))
        
        # 将新找到的边界点去重并排除前一层次的边界节点
        new_boundary_nodes = list(set(new_boundary_nodes) - set(boundary_nodes))
        
        final_boundary_nodes.append(new_boundary_nodes)
        
        # 将新边界点合并到当前的边界节点列表中
        boundary_nodes.extend(new_boundary_nodes)
        boundary_nodes = list(set(boundary_nodes))
        
        if not new_boundary_nodes:
        # 如果new_boundary_nodes为空，跳出for循环
            break

    #core_nodes = list(set(G.nodes()) - set(boundary_nodes))

    #num_nodes = G.number_of_nodes()
    # 最终得到不同级别的边界节点，存储在 final_boundary_nodes 中
    
    #以下是画图的部分
    plt.figure(figsize=(10, 10))

    pos = {node: (test_coord[node, 0], test_coord[node, 1]) for node in G.nodes()}

    # 定义颜色映射
    colors = cm.rainbow(np.linspace(0, 1, len(final_boundary_nodes)))

    # 可视化不同级别的边界节点
    for i, level_nodes in enumerate(final_boundary_nodes):
        nodes = nx.draw_networkx_nodes(G, pos, nodelist=level_nodes, node_color=colors[i], node_size=50, label=f'Level {i+1}')
        edges = nx.draw_networkx_edges(G, pos, edgelist=G.edges(), width=0.5, edge_color='gray')

   # nodes = nx.draw_networkx_nodes(G, pos, nodelist=core_nodes, node_color='yellow', node_size=50, label= 'Core Nodes')
    # 添加图例
    plt.legend(loc='upper right')
    plt.title("Boundary Nodes Visualization")
    plt.savefig(r'E:\01-TLS\02-CODEX\14-Delaunay\01-Plot\{}.png'.format(tls_name))
    plt.show()
    
    # 将final_boundary_nodes中的不同级别的边界节点与cellid进行匹配
    boundary_data = []
    for i, level_nodes in enumerate(final_boundary_nodes):
        for node in level_nodes:
            cell_id = cellid[node][0]  # 获取对应节点的cellid
            boundary_data.append({'Level': i + 1, 'CellID': cell_id})

    # 创建一个包含边界信息的DataFrame
    boundary_df = pd.DataFrame(boundary_data)
    final_boundary = pd.concat([final_boundary,boundary_df])
                    
final_boundary.to_csv(r'E:\01-TLS\02-CODEX\14-Delaunay\02-Delaunay Out\TLS Spatial Layer labels.csv')















