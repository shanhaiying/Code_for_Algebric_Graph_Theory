import networkx as nx
import numpy as np
import sympy as sp
import pandas as pd
from sympy import Matrix
from IPython.display import display, Math
from networkx.algorithms.isomorphism import GraphMatcher
import matplotlib.pyplot as plt

def eigenpair(M, spe=False, order=0):
    """
    计算实对称阵的第 k 大特征值及其对应的特征向量
    
    参数：
    M (numpy.ndarray): 实对称矩阵
    spe (bool): 是否返回所有特征值
    order (int or str): 需要返回的特征值及其对应特征向量的顺序，0 表示最大，1 表示第二大，依此类推；
                        "all" 表示返回所有特征值和特征向量
    返回：
    tuple or list: 第 k 大特征值及其对应的特征向量，或所有特征值和特征向量
    """
    # 确保矩阵是对称的
    if not np.allclose(M, M.T):
        return "A symmetric matrix is required."
    
    # 计算特征值和特征向量
    w, v = np.linalg.eigh(M)
    
    # 处理特征向量和特征值
    v = np.around(v, decimals=12)  # 类似于 zero_at(1e-12)
    w = np.around(w, decimals=12)  # 类似于 zero_at(1e-12)
    
    # 排序特征值及其对应的特征向量
    idx = np.argsort(w)[::-1]
    w = w[idx]
    v = v[:, idx]
    
    if spe:
        return w
    if order == "all":
        return w, v
    else:
        # 获取第 order 大的特征值及其对应的特征向量
        rec = v[:, order]
        if rec[0] != 0:
            rec = np.sign(rec[0]) * rec
        return w[order], rec

def plot_graph_with_labels(G, eigenvector=('L', 0), **kwargs):
    # 解析 kwargs 以设置绘图参数
    figsize = kwargs.get('figsize', (8, 6))
    scale = kwargs.get('scale', 1)
    vertex_size = kwargs.get('vertex_size', 700)
    vertex_labels = kwargs.get('vertex_labels', True)
    vertex_color = kwargs.get('vertex_color', "skyblue")
    title = kwargs.get('title', "Graph")

    # 根据 eigenvector 参数选择计算拉普拉斯矩阵或邻接矩阵的特征向量
    if eigenvector[0] == 'L':
        matrix = nx.laplacian_matrix(G)
    elif eigenvector[0] == 'A':
        matrix = nx.adjacency_matrix(G)
    else:
        raise ValueError("eigenvector parameter must be ('L', 0) for Laplacian or ('A', 0) for adjacency matrix")

    # 计算特征值和特征向量
    spectral_radius, eigenvector = eigenpair(matrix.toarray())
    eigenvector *= scale

    # 设置顶点标签
    labels = {i: f'{eigenvector[i]:.2f}' for i in range(len(eigenvector))}

    # 生成布局并绘制图
    pos = nx.spring_layout(G)
    plt.figure(figsize=figsize)
    nx.draw(G, pos, with_labels=vertex_labels, node_color=vertex_color, node_size=vertex_size, 
            font_size=12, font_weight='bold', edge_color='gray')

    # 在节点外部显示特征向量分量
    for node, (x, y) in pos.items():
        plt.text(x, y + 0.1, f'{eigenvector[node]:.2f}', horizontalalignment='center', fontsize=12, fontweight='bold')

    plt.title(f"{title} (Spectral Radius: {spectral_radius:.2f})")
    plt.show()
    nx.Graph.eigenpair = eigenpair_for_graph
