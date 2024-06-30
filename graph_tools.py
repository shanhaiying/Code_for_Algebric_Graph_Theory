import networkx as nx
import matplotlib.pyplot as plt

def plot_graph_with_labels(G, labels, **kwargs):
    """
    绘制图并在节点附近显示给定的向量标签。

    参数:
    G (networkx.Graph): 要绘制的图
    labels (list or array): 节点对应的向量
    kwargs: 其他绘图参数
        figsize (tuple): 图的尺寸
        vertex_size (int): 节点大小
        vertex_labels (bool): 是否显示节点标签
        vertex_color (str): 节点颜色
        title (str): 图的标题
    """
    # 解析 kwargs 以设置绘图参数
    figsize = kwargs.get('figsize', (8, 6))
    vertex_size = kwargs.get('vertex_size', 700)
    vertex_labels = kwargs.get('vertex_labels', True)
    vertex_color = kwargs.get('vertex_color', "skyblue")
    title = kwargs.get('title', "Graph")
    
    # 设置顶点标签
    labels = {i: f'{labels[i]}' for i in range(len(labels))}
    
    # 生成布局并绘制图
    pos = nx.spring_layout(G)
    plt.figure(figsize=figsize)
    nx.draw(G, pos, with_labels=vertex_labels, node_color=vertex_color, node_size=vertex_size, 
            font_size=12, font_weight='bold', edge_color='gray')
    
    # 在节点外部显示向量分量
    for node, (x, y) in pos.items():
        plt.text(x, y + 0.1, f'{labels[node]}', horizontalalignment='center', fontsize=12, fontweight='bold')
    
    plt.title(title)
    plt.show() 

# 扩展 NetworkX 的 Graph 类
nx.Graph.plot_graph_with_labels = plot_graph_with_labels
