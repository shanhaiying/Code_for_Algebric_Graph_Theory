# GraphDisplay.py

import matplotlib.pyplot as plt
import networkx as nx
from functools import wraps
from io import BytesIO
import base64

# 定义图形显示配置类
class GraphDisplayOptions:
    def __init__(self):
        self.figsize = (3, 3)
        self.node_size = 90
        self.edge_color = 'gray'
        self.node_color = 'skyblue'
        self.font_size = 8
        self.font_weight = 'bold'

# 初始化全局配置实例
display_options = GraphDisplayOptions()

# 函数：将NetworkX图转换为Base64编码的图像
def plot_to_base64(graph, options=display_options):
    """Convert a NetworkX graph plot to a base64 encoded image."""
    buffer = BytesIO()
    pos = nx.spring_layout(graph)
    plt.figure(figsize=options.figsize)
    nx.draw(graph, pos, 
            with_labels=True, 
            node_size=options.node_size, 
            node_color=options.node_color, 
            edge_color=options.edge_color, 
            font_size=options.font_size, 
            font_weight=options.font_weight)
    plt.axis('off')  
    plt.savefig(buffer, format='png', bbox_inches='tight', pad_inches=0)
    plt.close()
    buffer.seek(0)
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    return f'data:image/png;base64,{img_base64}'

# 装饰器：自定义NetworkX图的_repr_html_方法
def custom_repr_decorator(func):
    @wraps(func)
    def wrapper(self):
        img_tag = plot_to_base64(self)
        return f'<img src="{img_tag}" alt="Graph Visualization">'
    return wrapper

# 全局标志：控制是否使用自定义_repr_html_
use_custom_repr_html = True
if use_custom_repr_html:
    original_repr_html = getattr(nx.Graph, '_repr_html_', None)
    nx.Graph._repr_html_ = custom_repr_decorator(original_repr_html) if original_repr_html else custom_repr_decorator(lambda x: "")
else:
    nx.Graph._repr_html_ = original_repr_html