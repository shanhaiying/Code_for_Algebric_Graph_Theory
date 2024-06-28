import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import ipywidgets as widgets

from numpy import array
# 设置环境变量
import os
os.environ['OPENAI_API_KEY'] = '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'




# from ipywidgets import interact, interactive, fixed
# from IPython.display import clear_output, display, HTML, Image, Javascript,Math,SVG
load(DOT_SAGE+"/MyTools.sage")
load(DOT_SAGE+"/MyPolynomial.sage")
load(DOT_SAGE+"/GraphTheory.sage")
load(DOT_SAGE+"/myCombinatorics.sage")
load(DOT_SAGE+"/LinearAlgebra.sage")
load(DOT_SAGE+"/set_latex.sage")
load(DOT_SAGE+"/Tensor.sage")
load(DOT_SAGE+"/HyperMatrix.sage")
load(DOT_SAGE+"/MyUtility.sage")
load(DOT_SAGE+"/hypergraph_v2.sage")
load(DOT_SAGE+"/NameGraph.sage")

# 加载 Jupyter AI Magics 扩展
from IPython import get_ipython
ipython = get_ipython()

if 'jupyter_ai_magics' not in ipython.extension_manager.loaded:
    ipython.extension_manager.load_extension('jupyter_ai_magics')