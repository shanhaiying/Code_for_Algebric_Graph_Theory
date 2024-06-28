#设置matplotlib的中文字体支持，
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Microsoft YaHei']})

from sage.misc.sage_ostools import *
os.environ["PATH"]="/Users/haiyingshan/.sage/local/bin:/usr/local/bin:"+os.environ["PATH"]

latex.engine("xelatex")
latex.extra_preamble(r"\usetikzlibrary{hobby}")



