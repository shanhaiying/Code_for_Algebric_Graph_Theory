import functools
from itertools import groupby
from sage.misc.html import HtmlFragment
# from sage.misc.temporary_file import tmp_filename, tmp_dir, delete_tmpfiles
from sage.misc.latex import _run_latex_,_Latex_prefs, png
from IPython.display import Image, display
import base64,shutil,os

def Latex2Image(Object):
    engine = _Latex_prefs._option["engine"]
    file_name = tmp_filename() + ".png"
    png(Object, file_name, debug = False, engine = engine)
    return Image(filename=file_name)



def Graph2Image(Object):
    from IPython.display import Image
    engine = _Latex_prefs._option["engine"]
    tempfilename = tmp_filename()
    png_file = tempfilename + ".png"
    tex_file = tempfilename + ".tex"
    latex_code = r'''
\documentclass[tikz,border=2pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{hobby}
\begin{document}
%s
\end{document}''' % (Object.Latex())
    open(tex_file, 'w').write(latex_code)
    e = _run_latex_(tex_file, png=True)
    if e.find("Error") == -1:
        return Image(filename=png_file)
    else:
        raise RuntimeError("编译错误，请检查TeX代码")



# Hypergraph._repr_html_=myreprhtml

def my_repr_html(G,**kwds):
    from sage.repl.rich_output import get_display_manager
    import base64
    dm = get_display_manager()
    css=""
    if "css" in kwds:
        css=kwds.pop("css")
    # plt=G.plot(**kwds)
    out=G._rich_repr_(dm,**kwds)
    if not isinstance(out, sage.repl.rich_output.output_graphics.OutputImagePng):
        plt = G.plot(**kwds)
        out = plt._rich_repr_(dm)
    #     为了兼容Python2和 Python3
    if "encodebytes" in dir(base64):    
        pngstr=str(base64.encodebytes(out.png.get()), 'utf-8')
    else:
        pngstr=out.png.get_str().encode('base64')
    return html('<img src="data:image/png;base64,{0}" {1}>'.format(pngstr,css))
Graph._repr_html_=my_repr_html
Graphics._repr_html_=my_repr_html
DiGraph._repr_html_=my_repr_html



def reshape(L,dims,flat=True,fil=0):
    u"""
    当数据不足时, 填充fil的值
    """
    if len(L)==0:
        l=[]
        for i in dims[::-1]:
            l=[l*i]
        return l[0]
    l=L
    if flat:
        l=flatten(L)
    l=l+[fil]*(mul(dims)-len(l))
    for i in dims[::-1]:
        l=map(list,zip(*[iter(l)]*i))
    return list(l)[0]


def Reshape(tb, **kwds):
    nc = kwds.pop('nc', len(tb))
    ntb = []
    TB = reshape(tb, [ceil(len(tb)/nc), nc], flat=False, fil="")
    tb = [flatten(t) for t in TB]
    return tb


# def Html_table_row(self, file, row, header=False, **kwds):
#     from sage.plot.all import Graphics
#     from sage.misc.latex import latex
#     from sage.misc.html import math_parse
#     from sage.graphs.generic_graph import GenericGraph
#     import types

#     if isinstance(row, types.GeneratorType):
#         row = list(row)
#         print("1")
#     elif not isinstance(row, (list, tuple)):
#         row = [row]
#         print("2")

#     column_tag = "<th>%s</th>\n" if header else "<td>%s</td>\n"

#     if self._options['header_column']:
#         first_column_tag = '<th class="ch">%s</th>\n' if header else '<td class="ch">%s</td>\n'
#         print("3")
#     else:
#         first_column_tag = column_tag
#         print("4")

#     # First entry of row:
#     entry = row[0]
#     if isinstance(entry, (HyperG, GenericGraph, Graphics)):
#         print("5")
#         # file.write(first_column_tag % math_parse(entry))
#         file.write(first_column_tag % entry._repr_html_())
        
#     elif isinstance(entry, str):
#         print("6")
#         if bool(re.search('_|^|\\\\', entry)) and not bool(re.search('$', entry)):
#             entry = "$"+entry+"$"
#             print("7")
#         file.write(first_column_tag % math_parse(entry))
#         print("8")
#     else:
#         file.write(first_column_tag %
#                     ('<script type="math/tex">%s</script>' % latex(entry)))
#         print("9")

#     # Other entries:
#     for column in range(1, len(row)):
#         if isinstance(row[column],  (HyperG, GenericGraph, Graphics)):
#             file.write(column_tag % row[column]._repr_html_())
#             print("10")
#         elif isinstance(row[column], str):
#             if bool(re.search('_|^|\\\\', row[column])) and not bool(re.search('$', row[column])):
#                 row[column] = "$"+row[column]+"$"
#                 print(11)
#             file.write(column_tag % math_parse(row[column]))
#             print("12")
#         else:
#             print("13")
#             file.write(column_tag % (
#                 '<script type="math/tex">%s</script>' % latex(row[column])))


# table._html_table_row = Html_table_row


def Table(tb, **kwds):
    from pandas.io.formats.style import Styler
    def Proc(it):
        if isinstance(
                it, (HyperG, OrderedTree, GenericGraph, Graphics)):
            return it.plot()._repr_html_()
        elif isinstance(it, Styler):
            return it.render()
        else:
            return it
    islist=mul([isinstance(it,list) for it in tb])
    if ('nc' not in kwds) and islist==1:
        Tb = [[Proc(it) for it in row] for row in tb]
        return table(Tb, **kwds)
    else:
        nc = kwds.pop('nc', len(tb))
        Tb = flatten(tb,max_level=1)
        tb = [Proc(it) for it in Tb]
        TB = reshape(tb, [ceil(len(tb)/nc), nc], flat=False, fil="")
        return table(TB, **kwds)

def viewgraph(self, h=400, w=800):
    import random
    import string
    from sage.graphs.graph_plot_js import gen_html_code
    self.name(''.join(random.choice(string.lowercase) for x in range(6)))
    filename = gen_html_code(self)
    f = file(filename)
    data = f.readlines()
    f.close()
    out = file(self.name()+".html", 'w')
    out.writelines(data)
    out.close()
    return html('<iframe seamless frameborder="0" height="%d" width="%d" src="%s"></iframe>' % (h, w, "./"+self.name()+".html"))


Graph.viewgraph = viewgraph


def htmlview(Ts, G=True, inv=None, html=False, **options):
    import numpy as np
    A = np.array(Ts)
    gL,EnL=[list(A[:,i]) for i in [0,1]]
    if inv==None:
        inv="不变量"
    if G:
        GL=[Graphics2Data(G,img=False,**options) for G in gL]
    else:
        GL=[g.name() for g in gL]
    Tab=table([range(len(GL)),GL,EnL],header_column=["序号","图",inv]).transpose()
    if not html:
        pretty_print(sage.misc.html.html(Tab))
    else:
        return Tab
graphs_list.htmlview=htmlview


def showGA(GL,rc,fs,layout="spring"):
    GA=graphics_array([plot(G,layout=layout,vertex_size=20,graph_border=True, vertex_labels=False) for G in GL],rc[0],rc[1])
    show(GA,figsize=fs)


def latex_draw(GList,columns=4):
    import random
    import string
    s=string.lowercase+string.digits
    namelist=[]
    for i in range(len(GList)):
        filename=''.join(random.sample(s,6))
        namelist.append(filename)
        latex.eval('''
          '''+latex(GList[i]), locals(), filename=DATA+filename)
    imglst=["<img src=data/%s.png></img>"%fn for fn in namelist]
    pretty_print(table([imglst[i:i+columns] for i in xrange(0,len(imglst),columns)]))
    sleep(0.5)  #等待图像完全显示后将图像文件删除!
    #res=os.system("rm data/*")







def Graphics2Data(P,img=True,**kwds):
    from IPython.display import Image
    import PIL as pil
    imgdata = io.StringIO()
    # if type(P) in [DiGraph, Graph, sage.graphs.graph_plot.GraphPlot]:

    if isinstance(P,(Hypergraph)):
        if "html" in dir(P):
            return P.html
        else:
            labels=True
            if "labels" in  kwds:
                labels=kwds["labels"]
        return P.plot(labels=labels)
    elif isinstance(P,(Graph,DiGraph, sage.graphs.graph_plot.GraphPlot)):
        P=P.plot(**kwds)

    if type(P)==pil.Image.Image:
        P.save(imgdata, format='png')
        data=imgdata.getvalue()
    else:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        options = dict()
        try:
            options.update(P.SHOW_OPTIONS)
            options.update(P._extra_kwds)
            options.update(kwds)
            options.pop('layout')
            dpi = options.pop('dpi')
            transparent = options.pop('transparent')
            fig_tight = options.pop('fig_tight')
        except:
            pass
        fig_tight = options.pop('fig_tight')
        transparent = options.pop('transparent')
        dpi = options.pop('dpi')
        figure = P.matplotlib(**options)
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        figure.set_canvas(FigureCanvasAgg(figure))
        figure.tight_layout()
        opts = dict(dpi=dpi, transparent=transparent)
        figure.savefig(imgdata, format='png',**opts)
        data=imgdata.getvalue()
    if img:
        return Image(data=data)
    else:
        return '<img class="tongji_img", src="data:image/png;base64,{0}">'.format(data.encode('base64').replace('\n', ''))



def graph_to_html(self,**kwds):
    import StringIO
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    imgdata = StringIO.StringIO()
    P=self.plot(**kwds)
    options = dict()
    options.update(P.SHOW_OPTIONS)
    options.update(P._extra_kwds)
    dpi = options.pop('dpi')
    transparent = options.pop('transparent')
    fig_tight = options.pop('fig_tight')
    figure = P.matplotlib(**options)
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    figure.set_canvas(FigureCanvasAgg(figure))
    figure.tight_layout()
    opts = dict(dpi=dpi, transparent=transparent)
    figure.savefig(imgdata, format='png', **opts)
    return '<img src="data:image/png;base64,{0}">'.format(imgdata.getvalue().encode('base64').replace('\n', ''))
Graph.graph_to_html=graph_to_html
DiGraph.graph_to_html=graph_to_html



def graphs_list_save(GL,enum=True,**kwds):
    import StringIO
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    body=[]
    for i in range(len(GL)):
        imgdata = StringIO.StringIO()
        P=GL[i].plot(**kwds)
        options = dict()
        options.update(P.SHOW_OPTIONS)
        options.update(P._extra_kwds)
        dpi = options.pop('dpi')
        transparent = options.pop('transparent')
        fig_tight = options.pop('fig_tight')
        figure = P.matplotlib(**options)
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        figure.set_canvas(FigureCanvasAgg(figure))
        figure.tight_layout()
        opts = dict(dpi=dpi, transparent=transparent)
        figure.savefig(imgdata, format='png', **opts)
        body.append(['<img src="data:image/png;base64,{0}">'.format(imgdata.getvalue().encode('base64').replace('\n', ''))])
    return body

# def GraphList2TeX(gL,nc,filename="GraphList"):
#     GL=copy(gL)

#     s=r'''
#     % !TEX program = xelatex
#     \documentclass{ctexart}
#     \usepackage{subcaption}
#     \usepackage{tikz}
#     \usepackage{subcaption}
#     \begin{document}

#     \begin{figure}
#     \centering
#     '''
#     nc=3
#     for i in range(len(gL)):
#         G=GL.pop()
#         s+=r'''
#     \begin{subfigure}[b]{'''+str(round(1/3,2))+r'''\textwidth}
#     \centering'''
#         s+=G._latex_()
#         s+=r'''\caption{$T_{1}$}
#     \end{subfigure}\hfill'''

#         if i+1 % nc==0:
#             s+=r"\\[3mm]"+"\n"

#     s+=r'''
#     \caption{SageMath生成图}
#     \end{figure}

#     \end{document}
#     '''

#     f=open(filename+".tex","w")
#     f.write(s)
#     f.close()


# def Graph2TeX(G,filename="Graph"):   
#     pos=G.Pos()
#     tex=r'''
#     % !TEX program = xelatex
#     \documentclass{ctexart}
#     \usepackage{subcaption}
#     \usepackage{tikz}
#     \usepackage{subcaption}
#     \begin{document}

#     \begin{figure}
#     \centering
#     '''
#     tex += r'''
#     \begin{tikzpicture}[scale=3,
#     Node/.style={fill,circle,scale=.5}
#     ]'''
#     for i in G.vertices():
#         tex += r'''
#     \coordinate (v{0}) at {1};'''.format(i, tuple(map(lambda p:round(p,2),pos[i])))
#     for v in G.vertices():
#         label = "label={90:$"+latex(v)+"$}"
#         tex += r'''
#     \draw node[Node,label={90:$v_{%s}$}] at (v%s){{}};'''%(v,v)
#     #     tex += str(v)+r'''}$} ] at (v{0}){{}};'''.format(v)

#     for u,v in edges:
#          tex += r"""
#     \draw[color=black] (v{0})-- (v{1});""".format(u,v)

#     tex += r'''
#     \end{tikzpicture}
#     '''  
#     tex+=r'''
#     \caption{SageMath生成图}
#     \end{figure}

#     \end{document}
#     '''
#     f = open(filename+".tex", "w")
#     f.write(tex)
#     f.close()
#     return None
        
# Graph.Graph2TeX=Graph2TeX

#返回变量名
def variable_name(p):
    for name, value in globals().items():
        if id(value) == id(p):
            return name


#对列表元素按 method 的返回值进行分组.
def GroupBy(GL, method):
    Grps = groupby(sorted(GL, key=method), method)
    return dict([(key, list(group)) for key,group in Grps])


import numpy as np
# from scipy import interpolate
class TCBSpline():
	def __init__(self):
		self.c = 0.7
		self.b = 0
		self.t = 0.7
		self.ControlPoints = []
		self.subpoints = []

	@staticmethod
	def _calc_tangents_kochanek_bartel(control_points, t, c, b):
		
		tans = []
		tand = []

		cona = (1 - t) * (1 + b) * (1 - c) * 0.5
		conb = (1 - t) * (1 - b) * (1 + c) * 0.5
		conc = (1 - t) * (1 + b) * (1 + c) * 0.5
		cond = (1 - t) * (1 - b) * (1 - c) * 0.5

		for i in range(1, len(control_points) - 1):
			pa = control_points[i - 1]
			pb = control_points[i]
			pc = control_points[i + 1]

			x1 = pb[0] - pa[0]
			y1 = pb[1] - pa[1]
			# z1 = pb[2] - pa[2]
			x2 = pc[0] - pb[0]
			y2 = pc[1] - pb[1]
			# z2 = pc[2] - pb[2]

			tans.append((cona * x1 + conb * x2, cona * y1 + conb * y2))  # cona*z1+conb*z2
			tand.append((conc * x1 + cond * x2, conc * y1 + cond * y2))  # conc*z1+cond*z2

		return tans, tand

	@staticmethod
	def _calc_tangents_catmull_rom(control_points):

		tans = []
		for i in range(1, len(control_points) - 1):
			pa = control_points[i - 1]
			pb = control_points[i]
			pc = control_points[i + 1]

			x1 = pb[0] - pa[0]
			y1 = pb[1] - pa[1]
			# z1 = pb[2] - pa[2]
			x2 = pc[0] - pb[0]
			y2 = pc[1] - pb[1]
			# z2 = pc[2] - pb[2]

			tans.append((0.5 * (x1 + x2), 0.5 * (y1 + y2)))  # 0.5 * (z1 + z2)

		return tans

	@staticmethod
	def interpolate(control_points, t, c, b, closed= True, t_inc = 0.01):
		if closed:
			control_points = [control_points[-1]] + control_points + [control_points[0]]
		else:
			control_points = [control_points[0]] + control_points + [control_points[-1]]

		tans, tand = TCBSpline._calc_tangents_kochanek_bartel(control_points, t, c, b)
		# tans = tand = TCBSpline._calc_tangents_catmull_rom(control_points)

		if closed:
			control_points.append(control_points[2])
			tans.append(tans[0])
			tand.append(tand[0])

		final_lines = []
		for i in range(1, len(control_points) - 2):
			p0 = control_points[i]
			p1 = control_points[i + 1]
			m0 = tand[i-1]
			m1 = tans[i]
			# interpolate curve from p0 to p1
			final_lines.append((p0[0], p0[1]))
			t_iter = t_inc
			while t_iter < 1.0:
				t_iter_2 = t_iter ** 2
				t_iter_3 = t_iter ** 3

				h00 = 2*t_iter_3 - 3*t_iter_2 + 1
				h10 = 1*t_iter_3 - 2*t_iter_2 + t_iter
				h01 = -2*t_iter_3 + 3*t_iter_2
				h11 = 1*t_iter_3 - 1*t_iter_2
				px = h00*p0[0] + h10*m0[0] + h01*p1[0] + h11*m1[0]
				py = h00*p0[1] + h10*m0[1] + h01*p1[1] + h11*m1[1]
				#pz = h00*p0[2] + h10*m0[2] + h01*p1[2] + h11*m1[2]

				final_lines.append((px, py))
				t_iter += t_inc

			final_lines.append((p1[0], p1[1]))

		return final_lines





#将列表显示为每行nc个项的图表
# def reshape(L,nc):
#     nL=L+["</pre>"]*(bool(len(L)%nc<>0)*nc-(len(L)%nc))
#     return [nL[i:i+nc] for i in range(0,len(nL),nc)]


# def graphs_list_show(GL,enum=True,T=True,html=False,**options):
#     Tb=graphs_list_save(GL,enum=enum,**options)

#     Tab=table(Tb)
#     if T:
#         Tab=Tab.transpose()
#     if html:
#         return Tab._html_()
#     else:
#         pretty_print(sage.misc.html.html(Tab))



# def graphs_list_show(GL,enum="no",colnum=1,**options):
#     n=len(GL)
#     if enum=="no":
#         return table(reshape([table([graph_to_html(GL[i],**options)],align="center")._html_() for i in range(n)],colnum))
#     if enum=="num":
#         nbc=range(n)
#     if enum=="name":
#         nbc=[G.name() for G in GL]
#     return table(reshape([table([nbc[i],graph_to_html(GL[i],**options)],align="center",frame=True)._html_() for i in range(n)],colnum))


def graphs_list_show(GL,enum="no",colnum=1,**options):
    n=len(GL)
    if enum=="no":
        return table(reshape([graph_to_html(GL[i],**options) for i in range(n)],colnum))
    if enum=="num":
        nbc=range(n)
    if enum=="name":
        nbc=[G.name() for G in GL]
    return table(reshape(flatten([[nbc[i],graph_to_html(GL[i],**options)] for i in range(n)]),colnum*2))


def base64ToPng(H,filepath):
    import base64
    from PIL import Image
    from io import BytesIO
    data=str(H.html)[33:-3]
    im = Image.open(BytesIO(base64.b64decode(data)))
    im.save('%s.png'%filepath, 'PNG')


def html_add(H1, H2):
    return HtmlFragment(str(H1)+"</br>"+str(H2))


HtmlFragment.__add__ = html_add


def table_add(T1, T2):
    return table(list(T1._rows)+list(T2._rows))


table.__add__ = table_add
table.T=table.transpose

import inspect

def getVarname(x, nsp=locals()):
    """
    获取变量对应的变量名
    
    参数：
    x: 变量，用于查找对应的变量名
    nsp: 命名空间字典，默认为 locals()，存储变量名和变量值的映射关系
    
    返回值：
    变量 x 对应的变量名
    
    注意：
    - 函数默认使用 locals() 函数获取当前命名空间的局部变量，但也可以通过显式传递 nsp 参数来指定其他命名空间字典。
    - 该函数依赖于变量的内存地址来进行匹配，因此只能在同一命名空间中进行操作。
    """
    return {id(var): name for name, var in nsp.items()}[id(x)]


def getVarname1(x):
    """
    获取变量对应的变量名
    
    参数：
    x: 变量，用于查找对应的变量名
    
    返回值：
    变量 x 对应的变量名
    """
    caller_frame = inspect.currentframe().f_back
    nsp = caller_frame.f_locals
    return {id(var): name for name, var in nsp.items()}[id(x)]


#import pandas as pd

def dict2htm(dt,head=None):
    """
    将字典转换为 HTML 表格
    
    参数：
    dt: 字典，要转换为表格的数据
    
    返回值：
    字典数据转换后的 HTML 表格
    
    注意：
    - 函数依赖于 getVarname 函数，该函数用于获取变量名。
    """
    # 将字典转换为包含单列的 DataFrame
    ht = pd.DataFrame([dt]).transpose()
    
    # 设置列名为变量名
    if head==None:
         head=getVarname(dt)   
    ht.columns = [head]
    
    # 将 DataFrame 转换为 HTML 表格并返回
    return ht.to_html()




@functools.lru_cache(maxsize=None)
def seq_iter(n, f, x0):
    """
    {x_i}:x_{i+1}=f(x_i), x_0=x0
    返回数列的第n项  
    
    """
    if n < 1:
        return x0
    nx = f(x0)
    return seq_iter(n-1, f, nx)
