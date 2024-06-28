# -*- coding: utf-8 -*-

#########################################################################################
#      Here is  an enhancement to the part of  graph theory within sage                 #
#                                                                                       #
#       Copyright (C) 2014~2023 HaiyingShan <shan_haiying@tongji.edu.cn>                #
#                             Last Modified: 2023-10-19                                  #
#########################################################################################

# from matplotlib.cbook import flatten
from six import iteritems
from sage.graphs.graph_plot import GraphPlot
import re
from sage.structure.element import is_Matrix
from sage.graphs.generic_graph import GenericGraph


DEFAULT_HyperGraphPLOT_OPTIONS = {
    "vertex_size": 20,
    "vertex_labels": True,
    "edge_style": 'solid',
    "edge_thickness": 1,
    "edge_color": 'black',
    "edge_colors": None,
    "edge_labels": False,
    "iterations": 50,
    "color_by_label": False,
    "graph_border": False,
    "facecolor": 'green',
    "edgecolor": 'green',
    "partition": None,
    "normalize": True,
    "dist": .075,
    "vlabelshift": [90, 0.05],
    # "vertex_label_dist": 0.05,
    # "vertex_label_orient": 90,
    'edgeplot': "TCB-TCBSplines",
    "max_dist": 1.5,
    "edge_labels_background": "white"
}


class edge(tuple):
    def __init__(self, e):
        self.__name__ = "edge"
        self.weight = 1
        self.vertices = list(e)
        self.order = len(e)
        self.index = -1
        self.label = "e%s" % self.index
        self.G = None

    def __str__(self):
        return tuple(self.vertices).__str__()

    def issuperset(self, sub):
        return set(self.vertices).issuperset(set(sub))

    def __sub__(self, otherset):
        return edge([i for i in self if i not in otherset])

    def intersection(self, other):
        return edge(Set(self).intersection(Set(other)))

    def __add__(self, other):
        return edge(set(Set(self)+(Set(other))))

    def __le__(self, other):
        return self.__str__() <= str(other)

    def __ge__(self, other):
        return self.__str__() >= str(other)

    def __lt__(self, other):
        return self.__str__() < str(other)

    def __gt__(self, other):
        return self.__str__() > str(other)

    def __and__(self, other):
        return edge(set(Set(self) & (Set(other))))

    def __getitem__(self, key):
        return self.vertices.__getitem__(key)


class HyperG(GenericGraph):
    def __init__(self, *args, **kwds):
        self.__name__ = "Hypergraph"
        if len(args) == 0:
            data = IncidenceStructure([])
        else:
            data = args[0]
        if isinstance(data, GenericGraph) and not isinstance(data, HyperG):
            self.Vertices = data.vertices()
            self.Edges = list(map(edge, data.edges(labels=False)))

        elif isinstance(data, HyperG):
            self.Vertices = data.Vertices
            self.Edges = data.Edges

        elif is_Matrix(data):
            data = IncidenceStructure(data)
            self.Vertices = data.ground_set()
            self.Edges = list(map(edge, data.blocks()))

        elif isinstance(data, IncidenceStructure):

            self.Vertices = data.ground_set()
            self.Edges = list(map(edge, data.blocks()))

        elif isinstance(data, Integer):
            self.Vertices =range(data)
            self.Edges=[]

        elif hasattr(data, "__getitem__") and isinstance(data[0], edge):
            data = IncidenceStructure(data)
            self.Vertices = data.ground_set()
            self.Edges = list(map(edge, data.blocks()))

        elif hasattr(data, "__iter__"):
            data = IncidenceStructure(data)
            self.Vertices = data.ground_set()
            self.Edges = list(map(edge, data.blocks()))

        else:
            raise TypeError("data type error!")
        

        self.options = deepcopy(DEFAULT_HyperGraphPLOT_OPTIONS)
        self.labels = True

        self.Data = dict({})
        self.drawn = False
        self._alpha = 0
        self.alpha = 0
        self._incidence_graph = self.incidence_graph()
        self.Degrees = self.degrees()
        self.vLabels = {v: v for v in self.Vertices}
        
        if self.Edges != []:
            self.rank = max(map(len, self.Edges))
            self.corank = min(map(len, self.Edges))
            # self.reIndex()
            self.Pos(reset=True)
            self._edge_coloring = self.edge_coloring()
            
        else:
            
            if "V" in kwds:
                self.Vertices = kwds.pop('V')

            else:
                if not hasattr(self,"Vertices"):
                    self.Vertices = []
            self._pos = {v: vector([v, 0]) for v in self.Vertices}
            self.rank = 0
            self.corank = 0
            

        self._vLabPos = {}
        self.reIndex()

    def __repr__(self):
        return self.name()

    def __copy__(self):
        G = HyperG()
        keys = ['__name__', 'Vertices', 'Edges', 'options', 'labels', 'Data', '_alpha', 'alpha', 'rank',
                'corank', '_incidence_graph', '_pos', '_edge_coloring', 'vLabels', '_vLabPos']
        for key in keys:
            if key in self.__dict__:
                setattr(G, key, copy(self.__dict__[key]))
        return G

    def size(self):
        return len(self.Edges)

    def order(self):
        return len(self.Vertices)

    def vertices(self):
        return self.Vertices

    def reIndex(self):
        for i in range(self.size()):
            e = self.Edges[i]
            e.index = i
            e.G = self
            e.label = "e%s" % i

    def incidence_matrix(self):
        n = self.order()
        m = self.size()
        V = self.Vertices
        return matrix(n, m, lambda i, j: V[i] in self.Edges[j])

    def adjacency_matrix(self):
        return matrix(self.order(), lambda i, j: len(self.incidence_edges([i, j]))*(i != j))

    def distance_matrix(self):
        return Graph(self.adjacency_matrix()).distance_matrix()

    def distance(self, u, v):
        return Graph(self.adjacency_matrix()).distance(u, v)

    def all_paths(self, u, v):
        IG = self.incidence_graph()
        return IG.all_paths(u, v)

    def distance1(self, u, v):
        P = self.all_paths(u, v)
        return min([len(p[1::2]) for p in P])

    def degree(self, v):
        return self.incidence_graph().degree(v)

    def degrees(self):
        return {v: self.degree(v) for v in self.Vertices}

    def hyperedges(self):
        return self.Edges

    def intersection_graph(self, sizes=None, *arg):
        if sizes is None:
            sizes = PositiveIntegers()
        elif sizes in PositiveIntegers():
            sizes = (sizes,)

        def igedge(e, f):
            len_is = len(e & f)
            return len_is if e != f and len_is in sizes else 0
        return Graph([self.Edges, igedge], weighted=True)

    def incidence_graph(self):
        BG = BipartiteGraph()
        BG.add_vertices(self.Vertices, left=True)
        BG.add_vertices(self.Edges, right=True)
        BG.add_edges([(v, e) for v in BG.left for e in BG.right if v in e])
        return Graph(BG)

    def automorphism_group(self):
        from sage.groups.perm_gps.permgroup import PermutationGroup
        G = self.incidence_graph()
        part = [self.Vertices, self.Edges]
        ag = G.automorphism_group(partition=part)
        gens = [[tuple(c) for c in g.cycle_tuples() if c[0] in self.Vertices]
                for g in ag.gens()]
        return PermutationGroup(gens=gens, domain=self.Vertices)

    def is_connected(self):
        return self.incidence_graph().is_connected()

    def degree_sequence(self):
        return list(self.degrees().values())

    def __iter__(self):
        for b in self.Vertices:
            yield b

    def __getitem__(self, key):
        return self.incidence_edges(key)

    def is_uniform(self):
        return self.rank == self.corank

    def uniform(self):
        assert self.is_uniform() == True,  "This hypergraph is not uniform"
        return self.rank

    def incidence_edges(self, v):
        """
        Returns a list of edges incident to the vertices in v.

        Args:
            self: The graph object.
            v: A vertex or a collection of vertices.

        Returns:
            A list of edges incident to the vertices in v.
        """
        if not hasattr(v, "__iter__") or isinstance(v, (str, dict)):
            v = [v]
        v = Set(v)
        return [e for e in self.Edges if Set(e).issuperset(v)]

    def edge_coloring(self):
        return {}
        # return self.intersection_graph().coloring(algorithm="MILP")
   

    def _spring_layout(self):
        _ = self._incidence_graph.plot(iterations=50000, save_pos=True)

        return {k: vector([round(x, 3), round(y, 3)]) for k, (x, y) in self._incidence_graph.get_pos().items()}

    def Pos(self, reset=False):
        if not hasattr(self,"_pos"):
            reset=True
        elif self._pos=={}:
            reset=True
        if reset:
            pos = self.incidence_graph().layout_spring(iterations=500)
            self._pos = {k: vector([round(x, 3), round(y, 3)])
                         for k, (x, y) in pos.items()}
        for e in self.Edges:
            e.center = Mean([self._pos[v] for v in e])
        return self._pos

    # def set_vert_lab(self):
    #     deg, dist = self.options["vlabelshift"]
    #     for e in self.Edges:
    #         for v in e:
    #             if self.degree(v) == 1:
    #                 vc = self._pos[v]-self._pos[e]
    #                 arg = round((vc*vector([1, CC.0])).arg()*180/pi, 3)
    #                 self._vLabelShift[v] = [arg, dist]

    def setLposDeg1(self):
        """[summary]
        利用TCB样条线算法(又名Kochanek-Bartels样条线),
        计算超边中一度点标号位置.
        """
        deg, dist = self.options["vlabelshift"]
        vLvec = {}
        for e in self.Edges:
            if len(e) < 3:
                continue
            se = sorted(e, key=lambda v: arctan2(*(self._pos[v]-self._pos[e])))
            tans, tand = TCBSpline._calc_tangents_kochanek_bartel(
                [self._pos[v] for v in se]*2, 0.5, 0.5, 0.05)
            tans = list(map(vector, tans))
            tand = list(map(vector, tand))
            TM = matrix([((tans[k].normalized()-tand[k].normalized())*diagonal_matrix([1, -1])).normalized()
                         for k in range(e.order)])
            vLvec.update(dict(zip(se, TM.rows())))
        vLdeg1 = {v: -dist*vLvec[v]
                  for v in vLvec.keys() if self.degree(v) == 1}
        self._vLabPos.update(vLdeg1)

    def setLposDeg1(self):
        deg, dist = self.options["vlabelshift"]
        vLvec = {}
        z=exp(CC.0*deg*pi/180)
        x, y = z.real(), z.imag()
        R = matrix(RR, 2, 2, [x, y, -y, x])
        Pos = {v: vector(self._pos[v])*R for v in self}
        # self._pos = {k: Pos[k] for k in self}

    def setLabPos(self):
        try:
            deg, dist = self.options["vlabelshift"]
            c = dist*exp(deg*pi/180*CC.0)
            vlp = vector(RR, [c.real(), c.imag()])
            self._vLabPos = {v: vlp for v in self.Vertices}
            self.setLposDeg1()

        except KeyError:
            raise KeyError(
                "you have not specified positions for all the vertices")

    def delete_vertices(self, S,replace=True):
        V=self.Vertices
        nV=Set(V)-Set(S)
        E = self.Edges
        for v in S:
            E = Set(E)-Set(self.incidence_edges(v))
        IS=IncidenceStructure(points=nV,blocks=E)
        H=HyperG(IS)
        if replace:
            self.__init__(IS)
            return self
        return H

    def delete_edges(self, E,replace=True):
        V=self.Vertices
        Edges = Set(self.Edges)-Set(map(edge, E))
        IS=IncidenceStructure(points=V,blocks=Edges)
        H=HyperG(IS)
        if replace:
            self.__init__(IS)
            return self
        return H

    def add_edges(self, E, replace=True):
        Edges = Set(self.Edges)+Set(map(edge, E))
        V=Set(self.Vertices)+Set(flatten(Edges))
        IS=IncidenceStructure(points=V,blocks=Edges)
        H=HyperG(IS)
        if replace:
            self.__init__(IS)
            return self
        return H

    def disjoint_union(self, G2, relabel=False):
        return graphs_list.disjoint_union([self, G2], relabel=relabel)

    def __add__(self, other):
        if isinstance(other, HyperG):
            return self.disjoint_union(other, relabel=False)
        raise TypeError("adding {} and {} is not defined".format(
            type(self), type(other)))

    def orbit_subdivision(self):
        GP = self.automorphism_group()
        return GP.orbits()

        return newH

    def name(self, *args):
        if len(args) == 0:
            return self.__name__
        else:
            self.__name__ = args[0]

    def add_star(self, R, s):
        '''
        在超图$G$中添加s条包含点集R 的新边       

        '''
        R = set(R)
        k = self.rank-len(R)
        n = self.order()
        E = [set(ae) | R for ae in reshape(range(n, n+k*s), [s, k])]
        m = self.size()
        Edges = []
        for e in E:
            e = edge(e)
            e.index = m
            m = m+1
            Edges.append(e)
        Edges = Set(Edges)+Set(self.Edges)

        return HyperG(Edges)

    def add_stars(self, rootsk):
        H1 = self
        if isinstance(rootsk, dict):
            for it in [(k, rootsk[k]) for k in rootsk.keys()]:
                if not isinstance(it[0], (list, tuple)):
                    it = ([it[0]], it[1])
                H1 = H1.add_star(*it)
            return H1
        else:
            for k in rootsk:
                H1 = H1.add_star({k[0]}, k[1])
            return H1

    def edge_weight(self, e, X):
        r"""
        向量X所诱导的边权 $\prod_{v \in e}X[i]$        
        """
        w = 1
        V = self.vertices()
        index = dict(zip(V, range(self.order())))
        for i in e:
            w = w*X[index[i]]
        return w

    def charpoly_system(self, vari='x', D=False):
        X = [var(vari + str(i), domain=RR) for i in range(self.order())]
        E = self.hyperedges()
        V = self.vertices()
        index = dict(zip(V, range(self.order())))
        l = var('lambda_', domain=RR)
        equs = []
        k = self.uniform()
        for v in V:
            equs.append((l-D*self.degree(v))*X[index[v]] ^ (k-1) == - add(
                [self.edge_weight(e, X)/X[index[v]] for e in self.incidence_edges(v)]))
        return equs

    def eigen_equation(self, Type="A", l=None):
        X = var_list(x, self.order())
        V = self.vertices()
        r = self.uniform()
        l = var('lambda_', domain=RR)
        Va = vector(
            [add([mul(X[u] for u in set(e)-{v}) for e in self.incidence_edges(v)]) for v in V])
        Vd = vector([self.degree(v)*X[v] ^ (r-1) for v in V])
        Vl = vector([l*X[v] ^ (r-1) for v in self.vertices()])
        EV = Vl-(Type != "A")*Vd-(-1) ^ (Type == "L")*Va
        return EV.list()

    def reduced_eigenequs(self, Type="A"):
        GP = self.automorphism_group()
        orb = GP.orbits()
        Y = var_list("y", len(orb))
        M = self.eigen_equation(Type)
        X = var_list("x", self.order())
        reduced_eigeq = [M[o[0]] for o in orb]
        Sub = dict()
        for k in range(len(orb)):
            Sub.update({X[i]: Y[k] for i in orb[k]})
        red_eigeqs = [f.subs(Sub) for f in reduced_eigeq]
        return red_eigeqs

    def eigenVfun(self, l=None, X=None):
        if X is None:
            X = [var('x' + str(i), domain=RR) for i in range(self.order())]
        if l is None:
            l = var('lambda_', domain=RR)

        d = self.uniform()
        Ve = vector(map(lambda it: it ^ (d-1), X)).column()
        E = self.hyperedges()
        V = self.vertices()
        index = dict(zip(V, range(self.order())))
        return l*Ve-vector([self.degree(v)*X[index[v]] ^ (d-1) - add(mul(X[index[u]] for u in set(e)-{v}) for e in self.incidence_edges(v)) for v in V]).column()

    def Alpha_Spectral_Radius(self, alpha=0,edge_weighted=False, vertex_weighted=False):
        self.alpha = alpha
        if "alpha" not in self.Data:
            self.Data.update({"alpha": {}})
        if alpha not in self.Data["alpha"]:
            self.Data["alpha"].update({alpha: {}})
        Alpha = self.Data["alpha"][alpha]
        if "EP" not in Alpha:
            EP = self._Alpha_Spectral_Radius(alpha,edge_weighted=edge_weighted, vertex_weighted=vertex_weighted)
            Alpha.update({"EP": EP})
        return Alpha["EP"]

    def Alpha(self):
        return self.Alpha_Spectral_Radius()[0] ^ -self.rank

    def WeightIncidenceMatrix(self, alpha=0, reduced=False, dig=7,edge_weighted=False, vertex_weighted=False):
        self.alpha = alpha
        if "alpha" not in self.Data:
            self.Data.update({"alpha": {}})
        if alpha not in self.Data["alpha"]:
            self.Data["alpha"].update({alpha: {}})
        Alpha = self.Data["alpha"][alpha]
        if "WB" not in Alpha:
            WB = self._WeightIncidenceMatrix(
                alpha, dig=dig,edge_weighted=edge_weighted, vertex_weighted=vertex_weighted)
            Alpha.update({"WB": WB})
        if reduced:
            VD2 = [v for v in self.Vertices if self.Degrees[v] > 1]
            return Alpha["WB"][VD2, :]
        return Alpha["WB"]

    def Associated_Weight_Bipartite_Graph(self, alpha=0, W=True,edge_weighted=False,vertex_weighted=False):
        self.alpha = alpha
        if "alpha" not in self.Data:
            self.Data.update({"alpha": {}})
        if alpha not in self.Data["alpha"]:
            self.Data["alpha"].update({alpha: {}})
        Alpha = self.Data["alpha"][alpha]
        if "AWBG" not in Alpha:
            AWBG = self._Associated_Weight_Bipartite_Graph(alpha, W=W,edge_weighted=edge_weighted, vertex_weighted=vertex_weighted)
            Alpha.update({"AWBG": AWBG})
        return Alpha["AWBG"]

    def edges_reindex(self):
        Edges = sorted(self.Edges)
        for i in range(self.size()):
            Edges[i].index = i

    def Relabel(self, L=None, inv=False):
        G = copy(self)
        n = self.order()
        if type(L) == list:
            oL = L+list(Set(range(n))-Set(L))
            if inv:
                nL = range(n)[::-1]
            else:
                nL = range(n)
            rebdict = dict(zip(oL, nL))
        if type(L) == dict:
            K = L.keys()
            V = L.values()
            lK = Set(K+V)-Set(K)
            lV = Set(K+V)-Set(V)
            otV = Set(range(n))-Set(K+V)
            L.update(dict(zip(lK, lV)+zip(otV, otV)))
            rebdict = L
        if L == None:
            rebdict = {v: i for i, v in enumerate(
                sorted(self.Vertices, key=str))}
        return HyperG([edge([rebdict[v] for v in e]) for e in self.Edges])

    def graphplot(self, **options):
        return HyperGraphPlot(self, options=options)

    def plot(self, **options):
        if self.size() == 0:
            return plot(self.incidence_graph())
        redraw = options.setdefault('redraw', False)
        normalize = options.setdefault('normalize', True)
        if redraw:
            self.drawn = False
            if normalize:
                Vals = self.Valleys()
                deg = Vals[0]
                self.Rotate(deg=deg)

        if not (hasattr(self, "_graphics") and self.drawn):
            # if not self.drawn:
            #     Vals = self.Valleys()
            #     deg = Vals[0]
            #     self.Rotate(deg=deg)
            # self.setLabPos()
            self._graphics = self.graphplot(**options).plot(**options)
            self.drawn = True
        return self._graphics

    def show(self, method="matplotlib", **kwds):
        plot_kwds = {k: kwds.pop(k)
                     for k in DEFAULT_HyperGraphPLOT_OPTIONS if k in kwds}
        return self.graphplot(**plot_kwds).show(**kwds)

    def _rich_repr_(self, display_manager, **kwds):
        tp = display_manager.types
        if self.Vertices == []:
            return tp.OutputPlainText("This is an Empty HyperGraph")
        prefs = display_manager.preferences
        oldprefs = display_manager.preferences.text
        display_manager.preferences.text = None
        is_small = (0 < self.order() < 60)
        can_plot = (prefs.supplemental_plot != 'never')
        plot_graph = can_plot and (
            prefs.supplemental_plot == 'always' or is_small)

        # Under certain circumstances we display the plot as graphics
        # if plot_graph:
        plot_kwds = dict(kwds)
        # plot_kwds.setdefault('title', self.name())
        plot_kwds.setdefault('fontsize', '20')
        plot_kwds.setdefault('fontweight', 'bold')
        output = self.plot()._rich_repr_(display_manager)
        # output=self.show()
        return output
        #     if output is not None:
        #         return output
        # print(plot_graph)
        # # create text for non-graphical output
        # if can_plot:
        #     text = '{0} (use the .plot() method to plot)'.format(repr(self))
        # else:
        #     text = repr(self)

        # if (prefs.text == 'latex' and tp.OutputLatex in display_manager.supported_output()):
        #     print(tp.OutputLatex in display_manager.supported_output())
        #     return tp.OutputLatex(r'\text{{{0}}}'.format(text))
        # return tp.OutputPlainText(text)

    def Latex(self):

        from sage.rings.integer import Integer
        from sage.functions.trig import arctan2
        if hasattr(self, "vLabels"):
            vLab = self.vLabels
        else:
            vLab = {v: str(v) for v in self.vertices()}

        # from warnings import warn
        # warn("\nThe hypergraph is drawn as a set of closed curves. The curve "
        #      "representing a set S go **THROUGH** the vertices contained "
        #      "in S.\n A vertex which is encircled by a curve but is not located "
        #      "on its boundary is **NOT** included in the corresponding set.\n"
        #      "\n"
        #      "The colors are picked for readability and have no other meaning.")

        latex.add_package_to_preamble_if_available("tikz")

        if not latex.has_file("tikz.sty"):
            raise RuntimeError(
                "You must have TikZ installed in order to draw a hypergraph.")

        latex.add_to_preamble(r"\usetikzlibrary{hobby}")

        pos = self._pos
        tex = r'''
\begin{tikzpicture}[scale=3,
Edge/.style={line width=1mm,opacity = .6,line cap=round,line join=round},    
HyperEdge/.style={smooth cycle,tension=1},
Node/.style={fill,circle,scale=.2}
]'''

        # colors = ["black", "red", "green", "blue",
        #           "cyan", "magenta", "yellow", "pink", "brown"]
        # colored_sets = [(s, i)
        #                 for i, S in enumerate(self.edge_coloring()) for s in S]

        for i in self.Vertices:
            tex += r'''
\coordinate (v{0}) at {1};'''.format(i, pos[i])
        E=self.Edges
        for e in E:
            s=list(e)
            if len(s) == 2:
                tex += r"""
\draw[color={0}](v{1})--(v{2});""".format(e.color,e[0], e[1])
                continue
            epos = [vector(self._pos[v]) for v in e]
            ecent = Mean(epos)
            se = sorted(e, key=lambda v: arctan2(*(vector(self._pos[v])-ecent)))
            coordinates = "".join(map(lambda k: "(v{0}) ".format(k), se))
            tex += r'''
\draw[color={0}] plot [HyperEdge] coordinates {{{1}}};'''.format(e.color, coordinates)
        for v in self.Vertices:
            label = ",label={90:$"+latex(v)+"$}" if self.labels else ""
            tex += r'''
\draw node[Node{0}] at (v{1}){{}};'''.format(label, v)
        tex += r'''
\end{tikzpicture}
'''
        return tex

    def _latex_(self):
        return self.Latex()


class HyperGraphPlot(GraphPlot):
    # class HyperGraphPlot():
    def __init__(self, hypergraph, options):
        for k, value in iteritems(DEFAULT_HyperGraphPLOT_OPTIONS):
            if k not in options:
                options[k] = value
        self._plot_components = {}
        # self.vert_lab_pos = {}
        self._nodelist = hypergraph.vertices()
        self._graph = hypergraph
        self._options = options  # contains both plot and show options
        self._options.update(self._graph.options)
        # self.set_pos()

        self._plot_components['edges'] = []
        self.set_vertices()
        self.set_edges()

    def _repr_(self):
        return "GraphPlot object for %s" % self._graph

    # def set_pos(self):
    #     if self._graph._pos == {}:
    #         self._graph._pos = self._graph.Pos()
    #     for e in self._graph.Edges:
    #         epos = [self._graph._pos[v] for v in e]
    #         e.center = mean(epos)
    #     self._pos = {k: self._graph._pos[k] for k in self._graph.vertices()}

    def set_vertices(self, **vertex_options):
        from sage.misc.superseded import deprecation

        # Handle base vertex options
        voptions = {}

        for arg in vertex_options:
            self._options[arg] = vertex_options[arg]

        # First set defaults for styles
        voptions['markersize'] = 30
        voptions['marker'] = "."
        voptions['zorder'] = 7
        voptions['facecolor'] = '#000000'
        voptions['edgecolor'] = '#000000'

        self._plot_components['vertices'] = scatter_plot(list([self._graph._pos[v] for v in self._graph.Vertices]),
                                                         clip=False, **voptions)

        self._plot_components['vertex_labels'] = []

        self._graph.setLabPos()
        if self._graph.labels:
            for v in self._graph.Vertices:
                self._plot_components['vertex_labels'].append(text("${0}$".format(self._graph.vLabels[v]),
                                                                   vector(self._graph._pos[v]) + vector(self._graph._vLabPos[v]), rgbcolor=(0, 0, 0), zorder=8, fontsize="x-large"))

        self._plot_components.update({'others': []})

    def set_edges(self, **edge_options):
        Edges = self._graph.Edges
        if Edges == []:
            pass
        for arg in edge_options:
            self._options[arg] = edge_options[arg]  
        algorithm = self._graph.options["edgeplot"]
        if 'edge_colors' in edge_options:
            self._options['color_by_label'] = False
        if self._options['edge_labels_background'] == "transparent":
            self._options['edge_labels_background'] = "None"
        
        # Handle base edge options: thickness, linestyle
        eoptions = {}
        if 'edge_style' in self._options:
            from sage.plot.misc import get_matplotlib_linestyle
            eoptions['linestyle'] = get_matplotlib_linestyle(
                self._options['edge_style'],
                return_type='long')
        if 'edge_thickness' in self._options:
            eoptions['thickness'] = self._options['edge_thickness']

        # Set labels param to add labels on the fly
        labels = False
        if self._options['edge_labels']:
            labels = True
            self._plot_components['edge_labels'] = []

        # Make dict collection of all edges (keep label and edge color)
        # Edge_Coloring = dict(flatten([[(e, i) for e in eds] for i, eds in enumerate(
        #     self._graph._edge_coloring)], max_level=1))
        colors = ["black", "red", "green", "blue",
                  "cyan", "magenta", "yellow", "pink", "brown"]
        for e in Edges:
            # e.color="black"
            e.color = colors[e.index % len(colors)]
            # e.color = colors[Edge_Coloring[e] % len(colors)]

        
        for e in self._graph.Edges:
            epos = [vector(self._graph._pos[v]) for v in e]

            ecent = Mean(epos)
            sepos = sorted(epos, key=lambda p: arctan2(*(p-ecent)))
            if algorithm == "TCB-TCBSplines":
                Pts = TCBSpline.interpolate(sepos, 0.5, 0.5, 0, True, 0.05)
            else:
                tck, u = interpolate.splprep(
                    matrix(sepos)[range(len(e))+[0], :].columns(), s=0, per=True)
                xi, yi = interpolate.splev(np.linspace(0, 1, 100), tck)
                Pts = zip(xi, yi)

            self._plot_components['edges'].append(
                polygon2d(Pts, fill=False, color=e.color))

    def show(self, **kwds):
        from sage.graphs.graph_plot import DEFAULT_SHOW_OPTIONS
        for k, value in iteritems(DEFAULT_SHOW_OPTIONS):
            if k not in kwds:
                kwds[k] = value

        self.plot().show(**kwds)

    def plot(self, **kwds):
        # self.vLabelShift()
        G = Graphics()
        options = self._options.copy()
        options.update(kwds)
        G._set_extra_kwds(Graphics._extract_kwds_for_show(options))

        # Check the arguments
        # for o in options:
        # if o not in graphplot_options and o not in G._extra_kwds:
        #     raise ValueError("Invalid input '{}={}'".format(o, options[o]))

        for comp in self._plot_components.values():
            if not isinstance(comp, list):
                G += comp
            else:
                for item in comp:
                    G += item

        if self._options['graph_border']:
            xmin = G.xmin()
            xmax = G.xmax()
            ymin = G.ymin()
            ymax = G.ymax()
            dx = (xmax - xmin) / 10.0
            dy = (ymax - ymin) / 10.0
            border = (line([(xmin - dx, ymin - dy), (xmin - dx, ymax + dy),
                            (xmax + dx, ymax + dy), (xmax + dx, ymin - dy),
                            (xmin - dx, ymin - dy)], thickness=1.3))
            border.axes_range(xmin=(xmin - dx), xmax=(xmax + dx),
                              ymin=(ymin - dy), ymax=(ymax + dy))
            G += border
        G.set_aspect_ratio(1)
        G.axes(False)
        return G


# def Alpha_Spectral_Radius(H, alpha=0):
#     if H.Vertices != list(range(H.order())):
#         raise TypeError("计算特征值需要图的顶点从零按自然顺序编号!")
#     r = H.uniform()
#     n = H.order()
#     E = H.hyperedges()
#     D = H.degrees()
#     sg = {"A": 0, "Q": 1}
#     eigval_old = 2*max(D.values())
#     X = random_vector(RR, n, 0.5, 1)
#     go_on = True
#     while go_on:
#         Y = vector([(((D[v])*X[v] ^ (r-1))*alpha + (1-alpha)*add([mul([X[i]
#                                                                        for i in set(e)-{v}]) for e in H.incidence_edges(v)])) for v in H.vertices()])
#         eigval = max([Y[i]/(X[i] ^ (r-1))for i in range(n)])
#         X1 = (Y+X.apply_map(lambda i: operator.pow(i, r-1))
#               ).apply_map(lambda i: operator.pow(i, 1/(r-1)))/add(X)
# #         if (X-X1).norm()<1e-15:
#         if abs(eigval-eigval_old) < 1e-18:
#             go_on = False
#         X = X1
#         eigval_old = eigval
#     return eigval, X


# HyperG._Alpha_Spectral_Radius = alpha_spectral_radius


def alpha_spectral_radius(G, alpha=0,edge_weighted=False,vertex_weighted=False):
    n, m, k = G.order(), G.size(), G.rank
    if edge_weighted:
        if not hasattr(G,"Edge_Weigth"):
            G.Edge_Weigth=[1]*m
    else:
        G.Edge_Weigth=[1]*m 
    edge_weight=vector(G.Edge_Weigth)   
    IM = G.incidence_matrix()
    VD = {v: IM.rows()[v].dict().keys() for v in range(n)}
    ED = {e: IM.columns()[e].dict().keys() for e in range(m)}
    #指定点权是为了利用超图的商图来计算原超图的alpha谱半径
    if vertex_weighted:
        D=G.Vertex_Weigth
    else:
        D = IM*edge_weight   
    k1, a1 = 1/(k-1), 1-alpha
    eigval_old = 2*max(D)
    #初始近似谱半径
    X = random_vector(RR, n, 0.5, 1)
    #初始近似Perron向量
    go_on = True
    while go_on:
        vE = [edge_weight[e]*mul([X[v] for v in ED[e]]) for e in range(m)]
        Y0 = [alpha*D[v]+a1*X[v] ^ (-k)*add([vE[u] for u in VD[v]])
              for v in range(n)]
        #计算A_{alpha}(H)x=lambda x
        eigval = Y0[0]
        if abs(eigval-eigval_old) < 1e-12:
            go_on = False
        X = [(X[v]*(Y0[v]+1) ^ k1) for v in range(n)]
        # 对迭代向量进行规范化, 使其1-范数为1.
        nX = add(X)
        X = [t/nX for t in X]  
        eigval_old = eigval
    return eigval, vector(X)


HyperG._Alpha_Spectral_Radius = alpha_spectral_radius
HyperG.alpha_spectral_radius = alpha_spectral_radius



def WeigthedIncidenceMatrixForQuotientGraph(G, a=0):
    k = G.uniform()
    Gr=G.incidence_graph().automorphism_group()
    EObs=Set([Set(Gr.orbit(e)) for e in G.Edges])
    EObs=list(EObs)
    connected=False
    C=cartesian_product(EObs)
    C=list(C)
    while not connected:
        c=C.pop()
        if HyperG(c).is_connected():
            connected=True
    G1=HyperG(c)
    V=G1.vertices()
    G2=G1.Relabel(V)
    n=G2.order()
    m=G2.size()
    E1=G2.Edges
    def WBe(e):
        for OE in EObs:
            if tuple([V[v] for v in e]) in OE:
                OE=flatten(OE)
                return [OE.count(V[v]) for v in range(G2.order())]
    WB=matrix([WBe(e) for e in E1]).T
    G2.Vertex_Weigth=add(WB.T)
    B=G2.incidence_matrix()
    BD = B.dict().keys()
    G2.Edge_Weigth=[RR(mul([el for el in c if el!=0])^(1/k)) for c in WB.columns()]
    EWM=diagonal_matrix(RR,[mul([el for el in c if el!=0])^(1/k) for c in WB.columns()])
    radius, X = G2.alpha_spectral_radius(alpha=a,edge_weighted=True,vertex_weighted=True)
    WE=[mul([X[v] for v in e]) for e in E1]
    PV = [X[i] ^ k for i in range(G2.order())]
    WM1 = matrix(RR,n, m, {key: round_with_eps(radius^-1*WE[key[1]]/PV[key[0]]) for key in BD})
    WIB=matrix(RDF,(1-a)*WM1*EWM+a*radius^-1*WB)
    return matrix_round(WIB),WM1.change_ring(RDF),EWM.change_ring(RDF),WB
HyperG.WIM4QG=WeigthedIncidenceMatrixForQuotientGraph




def charpoly_system(G):
    n = G.order()
    R = PolynomialRing(ZZ, 'x', n)
    B = G.incidence_matrix()
    k = G.uniform()
    E = (k-1)*identity_matrix(SR, n)
    equs = []
    for i in G.vertices():
        Ei = G[i]
        Bi = B[:, [e.index for e in Ei]]
        Bi[i] = 0
        equs.append(add([R.monomial(*c) for c in Bi.columns()]))
    return equs


HyperG.charpoly_system = charpoly_system


def comp_charpoly(M):
    dG = DiGraph(M)
    KS = dG.strongly_connected_components()
    Ks = [K for K in KS if len(K) > 1]
    ndeg = len(KS)-len(Ks)
    if Ks == []:
        return 1
    f = mul([M[K, K].charpoly() for K in Ks])
    return f*f.parent()("x") ^ ndeg


def charpoly_H(G):
    n = G.order()
    m = G.uniform()
    d = (m-1)*n - n + 1
    dlist = [m-1]*n
    R = PolynomialRing(ZZ, 'x', n)
    flist = G.charpoly_system()
    xlist = R.gens()

    mons = IntegerVectors(d, n).list()
    mons_idx = {str(mon): idx for idx, mon in enumerate(mons)}
    mons_num = len(mons)
    mons_to_keep = []
    newflist = []
    exp_coef = [[f.exponents(), f.coefficients()] for f in flist]
    numer_matrix = zero_matrix(ZZ, mons_num, sparse=True)
    def is_reduced(mon): return sum(map(lambda x: x >= m-1, mon)) == 1
    for j, mon in enumerate(mons):
        # if monomial is not reduced, then we keep it in the
        # denominator matrix:
        if not is_reduced(mon):
            mons_to_keep.append(j)
        si_mon = R._macaulay_resultant_getS(mon, dlist)
        # Monomial is in S_i under the partition, now we reduce
        # the i'th degree of the monomial
        new_mon = list(mon)
        new_mon[si_mon] -= dlist[si_mon]
        new_f = [[[g[k] + new_mon[k] for k in range(n)]
                  for g in exp_coef[si_mon][0]], exp_coef[si_mon][1]]
        for i, mon in enumerate(new_f[0]):
            k = mons_idx[str(mon)]
            numer_matrix[j, k] = new_f[1][i]
    poly_num = comp_charpoly(numer_matrix)
    poly_denom = 1
    if mons_to_keep != []:
        denom_matrix = numer_matrix.matrix_from_rows_and_columns(
            mons_to_keep, mons_to_keep)
        poly_denom = comp_charpoly(denom_matrix)
    return poly_num.quo_rem(poly_denom)[0].factor()


HyperG.charpoly = charpoly_H



def Matching_Polynomial(H, a=0):
    if isinstance(H, Graph):
        H = HyperG(H)
    E = H.Edges
    V = H.Vertices
    n = H.order()
    k = H.uniform()
    M = Combinations(E)
    Matching = [m for m in M if Set(flatten(m)).cardinality() == len(m)*k]
    D = H.degree_sequence()
    return (1-a) ^ n*add([(-1)^len(m)*mul((x-a*D[v])/(1-a) for v in Set(V)-Set(flatten(m))) for m in Matching])


HyperG.Matching_Polynomial = Matching_Polynomial
Graph.Matching_Polynomial = Matching_Polynomial




def cartesian_product_H(H, G):
    V = H.vertices()
    U = G.vertices()
    E = H.Edges
    F = G.Edges
    nE = [list(cartesian_product([[v], e])) for v in V for e in F] + \
        [list(cartesian_product([e, [v]])) for v in U for e in E]
    D = {p: i for i, p in enumerate(cartesian_product([U, V]))}
    nnE = [[D[tuple(v)] for v in Set(map(tuple, e))] for e in nE]
    return HyperG(nnE)


HyperG.cartesian_product = cartesian_product_H




def Quotient_Graph(G):
    """返回超图的边自同构群诱导划分的商图

    Args:
        G (HyperG): 超图

    Returns:
        HyperG: 带有边权和点权的商图。
    """    
    Gr=G.incidence_graph().automorphism_group()
    EObs=Set([Set(Gr.orbit(e)) for e in G.Edges])
    k=G.uniform()
    connected=False
    C=cartesian_product(EObs)
    C=list(C)
    while not connected:
        c=C.pop()
        if HyperG(c).is_connected():
            connected=True
    #c中的元素由边自同构群诱导划分中的代表元（边）构成， 即每个划分块在c中仅仅有一条边。
    #且c导出的边超图是一个连通图。
    
    EW={}
    for EO in EObs:
        for e in c:
            if e in EO:
                EW.update({e:mul([flatten(EO).count(v)^(1/k) for v in e])})
    #计算边权
    
    
    G1=HyperG(c)
    V=G1.vertices()
    VW=[G.degree(v) for v in V]
    G2=G1.Relabel(V)
    E=G2.Edges
    G2.Edge_Weigth=[RR(EW[tuple(V[v] for v in e)]) for e in E]
    for e in E:
        e.weight=RR(EW[tuple(V[v] for v in e)])
    G2.Vertex_Weigth=VW
    return G2

HyperG.Quotient_Graph=Quotient_Graph



def split_vertex(H,v,n):
    """分裂超图所有含v的超边成n个不同的边

    Args:
        H (HyperG): 超图H
        v (int): 点v
        n (int): 分裂次数

    Returns:
        HyperG: 所得超图
    """    
    Ev=H.incidence_edges(v)
    remE=[e for e in H.Edges if e not in Ev]
    for e in Ev:
        for i in range(n):
            f=copy(list(e))
            f[f.index(v)]=f"{H.Index(e)[0]}{v}_{i}"
            remE.append(f)
    G=HyperG(remE)
    return G.Relabel()

HyperG.split_vertex=split_vertex

# def WeightIncidenceMatrix(self, a=0, rational=True, D=False, reduced=True, max_e=1E-10, dig=4):
#     k = self.uniform()
#     B = self.incidence_matrix()
#     n = self.order()
#     m = self.size()
#     E = self.Edges
#     radius, X = self.Alpha_Spectral_Radius(a)

#     def f(el):
#         res = el.nearby_rational(max_error=max_e)
#         if res.denominator() > 10:
#             return SR(round(el, dig))
#         else:
#             return QQ(res)
#     WB = matrix(n, m, lambda v, e: (
#         v in E[e])*(a+(1-a)*mul([X[i] for i in E[e]])/X[v] ^ k)*radius ^ (-1))
#     B1 = WB.apply_map(f)
#     Dic = B1.dict()
#     redD = {k: Dic[k] for k in Dic if Dic[k] != 1}
#     if not D:
#         if rational:
#             if reduced:
#                 Rm = matrix(redD)
#                 Row = sorted(Set(list(matrix(redD.keys()).T[0])))
#                 return block_matrix(1, 2, [matrix(Row).T, Rm[Row, :]])
#             else:
#                 return B1
#         else:
#             return WB

#     if D:
#         if reduced:
#             return redD
#         else:
#             return Dic


# HyperG._WeightIncidenceMatrix = WeightIncidenceMatrix

def WeightIncidenceMatrix(self, alpha=0, max_e=1E-10, dig=4,edge_weighted=False,vertex_weighted=False):
    """根据超图H的alpha-张量的Perron向量构造H的赋权点边关联阵

    Args:
        a (int, optional): alpha值. Defaults to 0.
        max_e (_type_, optional):误差精度. Defaults to 1E-10.
        dig (int, optional): 保留小数位数. Defaults to 4.
        edge_weighted (bool, optional): 是否为边赋权超图. Defaults to False.

    Returns:
        _type_:返回赋权点边关联阵
    """    
    k = self.uniform()
    BD = self.incidence_matrix().dict().keys()
    n = self.order()
    m = self.size()
    E = self.Edges
    V = self.Vertices
    radius, X = self.Alpha_Spectral_Radius(alpha=alpha,edge_weighted=edge_weighted, vertex_weighted=vertex_weighted)

    def f(el):
        el=(radius^-1)*(alpha+(1-alpha)*el)
        res = el.nearby_rational(max_error=max_e)
        if res.denominator() > 10:
            return SR(round(el, dig))
        else:
            return QQ(res)
    WE = [mul([X[V.index(i)] for i in E[j]]) for j in range(m)]
    PV = [X[i] ^ k for i in range(n)]
    WM = matrix(n, m, {key: self.Edge_Weigth[key[1]]*f(WE[key[1]]/PV[key[0]]) for key in BD})
    return WM


HyperG._WeightIncidenceMatrix = WeightIncidenceMatrix




def AdjacencyTensor(G,edge_weighted=False):
    rank = G.uniform()
    order = G.order()
    A=Tensor(shape=[order]*rank)
    if G.order() < max(G.vertices())+1:
        raise TypeError("vertex number error!")
    HE = G.Edges
    if not hasattr(G,"Edge_Weigth"):
        G.Edge_Weigth=[1]*G.size()
    for i in HE:
        for el in Permutations(i):
            A[*el] = i.weight if edge_weighted else 1
    return Tensor(((1/factorial(rank-1))*A).array, symmetric=True)


HyperG.AdjacencyTensor = AdjacencyTensor



def Index(G, e):
    """返回超图G中包含点集e的所有超边的指标(编号),
    方便访问G的赋权点边关联阵中的元素

    Args:
        G (HyperG): 超图G
        e (list): 顶点列表

    Returns:
        list: 超边的指标
    Example:

    """    
    return [eg.index for eg in G.incidence_edges(e)]


HyperG.Index = Index

HyperG.WIM = WeightIncidenceMatrix


def WIe(G, v, e, dig=4,alpha=0,edge_weighted=False):
    """返回连通超图H的与alpha-张量的Perron向量相应的赋权点边关联阵中对应点v, 包含含点集e的唯一边的元素.

    Args:
        G (HyperG): k-连通超图H
        v (int): 点v的编号
        e (list): 顶点列表
        dig (int, optional): 精度到小数点后位数. Defaults to 4.
        a (int, optional): alpha值. Defaults to 0.

    Returns:
        Real number: 元素值
    """    
    return G.WIM(dig=dig,alpha=alpha,edge_weighted=edge_weighted)[v, G.Index(e)][0].list()[0]


HyperG.WIe = WIe


def Associated_Weight_Bipartite_Graph(self, alpha=None, W=True,edge_weighted=False):
    if alpha == None:
        if hasattr(self, "alpha"):
            alpha = self.alpha
        else:
            alpha = 0
            self.alpha = alpha
    if W:
        redD = self.WeightIncidenceMatrix(alpha,edge_weighted=edge_weighted).dict()
    else:
        redD = self.IncidenceMatrix().dict()
    M = matrix(redD.keys())
    V, E = map(Set, map(list, M.T))
    VB = [r"$v_%s$" % v for v in V]+[r"$e_%s$" % e for e in E]
    EB = [(r"$v_%s$" % v, "e%s" % e, redD[(v, e)]) for (v, e) in redD]
    G = Graph(EB)
    G.weighted(True)
    bG = BipartiteGraph(G)
    bG.H = self
    return bG


HyperG._Associated_Weight_Bipartite_Graph = Associated_Weight_Bipartite_Graph


def myreprhtml(G, **kwds):
    if G.Edges == []:
        return html("This is an empty Hypergraph")
    if "labels" in kwds:
        labels = kwds.pop('labels')
        G.labels = labels
    output = G.Normalize_Pos().plot()._repr_html_()
    G.drawn = True
    return output


HyperG._repr_html_ = myreprhtml


def getPath(H, u, v, t):
    u'''
    返回始边为含u,v的边，终边为含t的边的沿着u,v方向的内部路，
    '''
    e = [e for e in H.Edges if e.issuperset([u, v])][0]
    EL = []
    while (t not in e):
        EL = EL+[e]
        e = [e for e in H.incidence_edges(v) if u not in e][0]
        other_v = {i for i in e if H.degree(i) == 2}-{v}
        if len(other_v) == 1:
            u, v = v, list(other_v)[0]
    EL.append(e)
    return EL


HyperG.getPath = getPath


def list2dtopanda(TB, rows, cols):
    import pandas as pd
    DF = pd.DataFrame(TB, columns=cols, index=rows)
    styler = DF.style
    DF1 = styler.set_table_styles(
        [dict(selector="th", props=[("text-align", "right"), ('background-color',
                                                              'lightgray'), ('border-style', 'solid'), ('border-width', '1px')])]
    )
    return DF1


def get_components(H, u, v, t,edge_weighted=False):
    Path = HyperG(H.getPath(u, v, t))
    Eind = {e: e.index for e in H.Edges}
    Ev = list(Graph([[v for v in e if v in {u, t} or Path.degree(
        v) != 1] for e in Path.Edges]).depth_first_search(u))
    sE = [Ev[i:i+2] for i in range(Path.size())]
    Ei = []
    for f in sE:
        for e in Path.Edges:
            if e.issuperset(f):
                Ei.append(Eind[e])

    alpha = H.alpha
    B = H.WeightIncidenceMatrix(H.alpha,edge_weighted=edge_weighted)
    subB = B[Ev, Ei]
    # import pandas as pd
    TB = list(map(list, subB))
    # DF = pd.DataFrame(TB, columns=Ei, index=Ev)
    # styler = DF.style
    # DF1 = styler.set_table_styles(
    #     [dict(selector="th", props=[("text-align", "right"), ('background-color',
    #                                                           'lightgray'), ('border-style', 'solid'), ('border-width', '1px')])]
    # )
    return list2dtopanda(TB, Ev, Ei), Path


HyperG.incident_matrix_components = get_components

def spanning_subgraph(G,E):
    """返回超图G的边支持子图

    Args:
        G (HyperG): 超图G
        E (list): 边子集或边指标子集

    Returns:
        HyperG: 子图
    """    
    if isinstance(E[0],Integer):
        indices=sorted(E)
    if isinstance(E[0],(tuple,list)):
        indices=[]
        for e in E:
            indices=indices+G.Index(e)
        indices=sorted(indices)
    B=G.incidence_matrix()[:,indices]
    return HyperG(B)
    
    
HyperG.spanning_subgraph=spanning_subgraph  

def induced_subgraph(G,U):
    """返回超图G的点导出子图

    Args:
        G (HyperG): 超图G
        U (list): 点子集

    Returns:
        HyperG: 子图
    """
    U=Set(U)
    B=[e for e in G.Edges if Set(e).issubset(U)]
    
    return HyperG(B)
    
    
HyperG.induced_subgraph=induced_subgraph 


def TotalGraph(H):
    redD = H.incidence_matrix().dict()
    M = matrix(redD.keys())
    V, E = map(Set, map(list, M.T))
    VB = ["v%s" % v for v in V]+["e%s" % e for e in E]
    EB = [("v%s" % v, "e%s" % e, redD[(v, e)]) for (v, e) in redD]
    G = Graph(EB)
    G.weighted(True)
    return G


HyperG.TotalGraph = TotalGraph


def CyclomaticNumber(H):
    TG = H.TotalGraph()
    return TG.size()-TG.order()+TG.connected_components_number()


HyperG.CyclomaticNumber = CyclomaticNumber

def WeightedLineGraph(H):
    V=H.edges()
    n=len(V)
    A=matrix(n,lambda i,j:(i!=j)*(V[i].intersection(V[j])).cardinality())
    return Graph(A,weighted=true)
    
HyperG.WeightedLineGraph = WeightedLineGraph

def adjacency_tensor(G,edge_weighted=False):
    order = G.uniform()
    dims = G.order()
    if G.order() < max(G.vertices())+1:
        raise TypeError("vertex number error!")
    HE = G.Edges
    if not hasattr(G,"Edge_Weigth"):
        G.Edge_Weigth=[1]*G.size()
    Ts = HM(tuple([dims]*order))
    for i in HE:
        for el in Permutations(i):
            Ts[list(el)] = i.weight if edge_weighted else 1
    return (1/factorial(order-1))*Ts


HyperG.adjacency_tensor = adjacency_tensor


def laplacian_tensor(H):
    A = H.adjacency_tensor()
    D = A.degree_tensor()
    return D-A


HyperG.laplacian_tensor = laplacian_tensor


def signless_laplacian_tensor(H):
    A = H.adjacency_tensor()
    D = A.degree_tensor()
    return D+A


HyperG.signless_laplacian_tensor = signless_laplacian_tensor


def alpha_tensor(H, alpha):
    A = H.adjacency_tensor()
    D = A.degree_tensor()
    return alpha*D+(1-alpha)*A


HyperG._alpha_tensor = alpha_tensor


def append_edges(H):
    IG = H.intersection_graph()
    return [e for e in IG.vertices(sort=True) if IG.degree(e) == 1]


HyperG.append_edges = append_edges


# 添加一条首边包含点集v的长为s的悬挂路
def add_path(H, v, s, rank=0):
    order = H.order()
    size = H.size()
    r = max(H.rank, rank)
    if type(v) != list:
        v = [v]
    if s == 0:
        return H
    e = edge(v+[order, .., order+r-len(v)-1])
    e.index = size
    E = [e]
    for i in range(s-1):
        e = edge([E[-1][-1], .., E[-1][-1]+r-1])
        e.index = size+i+1
        E.append(e)
    return H.add_edges(E)


HyperG.add_path = add_path


def add_paths(H, D, rank=0,copy=True):
    G=HyperG(H)
    if not isinstance(D, dict):
        D = dict(D)
    for j in D:
        if not isinstance(D[j], list):
            D[j] = [D[j]]
        for i in D[j]:
            G = HyperG(G.add_path(j, i, rank=rank))
    return G
HyperG.add_paths = add_paths


# def blowup_edge(self, k):
#     if isinstance(self, Graph):
#         self = HyperG(self)
#     nE = []
#     E = self.Edges
#     nV = self.Vertices
#     for i in range(self.size()):
#         e = E[i]
#         ne = edge(["v({1},{0})".format(j, i) for j in range(k)])
#         nE.append(e+ne)
#     nV = sorted(Set(flatten(nE)), key=str)
#     rebdict = {v: i for i, v in enumerate(
#         sorted(nV, key=str))}
#     NE = [edge([rebdict[v] for v in e]) for e in nE]
#     nG = HyperG(NE)
#     pos = nG.incidence_graph().layout_spring(iterations=500)
#     nG._pos = {k: vector([round(x, 3), round(y, 3)])
#                for k, (x, y) in pos.items()}
#     return nG


def blowup_edge(self, k):
    if isinstance(self, Graph):
        self = HyperG(self)
    n, m = self.order(), self.size()
    pE = reshape([n..n+m*k-1], dims=[m, k])
    return HyperG([self.Edges[i]+pE[i] for i in range(m)])


HyperG.blowup_edge = blowup_edge
Graph.blowup_edge = blowup_edge

def BlowSomeEdges(H,blowedges):
    """
    Generate a new graph by deleting edges and adding new ones.

    Parameters:
    - H: the input graph
    - blowedges: a dictionary of edges to be deleted and the number of nodes to be added

    Returns:
    - A new graph with the specified edges deleted and new edges added
    """
    order=H.order()
    oldEs=list(blowedges.keys())
    newEs=[]
    for e in blowedges:
        norder=order+blowedges[e]
        newEs.append(list(e)+list(range(order,norder)))
        order=norder 
    H=H.delete_edges(oldEs)
    return H.add_edges(newEs)

HyperG.BlowSomeEdges = BlowSomeEdges
Graph.BlowSomeEdges = BlowSomeEdges

# def blowup_vertex(G, d):
#     H = HyperG(G)
#     V = H.Vertices
#     repl = {v: [(v, i) for i in range(d)] for v in V}
#     nV = flatten(list(repl.values()), max_level=1)
#     nE = list(
#         map(edge, [flatten([repl[v] for v in e], max_level=1) for e in H.Edges]))
#     L = {v: i for i, v in enumerate(nV)}
#     nE = [edge([L[v] for v in e]) for e in nE]
#     return HyperG(IncidenceStructure(nE))


# Graph.blowup_vertex = blowup_vertex
# HyperG.blowup_vertex = blowup_vertex

def blowup_vertex(G, d):
    """
    Function to blow up vertices in a hypergraph.

    Parameters:
    - G: Hypergraph to blow up for vertices.
    - d: Dictionary or integer representing the number of replicas for each vertex.

    Returns:
    - H: Hypergraph after blowing up for vertices.

    Raises:
    - TypeError: If the input `d` is not a dictionary or an integer.

    Example usage:
    >>> G = HyperG(...)
    >>> d = {1: 3, 2: 2, 3: 1}
    >>> H = blowup_vertex(G, d)
    """

    # Convert G to a HyperG object
    H = HyperG(G)

    # Get the vertices of H
    V = H.Vertices

    # Check if d is a dictionary or an integer
    if not isinstance(d, dict):
        # If d is an integer, set the same number of replicas for all vertices
        n = d
        d = {v: d for v in V}

    # Check if d is a dictionary
    if isinstance(d, dict):
        # Create a replacement dictionary for each vertex
        repl = {}
        for key in d:
            repl.update({key: ["%d%d" % (key, i) for i in range(d[key])]})

    # Create a new edge list with replicas for each vertex
    nE = [flatten([repl[v] if v in d else v for v in e]) for e in H.Edges]

    # Create a new HyperG object with the modified edge list
    H = HyperG(IncidenceStructure(nE))

    # Relabel the vertices in H
    return H.Relabel()


Graph.blowup_vertex = blowup_vertex
HyperG.blowup_vertex = blowup_vertex


def Disjoint_Union(GS, relabel=False):

    nEs = []
    p = 0
    gn = 1
    isHyperG=isinstance(GS[0],HyperG)
    for H in GS:
        E = H.hyperedges() if isHyperG else H.edges(sort=False)
        nE = list(map(edge, [{"%s:%s" % (gn, v) for v in e} for e in E]))
        gn = gn+1
        for e in nE:
            e.index = p
            e.label = "e_{%s,%s}" % (gn, p)
            p = p+1
        nEs = nEs+nE
    nH = HyperG(nEs) if isHyperG else Graph(nEs)
    nH.vLabels = {v: "v_{{{0},{1}}}".format(*re.split(':', v)) for v in nH}

    if relabel:
        return nH.Relabel()
    else:
        return nH


graphs_list.disjoint_union = Disjoint_Union


def IsIsomorphic(G, H):
    """
    Check if two hypergraphs G and H are isomorphic.

    Parameters:
    - G: Hypergraph object representing the first hypergraph.
    - H: Hypergraph object representing the second hypergraph.

    Returns:
    - True if G and H are isomorphic, False otherwise.
    """
    return Hypergraph(G.Edges).is_isomorphic(Hypergraph(H.Edges))

# Add the IsIsomorphic function as a method to the Hypergraph class
HyperG.is_isomorphic = IsIsomorphic




# def merge_vertices(H, VP, relabel=False, nv="x"):
#     u'''
#     分组合并
#     VP为图 H的点子集组: VP=[V1,V2,...,Vk]
#     该函数依次将点子集 V1,...,Vk合并为点 x1,x2,...,xk
#     '''
#     E = H.hyperedges()
#     if len(VP) == 2 and not mul([isinstance(vp, list) for vp in VP]):
#         VP = [flatten(list(VP))]
#     nE = []
#     for e in E:
#         for i in range(len(VP)):
#             oe = edge(VP[i])
#             e1 = e-oe
#             if e != e1:
#                 e = e1+edge(["%s%d" % (nv, i)])
#         nE.append(e)
#     nG = HyperG(nE)
#     if relabel:
#         return nG.Relabel()
#     else:
#         return nG


# HyperG.merge_vertices = merge_vertices


def merge_vertices(H, VGs, relabel=False, nv="w_"):
    u'''
    分组合并
    VP为图 H的点子集组: VGs=[V1,V2,...,Vk]
    该函数依次将点子集 V1,...,Vk合并为点 x1,x2,...,xk
    '''
    nVGs = []
    for VG in VGs:
        if isinstance(VG, str):
            VG = VG.split()
        nVGs.append([Integer(it) if isinstance(it, str)
                     and str.isdigit(it) else it for it in VG])
#     assert Set(H.Vertices).issuperset(Set(add(nVGs,[]))), "待合并的点中存在不少图中的顶点"

    def No(nVGs, v):
        no = -1
        i = 0
        for i in range(len(nVGs)):
            if v in nVGs[i]:
                no = i
        if no == -1:
            return v
        else:
            return nv+str(no)
    DV = {v: No(nVGs, v) for v in H.Vertices}
    nG = HyperG([[DV[v] for v in e] for e in H.Edges])

    if relabel:
        return nG.Relabel()
    else:
        return nG


HyperG.merge_vertices = merge_vertices


# def Merge_Vertices(H, VG, relabel=False, nv="w"):
#     oe = edge(VG)
#     H = HyperG(H)
#     E = H.Edges
#     cE = flatten([[e for e in H.incidence_edges(v)] for v in oe], max_level=1)
#     nE = flatten([[edge([nv])+e-edge(v if isinstance(v, list) else [v])
#                    for e in H.incidence_edges(v)] for v in VG], max_level=1)
#     nG = HyperG(Set(E)+Set(nE)-Set(cE))
#     if relabel:
#         return nG.Relabel()
#     else:
#         return nG


# HyperG.Merge_Vertices = Merge_Vertices


# def Add_innerPath(H, vp, l, k=None, relabel=False, nv="w"):
#     G = copy(H)
#     if k == None:
#         k = G.rank
#     # if l==-1:

#     if l == 0:
#         G = G.Merge_Vertices(vp)
#     else:
#         if l == 1:
#             oe = edge(flatten(vp))
#             if oe.order < k:
#                 oe = oe+edge(["%s%s" % (nv, i) for i in range(k-oe.order)])
#             nE = [oe]
#         else:
#             vsnP = flatten([vp[0]]+["%s%s" % (nv, i)
#                                     for i in range(l*(k-1)+1-len(flatten(vp)))]+[vp[1]])
#             nE = list(map(edge, [vsnP[i:i+k] for i in [0, k-1..len(vsnP)-k]]))
#         G = HyperG(nE+G.Edges)
#     if relabel:
#         return G.Relabel()
#     return G


# HyperG.Add_innerPath = Add_innerPath


def Add_innerPath(H, VLP, l, k=None, relabel=False, nv="w"):
    G = copy(H)
    assert len(VLP) == 2, "Bad Group!"
    VF, VS = VLP
    k = G.rank
    if type(VF) != list:
        VF = [VF]
    if type(VS) != list:
        VS = [VS]
    ne = flatten(VLP)
    if k == None:
        k = G.rank
    if l == 0:
        G = G.merge_vertices(ne)
    if l == 1:
        nE = [ne+[nv+"%s" % i for i in range(max(H.rank-len(ne), 0))]]
    if l > 1:
        P = HyperPath(l, k)
        E = P.Edges
        aE = P.append_edges()
        ne1, ne2 = [sorted(e, key=lambda v:P.degree(v), reverse=True)
                    for e in aE]
        s, t = [e.index for e in aE]
        E[s], E[t] = ne1, ne2
        nE = [[nv+"%s" % v for v in e] for e in E]
        e1, e2 = nE[s], nE[t]
        e1[max(len(e1)-len(VF), 1):] = VF
        e2[max(len(e2)-len(VS), 1):] = VS
    G = HyperG(nE+G.Edges)
    if relabel:
        return G.Relabel()
    return G


HyperG.Add_innerPath = Add_innerPath


# def Add_innerPaths(H, D, relabel=False, nv="w"):
#     G = copy(H)
#     i = 0
#     for item in D.items():
#         vp, l = item
#         G = G.Add_innerPath(vp, l, nv="%s%s" % (nv, i))
#         i = i+1
#     if relabel:
#         return G.Relabel()
#     return G


# HyperG.Add_innerPaths = Add_innerPaths


def Add_innerPaths(H, D, relabel=False, nv="w"):
    """对超图进行添加内部路操作

    Args:
        H (HyperG): 超图
        D (dict): 字典，用来描述如何对超图添加内部路
        键为顶点有序对，值为点对间的内部路中2度点数的列表， 比如： D={(u,v):[2,1,3]}
        relabel (重标选项), optional): 是否对新的内部路中的顶点进行重新标号 Defaults to False.
        nv (str, optional): [新点标号字符串头]. Defaults to "w".

    Returns:
        [type]: [增加内部路后的超图，原超图不变]
    """    
    G = copy(H)
    i = 0
    for item in D.items():
        vp, l = item
        if not isinstance(l, list):
            G = G.Add_innerPath(vp, l, nv="%s%s" % (nv, i))
        else:
            for key in range(len(l)):
                G = G.Add_innerPath(vp, l[key], nv="%s%s%s" % (nv, i, key))
        i = i+1
    if relabel:
        return G.Relabel()
    return G


HyperG.Add_innerPaths = Add_innerPaths


def subdivide_edge(G, e, k):
    IE = G.incidence_edges(e)
    assert len(IE) == 1, "需要确定唯一边"
    e = IE[0]
    V2 = [v for v in e if G.degree(v) > 1]
    assert len(V2) < 3, "需要所选边中至多有两个度数大于 1 的点"
    if len(V2) == 1:
        u = [v for v in e if G.degree(v) > 1]
        return G.add_path(u, k-1)
    else:
        G1 = G.Add_innerPath(V2[:2], k)
        return G1.delete_edges(IE)


HyperG.subdivide_edge = subdivide_edge


def vertice_releasing(G, i, L):
    E = copy(G.Edges)
    e = E.pop(i)
    el = list(e)
    for v in L:
        el[el.index(v)] = "w_{%d}" % (v)
    E.insert(e.index, edge(el))
    return HyperG(E)


HyperG.vertice_releasing = vertice_releasing


def LinkG(H, u, w, G, n1, n2):
    """返回在超图H的点与G的顶点间用内部路连接所得超图

    Args:
        H (HyperG): 超图H
        u (int, list): 顶点或点对
        w (int): 超图G的根点
        G (HyperG): 超图G
        n1 (int): 内部路长n1
        n2 (int): 内部路长n2

    Returns:
        HyperG: 所得超图
    """    
    G1 = copy(H)
    if isinstance(u, (list, tuple)):
        u1, u2 = u
        H1 = G1.add_paths({u1: n1, u2: n2})
    else:
        H1 = G1.add_paths({1: 3, 10: 4})
    IG = H1.intersection_graph()
    Es = [e for e in IG.vertices() if IG.degree(e) == 1]
    Vs = map(max, sorted(Es, key=lambda e: e.index, reverse=True)[:2])
    DvG = {v: (w, G) for v in Vs}
    return H1.rooted_product(DvG)


HyperG.LinkG = LinkG


# def Coalescence_Hypergraphs(G, H, VLP, relabel=True):
#     assert len(VLP) == 2, "Bad Group!"
#     nH = G+H
#     VF, VS = VLP
#     if type(VF) != list:
#         VF = [VF]
#     if type(VS) != list:
#         VS = [VS]
#     VL = nH.vertices()
#     nVL = {"x%s" % i: "w_%s" % i for i in range(len(VF))}
#     vLabels = {v: "v_{{{0},{1}}}".format(
#         *re.split('X', v[1:])) for v in nH.vertices()}

#     vLabels.update(nVL)
#     assert len(VF) == len(
#         VS),  u"The number of groups of two graphs should be same"

#     vGs = []
#     for i in range(len(VF)):
#         if not isinstance(VF[i], list):
#             VF[i] = [VF[i]]
#         nt1 = ["v%sX%s" % (1, t) for t in VF[i]]
#         if not isinstance(VS[i], list):
#             VS[i] = [VS[i]]
#         nt2 = ["v%sX%s" % (2, t) for t in VS[i]]
#         vGs.append(nt1+nt2)

#     for i in range(len(vGs)):
#         vG = vGs[i]
#         nH = nH.Merge_Vertices(vG, relabel=False, nv="x%s" % (i))
#     nVS = Set(nH.vertices())
#     nVx = sorted(nVS-Set(VL))
#     nVD = {v: vLabels[v] for v in nH.vertices() if v not in nVx}
#     nVD2 = {nVx[i]: "w_%s" % i for i in range(len(nVx))}
#     nVD.update(nVD2)
#     if relabel:
#         nH = nH.Relabel(nVD)
#     return nH


# HyperG.coalescence_graphs = Coalescence_Hypergraphs


def Coalescence_Hypergraphs(G, H, VLP, relabel=False):
    """
        VF=[V1,V2,...,Vk], V1,...,Vk为 G 的点子集族
        VS=[U1,U2,...,Uk], U1,...,Uk 为 H 的点子集族
        VLP=[VF,VS]
        该函数依次将点子集 V1+U1,...,Vk+Uk合并为点 w1,w2,...,wk

    Args:
        G (_type_): _description_
        H (_type_): _description_
        VLP (_type_): _description_
        relabel (bool, optional): 是否对顶点重新编号. Defaults to False.

    Returns:
        HyperG: 返回所得超图
    """    

    assert len(VLP) == 2, "Bad Group!"
    VF, VS = VLP
    assert len(VF) == len(
        VS),  u"The number of groups of two graphs should be same"
    nVF, nVS = [], []

    for V in VF:
        if not isinstance(V, list):
            V = [V]
        nVF.append(["1:%s" % v for v in V])
    for V in VS:
        if not isinstance(V, list):
            V = [V]
        nVS.append(["2:%s" % v for v in V])
    VP = [flatten(t) for t in zip(nVF, nVS)]
    nH = G+H
    nH = nH.merge_vertices(VP, nv="w_")
    if relabel:
        nH = nH.Relabel()
    return nH


HyperG.coalescence_graphs = Coalescence_Hypergraphs


def weighted_incidence_graph(H):
    G = H.incidence_graph()
    G.weighted(True)
    Eg = G.edges(labels=False)
    edge_weighted = {e: H.WIe(e[1], e[0]) for e in Eg}
    G.set_edge_weight(edge_weighted)
    G1 = HyperG(G)
    G1 = G1.Normalize_Pos()
    pos = G1.Pos()
    G._pos = {e: pos[e] for e in G.Pos()}
    return G.graphplot(figsize=[5, 3], vertex_size=50, edge_labels=True, edge_labels_background="transparent", vertex_labels=False, vertex_colors={"red": H.Edges, "green": H.Vertices})


HyperG.weighted_incidence_graph = weighted_incidence_graph


def Rooted_Product(self, DvH, relabel=True):
    """
    超图的根积运算
    H=HyperCycle(3,3)
    G=VStarLike(2,2,2)
    G.rooted_product({3:(0,H),7:(0,H)})
    """
    rtV = list(DvH.keys())
    rtGs = DvH.values()
    HL = [t[1] for t in rtGs]
    rtHs = [t[0] for t in rtGs]
    vp = [["1:%s" % rtV[i], "%s:%s" % (i+2, rtHs[i])] for i in range(len(HL))]
    G1 = graphs_list.disjoint_union([self]+HL)
    G2 = G1.merge_vertices(vp)
    if relabel:
        G2 = G2.Relabel()
    return G2


HyperG.rooted_product = Rooted_Product


def Blow_Up(G, VS, relabel=False):
    """
    Performs a blow-up operation on a hypergraph.

    Parameters:
    - G: HyperG object
        The input hypergraph on which the blow-up operation is performed.
    - VS: list or dictionary
        The vertex set of the hypergraph `G` specifying the number of copies for each vertex.
    - relabel: bool, optional (default=False)
        Specifies whether to relabel the vertices of the resulting blown-up hypergraph.

    Returns:
    - H: HyperG object
        The blown-up hypergraph obtained after performing the blow-up operation.

    Notes:
    - The blow-up operation replaces each vertex `v` in the input hypergraph `G` with a set of vertices,
      where the number of copies for each vertex is specified by the `VS` parameter.
    - The resulting blown-up hypergraph `H` is a hypergraph.

    Example:
    >>> G = HyperG([(0, 1), (1, 2)])
    >>> VS = {1:2, 2:3]
    >>> H = G.Blow_Up(VS)
    >>> H.Vertices
        ['2_2', '0', '2_1', '2_0', '1_0', '1_1']
    >>> H.Edges
        [('2_2', '1_0'), ('2_2', '1_1'), ('0', '1_0'), ('0', '1_1'), ('2_1', '1_0'), ('2_1', '1_1'), ('2_0', '1_0'), ('2_0', '1_1')]
    """
    E = G.Edges
    nE = []
    for e in E:
        nE = nE + cartesian_product([[f"{v}_{u}" for u in range(VS[v])] if v in VS else [f"{v}"] for v in e]).list()
    H = HyperG(nE)
    return H.Relabel() if relabel else H

HyperG.Blow_Up = Blow_Up



Mean=lambda L:add(L)/len(L)

def Rotate_Shift(H, Vp, degree=0, orgin=0, ndigits=3):
    """对图或超图H的图形表示进行整体围绕点 orgin进行旋转, 使得点对或点集对的重心连线与水平方向成夹角degee度

    Args:
        H ([Graph, HyperG]): 待处理的图或超图
        Vp ([type]): [description]
        degree (int, optional): 水平角度. Defaults to 0.
        orgin (int, optional): [旋转中心的顶点]. Defaults to 0.
        ndigits (int, optional): [坐标计算精度]. Defaults to 3.

    Returns:
        [Graph, HyperG]: 修正顶点坐标的图或超图.
    """
    
    a = degree*pi/180
    R = matrix(RR, 2, 2, [cos(a), sin(a), -sin(a), cos(a)])
    V = H.vertices()
    if isinstance(H, Graph):
        pos = H.Pos()
    else:
        pos = H._pos
    vc = operator.sub(*[Mean([vector(pos[v])
                              for v in flatten([P])]) for P in Vp])
    ct = Mean([pos[i] for i in V])
    nvc = vc.normalized()
    T = nvc[0]*identity_matrix(2)+nvc[1]*matrix(2, [0, 1, -1, 0])
    npos = {i: R*T*(pos[i]-ct)+ct for i in pos}
    opos = npos[orgin]
    H._pos = {i: vector(map(lambda p: round(p, ndigits),
                            npos[i]-opos)) for i in npos}
    if isinstance(H, Graph):
        return H.plot()
    else:
        return H.plot(redraw=True)


HyperG.Rotate_Shift = Rotate_Shift


def Valleys(nH):
    import numpy as np
    import peakutils
    if not nH._pos:
        nH.Pos()
    PM = matrix([nH._pos[v] for v in nH.vertices()])
    def gap(X): return max(X)-min(X)
    Area = +oo
    D = 0
    LAD = []
    xData = [0, 0.5, .., 180]
    for deg in xData:
        c = exp(CC.0*deg*pi/180)
        x, y = c.real(), c.imag()
        R = matrix(RR, 2, 2, [x, y, -y, x])
        dx, dy = list(map(gap, (PM*R).columns()))
        area = dx*dy ^ 2
        LAD.append([deg, area])
    MT = matrix(LAD)
    cb = np.array(-MT.T[1])
    indexes = list(peakutils.indexes(cb, thres=0.01/max(cb), min_dist=50))
    return matrix(sorted(MT[indexes], key=lambda x: x[1], reverse=False)).T[0]


HyperG.Valleys = Valleys
Graph.Valleys = Valleys


def Rotate(nH, deg, digits=3, scale=1):
    PM = matrix([nH._pos[v] for v in nH.vertices()])
    c = exp(CC.0*deg*pi/180)
    x, y = c.real(), c.imag()

    R = scale*matrix(RR, 2, 2, [x, y, -y, x])
    Pos = {v: nH._pos[v]*R for v in nH._pos}
    ct = Mean([Pos[i] for i in nH.vertices()])
    nH._pos = {
        v: vector(SR, map(lambda x: round(x, digits) if abs(x-round(x, 0)) > 0.1 ^ (deg) else int(round(x, 0)), Pos[v]-ct)) for v in Pos}
    nH.drawn = False
    return nH


HyperG.Rotate = Rotate
Graph.Rotate = Rotate


def Normalize_Pos(nH, digits=3, rank=0, scale=1):
    Vals = nH.Valleys()
    deg = Vals[rank % len(Vals)]
    H = copy(nH)
    H.Rotate(deg=deg, digits=digits, scale=scale)
    return H


HyperG.Normalize_Pos = Normalize_Pos
Graph.Normalize_Pos = Normalize_Pos


def Graph2TeX(G):
    tempfilename = tmp_filename()
    tex_file = tempfilename + ".tex"
    tex = r'''
% !TEX program = xelatex
\documentclass{ctexart}
\usepackage{subcaption}
\usepackage{tikz}
\usetikzlibrary{hobby}
\begin{document}

\begin{figure}
\centering
'''
    tex += G.Latex()
    tex += r'''
\caption{%s}
\end{figure}

\end{document}
''' % G.name()
    f = open(tex_file, "w")
    f.write(tex)
    f.close()
    os.system('open "%s"' % (tex_file))


HyperG.Graph2TeX = Graph2TeX


# def HyperPath(e, k=2):
#     return HyperG(graphs.PathGraph(e+1)).blowup_edge(k-2)


def HyperPath(l, k=2):
    return HyperG(block_matrix(1, 2, [matrix([[0..l-1], [1..l]]).T,
                                      matrix(reshape([l+1..l*(k-1)], dims=[l, k-2]))]).rows())


def HyperCycle(e, k=2):
    """
    根据给定的参数生成一个超图。

    参数:
    - e: int，环的大小。
    - k: int，可选参数，默认值为2。边的膨胀数。

    返回值:
    - HyperG，生成的超图对象。
    """
    if e==2 and k>2:
        return HyperG([[0..k-1],[0,k-1]+[k..2*k-3]])
    return HyperG(graphs.CycleGraph(e)).blowup_edge(k-2)


def HyperStar(e, k=2):
    return HyperG(graphs.StarGraph(e)).blowup_edge(k-2)


def Bicyclic_C1(p, k=3, relabel=True):
    G = HyperG(2*graphs.PathGraph(2))
    G = G.blowup_edge(k-2)
    ME = matrix(G.hyperedges())
    vPs = ME[:, :2].rows()+ME[:, 2].T.rows()
    D = dict(zip(vPs, p))
    G = G.Add_innerPaths(D)
    if relabel:
        G = G.Relabel()
    G.Normalize_Pos()
    G.name(LatexExpr(r"$C_1(%s,%s;%s)$" % tuple(p)))
    return G


def Bicyclic_C2(p, k=3, relabel=True):
    p = tuple(p)
    G = HyperG(2*graphs.PathGraph(2))
    G = G.blowup_edge(k-2)
    D = dict(zip(matrix(G.hyperedges()).columns(), p))
    G = G.Add_innerPaths(D)
    if relabel:
        G = G.Relabel()
    G.name(LatexExpr(r"$\Theta_2(%s,%s,%s)$" % p))
    return G


def Bicyclic_C3(p, k=4, relabel=True):
    G = HyperG(graphs.PathGraph(2))
    G = G.blowup_edge(k-2)
    D = dict(zip([(0, 1), (2, 3)], p))
    G = G.Add_innerPaths(D)
    if relabel:
        G = G.Relabel()
    G.name(LatexExpr(r"C_3(%s,%s)" % tuple(p)))
    return G


# 生成 C1 型超图

def Bicyclic_C1s(n, k=3):
    GL = []
    P = []
    tpart = [list(map(lambda i:i-1, c))
             for c in Compositions(n+2, inner=[3, 2], length=2)]
    for t in tpart:
        for p in Partitions(t[0], length=2, min_part=1):
            P.append(p+[t[1]])
    for p in P:
        H = Bicyclic_C1(p, k)
        GL.append(H)
    return GL


# 生成 C2 型超图


def Bicyclic_C2s(n, k):
    GL = []
    P = list(Partitions(n, length=3))
    for p in P:
        H = Bicyclic_C2(p, k)
        GL.append(H)
    return GL


def VStarLike(*args, rank=3):
    """
    设$H$是由顶点$v$上分别长出长为$n_1,n_2,\cdots,n_s$的悬挂路所得的$k$一致超图, 称$H$点心像星超树, 
    记做$VStarlike^{[k]}(n_1,n_2,\cdots, n_s)$

    Args:
        rank (int, optional): k一致. Defaults to 3.

    Returns:
        HyperG: 点心像星超树
    """    
    G = HyperG(V=[0])
    return G.add_paths({0: list(args)}, rank=rank)


def EStarLike(*args,rank=None):
    """
    设$H$是由超边$e$的$s$个不同的 1 度点分别长出长为$n_1,n_2,\cdots,n_s$的悬挂路所得的$k$一致超图, 称$H$边心像星超树, 
    记做$EStarlike^{[k]}(n_1,n_2,\cdots, n_s)$

    Args:
        rank (int, optional): k一致. 默认值:长出的悬挂边数.

    Returns:
        HyperG: 边心像星超树
    """    
    m = len(args)
    if rank==None:
        rank=m
    if rank<m:
        rank=m
    G = HyperPath(1, rank)
    D = {i: args[i] for i in range(m)}
    return G.add_paths(D)


def HyperTree1(l):
    G = HyperG(4*graphs.PathGraph(2), k=4)
    G = G.blowup_edge(2)
    H = G.Add_innerPath([(0, 2), (4, 6)], l)
    return H.Relabel()


# def H2222(l,k=3):
#     H=graphs_list.disjoint_union([HyperPath(2,k)]*4)
#     return H.Add_innerPath((["1:0","2:0"],["3:0","4:0"]),l,relabel=True)


def Gee(*args, k=3):
    """   
     由$\mathbb{P}_m$分别在悬挂点$u_1,u_2,...,u_s;v_1,v_2,...,v_t$
    依次长为$n_1,n_2,...,n_s;n'_1,n'_1,...,n'_t$ 个悬挂路所得图, 其最大度为2.

    Args:
        k (int, optional):k一致. Defaults to 3.

    Returns:
        HyperG: 超图
    Examples:
        #在长为9的loose超路的各在两个悬挂边的两个悬挂点上长悬挂边所得4一致超图
        >>> Gee(9,[1,1],[1,1],k=4)
    """    
    m = args[0]
    if m==1 and len(args[1]+args[2])>k:
        k=len(args[1]+args[2])
    P = HyperPath(m, k)
    E = P.Edges
    e1, e2 = E[0], E[-1]
    e1d1s = [e1[0]]+e1[2:]
    if m==1:
        e1d1s=e1
        assert len(e1d1s) >= len(args[1]+args[2]), "粘接点数不足, 请增大边秩k"
    e2d1s = [v for v in e2 if P.degree(v) == 1]

    if isinstance(args[1], list):
        if m==1:
            e2d1s=e2d1s[::-1]
            
        assert len(e1d1s) >= len(args[1]) and len(e2d1s) >= len(args[2]), "粘接点数不足, 请增大边秩k"
        D = dict(list(zip(e1d1s, args[1]))+list(zip(e2d1s, args[2])))
    else:
        n1, n2, n3, n4 = args[1:]
        u, v = e1[0], e1[2]
        s, t = e2[-2:]
        D = {u: n1, v: n2, s: n3, t: n4}
    return P.add_paths(D)



def Gev(*args, k=3):
    r"""
    由$\mathbb{P}_m$分别在悬挂点$u_1,u_2,...,u_s;v_1,v_2,...,v_t$依次长为$n_1,n_2,...,n_s;n'_1,n'_1,...,n'_t$ 个悬挂路所得图, 其最大度为2.
    """
    m = args[0]
    P = HyperPath(m, k)
    E = P.Edges
    e1, e2 = E[0], E[-1]
    e1d1s = [e1[0]]+e1[2:]
    if m==1:
        e1d1s=e1
        assert len(e1d1s)-1 >= len(args[1]), "粘接点数不足, 请增大边秩k"
    e2d1s = [v for v in e2 if P.degree(v) == 1]
    if isinstance(args[1], list):
        assert len(e1d1s) >= len(args[1]), "粘接点数不足, 请增大边秩k"
        D1 = dict(zip(e1d1s[:len(args[1])], args[1]))
        if m>1:
            D2 = dict(zip(e2d1s[:1], [flatten(args[2])]))
        else:
            D2 = dict(zip(e2d1s[::-1][:1], [flatten(args[2])]))
        D = dict(list(D1.items())+list(D2.items()))
    else:
        n1, n2, n3, n4 = args[1:]
        u, v = e1[0], e1[2]
        s, t = e2[-2:]
        D = {u: n1, v: n2, s: [n3, n4]}
    return P.add_paths(D)


def Gvv(m, VL, VR, k=3):
    r"""
    由$\mathbb{P}_m$分别在悬挂点$u_1;v_1$依次长为$VL;VR$ 个悬挂路所得图, 其最大度为2.
    """
    P = HyperPath(m, k)
    E = P.Edges
    e1, e2 = E[0], E[-1]
    u, v = e1[0], e1[2]
    s, t = e2[-2:]
    return P.add_paths({u: VL, s: VR})

# def Gvv(m, n1, n2, n3, n4,k=3):
#     r"""
#     由$\mathbb{P}_m$分别在悬挂点$u_1;v_1$依次长为$n_1,n_2;n_3,n_4$ 个悬挂路所得图, 其最大度为2.
#     """
#     P = HyperPath(m, k)
#     E = P.Edges
#     e1, e2 = E[0], E[-1]
#     u, v = e1[0], e1[2]
#     s, t = e2[-2:]
#     return P.add_paths({u: [n1,n2], s: [n3, n4]})


def GWG(m, De={}, Dv={}, k=3):
    """由$\mathbb{P}_m$分别在它的各点长若干悬挂路所得图G.

    Args:
        m (intege): 主干路的长
        De (dict, optional): 边上长悬挂路的分布字典. Defaults to {}.
        Dv (dict, optional): 非悬挂点上长悬挂路的分布字典. Defaults to {}.
        k (int, optional): k一致 Defaults to 3.

    Returns:
        HyperG: G
    Examples:
        >>> G = GWG(5,{0:[1,1],2:[3]},{3:[2,2]})
    在$P_5$的第0条边的两个悬挂点上长悬挂边, 在第2条边的1个悬挂点上长为3悬挂路, 在主干路第3个非悬挂点上长两个长为2悬挂路所得图G.
    """    

    for key in De:
        if not isinstance(De[key], list):
            De[key] = [De[key]]
        le = len(De[key])+(1 if key in {0, -1, m-1} else 2)
        if le > k:
            k = le
    if m == 1:
        k = k-1
    P = HyperPath(m, k)
    E = P.Edges
    V1ds = {e: [v for v in e if P.degree(v) == 1] for e in E}
    V2ds = {e: [v for v in e if P.degree(v) > 1] for e in E}
    D = {}
    for key in De:
        D.update(dict(zip(V1ds[E[key]][:len(De[key])], De[key])))
    # D.update(dict({V2ds[E[key]][0]: flatten(Dv[key]) for key in Dv})
    D.update(dict({key: flatten(Dv[key]) for key in Dv}))
    return P.add_paths(D)


def Ca(m, a=1, k=3):
    """在长为m的超圈的1度点上接长为a的悬挂路所得超图

    Args:
        m (int): 圈长
        a (int, optional): 悬挂路长. Defaults to 1.
        k (int, optional): k一致. Defaults to 3.

    Returns:
        HyperG: 所得超图
    """    
    H = HyperCycle(m, k=k)
    n = H.order()
    H = H.add_paths({n-2: a})
    H.name("$Ca(%s)$" % m)
    return H


def Cva(m, a=1, k=3):
    """在长为m的超圈的2度点上接长为a的悬挂路所得超图
    Args:
        m (int): 圈长
        a (int, optional): 悬挂路长. Defaults to 1.
        k (int, optional): k一致. Defaults to 3.
    Returns:
        HyperG: 所得超图
    """


    H = HyperCycle(m, k=k)
    n = H.order()
    H = H.add_paths({0: a})
    H.name("$Cva(%s)$" % m)
    return H


def U1(n, i):
    H = 2*HyperPath(2, 3)
    H = H.Relabel()
    H = H.Add_innerPaths({(0, 5): ceil(i), (3, 8): floor(n-i-4)})
    return H.Relabel()


def U2(n):
    G = HyperPath(3, 3)
    G = G.Add_innerPaths({(0, 4): n-3})
    return G.Relabel()


def CompleteHyperG(n, r):
    """
    This function calculates the Complete r Uniform Hypergraph of order n.
    
    Parameters:
        - n (int): The order of the hypergraph.
        - r (int): The uniform of the hypergraph.
    
    Returns:
        - HyperG: The Complete r Uniform Hypergraph of order n.
    """
    return HyperG(hypergraphs.CompleteUniform(n, r))