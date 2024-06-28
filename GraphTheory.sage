# -*- coding: utf-8 -*-

#########################################################################################
#      Here is  an enhancement to the part of  graph theory within sage                 #
#                                                                                       #
#       Copyright (C) 2014~2023 HaiyingShan <shan_haiying@tongji.edu.cn>                #
#                             Last Modified: 2024-10-13                                  #
#########################################################################################

from sage.matrix.constructor import matrix
import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import scipy.integrate
import sage.graphs.graph_generators
import sage.graphs.graph_plot
from scipy.linalg import eigh as largest_eigh
# from scipy.sparse.linalg.eigen.arpack import eigsh as largest_eigsh
from PIL import Image
import io
import base64


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



def Index(Edges):
    return [e.index for e in Edges]

# def attach_Edges_Class(G, label=False):
#     if type(G) == Graph and not label:
#         E = G.edges(labels=label)
#     else:
#         E = G.edges()
#     G.Edges = []
#     k = 2 if type(G) == Graph else G.uniform()
#     for i in range(G.size()):
#         e = edge(E[i])
#         e.G = G
#         e.order = k
#         e.vertice = sorted(e[:k])
#         e.index = i
#         G.Edges.append(e)


# Graph.attach_Edges_Class = attach_Edges_Class
# Hypergraph.attach_Edges_Class = attach_Edges_Class


# def Edge2Dict(G):
#     if hasattr(G, "Edges"):
#         G.attach_Edges_Class()
#     G.D = {tuple(e.vertice): e for e in G.Edges}


# Hypergraph.Edge2Dict = Edge2Dict
# Graph.Edge2Dict = Edge2Dict


# def RichRepr(self, display_manager, **kwds):
#        tp = display_manager.types
#         if self.Vertices == []:
#             return tp.OutputPlainText("This is an Empty HyperGraph")
#         prefs = display_manager.preferences
#         oldprefs = display_manager.preferences.text
#         display_manager.preferences.text = None
#         is_small = (0 < self.order() < 60)
#         can_plot = (prefs.supplemental_plot != 'never')
#         plot_graph = can_plot and (
#             prefs.supplemental_plot == 'always' or is_small)

#         # Under certain circumstances we display the plot as graphics
#         # if plot_graph:
#         plot_kwds = dict(kwds)
#         # plot_kwds.setdefault('title', self.name())
#         plot_kwds.setdefault('fontsize', '20')
#         plot_kwds.setdefault('fontweight', 'bold')
#         output = self.plot()._rich_repr_(display_manager)
#         # output=self.show()
#         return output

# Graph._rich_repr_=_rich_repr_
# =======================图的矩阵表示与谱:开始==================================


def degree_matrix(self):
    A = self.am()
    D = diagonal_matrix((A*ones_matrix(self.order(), 1)).list())
    return D


Graph.degree_matrix = degree_matrix


def signless_laplacian_matrix(self):
    return self.degree_matrix()+self.adjacency_matrix()


Graph.signless_laplacian_matrix = signless_laplacian_matrix


def is_transmission_regular(G):
    return len(Set(add(G.distance_matrix()).list())) == 1


Graph.is_transmission_regular = is_transmission_regular


def IsSeparable(G):
    f=G.charpoly()
    return (f.gcd(f.diff(x))).degree()==0
Graph.IsSeparable=IsSeparable


def alpha_matrix(G, alpha):
    A = G.adjacency_matrix()
    n = G.order()
    D = diagonal_matrix(A*vector([1]*n))
    return alpha*D+(1-alpha)*A


Graph.alpha_matrix = alpha_matrix
DiGraph.alpha_matrix = alpha_matrix

def Distance_alpha_matrix(G, alpha):
    D=G.distance_matrix()
    T= diagonal_matrix( add(G.distance_matrix()))
    return alpha*T+(1-alpha)*D


Graph.Distance_alpha_matrix = Distance_alpha_matrix



def ABC_matrix(G):
    return matrix(G.order(), lambda i, j: G.has_edge(i, j)*sqrt((G.degree(i)+G.degree(j)-2)/(G.degree(i)*G.degree(j))))


Graph.ABC_matrix = ABC_matrix


def H_Spectra_Radius(Q):
    EI = Q.eigenvalues()
    return max([e for e in EI if e in RR])

# 幂敛法计算非负矩阵的谱半径(对非对称阵，算法不知是否一定收敛，要查文献)


def spectral_radius(A, niter=2000):
    tol = 1e-18
    n = A.dimensions()[0]
    eigvec = random_vector(RR, n, 0.5, 1).normalized(2)
    dc = A.dict()
    eigval_old = add([mul([eigvec[t] for t in key])*dc[key]
                      for key in dc.keys()])
    i = 0
    while 1:
        i = i+1
        vlst = []
        for v in range(n):
            row = A[v].dict()
            vlst.append(add([eigvec[t]*row[t] for t in row.keys()]))
        eigvec1 = vector(vlst).normalized(2)
        eigval_new = add([mul([eigvec1[t] for t in key])*dc[key]
                          for key in dc.keys()])
        if (abs(eigval_new-eigval_old)/eigval_new) < tol or i > niter:
            return eigval_new, eigvec1
        eigvec = eigvec1
        eigval_old = eigval_new


def Spectral_Radius_digraph(G, a, niter=2000):
    Q = G.alpha_matrix(a)
    return spectral_radius(Q)


DiGraph.Spectral_Radius = Spectral_Radius_digraph


def RandomCactusGraph(lc):
    G = DiGraph(1)
    for c in lc:
        v = randint(0, G.order()-1)
        G = G.add_innerPaths({v: c-1})
    return G


sage.graphs.digraph_generators.digraphs.RandomCactusGraph = RandomCactusGraph


def normalized_laplace_matrix(self):
    M = self.laplacian_matrix(normalized=True)
    return M


Graph.normalized_laplace_matrix = normalized_laplace_matrix

# Definition of Distance Laplacian


def distLap(G):
    n = G.order()
    ones = ones_matrix(n, 1)
    dg = G.distance_matrix()
    trg = dg*ones
    DL = diagonal_matrix(vector(trg))-dg
    return DL


def eigenpair(M, spe=False, order=0):
    M = matrix(CC, M)
    if not M.is_hermitian():
        return "A hermitian matrix is required."
    w, v = np.linalg.eigh(M)
    if spe:
        return list(map(lambda n: round(n) if abs(n-round(n)) < 0.0000000001 else n, sorted(list(vector(w).zero_at(1e-12)), reverse=True)))
    v = matrix(v).zero_at(1e-12)
    w = diagonal_matrix(w).dense_matrix().zero_at(1e-12).diagonal()
    enulst = dict([(i, j) for i, j in enumerate(w)])
    ret = sorted(enulst.items(), key=lambda x: x[1], reverse=True)
    Od = [i[0] for i in ret]
    if order == "all":
        return diagonal_matrix(w), matrix(v)
    else:
        rec = v.columns()[Od[order]]
        if rec[0] != 0:
            rec = sgn(rec[0])*rec
        return w[Od[order]], rec

# def eigenpair(M, spe=False, order=0):
#     M = matrix(RR, M)
#     if not M.is_symmetric():
#         return "A symmetric matrix is required."
#     w, v = np.linalg.eigh(M)
#     w = w[::-1]
#     v = v[::-1]
#     if spe:
#         return w
#     if order == "all":
#         return diagonal_matrix(w), matrix(v)
#     else:
#         return w[order], vector(v[order])


# def hermitian_adjacency_matrix(H):
#     G = H.to_undirected()
#     A1 = G.adjacency_matrix()
#     A = H.adjacency_matrix()
#     B = A1-A
#     return A+(I-1)*B.T-B*I


def hermitian_adjacency_matrix(H):
    n = H.order()
    sparse = True
    if n <= 256 or H.density() > 0.05:
        sparse = False
    I = CC(sqrt(-1))
    D = {}
    for i, j, k in H.edge_iterator():
        D[i, j] = I
        if (j, i) in D:
            D[i, j] = 1
            D[j, i] = 1
        else:
            D[j, i] = -I
    return matrix(CC, n, n, D, sparse=sparse)


DiGraph.hermitian_adjacency_matrix = hermitian_adjacency_matrix


def charpolys_for_graph(self, type="L"):
    r"""[summary]

    Args:
        type (str, optional): [description]. Defaults to "L".

    Returns:
        [type]: [description]
    """    
    r"""
        Return the characteristic polynomial of the adjacency matrix of the
        (di)graph.

        Let `G` be a (simple) graph with adjacency matrix `A`. Let `I` be the
        identity matrix of dimensions the same as `A`. The characteristic
        polynomial of `G` is defined as the determinant `\det(xI - A)`.

        .. NOTE::

            ``characteristic_polynomial`` and ``charpoly`` are aliases and thus
            provide exactly the same method.

        INPUT:

        - ``x`` -- (default: ``'x'``); the variable of the characteristic
          polynomial

        - ``laplacian`` -- boolean (default: ``False``); if ``True``, use the
          Laplacian matrix

        .. SEEALSO::

            - :meth:`kirchhoff_matrix`

            - :meth:`laplacian_matrix`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.characteristic_polynomial()
            x^10 - 15*x^8 + 75*x^6 - 24*x^5 - 165*x^4 + 120*x^3 + 120*x^2 - 160*x + 48
            sage: P.charpoly()
            x^10 - 15*x^8 + 75*x^6 - 24*x^5 - 165*x^4 + 120*x^3 + 120*x^2 - 160*x + 48
            sage: P.characteristic_polynomial(laplacian=True)
            x^10 - 30*x^9 + 390*x^8 - 2880*x^7 + 13305*x^6 -
            39882*x^5 + 77640*x^4 - 94800*x^3 + 66000*x^2 - 20000*x
        """
    M = self.get_Matrix(type)
    return M.charpoly()


Graph.charpolys_for_graph = charpolys_for_graph
Graph.charpolys = charpolys_for_graph
DiGraph.charpolys = charpolys_for_graph


def Weighted_Adjancency_Matrix(H,f):
    n=H.order()
    return matrix(n,lambda i,j:(i in H[j])*f(H.degree(i),H.degree(j)))
Graph.Weighted_Adjancency_Matrix=Weighted_Adjancency_Matrix





def eigenpairs(self, type="L", order=0, spe=False):
    """
    计算图的各类特征值，谱及特征向量

    例:计算Petersen图的邻接谱
    ::
        sage: G=graphs.PetersenGraph()
        sage: G.eigenpairs(type="A",spe=True)
        [3, 1, 1, 1, 1, 1, -2, -2, -2, -2]
    """
    M = self.get_Matrix(type)
    return eigenpair(M, spe=spe, order=order)


Graph.eigenpairs = eigenpairs
DiGraph.eigenpairs = eigenpairs


def Spectral_Radius(self, type="L", order=0, spe=False):
    """
    计算图的各类特征值，谱及特征向量

    例:计算Petersen图的邻接谱
    ::
        sage: G=graphs.PetersenGraph()
        sage: G.eigenpairs(type="A",spe=True)
        [3, 1, 1, 1, 1, 1, -2, -2, -2, -2]
    """
    M = self.get_Matrix(type)
    return largest_eigsh(matrix(RR, M).numpy(), 1)


Graph.Spectral_Radius = Spectral_Radius


def plot_by_eigenvector(G, cols, **options):
    '''
    根据图的特征向量来画图

    INPUT:

    - ``**options`` -- 输入特征向量的类型，如A,L,Q, W 等

    - ``cols`` -- 用第几个向量来画图， 如（-1,-2）等

    例:

    以laplacian矩阵的第二，三小特征向量画图::

        sage: G.plot_by_eigenvector((-2,-3),type="L")
        ...

    '''
    T = options.pop("type", "L")
    D, V = G.eigenpairs(type=T, order="all")
    pos = dict(zip(G.vertices(), zip(*[V.T[i] for i in cols])))
    if len(cols) == 2:
        G.set_pos(pos)
        G.relabel()
        return G.plot(**options)
    else:
        return G.plot3d(pos3d=pos, **options)


Graph.plot_by_eigenvector = plot_by_eigenvector


# Spectral decomposition
def Spectral_decomposition(G):
    esp = G.am().eigenspaces_left()
    m = len(esp)
    evs = [sp[0] for sp in esp]
    P = block_matrix(
        1, m, [matrix(QQbar, sp[1].matrix().T).QR(full=False)[0] for sp in esp])
    return evs, P


Graph.Spectral_decomposition = Spectral_decomposition


def angle(G):
    ei, sp = Spectral_decomposition(G)
    return ei, matrix([[r.norm() for r in sp.subdivision(0, i)] for i in range(len(ei))])


Graph.angle = angle


def sort_vertice_by_perron_vector(G, T="A"):
    V = range(G.order())
    osub = orbit_subdivision(G)
    r, v = G.eigenpairs(type=T, order=0)
    sos = sorted([(v[t[0]], t) for t in osub], reverse=True)
    spt = [t[1] for t in sos]
    bk = map(len, spt)
    H = G.relabel(dict(zip(flatten(spt), V)), inplace=False)
    LB = subd_list(V, bk)
    D = dict(zip(rainbow(len(LB)), LB))
    return H, H.plot(vertex_colors=D)


Graph.sort_vertice_by_perron_vector = sort_vertice_by_perron_vector


def Add_path(self, vertices, copy=False):
    if copy:
        H = copy(self)
    else:
        H = self
    if not vertices:
        return
    vert1 = vertices[0]
    for v in vertices[1:]:
        H.add_edge(vert1, v)
        vert1 = v
    H.set_pos(H.layout_default())
    if copy:
        return H


Graph.Add_path = Add_path


def AddHPath(G,r,h,H,t,l):
    HP = graphs_list.union_graphs([H]*l)
    for i in range(l-1):
        HP.merge_vertices([(i,t),(i+1,h)])
    G=G.sew_graphs(HP,{r:[(0,h)]},relabel=False)
    return G.Relabel()
Graph.AddHPath=AddHPath


def AddHPaths(G,h,H,t,D):
    """在图G中的某些点上添加H边路

    Args:
        G ([Graph]): [被添加H边路的图]
        h ([vertex]): [H的首点]
        H ([Graph]): [做为边的图H]
        t ([type]): [H的尾边]
        D ([dict]): [顶点与路长列表构成的字典]

    Returns:
        [Graph]: [所得图]
    """    
    for j in D:
        for i in D[j]:
            G=G.AddHPath(j,h,H,t,i)
    return G
Graph.AddHPaths=AddHPaths




def sorted_label(H, v, rev=False):
    u"""
    根据向量v中数量大小对图的顶点重新编号
    rev=False时， 从小到大排，否则从大到小排
    """
    repl = dict(zip([t[0] for t in sorted(enumerate(v/v[0]),
                                          lambda x, y: cmp(x[1], y[1]), reverse=rev)], range(len(v))))
    G = H.relabel(repl, inplace=False)
    return G


Graph.sorted_label = sorted_label


def get_appand_paths(H):
    """
    生成图的所有悬挂路的列表
    """
    paths = []
    V1 = [v for v in H.vertices() if H.degree(v) == 1]
    for v in V1:
        con = True
        appand_path = [v]
        while con:
            N = H[v]
            if len(N) <= 2:
                N1 = list(Set(N)-Set(appand_path))
                appand_path = appand_path+N1
                v = N1[0]
            else:
                con = False
        paths.append(appand_path)
    return paths


Graph.get_appand_paths = get_appand_paths

# 计算商矩阵，


def Quotient_Matrix(A, Bk=[], sym=false, **kwds):
    """
    计算矩阵关于给定划分的商矩阵。

    参数：
    A (matrix)：输入的矩阵。
    Bk (list, optional)：用于划分矩阵的子块列表。默认为空列表。
    sym (bool, optional)：是否对称化输出的矩阵。默认为False。
    **kwds：其他关键字参数，用于传递给Subdivide函数。

    返回值：
    Q (matrix)：商空间矩阵。
    B (list)：划分后的子块列表。
    D (matrix)：对角矩阵，用于对称化输出的矩阵。仅在sym为True时返回。

    """
    if Bk == []:
        diag = A.diagonal()
        diagK = groupbyinv(range(len(diag)), lambda i: diag[i])[0]
        P1 = SetPartition(Graph(A).automorphism_group(
            edge_labels=True).orbits())
        P2 = SetPartition(diagK.values())
        Bk = list(map(list, P1*P2))
    r = len(Bk)
    B = Subdivide(A, Bk, **kwds)
    ns = [(len(Bk[i]) if type(Bk[i]) == list else Bk[i]) for i in range(r)]
    Q = matrix([[add(add(B.subdivision(i, j)))/ns[i]
                 for j in range(r)] for i in range(r)])
    if sym:
        D = diagonal_matrix(list(map(sqrt, ns)))
        Q = D*Q*D ^ -1
        Q=1/2*(Q+Q.T)
        return Q, B, D ^ -1
    return Q, B


def Count_Subgraph(G, H):
    m = (H).subgraph_search_count(H)
    return (G).subgraph_search_count(H)/m


Graph.Count_Subgraph = Count_Subgraph


def HammingGraph(d, q):
    return graphs.CompleteGraph(q) ^ d


sage.graphs.graph_generators.graphs.HammingGraph = HammingGraph


def polynomials_for_addstarsgraph(self, D, type="L"):
    k = 1
    M = self.get_Matrix(type)
    n = self.order()
    DM = matrix(n, {(i, i): D[i] for i in D})
    return ((x * identity_matrix(n) - (k+(1/(x-k)))*DM-M).det()*((x-k) ^ DM.trace())).factor()


Graph.polynomials_for_addstarsgraph = polynomials_for_addstarsgraph


# 对图序列Gs按不变量Inv(G到RR的一个函数)排序(升序)，其中图格式可为图的实例，也可为图的graph6_string格式，重复图被滤掉
def sortGraphbyInvariant(Gs, Invs, rev=False, ndigits=5):
    if type(Invs) == type([]):
        if type(Gs[0]) != type([]):
            DT = [[T]+[inv(T) for inv in Invs] for T in Gs]
        else:
            DT = [[T]+[inv(T[0]) for inv in Invs] for T in Gs]
    else:
        if type(Gs[0]) != type([]):
            DT = [[T]+[Invs(T)] for T in Gs]
        else:
            DT = [[T]+[Invs(T[0])] for T in Gs]
    ret = sorted(DT, key=lambda x: x[1], reverse=rev)
    for dt in ret:
        dt[1] = round(dt[1], ndigits)
    return ret


graphs_list.sortGraphbyInvariant = sortGraphbyInvariant


def abs_polynomial(f):
    P = ZZ['x']
    return P(map(abs, P(f).coeffs()))


def rev_polynomial(f, n=None):
    P = ZZ['x']
    pf = P(f)
    m = pf.degree()
    if n == None or n < m:
        n = m
    return P((P(f).coeffs()+[0]*(n-pf.degree()))[::-1])


# 剖分加星能量比较,其中s为剖分次数, t为加星的边数, eh,eg分别为图H,G的被剖分边, vh,vg分别为图H,G的加星点(用数值计算, 一定要对被积函数进行化简，否则可能计算不正确）Ene
def De_subdivide_claw(H, eh, G, eg, s, vh, vg, t, wep=False, epsabs=1e-9, epsrel=1e-9):
    Hn = H.copy()
    Gn = G.copy()
    Hn.subdivide_edge(eh, s)
    Gn.subdivide_edge(eg, s)
    hs0 = Hn.absch_addclaws({vh: var("k")})
    gs1 = Gn.absch_addclaws({vg: var("k")})
    f = (hs0/gs1).simplify_full()
    rev, eps = numerical_integral(
        log(f.subs(k=t)), [0, +oo], eps_abs=epsabs, eps_rel=epsrel)
    if wep:
        return (N(rev/(2*pi), prec=100, digits=50), N(eps/(2*pi), prec=100, digits=50))
    else:
        return N(rev/(2*pi), prec=100, digits=50)


def HammingGraph(d, q):
    return graphs.CompleteGraph(q) ^ d


sage.graphs.graph_generators.graphs.HammingGraph = HammingGraph


def De_subdivide_edge(H, eh, G, eg, s, vh, vg, wep=False, epsabs=1e-9, epsrel=1e-9):
    Hn = H.copy()
    Gn = G.copy()
    Hn.subdivide_edge(eh, s)
    Gn.subdivide_edge(eg, s)
    H1 = Hn.Delete_Vertices([vh])
    G1 = Gn.Delete_Vertices([vg])
    h = H1.abscharpoly()
    g = G1.abscharpoly()
    f = h/g
    rev, eps = numerical_integral(
        log(f), [0, +oo], eps_abs=epsabs, eps_rel=epsrel)
    if wep:
        return (N(rev/(2*pi), prec=100, digits=50), N(eps/(2*pi), prec=100, digits=50))
    else:
        return N(rev/(2*pi), prec=100, digits=50)

# 二分图的能量差的计算


def Energy_Diff(H, G, wep=False, epsabs=1e-9, epsrel=1e-9):
    h = H.abscharpoly()
    g = G.abscharpoly()
    f = h/g
    rev, eps = numerical_integral(
        log(f), [0, +oo], eps_abs=epsabs, eps_rel=epsrel)
    if wep:
        return (N(rev/(2*pi), prec=100, digits=50), N(eps/(2*pi), prec=100, digits=50))
    else:
        return N(rev/(2*pi), prec=100, digits=50)

# 计算图或实多项式的能量


def Energy(self, poly=False, wep=False):
    if "matrix" in str(type(self)):
        return round_with_eps(add(matrix(RDF, self).singular_values()))
    if poly:
        if isinstance(self, Graph):
            f = self.charpoly()
            n = self.order()
    if "graph" in str(type(self)):
        return round_with_eps(add(map(abs, np.linalg.eigvalsh(self.am()))))
    if "expression" in str(type(self)):
        f = ZZ['x'](self)
        n = f.degree()
    if "polynomial" in str(type(self)):
        f = self
        n = f.degree()
    else:
        f = ZZ['x'](self)
        n = f.degree()
    L = rev_polynomial(f, n).coeffs()
    p = ZZ['x']([L[i]*(-I) ^ i*(i % 2 == 0) for i in range(len(L))])
    q = ZZ['x']([L[i]*(-1) ^ floor((i+1)/2)*(i % 2 == 1)
                 for i in range(len(L))])
    rev, eps = numerical_integral(1/(x ^ 2) * log(p ^ 2+q ^ 2), -Infinity, +
                                  Infinity, max_points=200, eps_abs=1e-12, eps_rel=1e-12, rule=4)
    if wep:
        return (N(rev/(2*pi), prec=100, digits=50), N(eps/(2*pi), prec=100, digits=50))
    else:
        return N(rev/(2*pi), prec=100, digits=50)


Graph.Energy = Energy


def reversed_polynomial(f, d):
    co = f.coeffs()
    if (d-len(co)+1) > 0:
        co = co+(d-len(co)+1)*[0]
    rec = list(reversed(co))
    return add([x ^ k*rec[k] for k in range(len(rec))])


def energy_polynomial(f, n):
    rf(x) = reversed_polynomial(f, n)
    assume(x, 'real')
    p(x) = real_part(rf(-i*x), hold=True).simplify()
    q(x) = imag_part(rf(-i*x), hold=True).simplify()
    return RR((1/(2*pi))*integral(1/(x ^ 2) * log(p(x) ^ 2+q(x) ^ 2), (x, -oo, +oo)))
# 这里的积分integral算法可能判断收敛性出错


def abspm(f):
    cl = f.dict()
    for c in cl:
        cl[c] = abs(cl[c])
    ZP.<x > = ZZ[]
    return ZP(cl)


def abscharpoly(self):
    """输入:图 G

    Returns:
        [多项式]: [图G的系数取绝对值的特征多项式]
    """
    f = self.charpoly()
    return abspm(f)


Graph.abscharpoly = abscharpoly


def absch_addclaws(self, S):
    S = dict(S)
    key = S.keys()
    Tw = add(S.values())
    f = 0*x
    for sub in Subsets(key):
        H = self.copy()
        H.delete_vertices(sub)
        f = f + mul([S[i] for i in sub])*x ^ (-len(sub))*abspm(H.charpoly())
    F = x ^ Tw*f
    F = F.expand()
    return F


Graph.absch_addclaws = absch_addclaws


def sum_connectivity(self):
    return add([(self.degree(e[0])+self.degree(e[1])) ^ -0.5 for e in self.edges()])


Graph.sum_connectivity = sum_connectivity


# =======================图的矩阵表示与谱:结束==================================

# =======================图的构造与变换:开始==================================

def Add_cycles(self, *args, **kwargs):
    """
    将圈添加到图中。

    参数:

    args: Union[int, List[int], Dict[int, List[int]]]，位置参数，表示要添加的圈。
    如果是一个整数，表示要添加一个长度为该整数的圈。
    如果是一个列表，表示要添加给定顶点的圈。
    如果是一个字典，表示要添加多个圈，键为起始顶点，值为包含顶点的列表。
    kwargs: bool，关键字参数，表示是否复制图对象。
    copy: 如果为True（默认），在原图的副本上进行操作。
    如果为False，直接在原图上进行操作。
    返回值:

    Graph，添加圈后的图对象
    """
    relabel = kwargs.pop('copy', True)
    if copy:
        H = self.copy()
    else:
        H = self
    if(len(args) == 1):
        D = args[0]
        if isinstance(args[0], dict):
            for j in D:
                for i in D[j]:
                    H.add_cycle([j]+list(range(H.order(), H.order()+i-1)))
            return H
        else:
            H.add_cycle(D)
            return H
    else:
        return H


Graph.Add_cycles = Add_cycles


def Subdivide_Edge(G,e,k):
    H=copy(G)
    H.subdivide_edge(e,k)
    return H
Graph.Subdivide_Edge=Subdivide_Edge


def partStrongProduct(G,sV):
    K2=graphs.PathGraph(2)
    sG=G.subgraph(sV)
    rG=G.Delete_Vertices(sV)
    pG=sG.strong_product(K2)
    rV=rG.vertices()
    pV=pG.vertices()
    nE=[e for e in cartesian_product([pV,rV]) if e[1] in G[e[0][0]]]+rG.edges(labels=False)+pG.edges(labels=False)
    nG=Graph(nE)
    nG.relabel()
    return nG
Graph.partStrongProduct=partStrongProduct

def BiRelationProduct(G,H,R):
    E=G.disjoint_union(H).edges(labels=False)
    Er=[((0,e[0]),(1,e[1])) for e in R]
    resG=Graph(E+Er)
    resG.relabel()
    return resG

Graph.BiRelationProduct=BiRelationProduct   


def Add_path(self, L, copy=True):
    if copy:
        H = self.copy()
    else:
        H = self
    for l in L:
        H.add_path(l)
    return H


DiGraph.Add_path = Add_path
Graph.Add_path = Add_path


def add_paths(self, D, copy=True):
    """"
    Returns the graph obtained from Graph by attaching some append paths to its vertices.
    EXAMPLES::

        sage: G=graphs.CycleGraph(5)
        sage: H=G.add_paths({0:[2,3],3:[3,5],4:[2,2]})
        sage: show(H,layout="spring")

    """
    D = dict(D)
    if copy:
        H = self.copy()
    else:
        H = self
    for j in D:
        for i in D[j]:
            H.add_path([j]+[H.order()..H.order()+i-1])
    return H+Graph()


Graph.add_paths = add_paths
Graph.Add_Paths = add_paths


def Delete_Vertices(self, vertices, copy=True):
    if copy:
        H = self.copy()
    else:
        H = self
    H.delete_vertices(vertices)
    return H


Graph.Delete_Vertices = Delete_Vertices
DiGraph.Delete_Vertices = Delete_Vertices


def Delete_Edges(self, edges, copy=True):
    if copy:
        H = self.copy()
    else:
        H = self
    H.delete_edges(edges)
    if copy:
        return H


Graph.Delete_Edges = Delete_Edges
DiGraph.Delete_Edges = Delete_Edges


def Add_Edges(self, edges, copy=True):
    if self.weighted:
        Edges = [list(e)+[1] if len(e) == 2 else list(e) for e in edges]
        edges = Edges
    if copy:
        H = self.copy()
    else:
        H = self
    H.add_edges(edges)
    if copy:
        return H


Graph.Add_Edges = Add_Edges
BipartiteGraph.Add_Edges = Add_Edges
DiGraph.Add_Edges = Add_Edges


def add_stars(self, D, copy=True):
    u"""
    Attach some append edges at some vetices of graph G. 中文

    EXAMPLES: Attach 2 append edges at vertex 1 and 3 append edges at vertex 3 of $P_5$.

    ::

        sage: G=graphs.PathGraph(5)
        sage: H=G.add_stars([[1,2],[3,3]])
        sage: show(H,layout="spring")

    """
    D = dict(D)
    if copy:
        H = self.copy()
    else:
        H = self
    label = [1] if H.weighted() else []
    for j in D:
        for i in range(D[j]):
            H.add_edge([j, H.order()]+label)
    return H


Graph.add_stars = add_stars


def seidel_switching(self, L):
    U = Set(L)
    G = copy(self)
    VG = Set(G.vertices())
    V = VG.difference(U)
    for u in U:
        for v in V:
            if G.has_edge(u, v):
                G.delete_edge(u, v)
            else:
                G.add_edge(u, v)
    return G


Graph.seidel_switching = seidel_switching


def seidel_matrix(self):
    return (self.complement()).adjacency_matrix()-self.adjacency_matrix()


Graph.seidel_matrix = seidel_matrix


# 对图进行标号
def Relabel(H, L=None, inv=False):
    """对图H 的顶点集 V进行重新标号
    Args:
        H ([Graph|DiGraph]): [要重标号的图]
        L ([列表|字典], optional): [参数若为列表则将 L[i]标号改为 i; 若为字典则将标号L[i]标为 i ]. Defaults to None.
        inv (bool, optional): [是否逆向标号]. Defaults to False.

    Returns:
        [type]: [重新标号的图]
    """
    if L == None:
        L = sorted(H.vertices(sort=True))
    G = copy(H)
    n = H.order()
    if type(L) == list:
        oL = L+list(Set(range(H.order()))-Set(L))
        if inv:
            nL = range(n)[::-1]
        else:
            nL = range(n)
        rebdict = dict(zip(oL, nL))
    if type(L) == dict:
        K = list(L.keys())
        V = list(L.values())
        lK = Set(K+V)-Set(K)
        lV = Set(K+V)-Set(V)
        otV = Set(range(n))-Set(K+V)
        L.update(dict(list(zip(lK, lV))+list(zip(otV, otV))))
        rebdict = L
    G.relabel(rebdict)

    return G


Graph.Relabel = Relabel
DiGraph.Relabel = Relabel
# 对图进行边赋权


def set_edge_weight(self, *args, **kwoptions):
    """
    set the weights for edges in graph.

    sage: G=graphs.CycleGraph(10)
    sage: edges_weights={(1,2):3,4:[(5,4),(3,2)],(7,8):2}
    sage: G.set_edge_weight(edges_weights,default_weight=1)
    sage: G.show(edge_labels=True)
    sage: G.weighted_adjacency_matrix()
    [0 1 0 0 0 0 0 0 0 1]
    [1 0 3 0 0 0 0 0 0 0]
    [0 3 0 4 0 0 0 0 0 0]
    [0 0 4 0 1 0 0 0 0 0]
    [0 0 0 1 0 4 0 0 0 0]
    [0 0 0 0 4 0 1 0 0 0]
    [0 0 0 0 0 1 0 1 0 0]
    [0 0 0 0 0 0 1 0 2 0]
    [0 0 0 0 0 0 0 2 0 1]
    [1 0 0 0 0 0 0 0 1 0]
    sage: H.weighted_adjacency_matrix()
    [0 1 0 1 1 0 1 1 0 0 0 1]
    [1 0 1 0 0 0 0 0 0 0 0 0]
    [0 1 0 1 0 0 0 0 0 0 0 0]
    [1 0 1 0 0 0 0 0 0 0 0 0]
    [1 0 0 0 0 1 0 0 0 0 0 0]
    [0 0 0 0 1 0 1 0 0 0 0 0]
    [1 0 0 0 0 1 0 0 0 0 0 0]
    [1 0 0 0 0 0 0 0 1 0 0 0]
    [0 0 0 0 0 0 0 1 0 1 0 0]
    [0 0 0 0 0 0 0 0 1 0 1 0]
    [0 0 0 0 0 0 0 0 0 1 0 1]
    [1 0 0 0 0 0 0 0 0 0 1 0]
    """
    ew = []
    if len(args) == 0:
        edges_weights = {}
    else:
        edges_weights = args[0]
    default_weight = kwoptions.pop("default_weight", 1)
    self.weighted(True)
    for i in edges_weights:
        if type(i) == type((0, 0)):
            self.set_edge_label(i[0], i[1], edges_weights[i])
            ew.append(i)
        else:
            ew = ew+edges_weights[i]
            for edges in edges_weights[i]:
                self.set_edge_label(edges[0], edges[1], i)
    for u, v, l in self.edges(sort=True):
        if l == None and (u, v) not in ew:
            self.set_edge_label(u, v, 1)


Graph.set_edge_weight = set_edge_weight
DiGraph.set_edge_weight = set_edge_weight


# 定义一个图的Path Tree
def path_tree(G, v):
    Paths = flatten([G.all_paths(v, u) for u in G.vertices()], max_level=1)
    A = matrix(len(Paths), lambda i, j: Paths[i] == Paths[j][:-1])
    return Graph(A+A.T)


Graph.path_tree = path_tree


def edge2graph(G, e, v1, v2, H):
    E = G.edges(labels=False)
    Vh = H.vertices()
    if (e not in E) or (v1 not in Vh) or (v2 not in Vh):
        print("Please see the usage of method")
        return
    for v in (v1, v2):
        Vh.remove(v)
    dic = dict([(i, "%d,%d:%d" % (e[0], e[1], i))
                for i in Vh]+[(v1, e[0]), (v2, e[1])])
    H1 = H.relabel(dic, inplace=False)
    E1 = H1.edges()
    G1 = G.Delete_Edges([e])
    G1 = G1.Add_Edges(E1)
    return G1


def replace_edges_by_graph(self, GL, Es, Vp, relabel=False):
    """
    return a graph which obtained from G by replacing some
    edges by graphs

    sage: Graph.replace_edges_by_graph=replace_edges_by_graph
    sage: GL=[graphs.CycleGraph(5)]*5
    sage: G=graphs.CycleGraph(5)
    sage: Es=G.edges()
    sage: Vp=[(0,2)]*5
    sage: G1=G.replace_edges_by_graph(GL,Es,Vp,relabel=True)
    sage: G1.show(layout="spring")

    """
    tm = zip(GL, Es, Vp)
    nedges = []
    for it in tm:
        H, e, vp = it
        V = H.vertices()
        dic = dict([(v, "%d,%d:%d" % (e[0], e[1], v))
                    for v in V if v not in vp]+zip(vp, e))
        H1 = H.relabel(dic, inplace=False)
        nedges = nedges+H1.edges()
    T1 = self.Delete_Edges(Es)
    T1.add_edges(nedges)
    if relabel:
        T1.relabel()
        return T1
    else:
        return T1


Graph.replace_edges_by_graph = replace_edges_by_graph

# 在图类中寻找同谱图


def findCospe(glist, *args, **kwargs):
    TP = kwargs.pop('type', "A")
    if len(args) == 0 or args[0] == None:
        tCharPolys = {}
        for H in glist:
            charPoly = H.charpolys_for_graph(type=TP)
#             st=H.graph6_string()
            if charPoly in tCharPolys:
                tCharPolys[charPoly].append(H)
            else:
                tCharPolys[charPoly] = [H]
        cosplt = [tCharPolys[k]
                  for k in tCharPolys.keys() if len(tCharPolys[k]) > 1]
        return cosplt
    else:
        G = args[0]
        chp = G.charpolys_for_graph(type=TP)
        for H in glist:
            charPoly = H.charpolys_for_graph(type=TP)
            if charPoly == chp:
                return H
        print("There exist no graphs copectral with G")


graphs_list.findCospe = findCospe
#
#
# #在图类中寻找同谱图
# def findCospe(glist,G,type="L"):
#     if G==None:
#         tCharPolys = {}
#         for H in glist:
#             charPoly = H.charpolys_for_graph(type=type)
#             st=H.graph6_string()
#             if charPoly in tCharPolys:
#                 tCharPolys[charPoly].append(st)
#             else:
#                 tCharPolys[charPoly] = [st]
#         cosplt=[tCharPolys[k] for k in tCharPolys.keys() if len(tCharPolys[k]) > 1]
#         return cosplt
#     else:
#         chp=G.charpolys_for_graph(type=type)
#         for H in glist:
#             charPoly = H.charpolys_for_graph(type=type)
#             st=H.graph6_string()
#             if charPoly==chp:
#                 return H
#         print("There exist no graphs copectral with G");
#
# graphs_list.findCospe=findCospe

# 广义的H联图


def H_join(self, GL, UL=None, relabel=True):
    """
    对图进行H-join操作，生成一个新的图。

    参数:
    - GL: List[Graph] 或 dict，要进行H-join操作的图列表。如果是一个列表，每个元素是一个图对象；如果是一个字典，键是图的索引，值是图对象。
    - UL: List[List[int]]，可选参数，默认值为None。图的顶点列表，表示要与H-join的图进行连接的顶点子集。
    - relabel: bool，可选参数，默认值为True。指定是否重新标记生成的图的顶点。

    返回值:
    - Graph，生成的新图对象。
    """
    if type(GL) == list:
        if self.order() != len(GL):
            print("图H的阶数要与图列表的长度相同")
            return
    if type(GL) == dict:
        gl = [Graph(1)]*self.order()
        for k in GL:
            gl[k] = GL[k]
        GL = gl
    V = self.vertices(sort=False)
    VL = []
    EL = []
    i = 0
    nUL = []
    for G in GL:
        VL.append([(i, v) for v in G.vertices(sort=False)])
        for (u, v) in G.edges(labels=False,sort=False):
            EL.append(((i, u), (i, v)))
        i = i+1

    if UL:
        i = 0
        for U in UL:
            nUL.append([(i, v) for v in U])
            i = i+1
    else:
        nUL = VL
    for (u, v) in self.edges(labels=False,sort=False):
        EL = EL+list(map(tuple, cartesian_product([nUL[u], nUL[v]])))
    G = Graph(0)
    G.add_vertices(flatten(VL, max_level=1))
    G.add_edges(EL)
    if relabel:
        G.relabel()
    return G


Graph.H_join = H_join

def Quo_Mat(G, type="A",sym=False):
    r"""返回图G在其轨道划分下的$\alpha$矩阵
    Args:
        G (Graph): [description]
        a (int, optional): [description]. Defaults to 0.
    Returns:
        [type]: [description]
    """
    M = G.get_Matrix(type)
    K = G.orbit_subdivision()
    return Quotient_Matrix(M, K, sym=sym)



Graph.Quotient_Matrix = Quo_Mat



def HJQM(H, dND):
    r"""计算H联图的alpha矩阵的商矩阵
    H.H_joion(GL), GL为替换H的各顶点的图序列, 可为字典格式${i:Gi|i \in V(H)}$

    Args:
        H (Graph)): H联图的基图
        dND (dict): GL中各图的(阶数,边数)的字典 （可含符号变量）
    Returns:
        matrix: H联图的alpha矩阵的商矩阵
    """
    a = var("alpha")
    n = H.order()
    N = ones_matrix(SR, n, 1)
    D = zero_matrix(SR, n, 1)
    for k in dND:
        N[k, 0] = dND[k][0]
        D[k, 0] = 2*dND[k][1]/dND[k][0]
    Qam = H.am()*diagonal_matrix(vector(N))
    Qd = 2*diagonal_matrix(D.list())+diagonal_matrix(add(Qam.T))
    Q = (1-a)*Qam+diagonal_matrix(D.list())+a*diagonal_matrix(add(Qam.T))
    return Q


Graph.HJQM = HJQM

# def add_innerPaths(self,u,v,L,relabel=True):
#     H=copy(self)
#     EL=[]
#     for p in range(len(L)):
#         L1=[u]+[(p,k) for k in range(L[p])]+[v]
#         EL=EL+[L1[i:i+2] for i in range(len(L1)-1)]
#     H=H.Add_Edges(EL,copy=True)
#     if relabel:
#         H.relabel()
#
#     return H


def add_innerPaths(G, *args, **kwargs):
    relabel = kwargs.pop('relabel', False)
    H = copy(G)
    EL = []
    if(len(args) == 1) and isinstance(args[0], dict):
        vp_pls = args[0]
        ct = 0
        for t in vp_pls.keys():
            ct = ct+1
            if isinstance(t, Integer):
                u, v = (t, t)
            elif isinstance(t, tuple):
                u, v = t
            data = vp_pls[t]
            if isinstance(data, Integer):
                L = [data]
            else:
                L = data
            for p in range(len(L)):
                L1 = [u]+[("%d_%d_%d" % (ct, p, k)) for k in range(L[p])]+[v]
                EL = EL+[L1[i:i+2] for i in range(len(L1)-1)]
    H = H.Add_Edges(EL, copy=True)
    H.set_pos(H.layout("spring"))
    if relabel:
        return H.Relabel(G.vertices(sort=True) +
                  sorted(Set(H.vertices(sort=False))-Set(G.vertices(sort=False))))
    else:
        return H


Graph.add_innerPaths = add_innerPaths
DiGraph.add_innerPaths = add_innerPaths


def join_graphs(GL, relabel=False):
    E = []
    V = []
    for k in range(len(GL)):
        for e in GL[k].complement().edges(labels=False):
            E = E+[((e[0], k), (e[1], k))]
        for v in GL[k].vertices():
            V = V+[(v, k)]
    H = Graph(0)
    H.add_vertices(V)
    H = H.complement()
    H.delete_edges(E)
    if relabel:
        H.relabel()
    return H


graphs_list.join_graphs = join_graphs


def part_join(G, H, sU=None, sV=None, relabel=True):
    
    E = []
    V = []
    GL = [G, H]
    for k in range(len(GL)):
        for e in GL[k].edges(labels=False,sort=False):
            E = E+[((k, e[0]), (k, e[1]))]
        for v in GL[k].vertices(sort=False):
            V = V+[(k, v)]
    if sU == None:
        sU = G.vertices(sort=False) 
    if sV == None:
        sV = H.vertices(sort=False)
    for u in sU:
        for v in sV:
            E = E+[((0, u), (1, v))]
    H = Graph(E)
    H.add_vertices(V)
    if relabel:
        H.relabel()
    return H


Graph.part_join = part_join


# 对图列表进行不交并运算
def union_graphs(GL, relabel=False):
    E = []
    V = []
    G = GL[0]
    weighted = G.weighted()
    if weighted:
        for G in GL:
            G.weight_edges()
    directed = G._directed
    Grp = DiGraph if directed else Graph
    for k in range(len(GL)):
        for e in GL[k].edges(sort=False):
            E = E+[((k, e[0]), (k, e[1]), e[2])]
        for v in GL[k].vertices(sort=False):
            V = V+[(k, v)]
    H = Grp(E)
    H.weighted(weighted)
    H.add_vertices(V)
    if relabel:
        H.relabel()
    return H


graphs_list.union_graphs = union_graphs

def graphs_union(GL,relabel=False):
        H=Graph(reduce(list.__add__,[[((i,e[0]),(i,e[1])) for e in GL[i].edges(sort=False,labels=False)] for i in range(len(GL))]))
        if relabel:
            H=H.Relabel()
        return H

graphs_list.graphs_union = graphs_union


def NumberOfWalks(G,Starts,Ends=None,nPass=[],Pass=[],l=1):
    """计算图G中起点在Starts中， 终点在Ends中长为l且经过Pass中的所有点， 不过nPass中的任意点的途径总数

    Args:
        G ([Graph]): 要考虑的图
        Starts ([list]): [起点列表]
        Ends ([list], optional): [终点列表]. Defaults to None.
        nPass (list, optional): [绕开点列表]. Defaults to [].
        Pass (list, optional): [经过点列表]. Defaults to [].
        l (int, optional): [途径长]. Defaults to 1.

    Returns:
        [int]: [途径数]
    """    
    if Ends==None:
        Ends=Starts
    A=G.am()
    B=copy(A)
    B[nPass,:]=0
    B[:,nPass]=0
    S=list(map(list,Set(Pass).subsets()))
    Res=0*B
    for s in S:
        T=copy(B)
        T[s,:]=0
        T[:,s]=0
        Res=Res+(-1)^len(s)*T^l
    return add(add(Res[Starts,Ends]))
Graph.NumberOfWalks=NumberOfWalks


#将图中的边替换为图
def edge2H(G,h,H,t):
    E=G.edges(sort=False,labels=False)
    GH=graphs_list.union_graphs([H]*G.size())
    GC=GH.connected_components_subgraphs()
    nE=[list((GC[i].Relabel({(i,h):E[i][0],(i,t):E[i][1]})).edges(sort=False,labels=False)) for i in range(len(E))]
    nG=Graph(reduce(list.__add__,nE))
    nG.relabel()
    return nG
Graph.edge2H=edge2H

def RandomBlockLikeGraph(H,h,t,m):
    G=H
    for _ in range(m):
        nG=graphs_list.union_graphs([G,H])
        u=randint(0,G.order()-1)
        nG.merge_vertices([(0,u),(1,h)])
        G=nG.Relabel()
    return G
Graph.RandomBlockLikeGraph=RandomBlockLikeGraph

# 旋转边变换
def Rotate_Edges(G, e):
    u, v = e[:2]
    E = G.edges()
    if Set(G[u]).intersection(Set(G[v])) != Set():
        return None
    else:
        E1 = [e for e in E if v in e and e[0] != u]
        E2 = []
        for e in E1:
            e = list(e)
            for i in range(2):
                if e[i] == v:
                    e[i] = u
                    E2.append(e)
        return G.Delete_Edges(E1).Add_Edges(E2)


Graph.Rotate_Edges = Rotate_Edges


def coalescence_graphs(GL, VL, relabel=True):
    G = graphs_list.union_graphs(GL)
    MV = [(i, VL[i]) for i in range(len(GL))]
    G.merge_vertices(MV)
    if relabel:
        G.relabel()
    return G


graphs_list.coalescence_graphs = coalescence_graphs


def rooted_product(self, DvH, relabel=True):
    rtV = [*(DvH.keys())]
    rtGs = DvH.values()
    HL = [t[1] for t in rtGs]
    rtHs = [t[0] for t in rtGs]
    G = graphs_list.union_graphs([self]+HL)
    vp = [[(0, rtV[i]), (i+1, rtHs[i])] for i in range(len(HL))]
    for p in vp:
        G.merge_vertices(p)
    if relabel:
        G.relabel()
    return G


Graph.rooted_product = rooted_product


def RootedProduct(G,H,v):
    return G.rooted_product({i:(v,H) for i in G.vertices()})
Graph.RootedProduct=RootedProduct   



def neps_graph(Gs, B):
    """
    Construct a Non-complete Expansion p-Sum (NEPS) of given graphs.
    :param graphs: A list of graphs G_1, ..., G_n.
    :param B: A set of n-tuples (binary vectors) which determines adjacency in the NEPS.
    :return: A graph representing the NEPS of the input graphs.
    """
    from itertools import product as pt
    
    # Create the Cartesian product of the vertex sets of the graphs
    vertices = list(pt(*(G.vertices() for G in Gs)))
    n = len(Gs)  # Number of graphs
    
    # Define the NEPS graph
    neps = Graph()
    neps.add_vertices(vertices)
    VC=Combinations(vertices,2)
    # Check pairs of vertices for adjacency
    for u,v in VC:
        # Check if there exists a beta in B such that u and v are adjacent according to beta
        if any(all((u[i] == v[i] if beta[i] == 0 else Gs[i].has_edge(u[i], v[i])) for i in range(n)) for beta in B):
            neps.add_edge(u, v)
    return neps


def Generilized_RootedProduct(G,H,W,r,B={(1,0),(0,1)}):
    if G.order()!=H.order():
        print("$G,H$必须是同阶的")
        return
    Gs=[H,W]
    neps=neps_graph(Gs, B)
    rV=[(v,r) for v in Gs[r]]
    rE=neps.subgraph(rV).edges()
    aE=G.Relabel({u:(u,r) for u in G.vertices()}).edges()
    neps1=neps.Delete_Edges(rE)
    neps1=neps1.Add_Edges(aE)
    return neps1
    
    





def PartProduct(self, sV, H, prod=Graph.cartesian_product):
    out = {i: [u for u in self.neighbors(i) if u not in sV] for i in sV}
    G = self.subgraph(sV)
    H0 = prod(G, H)
    V1 = H0.vertices()
    G1 = copy(self)
    G1.delete_vertices(sV)
    E1 = G1.edges(labels=False)
    E2 = []
    for u in V1:
        ol = out[u[0]]
        for v in ol:
            E2.append((u, v))
    H0.add_edges(E1+E2)
    return H0


Graph.PartProduct = PartProduct


# 计算图的幂（按卡氏积）
def Graph_Power(G, n):
    H = reduce(Graph.cartesian_product, [G]*n)
    H.relabel(map(tuple, map(flatten, H.vertices())))
    return H

def Join(G,H):
    if isinstance(H, Graph):
        return Graph.join(G, H, labels="integers")
    if isinstance(H, dict):
        return G.H_join(H)
Graph.__or__ = Join
Graph.__pow__ = Graph_Power
Graph.__invert__ = Graph.complement
Graph.__neg__ = Graph.complement



def RandBipartite(m, n):
    import random
    import string
    q = False
    while not q:
        G = BipartiteGraph(random_matrix(GF(2), m, n))
        q = G.is_connected()
    G.name(''.join(random.choice(string.lowercase) for x in range(6)))
    return G


def Partitioned_tensor_product(Mp, *args):
    L = []
    if not isinstance(Mp, list):
        L = [Mp]+list(args)
    else:
        L = Mp
    bt = [len(t)+_sage_const_1 for t in L[_sage_const_0].subdivisions()]
    Ablks, Bblks = [[m.subdivision(i, j) for i in range(
        bt[_sage_const_0]) for j in range(bt[_sage_const_1])] for m in L]
    return block_matrix(bt[_sage_const_0], bt[_sage_const_1], [Ablks[i].tensor_product(Bblks[i]) for i in range(len(Ablks))])


def Khatri_Rao_Product(A, B):
    bt1 = map(len, A.subdivisions())
    bt2 = map(len, B.subdivisions())
    if bt1 != bt2:
        print("输入的分块矩阵需要具有相同的行块数和列块数!")
        return
    else:
        return block_matrix(bt1[0]+1, bt1[1]+1, [A.subdivision(i, j).tensor_product(B.subdivision(i, j)) for i in [0..bt1[0]] for j in [0..bt1[1]]])


def Tracy_Singh_Product(A, B):
    bt1 = map(len, A.subdivisions())
    bt2 = map(len, B.subdivisions())
    R = [(i, j) for i in [0..bt1[0]] for j in [0..bt2[0]]]
    C = [(i, j) for i in [0..bt1[1]] for j in [0..bt2[1]]]
    return block_matrix(len(R), len(C), [A.subdivision(r[0], c[0]).tensor_product(B.subdivision(r[1], c[1])) for r in R for c in C])

# 对分块矩阵进行置换


def perm_block(A, R, C=None):
    if not C:
        C = R
    r, c = len(R), len(C)
    return block_matrix(r, c, [A.subdivision(R[i], C[j]) for i in range(len(R)) for j in range(len(C))])


# 判断图H,G拟序是否可比
def checksign1(f): return Set(map(lambda n: sgn(n), f.coefficients()))


def is_quasiorder_comparable(H, G):
    f = H.abscharpoly()-G.abscharpoly()
    return checksign1(f)


graphs_list.is_quasiorder_comparable = is_quasiorder_comparable


# 计算k剖分二分图的绝对值特征多项式， 其中h0为基图的绝对值特征多项式，h1为1剖分图的绝对值特征多项式
def gen_poly(k, h0, h1):
    P = ZZ['x']
    var("y")
    equ = y ^ 2-x*y-1 == 0
    rt = equ.solve(y, solution_dict=True)
    w1, w2 = (t[y] for t in rt)
    return P(((w1 ^ k*(w2*h0-h1)+w2 ^ k*(h1-w1*h0))/(w2-w1)).simplify_full())


def gen_poly1(k, h0, h1):
    P = ZZ['x']
    var("y")
    equ = y ^ 2-x*y-1 == 0
    rt = equ.solve(y, solution_dict=True)
    w1, w2 = (t[y] for t in rt)
    return P(((h1*(w2 ^ k-w1 ^ k)+h0*(w2 ^ (k-1)-w1 ^ (k-1)))/(w2-w1)).simplify_full())


def abscharpoly_subdivide_edge(H, eh, k):
    H1 = H.copy()
    H1.subdivide_edge(eh, 1)
    h0 = H.abscharpoly()
    h1 = H1.abscharpoly()
    return gen_poly1(k, h0, h1)


# 计算路的绝对值特征多项式

def abscharpoly_path(n):
    var("y")
    equ = y ^ 2-x*y-1 == 0
    rt = equ.solve(y, solution_dict=True)
    w1, w2 = (t[y] for t in rt)
    return ((w2 ^ (n+1)-w1 ^ (n+1))/(w2-w1)).simplify_full()


def abscharpoly_path(n):
    return add([binomial(n-k, k)*x ^ (n-2*k) for k in [0..floor(n/2)]])


# 在图H1,H2的顶点u,v间连一条长为k的路，其中H1,H2的顶点集都是自然数，返回的图中新内部路顶点为(1,1),...,(1,k-1);
# H1中的u点对就顶点为(0,u);H2中的v点对应顶点为(2,v)
# def linkGraph(H1,H2,u,v,k,relabel=true):
#     P2=graphs.PathGraph(k+1)
#     G=graphs_list.union_graphs([H1,P2,H2])
#     G.merge_vertices([(0,u),(1,0)])
#     G.merge_vertices([(2,v),(1,k)])
#     if relabel:
#         G.relabel()
#     return G

def sew_graphs(G, H, Vpairs, relabel=True):
    """
    将图G的一些顶点与H的一些点子集合并， 其中参数Vpairs是一个字典, 键是G中的点， 值是H中的点子集
    """
    G = graphs_list.union_graphs([G, H])
    for v in Vpairs:
        G.merge_vertices([(0, v)]+[(1, u) for u in Vpairs[v]])
    if relabel:
        G.relabel()
    return G

Graph.sew_graphs=sew_graphs
graphs_list.sew_graphs = sew_graphs


def linkGraph(H1, H2, u, v, k, relabel=true):
    if k not in NN:
        print("剖分次数应为非负整数！")
        return
    G = graphs_list.union_graphs([H1, H2])
    G = G.Add_Edges([((0, u), (1, v))], copy=True)
    if k == 0:
        G.merge_vertices([(0, u), (1, v)])
    else:
        G.subdivide_edge(((0, u), (1, v)), k-1)
    if relabel:
        G.relabel()
    return G


def G_starlike_H(G,u,K,v,H):
    """在图G的u点长若干条悬挂路， 每个悬挂路的一度点再与H的一个拷贝中的v点粘接

    Args:
        G ([Graph]): [图G]
        u ([int]): [顶点]
        K ([list]): [这些悬挂路的长]
        v ([int]): [顶点v]
        H ([Graph]): [图H]

    Returns:
        [Graph]: [所得图]
    """    
    H0=G.add_paths({u:K})
    H1=H*len(K)
    n=G.order()
    m=H.order()
    D={n+add(K[:i+1])-1:[m*i+v] for i in range(len(K))}
    return sew_graphs(H0,H1,D)
Graph.G_starlike_H=G_starlike_H


def caterpillar(d, D):
    return graphs.PathGraph(d+1).add_stars(D)


sage.graphs.graph_generators.graphs.caterpillar = caterpillar


def Random_Unicyclic_Graph(bk, c):
    Trees = [graphs.RandomTree(i) if i > 1 else Graph(1) for i in bk]
    G = graphs.CycleGraph(c)
    n = G.order()
    Gr = graphs_list.union_graphs([G]+Trees)
    Hr = copy(Gr)
    mergeparts = [((0, randint(0, n-1)), (k+1, 0)) for k in range(len(Trees))]
    DG = DiGraph(mergeparts)
    Ccomp = DG.connected_components()
    for vs in Ccomp:
        Gr.merge_vertices(vs)
    Gr.relabel()
    return Gr


sage.graphs.graph_generators.graphs.Random_Unicyclic_Graph = Random_Unicyclic_Graph


def RandomKCyclicGraph(n, k):
    """
    生成一个随机的n阶连通的k圈图。
    """
    if not 0 <= k <= binomial(n-1, 2):
        raise ValueError('参数k超出合理范围！')
    H = graphs.RandomTree(n)
    cE = list(H.complement().edges(sort=False))
    return H.Add_Edges([cE.pop(randint(0, len(cE)-1)) for _ in range(k)])


sage.graphs.graph_generators.graphs.RandomKCyclicGraph = RandomKCyclicGraph


def rand_cactus_graph(cn, ml):
    G = Graph(1)
    for i in [1..cn]:
        v = randint(0, G.order()-1)
        H = graphs.CycleGraph(randint(3, ml))
        u = randint(0, H.order()-1)
        G = graphs_list.coalescence_graphs([G, H], [v, u])
    return G


sage.graphs.graph_generators.graphs.rand_cactus_graph = rand_cactus_graph


# 对图的邻接矩阵进行赋高斯整数权
def weight_edges(G, *args):
    G.weighted(True)
    A = G.am()
    for e in G.edges():
        if e[2] == None:
            G.set_edge_label(e[0], e[1], 1)
            A[e[:2]] = 1
    if len(args) == 0:
        return A

    for k in D:
        A[k] = D[k]
        if G.isinstance(Graph):
            A[k[::-1]] = D[k]
        G.set_edge_label(k[0], k[1], D[k])
    return A


Graph.weight_edges = weight_edges
DiGraph.weight_edges = weight_edges

def get_edges_weight(G,e=None):
    E=G.edges()
    weight_edges={(e[0],e[1]):e[2] for e in E}
    if e==None:
        return weight_edges
    else:
        return weight_edges[e]
Graph.get_edges_weight=get_edges_weight


def complex_weight(G, E):
    D = dict(flatten([[[e, ZZ[I]([0, 1])], [e[::-1], -ZZ[I]([0, 1])]]
                      for e in E], max_level=1))
    return weight_edges(G, D)


Graph.complex_weight = complex_weight


def BaseGraph(G,relabel=True):
    u'''
    Return the BaseGraph of C-cyclic graph G   

    '''
    redG = copy(G)
    while min(redG.degree_sequence()) == 1:
        for v in redG.vertices(sort=False):
            if redG.degree(v) == 1:
                redG.delete_vertex(v)
    if relabel:
        redG=redG.Relabel()
    return redG


Graph.BaseGraph = BaseGraph


def slip_graph(self, H, d, i):
    """
    将图$G$上的一条悬挂路$P_d$上滑移图$H$到顶点$v_i$
    $G$接在$P_{d+1}$的顶点$v_1$上面。

    INPUT:

    -``H``--移动的子图
    -``d``--路的长度
    -``i``--移动$H$到$P_{d+1}$的顶点$v_{i+1}$上，$v_1$到$v_{i+1}$的距离恰为$i$

    OUTPUT:

    -:class `Graph`

    EXAMPLES:

        sage: G=graphs.CycleGraph(6)
        sage: H=graphs.StarGraph(5)
        sage: sG=G.slip_graph(H,10,3)

    """
    if i > d-1:
        return
    H1 = graphs_list.coalescence_graphs([graphs.PathGraph(d+1), self], [1, 0])
    H2 = graphs_list.coalescence_graphs([H1, H], [i+1, 0])
    H2.name(r"$\Gamma_{%d}(%d,%d)$" % (d, i, d-i-1))
    return H2


Graph.slip_graph = slip_graph


def pick_vertices(G, *args):
    '''
    返回图的最大度，最小度及度分布等信息

    INPUT:

    - ``args[0]=k`` -- 返回第k小度的信息

    - ``args=None`` -- 返回所有度的信息

    例:

    找出最大度及所有最大度点::

        sage: d,VL=G.pick_vertices(-1)
        sage: d,VL
        ...

    找出最小度及所有最小度点::

        sage: d,VL=G.pick_vertices(0)
        sage: d,VL
        ...


    返回所有点的度分布::
        sage: D=G.pick_vertices()
        sage: D
        ...
    '''
    D = dict()
    for v in G.vertex_iterator():
        D.setdefault(G.degree(v), []).append(v)
    keys = D.keys()
    if len(args) > 0:
        order = args[0]
        key = keys[order]
        return key, D[key]
    else:
        return D


Graph.pick_vertices = pick_vertices


# 图G关于顶点v的距离划分
def distance_partition(self, v):
    xvs = self.vertices()
    E = dict()
    for u in xvs:
        dj = self.distance(v, u)
        E.setdefault(dj, []).append(u)
    return E.values()


Graph.distance_partition = distance_partition


# 转化有序对列表为字典， 有序对的第一元为键， 第二元为值。
def list2dict(L):
    D = dict()
    for p in L:
        D.setdefault(p[0], []).append(p[1])
    return D


def corona_graphs(self, H):
    n = self.order()
    G0 = self.add_stars([[v, 1] for v in self.vertices(sort=False)])
    GL = [Graph(1)]*n+[H]*n
    return G0.H_join(GL)


Graph.corona_graphs = corona_graphs


def genreg(n, d, g=None):
    """
    Returns all connected d-regular graphs on n vertices (with girth g), or returns False on
    failure.

    Since every edge is incident to two vertices, n*d must be even.

    INPUT:

    -  ``n`` - number of vertices

    -  ``d`` - degree

    -  ``g`` - girth


    The function depend on M. Meringer's program genreg.
    see http://www.mathe2.uni-bayreuth.de/markus/reggraphs.html


    """
    args = (" %s %s" % (n, d))+(" %s" % g if g != None else "")
    cmd = DOT_SAGE+"genreg %s  -a" % args

    import subprocess
    import re
    try:
        wd = os.environ['TMPDIR']
        s = subprocess.check_output("ls", cwd=wd, stderr=subprocess.PIPE).decode()
        if ".asc" in s:
            temp = subprocess.check_output(
                ["/bin/sh", "-c", "rm *.asc"], cwd=wd, stderr=subprocess.PIPE).decode()
        res = subprocess.call(
            cmd, cwd=wd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        filepath = wd + \
            subprocess.check_output(
                ["/bin/sh", "-c", "ls *.asc"], cwd=wd, stderr=subprocess.PIPE).decode()[:-1]
        txt = open(filepath, "r").read()
        GL = []
        graphs = re.findall(r"Graph \d+:\n([\s\S]*?)Ordnung:", txt)

        for graph_str in graphs:
            # 解析邻接列表
            adj_list, tail_width = re.findall(r"([\s\S]*?)Taillenweite: (\d+)", graph_str)[0]
            # 将邻接列表的文本格式转换为字典格式
            adj_dict = {}
            for line in adj_list.strip().split("\n"):
                node, neighbors = line.split(":")
                node = int(node.strip())
                neighbors = list(map(int, neighbors.strip().split()))
                adj_dict[node] = neighbors
            G = Graph(adj_dict)
            GL.append(G.graph6_string())
        return GL
    except Exception:
        return False


sage.graphs.graph_generators.graphs.genreg = genreg

# =======================图的构造与变换:结束==================================

# 原Cayley图对群来说没有环，需要改定义


def cayley_graph(self, side="right", simple=False, elements=None, generators=None, connecting_set=None):
    r"""
    Return the Cayley graph for this finite semigroup.

    INPUT:

    - ``side`` -- "left", "right", or "twosided":
      the side on which the generators act (default:"right")
    - ``simple`` -- boolean (default:False):
      if True, returns a simple graph (no loops, no labels,
      no multiple edges)
    - ``generators`` -- a list, tuple, or family of elements
      of ``self`` (default: ``self.semigroup_generators()``)
    - ``connecting_set`` -- alias for ``generators``; deprecated
    - ``elements`` -- a list (or iterable) of elements of ``self``

    OUTPUT:

    - :class:`DiGraph`

    EXAMPLES:

    We start with the (right) Cayley graphs of some classical groups::

        sage: D4 = DihedralGroup(4); D4
        Dihedral group of order 8 as a permutation group
        sage: G = D4.cayley_graph()
        sage: show(G, color_by_label=True, edge_labels=True)
        sage: A5 = AlternatingGroup(5); A5
        Alternating group of order 5!/2 as a permutation group
        sage: G = A5.cayley_graph()
        sage: G.show3d(color_by_label=True, edge_size=0.01, edge_size2=0.02, vertex_size=0.03)
        sage: G.show3d(vertex_size=0.03, edge_size=0.01, edge_size2=0.02, vertex_colors={(1,1,1):G.vertices()}, bgcolor=(0,0,0), color_by_label=True, xres=700, yres=700, iterations=200) # long time (less than a minute)
        sage: G.num_edges()
        120

        sage: w = WeylGroup(['A',3])
        sage: d = w.cayley_graph(); d
        Digraph on 24 vertices
        sage: d.show3d(color_by_label=True, edge_size=0.01, vertex_size=0.03)

    Alternative generators may be specified::

        sage: G = A5.cayley_graph(generators=[A5.gens()[0]])
        sage: G.num_edges()
        60
        sage: g=PermutationGroup([(i+1,j+1) for i in range(5) for j in range(5) if j!=i])
        sage: g.cayley_graph(generators=[(1,2),(2,3)])
        Digraph on 120 vertices

    If ``elements`` is specified, then only the subgraph
    induced and those elements is returned. Here we use it to
    display the Cayley graph of the free monoid truncated on
    the elements of length at most 3::

        sage: M = Monoids().example(); M
        An example of a monoid: the free monoid generated by ('a', 'b', 'c', 'd')
        sage: elements = [ M.prod(w) for w in sum((list(Words(M.semigroup_generators(),k)) for k in range(4)),[]) ]
        sage: G = M.cayley_graph(elements = elements)
        sage: G.num_verts(), G.num_edges()
        (85, 84)
        sage: G.show3d(color_by_label=True, edge_size=0.001, vertex_size=0.01)

    We now illustrate the ``side`` and ``simple`` options on
    a semigroup::

        sage: S = FiniteSemigroups().example(alphabet=('a','b'))
        sage: g = S.cayley_graph(simple=True)
        sage: g.vertices()
        ['a', 'ab', 'b', 'ba']
        sage: g.edges()
        [('a', 'ab', None), ('b', 'ba', None)]

    ::

        sage: g = S.cayley_graph(side="left", simple=True)
        sage: g.vertices()
        ['a', 'ab', 'b', 'ba']
        sage: g.edges()
        [('a', 'ba', None), ('ab', 'ba', None), ('b', 'ab', None),
        ('ba', 'ab', None)]

    ::

        sage: g = S.cayley_graph(side="twosided", simple=True)
        sage: g.vertices()
        ['a', 'ab', 'b', 'ba']
        sage: g.edges()
        [('a', 'ab', None), ('a', 'ba', None), ('ab', 'ba', None),
        ('b', 'ab', None), ('b', 'ba', None), ('ba', 'ab', None)]

    ::

        sage: g = S.cayley_graph(side="twosided")
        sage: g.vertices()
        ['a', 'ab', 'b', 'ba']
        sage: g.edges()
        [('a', 'a', (0, 'left')), ('a', 'a', (0, 'right')), ('a', 'ab', (1, 'right')), ('a', 'ba', (1, 'left')), ('ab', 'ab', (0, 'left')), ('ab', 'ab', (0, 'right')), ('ab', 'ab', (1, 'right')), ('ab', 'ba', (1, 'left')), ('b', 'ab', (0, 'left')), ('b', 'b', (1, 'left')), ('b', 'b', (1, 'right')), ('b', 'ba', (0, 'right')), ('ba', 'ab', (0, 'left')), ('ba', 'ba', (0, 'right')), ('ba', 'ba', (1, 'left')), ('ba', 'ba', (1, 'right'))]

    ::

        sage: s1 = SymmetricGroup(1); s = s1.cayley_graph(); s.vertices()
        [()]

    TESTS::

        sage: SymmetricGroup(2).cayley_graph(side="both")
        Traceback (most recent call last):
        ...
        ValueError: option 'side' must be 'left', 'right' or 'twosided'

    .. TODO::

        - Add more options for constructing subgraphs of the
          Cayley graph, handling the standard use cases when
          exploring large/infinite semigroups (a predicate,
          generators of an ideal, a maximal length in term of the
          generators)

        - Specify good default layout/plot/latex options in the graph

        - Generalize to combinatorial modules with module generators / operators

    AUTHORS:

    - Bobby Moretti (2007-08-10)
    - Robert Miller (2008-05-01): editing
    - Nicolas M. Thiery (2008-12): extension to semigroups,
      ``side``, ``simple``, and ``elements`` options, ...
    """
    from sage.graphs.digraph import DiGraph
    from groups import Groups
    if not side in ["left", "right", "twosided"]:
        raise ValueError("option 'side' must be 'left', 'right' or 'twosided'")
    if elements is None:
        assert self.is_finite(), "elements should be specified for infinite semigroups"
        elements = list(self)
    elements_set = set(elements)
    if simple:
        result = DiGraph()
    else:
        result = DiGraph(multiedges=True, loops=True)
    result.add_vertices(elements)

    if connecting_set is not None:
        generators = connecting_set
    if generators is None:
        generators = self.semigroup_generators()
    if isinstance(generators, (list, tuple)):
        generators = dict((self(g), self(g)) for g in generators)
    left = (side == "left" or side == "twosided")
    right = (side == "right" or side == "twosided")

    def add_edge(source, target, label, side_label):
        """
        Skips edges whose targets are not in elements
        Return an appropriate edge given the options
        """
        if target not in elements_set:
            return
        if simple:
            result.add_edge([source, target])
        elif side == "twosided":
            result.add_edge([source, target, (label, side_label)])
        else:
            result.add_edge([source, target, label])
    for x in elements:
        for i in generators.keys():
            if left:
                add_edge(x, generators[i]*x, i, "left")
            if right:
                add_edge(x, x*generators[i], i, "right")
    return result


Semigroups.cayley_graph1 = cayley_graph


def mygraphmul(self, *args):
    """
    Returns the sum of a graph with itself n times.

    EXAMPLES::

        sage: G = graphs.CycleGraph(3)
        sage: H = G*3; H
        Cycle graph disjoint_union Cycle graph disjoint_union Cycle graph: Graph on 9 vertices
        sage: H.vertices()
        [0, 1, 2, 3, 4, 5, 6, 7, 8]
        sage: H = G*1; H
        Cycle graph: Graph on 3 vertices
    """
    n = args[0]
    if isinstance(n, (int, Integer)):
        if n < 0:
            raise TypeError(
                'multiplication of a graph and a negative integer is not defined')
        if n==0:
            return Graph(0)

        if n == 1:
            from copy import copy
            return copy(self)
        return sum([self]*(n-1), self)
    else:
        return reduce(Graph.cartesian_product, flatten([self]+list(args)))


Graph.__mul__ = mygraphmul
Graph.__rmul_ = mygraphmul


def mygraphsub(self, H):
    G = copy(self)
    return G.Delete_Edges(H.edges())


Graph.__sub__ = mygraphsub


def get_Matrix(self, Type):
    if Type == "L":
        return self.laplacian_matrix()
    if Type == "A":
        return self.am()
    if Type == "ABC":
        return self.ABC_matrix()
    if Type == "W":
        self.set_edge_weight()
        return self.weighted_adjacency_matrix()
    if isinstance(Type, type(Graph.order)):
        return self.Weighted_Adjancency_Matrix(Type)
    if Type == "Q":
        return self.signless_laplacian_matrix()
    if Type == "H":
        return self.hermitian_adjacency_matrix()
    if Type == "N":
        return self.normalized_laplace_matrix()
    if Type == "D":
        return self.distance_matrix()
    if isinstance(Type,(tuple,list)):
        if Type[0]=="DA":
            return self.Distance_alpha_matrix(Type[1])
    if Type == "DL":
        return self.distLap()
    if Type == "S":
        return self.Sombor()
    if Type in RR or Type in SR or Type.is_monomial():
        return self.alpha_matrix(Type)




Graph.get_Matrix = get_Matrix
DiGraph.get_Matrix = get_Matrix


def AddWeightedEdges(G, D):
    G=G.Add_Edges(D.keys())
    G.set_edge_weight(D)
    return G
Graph.AddWeightedEdges = AddWeightedEdges




def inertias_of_graph(self, type="A"):
    """
    计算图的各类矩阵的惯性指数
    返回（正惯性指数，负惯性指数，零度)

    例:计算Petersen图的邻接谱
    ::
        sage: G=graphs.PetersenGraph()
        sage: G.Signature(type="A")
        (6, 4, 0)
    """
    M=self.get_Matrix(type)
    coe = map(sgn, M.charpoly().coefficients())
    rk = M.rank()
    sn = [mul(coe[i:i+2]) for i in range(rk-1)]
    p = sn.count(-1)
    n = rk-p
    z = self.order()-rk
    return p, n, z


Graph.inertias_of_graph = inertias_of_graph
DiGraph.inertias_of_graph = inertias_of_graph


def signature(G, type="A"):
    inertias = G.inertias_of_graph(type=type)
    return inertias[0]-inertias[1]


Graph.signature = signature


def is_biregular_bipartite(G):
    if G.is_bipartite():
        G = BipartiteGraph(G)
    else:
        return False
    A = G.reduced_adjacency_matrix()
    return Set(list(add(A))).cardinality()+Set(list(add(A.T))).cardinality() == 2


Graph.is_biregular_bipartite = is_biregular_bipartite


def is_WF(G):
    
    return bool(mul([G.subgraph(G[v]).is_forest() for v in G]))


Graph.is_WF = is_WF


def fan_graph(n):
    G = Graph(1).join(graphs.PathGraph(n-1))
    G.relabel()
    return G


sage.graphs.graph_generators.graphs.fan_graph = fan_graph


def rose_graph(D):
    R = Graph(1)
    for d in D:
        n = R.order()
        R.add_cycle([0]+range(n, n+d-1))
    return R


sage.graphs.graph_generators.graphs.rose_graph = rose_graph


def random_oriente(G):
    return DiGraph([e[::(-1) ^ randint(0, 1)] for e in G.edges(labels=False)])


Graph.random_oriente = random_oriente


# 对矩阵按行型列型进行分块
def subdivide(N, rBk, cBk=None, bksize=True):
    M = copy(N)
    if cBk == None:
        cBk = rBk
    if bksize:
        rsd = Composition(rBk).partial_sums()[:-1]
        csd = Composition(cBk).partial_sums()[:-1]
        M.subdivide(rsd, csd)
    else:
        M.subdivide(rsd, csd)
    return M


def random_symmetric_matrix(n, x=0, y=1):
    B = matrix(n, lambda i, j: (i > j)*randint(x, y))
    return diagonal_matrix([randint(x, y) for _ in range(n)])+B+B.T


def subdivide1(N, rBk, cBk=None, bksize=True, quo=False, sym=False):
    M = copy(N)
    if cBk == None:
        cBk = rBk
    if bksize:
        rsd = Composition(rBk).partial_sums()[:-1]
        csd = Composition(cBk).partial_sums()[:-1]
        M.subdivide(rsd, csd)
        if quo:
            nrb = len(rBk)
            ncb = len(cBk)
            Q = matrix(nrb, ncb, lambda i, j: add(M.subdivision(i, j).list()))
            K = diagonal_matrix([sqrt(rb) for rb in cBk])
            if sym:
                return ~K ^ 2*Q
            return ~K*Q*~K
        else:
            return M
    else:
        nrb = len(rBk)
        ncb = len(cBk)
        B = block_matrix(len(rBk), len(cBk), [
                         M[r, c] for r in rBk for c in cBk])
        if quo:
            Q = matrix(nrb, ncb, lambda i, j: add(B.subdivision(i, j).list()))
            K = diagonal_matrix([sqrt(len(rb)) for rb in cBk])
            if sym:
                return ~K*Q*~K
            else:
                return ~K ^ 2*Q
        else:
            return B
    return M


def K_subdivide_graph(self, k, relabel=False):
    E = self.edges(labels=False)
    G = Graph(self.order())
    for i in range(len(E)):
        G.add_path([E[i][0]]+[(i, j) for j in range(k)]+[E[i][1]])
    if relabel:
        G.relabel()
    return G


Graph.K_subdivide_graph = K_subdivide_graph


def Kronecker_Sum(A, B):
    D = [A.dimensions(), B.dimensions()]
    for i in D:
        if i[0] != i[1]:
            print("请输入方阵")
            return
    return A.tensor_product(identity_matrix(D[1][0]))+(identity_matrix(D[0][0])).tensor_product(B)


# 开花，为构造广义线图
def add_blossom(self, d):
    H = self.copy()
    H.allow_multiple_edges(True)
    ne = []
    n = self.order()
    for i in d:
        for j in [1..d[i]]:
            ne = ne+[[i, n]]
            n = n+1
    H.add_edges(2*ne)
    return H


Graph.add_blossom = add_blossom


def generalized_line_graph(self):
    M = self.incidence_matrix()
    B = matrix(GF(2), M)
    return Graph(B.T*B)


Graph.generalized_line_graph = generalized_line_graph


def gline_graph_forbidden_subgraphs():
    sglfg = ['DIk', 'DFw', 'DNw', 'E?Fg', 'E?Bw', 'E?Fw', r'EC\w', 'E?^w', 'E?~w', 'E@~w', 'EF~w', 'E@QW', 'EHQW', 'E@`w', 'E@]o',
             'E@Rw', 'EAMw', 'E`]o', 'E@^W', 'E@Vw', 'EK]w', 'EBjw', 'EK\\w', 'E@^w', 'ELrw', 'EBnw', 'EJ]w', 'EK~w', 'EJnw', 'EL~w', 'EN~w']
    return [Graph(s) for s in sglfg]


sage.graphs.graph_generators.graphs.gline_graph_forbidden_subgraphs = gline_graph_forbidden_subgraphs


def is_gline_graph(g, certificate=False):
    from sage.graphs.graph import Graph
    g._scream_if_not_simple()
    for fg in sage.graphs.graph_generators.graphs.gline_graph_forbidden_subgraphs():
        h = g.subgraph_search(fg, induced=True)
        if h is not None:
            if certificate:
                return (False, h)
            else:
                return False
    if not certificate:
        return True
    if g.is_connected():
        R, isom = root_graph(g)
    else:
        R = Graph()
    for gg in g.connected_components_subgraphs():
        RR, _ = root_graph(gg)
        R = R + RR
    _, isom = g.is_isomorphic(R.line_graph(labels=False), certificate=True)
    return (True, R, isom)


Graph.is_gline_graph = is_gline_graph


def groupbyinv(GL, inv):
    u"""
    利用不变量inv(一个函数)对列表类对象GL进行分组;
    可以用来给定图类中的同谱图。
    T=graphs.trees(11)
    inv=lambda G: G.signless_laplacian_matrix().charpoly()
    res=groupbyinv(T,inv)
    graphs_list.show_graphs(flatten(res[1]))
    """
    groups = {}
    for G in GL:
        key = inv(G)
        if key in groups:
            groups[key].append(G)
        else:
            groups[key] = [G]
    GL = []
    for t in groups:
        if len(groups[t]) > 1:
            GL.append(groups[t])
    return groups, GL


def All_Paths(self, start, end, length="all",induced=False):
    if self.is_directed():
        iterator = self.neighbor_out_iterator
    else:
        iterator = self.neighbor_iterator
    if start == end:
        return [[start]]
    all_paths = []      # list of
    act_path = []       # the current path
    act_path_iter = []  # the neighbor/successor-iterators of the current path
    done = False
    s = start
    while not done:
        if s == end:      # if path completes, add to list
            if length == "all":
                all_paths.append(act_path+[s])
            else:
                if len(act_path+[s]) == (length+1):
                    if induced:
                        if (self.subgraph(act_path+[s])).size()==length:
                            all_paths.append(act_path+[s])
                    else: 
                        all_paths.append(act_path+[s])

        else:
            if s not in act_path:   # we want vertices just once in a path
                act_path.append(s)  # extend current path
                # save the state of the neighbor/successor-iterator of the current vertex
                act_path_iter.append(iterator(s))
        s = None
        while (s is None) and not done:
            try:
                # try to get the next neighbor/successor, ...
                s = next(act_path_iter[-1])
            except (StopIteration):         # ... if there is none ...
                act_path.pop()              # ... go one step back
                act_path_iter.pop()
            if len(act_path) == 0:            # there is no other vertex ...
                done = True                 # ... so we are done
    return all_paths


Graph.All_Paths = All_Paths


def orbit_subdivision(H):
    """
    返回图的轨道划分
    """
    GP = H.automorphism_group()
    return list(map(list,GP.orbits()))


Graph.orbit_subdivision = orbit_subdivision








#将图G的点子集U中的每个点用k*P_1代替
def Blowup(G,U,k):
    V1=Set(G.vertices()).intersection(Set(U))
    if isinstance(k,int):
        VP={v:Graph(k) for v in V1}
    if isinstance(k,(tuple,list)) and len(k)>=len(V1):
        VP={V1[i]:Graph(k[i]) for i in range(len(V1))}
    return G.H_join(VP,relabel=True)
Graph.Blowup=Blowup






def layout_Tree(self, tree_orientation="down", tree_root=None):
    n = self.order()
    vertices = self.vertices()
    if tree_root is None:
        root = self.center()[0]
    else:
        root = tree_root

    pos = {}

    # The children and parent of each vertex
    children = {root: sorted(self.neighbors(root))}
    parent = {u: root for u in children[root]}

    stack = [list(children[root])]
    stick = [root]

    obstruction = [0.0]*self.num_verts()

    if tree_orientation in ['down', 'left']:
        o = -1
    elif tree_orientation in ['up', 'right']:
        o = 1
    else:
        raise ValueError(
            'orientation should be "up", "down", "left" or "right"')

    def slide(v, dx):
        """
        shift the vertex v and its descendants to the right by dx

        Precondition: v and its descendents have already had their
        positions computed.
        """
        level = [v]
        while level:
            nextlevel = []
            for u in level:
                x, y = pos[u]
                x += dx
                obstruction[y] = max(x+1, obstruction[y])
                pos[u] = x, y
                nextlevel += children[u]

            level = nextlevel

    while stack:
        C = stack[-1]

        # If all the children of stick[-1] have been given a position
        if len(C) == 0:
            p = stick.pop()
            stack.pop()
            cp = children[p]
            y = o*len(stack)

            if len(cp) == 0:
                # If p has no children, we draw it at the leftmost position
                # which has not been forbidden
                x = obstruction[y]
                pos[p] = x, y
            else:
                # If p has children, we put v on a vertical line going
                # through the barycenter of its children
                x = sum([pos[c][0] for c in cp])/len(cp)
                pos[p] = x, y
                ox = obstruction[y]
                if x < ox:
                    slide(p, ox-x)
                    x = ox

            # If the vertex to the right of p has not children, we want it
            # at distance 1 from p
            obstruction[y] = x+1

        # Otherwise, we take one of the children and add it to the
        # stack. Note that this vertex is removed from the list C.
        else:
            t = C.pop()

            pt = parent[t]
            ct = [u for u in sorted(self.neighbors(t)) if u != pt]
            children[t] = ct

            for c in ct:
                parent[c] = t

            stack.append([c for c in ct])
            stick.append(t)

    if tree_orientation in ['right', 'left']:
        return {p: (py, px) for p, (px, py) in pos.iteritems()}

    return pos


Graph.layout_Tree = layout_Tree

# def GRelabel(G):
#     H=copy(G)
#     H.relabel()
#     return H
# Graph.Relabel=GRelabel


# ==================================下面几个函数是对图中的点及标号进行控制=================================

def ShiftPos(G,V=None, shift=[0,0], scale=[1, 1]):
    """对图的局部进行移动和放缩      

    Args:
        G ([Graph, HyperG]): 要移动的图。
        V ([list], optional): [要移动的局部顶点集]. Defaults to None.
        shift (list, optional): [偏移量]. Defaults to [0,0].
        scale (list, optional): [横向放缩比例， 纵向放出比例]. Defaults to [1, 1].
    """
    if "pos_matrix" not in dir(G):
        G.Pos()
    if V==None:
        V=G.vertices()
    nPos = G.pos_matrix+ones_matrix(G.order(), 1)*matrix(shift)
    nPos = nPos*diagonal_matrix(scale)
    G._pos.update({v: nPos[v] for v in V})
    G.VerticesLabel()
    G.LabelShift()
    G._nodelist = ["$%s$" % t for t in G.vertices()]


Graph.ShiftPos = ShiftPos
DiGraph.ShiftPos = ShiftPos


def Pos(G, pos=None):
    """初始化图G位置相关内部变量

    Args:
        G ([Graph]): [description]
        pos ([dict]], optional): [顶点位置词典]]. Defaults to None.

    Returns:
        [type]: [description]
    """    
    POS = {}
    if pos == None and G._pos ==None:
        G._pos = G.layout(layout="spring",save_pos=True)
    elif (pos != None):
        G._pos = pos
    pos = G._pos
    G.pos_matrix = matrix(RR, map(tuple, G._pos.values()))
    for i in pos:
        POS[i] = vector(RR, pos[i])
    G._POS = POS
    G._pos = POS
    return G._POS


Graph.Pos = Pos
DiGraph.Pos = Pos


def VerticesLabel(G, vlabel=None, name="v"):
    if vlabel == None:
        vlabel = {v: r"$v_{%d}$" % (v) for v in G.vertices()}
    G._VerticesLabel = vlabel
    return vlabel


Graph.VerticesLabel = VerticesLabel


def Init(G):
    """对图顶点位置和标号等信息进行重置

    Args:
        G ([Graph, HyperG]): 要处理的图
    """
    G.Pos(G.layout("spring"))
    G.VerticesLabel()
    G.LabelShift()
    G._nodelist = ["$%s$" % t for t in G.vertices()]


Graph.Init = Init




def local_spring(H, Pos):
    import networkx as nx
    G = nx.Graph()
    G.add_edges_from(H.edges(labels=False))
    pos = nx.spring_layout(G, pos=Pos, fixed=Pos.keys(), k=0.36, iterations=50)
    H._pos = pos


Graph.local_spring = local_spring

# def Sym_Vertices(H):
#     V = H.vertices()
#     if isinstance(H, Graph):
#         pos = H.Pos()
#     else:
#         pos = H._pos
#     absPos = {v: vector(map(abs, pos[v])) for v in pos}
#     D = matrix(H.order(), lambda i, j: (absPos[i]-absPos[j]).norm())
#     for d in sorted([v for v in sorted(Set(D.list())) if 0 < v < 0.5], reverse=True):
#         eqivG = Graph(matrix(H.order(), lambda i, j: (
#             absPos[i]-absPos[j]).norm() < d and i != j))
#         Bks = [bk for bk in eqivG.connected_components()]
#         if max(map(len, Bks)) == 2:
#             break
#     Bks = [bk for bk in Bks if len(bk) == 2]
#     for bk in Bks:
#         u, v = bk
#         pos[u][0] = -pos[v][0]
#         pos[u][1] = pos[v][1]

#     H._pos = pos
#     if isinstance(H, Graph):
#         return H.plot()
#     else:
#         return H.plot(redraw=True)


# Hypergraph.Sym_Vertices = Sym_Vertices
# Graph.Sym_Vertices = Sym_Vertices

# 将图H 的点集组 Vp的两个子集的中心旋转至水平,且两中心的连线中点为原点

def Rotate_Shift(H, Vp=None, degree=0, orgin=0, ndigits=3):
    """
    Rotate and shift the vertices of a graph or hypergraph.

    This function rotates a graph or hypergraph around an origin by a given angle and then shifts the vertices so that
    the specified origin is at the zero vector. This function works both for objects of the `Graph` and `HyperG` class.

    Parameters:
        H (Graph or HyperG): The graph or hypergraph to rotate.
        Vp (list): The pivot vertices for rotation.
        degree (int, optional): The angle of rotation in degrees. Default is 0.
        orgin (int, optional): The vertex to set to the zero vector after rotation. Default is 0.
        ndigits (int, optional): The number of decimal places for rounding the new positions. Default is 3.

    Returns:
        A plot of the rotated graph or hypergraph.

    Example:
        >>> G = Graph()
        >>> G.Rotate_Shift(Vp=[0, 1], degree=45)
    """

    # function implementation here    
    a = degree*pi/180
    R = matrix(RR, 2, 2, [cos(a), sin(a), -sin(a), cos(a)])
    V = H.vertices()
 
    if isinstance(H, Graph):
        pos = H.Pos()
    else:
        pos = H._pos
    if Vp!=None:
        vc = operator.sub(*[Mean([vector(pos[v])
                                  for v in flatten([P])]) for P in Vp])
        ct = Mean([pos[i] for i in V])
        nvc = vc.normalized()
        T = nvc[0]*identity_matrix(2)+nvc[1]*matrix(2, [0, 1, -1, 0])
        npos = {i: R*T*(pos[i]-ct)+ct for i in pos}
    else:
        npos=pos
    opos = npos[orgin]
    H._pos = {i: vector(map(lambda p: round(p, ndigits),
                            npos[i]-opos)) for i in npos}
    if isinstance(H, Graph):
        return H.plot()
    else:
        return H.plot(redraw=True)


# HyperG.Rotate_Shift = Rotate_Shift
Graph.Rotate_Shift = Rotate_Shift


def layout_grid(G):
    '''
    按纵坐标对点集进行分类,
    然后将它们的位置设为整点
    '''
    eqivG = Graph(matrix(G.order(), lambda i, j: abs(
        (G._pos[i]-G._pos[j])[1]) < 0.05 and i != j))
    Bks = eqivG.connected_components()

    L = [sorted(bk, key=lambda x:G._pos[x][::-1])
         for bk in Bks]
    L = sorted(L, key=lambda x:G._pos[x[0]][1])
    for i in range(len(L)):
        for j in range(len(L[i])):
            G._pos[L[i][j]] = vector([j, i])


Graph.layout_grid = layout_grid



# #对图G的顶点子集 VL 进行水平或竖直均匀分布
# def uniform_distribution(G,VL=None,Y=True):
#     D=G._pos
#     if Y:
#         cn=1
#     else:
#         cn=0
#     if VL==None:
#         VL=G.vertices()
#     YLst=sorted(Set(list(matrix(D.values()).columns()[cn])))
#     yn=(len(YLst)-1)
#     yh=YLst[-1]-YLst[0]
#     DY=dict(zip(YLst,[-yh/2*i+(1-i)*yh/2 for i in [0,1/yn,..,1]]))
#     if Y:
#         G._pos.update({k:vector([D[k][0],round(DY[D[k][1]],3)]) for k in VL})
#     else:
#         G._pos.update({k:vector([round(DY[D[k][0]],3),D[k][1],]) for k in VL})
#     return G
# Graph.uniform_distribution=uniform_distribution

# #使图的点子集VL的点水平或竖直吸附格线
# def SnapXY(H,VL=None,X=True,markunit=0.1):
#     G=copy(H)
#     D=G._pos
#     if VL==None:
#         VL=G.vertices()
#     if X:
#         G._pos={k:vector([round(D[k][0]/markunit)*markunit,D[k][1]]) for k in VL}
#     else:
#         G._pos={k:vector([D[k][0],round(D[k][1]/markunit)*markunit]) for k in VL}
#     return G
# Graph.SnapXY=SnapXY


#使得图的点集的中心或知道点 c 为坐标系原点
Mean=lambda L:add(L)/len(L)

def CenterPos(G,c=None):
    H=copy(G)
    ndigits=3
    c=None
    V = H.vertices()
    if isinstance(H, Graph):
        pos = H.Pos()
    else:
        pos = H._pos
    if c==None:
        center = Mean([pos[k] for k in pos])
    else:
        center=pos[c]
    pos = {k: pos[k]-center for k in pos}
    pos = {i: vector(map(lambda p: round(p, ndigits), pos[i])) for i in pos}
    return H


Graph.CenterPos = CenterPos


def LabelShift(G,vshift=None):
    """对图的顶点的标号位置进行平移。

    Args:
        G (Graph, HyperG): 待处理的图。
        vshift (dict, optional): 要移动顶点的标号的位置平移量字典. Defaults to None.

    Returns:
        [dict]: 各顶点标号位置字典
    """
    vpos=G.Pos()
    if vshift==None:
        vshift={v:vector([0,0]) for v in G.vertices()}
    if type(vshift)!=dict:
        vshift={v:vector([0.1,0.1]) for v in G.vertices()}
    G._vLabelPos={v:vpos[v]+vshift[v] for v in G.vertices()}
    return G._vLabelPos
Graph.LabelShift=LabelShift

# # 改变部分点集标号的位置
def vSet_Shift(G,U,shift):
    vpos=G.Pos()
    for v in U:
        G._vLabelPos[v]=vpos[v]+vector(shift)
Graph.vSet_Shift=vSet_Shift


def posTransformation(G,M=identity_matrix(2),v=vector([0,0]),U={}):
    V=G.vertices()
    U=Set(V) if U=={} else Set(V)&Set(U)
    subM=G.pos_matrix[[i for i in U],:]
    nM=(subM-matrix([v]*len(U)))*M
    G.pos_matrix[[i for i in U],:]=nM
    G._pos={v:G.pos_matrix[v] for v in V}
    G._POS=G._pos
    G.LabelShift()

Graph.posTransformation=posTransformation


# #自动调整顶点位置
def NormalizedPos(G, *args):
    #     G.Init()
    degree = 0
    VL = G.vertices()
    if len(args) == 1:
        VL = args[0]
    if len(args) == 2:
        degree = args[1]

    var('a, b, x')
    model(x) = a * x + b
    pos = G.Pos()
    sol = find_fit([pos[v] for v in VL], model)
    alpha = arctan(a.subs(sol[0]))+pi/2+degree
    R = matrix(2, 2, [cos(alpha), -sin(alpha), sin(alpha), cos(alpha)])
    ct = add(G._POS.values())/G.order()
    G.posTransformation(M=R, v=ct)
    vpos = G._pos
#     G._pos={k:matrix_round(vpos[k],ndigits=1,eps=1e-2) for k in vpos.keys()}
    G._pos = {k: vpos[k].apply_map(lambda x: round(x, 2)) for k in vpos.keys()}
    return G


Graph.NormalizedPos = NormalizedPos


# def _Normalize(H, Vp, ndigits=3):
#     V = H.vertices()
#     if isinstance(H, Graph):
#         pos = H.Pos()
#     else:
#         pos = H._pos
#     opos = Vp[0]
#     ct = Vp[0]
#     vt = Vp[1]-ct
#     x, y = vt.normalized()
#     R = matrix(RR, 2, 2, [x, y, -y, x])
#     npos = {i: R*(pos[i]-ct).apply_map(lambda p: round(p, 3)) for i in pos}
#     H._pos = npos
#     if isinstance(H, Graph):
#         return H.plot()
#     else:
#         return H.plot(redraw=True)


# Hypergraph._Normalize = _Normalize


# def Scale(H, ratio=1):
#     H._pos = {k: ratio*H._pos[k] for k in H._pos}
#     return H


# Hypergraph.Scale = Scale

Plot_Option = copy(sage.graphs.graph_plot.DEFAULT_PLOT_OPTIONS)
options = {
    'vertex_size': 20,
    'max_vertex_num_plot': 50,
    'vertex_labels': True,
    'vertex_color': "red",
    'edge_style': 'solid',
    "color_by_label": {-1: "green", 1: "magenta"},
    "edge_colors": {},
    'labelshift': (0, 0.05),
    'title': '$G$',
    'max_vertex_num_plot': 50,
    'save_pos': True,

}
Plot_Option.update(options)
Plot_Option['max_vertex_num_plot'] = 50





def GraphPlot(G, **options):
    edge_colors = options.pop('edge_colors', dict())
    color_by_label = options.pop('color_by_label', False)
    if color_by_label == True:
        colorbylabels = G._color_by_label()
    elif color_by_label != False:
        colorbylabels = G._color_by_label(color_by_label)
    else:
        colorbylabels = dict()

    nec = {k: map(lambda t: tuple(sorted(t)), [
                  it for it in v if G.has_edge(it)]) for k, v in edge_colors.items()}
    cled = flatten(nec.values(), max_level=1)
    nec1 = {k: map(lambda t: list(t)+[G.edge_label(*t)], v)
            for k, v in nec.items()}
    ncbylb = {k: [t for t in v if t[:2] not in cled]
              for k, v in colorbylabels.items()}
    cmdkeys = Set(nec1.keys()) & Set(ncbylb.keys())
    merg_eedge_color = dict()
    merg_eedge_color.update({k: nec1[k]+ncbylb[k] for k in cmdkeys})
    merg_eedge_color.update({k: nec1[k]
                             for k in nec1.keys() if k not in cmdkeys})
    merg_eedge_color.update({k: ncbylb[k]
                             for k in ncbylb.keys() if k not in cmdkeys})
    options.update({'edge_colors': merg_eedge_color})
    labelshift = options.pop('labelshift', (0, 0))
    G.vSet_Shift(range(G.order()), labelshift)
    gp = G.graphplot(**options)
    gp._plot_components.update({'others': []})
    gp._plot_components['vertex_labels'] = []
    vPos = G._POS
    vLb = G._nodelist
    vLp = G._vLabelPos
    for v in G.vertices():
        gp._plot_components['vertex_labels'].append(
            text(vLb[v], vLp[v], fontsize=12, color="blue", zorder=8))
    return gp


Graph.GraphPlot = GraphPlot
BipartiteGraph.GraphPlot = GraphPlot







def Plot(G,**options):
    return G.GraphPlot(**options).plot()
Graph.Plot=Plot
BipartiteGraph.Plot=Plot


def plot_graph_with_labels(G, *args, **kwargs):
    E1 = kwargs.pop('edge_colors', dict())
    figsize = kwargs.pop('figsize', [3, 3])
    egv = kwargs.pop('eigenvector', ("A", 0))
    scale=kwargs.pop('scale',1)
    vertexsize = kwargs.pop('vertex_size', 20)
    vertex_labels = kwargs.pop('vertex_labels', True)
    vertex_color = kwargs.pop('vertex_color', "red")
    color_by_label = kwargs.pop('color_by_label', True)
    edge_colors = kwargs.pop('edge_colors', dict())
    save_pos = kwargs.pop('save_pos', True)
    title = kwargs.pop('title', G.name())
    E = G.edges(labels=False)
    G.Init()
    if len(args) == 0:
        r, v1 = G.eigenpairs(type=egv[0], order=egv[1])
        v1=scale*v1
        G._nodelist = [("$v_{%s}$:%.2f" % t).rstrip('0').rstrip('.')
                       for t in list(enumerate(map(lambda x:round(x+0.0, 2), v1)))]
    else:
        vL = args[0]
        if vL == []:
            G._nodelist = ["$v_{%s}$" % t for t in G.vertices()]
        else:
            G._nodelist = vL
    options = {
        'figsize': figsize,
        'vertex_size': vertexsize,
        'vertex_labels': vertex_labels,
        'vertex_color': vertex_color,
        "color_by_label": color_by_label,
        "edge_colors": edge_colors,
        'labelshift': (0, 0.1),
        'title': title,
        'save_pos': save_pos
    }
    G = G.NormalizedPos()
    return G.Plot(**options)


Graph.plot_graph_with_labels = plot_graph_with_labels


def Set_Label_with_Eigenvector(G, *args, **kwargs):
    egv = kwargs.pop('type', ("A", 0))
    vertex_labels = kwargs.pop('vertex_labels', True)
    vertex_color = kwargs.pop('vertex_color', "red")
    save_pos = kwargs.pop('save_pos', True)
    title = kwargs.pop('title', G.name())
    E = G.edges(labels=False)
    G.Init()
    if len(args) == 0:
        r, v1 = G.eigenpairs(type=egv[0], order=egv[1])
        G._nodelist = [("$v_{%s}:%.2f" % t).rstrip('0').rstrip(
            '.')+"$" for t in list(enumerate(map(lambda x:round(x+0.0, 2), v1)))]
    else:
        vL = args[0]
        if vL == []:
            G._nodelist = ["$v_{%s}$" % t for t in G.vertices()]
        else:
            G._nodelist = vL
    G = G.NormalizedPos()
    return G


Graph.Set_Label_with_Eigenvector = Set_Label_with_Eigenvector


def longest_cycle(G):
    l=2
    maxH=None
    for e in G.edges():
        H=G.longest_path(*e[:2]).Add_Edges([e])
        if H.order()>l:
            maxH=H
    return H
Graph.longest_cycle=longest_cycle

def Find_all_hamiltonian_cycles(H):
    G = copy(H)
    v = G.vertices()[0]
    n = G.order()
    N = G[v]
    Cycles = []
    for u in N:
        cycle = G.all_paths(u, v)
        Cycles = Cycles+[c for c in cycle if len(c) == n]
        G.delete_edge(u, v)
    return Cycles


Graph.Find_all_hamiltonian_cycles = Find_all_hamiltonian_cycles


def Find_hamiltonian_cycle(H):
    G = copy(H)
    v = G.vertices()[0]
    n = G.order()
    N = G[v]
    Cycles = []
    for u in N:
        cycle = G.all_paths(u, v)
        for c in cycle:
            if len(c) == n:
                return c
        G.delete_edge(u, v)


Graph.Find_hamiltonian_cycle = Find_hamiltonian_cycle


def Find_Max_cycle(H):
    """找出图H的最长圈

    Args:
        H (Graph):图H

    Returns:
        Graph: 子图圈
    """    
    MC = [H.vertices()[0]]
    G = copy(H)
    while G.order()>max(2,len(MC)):
        v = G.vertices()[0]
        N = G[v]
        for u in N:
            cycle = G.all_paths(u, v)
            for c in cycle:
                if len(c) > len(MC):
                    MC= c
            G.delete_edge(u, v)
        G=G.Delete_Vertices([v])
    H=Graph(0)
    H.add_cycle(MC)
    return H
Graph.Find_Max_cycle = Find_Max_cycle



def plotbycycle(G):
    G = graphs.IcosahedralGraph()
    C = G.Find_hamiltonian_cycle()
    P = dict(zip(C, G.vertices()))
    G.relabel(P)
    return G.plot(layout="circular")


Graph.plotbycycle = plotbycycle


def Latex(G, Lb=True):
    """Generate the Latex Source code of Graph G.
    G {Graph} -- [Graph]
    """
    if Lb:
        def VLB(s): return "$v_{%s}$" % (s)
    else:
        def VLB(s): return s
    pos = G._pos

    G._nodelist = [VLB(str(t)) for t in G.vertices()]
    labels = G._nodelist
    V = G.vertices()
    D = {V[i]: i for i in range(G.order())}
    pos = G.Pos()
    if not '_nodelist' in dir(G):
        G._nodelist = ["$v_{%s}$" % t for t in range(G.order())]
    labels = G._nodelist
    edges = G.edges(labels=False)
    tex = r'''
\begin{tikzpicture}[scale=3,
Node/.style={fill,circle,scale=.5}
]'''
    for v in G.vertices():
        tex += r'''
\node[Node,label={{[label distance=0pt]0: {2}}}] ({0})  at {1} {{}};'''.format(v, (tuple(map(lambda v: round(v, 2), pos[v]))), labels[D[v]])
    for u, v in edges:
        tex += r"""
\draw[color=black]({0})--({1});""".format(u, v)
    tex += r'''
\end{tikzpicture}
'''
    return tex


Graph.Latex = Latex


def Graph2TeX(G):
    pos = G.Pos()
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
\caption{SageMath生成图}
\end{figure}

\end{document}
'''
    f = open(tex_file, "w")
    f.write(tex)
    f.close()
    os.system('open "%s"' % (tex_file))


Graph.Graph2TeX = Graph2TeX
# Hypergraph.Graph2TeX = Graph2TeX


def estrada_index(G):
    sp = eigenpair(G.weighted_adjacency_matrix(), spe=True)
    return add(map(exp, map(RR, sp)))


Graph.estrada_index = estrada_index


def AMdb(G):
    return matrix(G.order(), lambda i, j: (1/G.degree(i)+1/G.degree(j))*(i in G[j]))



Graph.AMdb = AMdb

def Sombor(G,k=2):
    return matrix(G.order(), lambda i, j: (G.degree(i)^k+G.degree(j)^k)*(i in G[j]))

Graph.Sombor=Sombor


def findDiCycles(H, *args):
    B, C = H.blocks_and_cut_vertices()
    maxlenC = max(map(len, B))
    L = []
    for lenc in [2..maxlenC]:
        L = L+list(H.subgraph_search_iterator(digraphs.Circuit(lenc)))
    Ecs = []
    for c in L:
        lenc = len(c)
        Ecs.append(Set([tuple([c[i], c[(i+1) % lenc]]) for i in range(lenc)]))
    Cycles = [DiGraph(map(tuple, s)) for s in Set(Ecs)]
    if len(args) > 0:
        sv = args[0]
        if isinstance(sv, Integer):
            sv = [sv]
        return [c for c in Cycles if Set(c.vertices()).issuperset(Set(sv))]
    return Cycles


DiGraph.findCycles = findDiCycles

# H.findCycles(VL) 寻找图H中包含点子集VL的所有圈
def findCycles(H, incVs=None, l=True):
    """
    查找输入的图 H 的所有圈。
    incVs：指定的一组顶点，只返回包含这些顶点的圈
    l：是否返回圈的深度优先遍历
    返回一个列表，包含所有符合条件的圈的深度优先遍历
    """

    # 获取图 H 的块和割点
    B, C = H.blocks_and_cut_vertices()
    # 找到图 H 的所有块中最大的圈的长度
    maxlenC = max(map(len, B))
    # 初始化一个列表 L，用于存储所有找到的圈
    L = []
    # 对所有可能的圈长进行循，从 3 到最大圈长
    for lenc in [3..maxlenC]:
        # 使用 subgraph_search_iterator 方法找到所有长度为 lenc 的圈
        # 并将它们添加到列表 L 中
        L = L+list(H.subgraph_search_iterator(graphs.CycleGraph(lenc)))
    # 初始化一个列表 Gs，用于存储所有找到的圈的信息
    Gs=[]
    # 初始化一个列表 unique_L，用于存储所有不同的圈
    unique_L=[]
    for G in L:
        # 如果 G 不同于列表 unique_L 中的任何一个圈
        if not any(G == H for H in unique_L):
            # 将 G 添加到列表 unique_L 中
            unique_L.append(G)
            # 将 G 和它的深度优先遍历添加到列表 Gs 中
            Gs.append([G,list(G.depth_first_search(G.vertices()[0]))])
 
    # 如果输入的 incVs 是一个整数，将它转换成只包含该整数的集合
    # 否则，将它转换成集合（如果它还不是集合的话）
    if isinstance(incVs, Integer):
        incVs=Set([incVs]) 
    else:
        incVs=Set(incVs) 
    # 返回所有包含 incVs 中所有顶点的圈的深度优先遍历的列表
    # 如果 l 为 True，则返回深度优先遍历，否则返回空列表
    return [gs[1*l] for gs in Gs if Set(gs[0].vertices()).issuperset(Set(incVs))]

# 将函数作为 Graph 类的方法添加
Graph.findCycles = findCycles

def gen_simple_cycles(H):
    # 将图 H 中的所有边表示为一个列表 E，并计算 H 的大小 m。
    E = list(H.edges(labels=False))
    m = H.size()
    # 创建一个 m 维二元有限域上的向量空间 R。
    R = FiniteField(2)^m
    # 使用 cycle_basis() 方法计算 H 的一个基圈组 Bs，
    # 并将每个圈表示为一个 0-1 向量，存储在列表 CBV 中。
    Bs = [[E.index(tuple(sorted(e[:2]))) for e in c] for c in H.cycle_basis(output='edge')]
    CBV = [R([e in B for e in range(m)]) for B in Bs]
    # 使用 subspace_with_basis() 方法生成 CBV 中向量构成的向量空间 SubR。
    SubR = R.subspace_with_basis(CBV)
    # 遍历 SubR 中除了零向量以外的所有向量 el，对应每个向量生成一个简单圈的子图 G。
    for el in SubR[1:]:
        G = Graph([E[i] for i in el.support()], format='list_of_edges')
        # 检查 G 是否是欧拉图，并且每个顶点的度数是否为 2。
        if G.is_eulerian() and max(G.degree(), default=0) == 2:
            # 如果满足条件，则输出 G 的欧拉回路中除了起点和终点以外的所有顶点，
            # 即生成的简单圈。
            yield G.eulerian_circuit(return_vertices=True)[1]
Graph.findCycles1 = gen_simple_cycles

def findCycles2(G, k=0, induced=False):
    # Find all simple cycles in the directed graph
    all_cycles = [c for c in G.to_directed().all_simple_cycles() if len(c) >= 4 and c[1] < c[-2]]

    if k == 0:
        return all_cycles
    else:
        # Filter cycles based on length
        k_cycles = [c for c in all_cycles if len(c) == k + 1]

        if not induced:
            return k_cycles
        else:
            # Filter cycles based on being induced
            return [c for c in k_cycles if G.subgraph(c[:-1]).size() == k]

# Add the function as a method to the Graph class
Graph.findCycles = findCycles2


# def findCycles(G):
#     Cs = []
#     E = G.edges(labels=False)
#     H = copy(G)
#     for e in E:
#         H = H.Delete_Edges([e])
#         Cs = Cs+H.all_paths(*e)
#     return map(Cycle, Cs)


# Graph.findCycles = findCycles


# 取矩阵A的余子阵
def CoSubmatrix(A, *args):
    R, C = [], []
    if len(args) == 0:
        return A
    if len(args) == 1:
        R, C = [args[0]]*2
    if len(args) > 1:
        R, C = args[:2]
    m, n = A.dimensions()
    coR = [i for i in range(m) if i not in R]
    coC = [j for j in range(n) if j not in C]
    return A[coR, coC]

# def plot_graph_with_vertexlabels(G, **kwargs):
#     E1 = kwargs.pop('edge_colors', dict())
#     E=G.edges()
#     E3=dict([])
#     for key in E1:
#         E3.update({key:[tuple(list(e)+[G.edge_label(*e[:2])]) for e in E1[key]]})
#     E0=[e for e in E if e[2]<0]
#     E00=[e[:2] for e in E0]+map(lambda t:tuple(sorted(t)), flatten(E1.values(),max_level=1))
#     E2=[e for e in E if e[:2] not in E00]
#     E3.update({"green":E0,"black":E2})
#     lsft=G.LabelShift()
#     G.vSet_Shift(range(G.order()),[0,0.05])
#     vLb=G.VerticesLabel([("$v_{%s}:%.2f"%t).rstrip('0').rstrip('.')+"$" for t in list(enumerate(map(lambda x:round(x+0.0,2),v1)))])
#     return G.Plot(vertex_size=20,vertex_color ="red",edge_colors=E3)
# Graph.plot_graph_with_vertexlabels=plot_graph_with_vertexlabels


# def plot_graph_with_vertexlabels(G):
#     lsft=G.LabelShift()
#     G.vSet_Shift(range(G.order()),[0,0.05])
#     vLb=G.VerticesLabel([("$v_{%s}:%.2f"%t).rstrip('0').rstrip('.')+"$" for t in list(enumerate(map(lambda x:round(x+0.0,2),v1)))])
#     d={1: "blue", -1: "green"}
#     return G.Plot(vertex_size=20,vertex_color ="red",color_by_label={-1:(0.2,1,0.8), 1:'black'})
# Graph.plot_graph_with_vertexlabels=plot_graph_with_vertexlabels


def Merge_And_Append_Edge(G, e):
    if e[1] not in G[e[0]]:
        return "%s is not in Graph G" % str(e)
    else:
        label = G.edge_label(*e[:2])
    H = copy(G)
    H.merge_vertices(e[:2])
    return H.Add_Edges([e+[label]])


Graph.Merge_And_Append_Edge = Merge_And_Append_Edge

def Generalized_Theta(K):
    G=DiGraph(2)
    K1=[k-1 for k in K[0]]
    K2=[k-1 for k in K[1]]    
    D={(0,1):K1,(1,0):K2}
    H=G.add_innerPaths(D,relabel=True)
    return H


def Graphs_WithDegreeSequence(D):
    u'''
    生成给定度序列的所有简单图

    效率较低， 可能是python的原因

    参考文献:
    Király Z. Recognizing graphic degree sequences and generating all realizations[J]. Tech. Rep. Egres TR-2011-11, 2012.

    '''
    import copy
    # graphy node class
    global g_ma_graph
    global g_output_file
    global GL

    class gnode:
        def __init__(self, id, value):
            self.id = id
            self.value = value

    def getnode(D):
        n = len(D)
        node = []
        for i in range(n):
            node.append(gnode(i, D[i]))
        return node
    # 3

    # binary tree node
    #
    class treenode:
        def __init__(self, parent, l=False, r=False, f=False):
            self.parent = parent
            self.Lisred = l
            self.Risred = r
            self.flag = f
            self.left = None
            self.right = None

    def createtree(n, parent):
        if n <= 0:
            return
        else:
            parent.left = treenode(parent)
            parent.right = treenode(parent)
            createtree(n-1, parent.left)
            createtree(n-1, parent.right)
    # get binary tree
    ######################################################

    def program_init(D):
        n = len(D)
        for i in D:
            if i < 0 or i > n-1:
                return False
        a = sum(D)
        if a % 2 == 1:
            return False
        else:
            return True
    # 判断参数是否正确

    def isgraphic(n, D):
        if program_init(D) and LRG(n, D):
            return True
        else:
            return False
    # 判断是否可图化

    def get_ma_graph(n):
        T = []
        ma_graph = []
        for i in range(n):
            for j in range(n):
                T.append(0)
            ma_graph.append(T)
            T = []
        return ma_graph
    # 构造空的邻接矩阵

    def graph_ouput():
        G = Graph(matrix(g_ma_graph))
        s = G.canonical_label().graph6_string()
        if s not in GL:
            GL.append(s)
    # 	g_output_file.writelines(str(g_ma_graph)+'\n')

    def get_aa(D, P, r, a, s, n):
        jieguo = []
        for i in range(n):
            if i < r:
                jieguo.append(D[i]-P[i])
            elif i > r and i <= r+s-a:
                jieguo.append(D[i]-1)
            else:
                jieguo.append(D[i])
        jieguo.sort(reverse=True)
        return jieguo
    #####################################################
    # LRG

    def LRG(n, D):
        w = n
        b = 0
        s = 0
        c = 0
        for k in range(1, n+1):
            b = b+D[k-1]
            c = c+w-1
            while (w > k and D[w-1] <= k):
                s = s+D[w-1]
                c = c-k
                w -= 1
            if(b > c+s):
                return False
            elif w == k:
                return True
    # LRG-extention

    def LRGE(n, r, D):
        w = n
        b = 0
        s = 0
        c = 0
        l1 = 0
        l2 = 0
        ww = [0 for i in range(n)]
        delta = [0 for i in range(n)]
        for k in range(1, n+1):
            b = b + D[k-1]
            c = c + w - 1
            while (w > k and D[w-1] <= k):
                s = s + D[w-1]
                c = c - k
                w -= 1
            ww[k-1] = w + 1
            delta[k-1] = c + s - b
            if k < r:
                if (delta[k-1] == 0 or (delta[k-1] == 1 and D[r-1] <= k)):
                    l1 = k
                if (delta[k-1] == 0 and D[r-1] <= k and l2 == 0):
                    l2 = ww[k-2]-1
            if delta[k-1] < 0:
                return False
            elif w == k:
                kstar = k
                return [True, max(l1, l2)]
    # 3

    def re_part(listD, n):
        if n == 0:
            graph_ouput()
        else:
            s = listD[n-1].value
            if s == 0:
                re_part(copy.deepcopy(listD[0:n-1]), n-1)
            else:
                btree = treenode(None)
                createtree(n-1, btree)
                r = 0  # current deepth
                a = 0  # size of P
                P = [0 for i in range(n-1)]
                D = [listD[i].value for i in range(n-1)]
                downward = True
                while True:
                    if downward == True:
                        if r < n-1:
                            if a < s:
                                btree.Lisred = True
                                if isgraphic(n-1, get_aa(D, P, r, a, s, n-1)) == True:
                                    btree.Risred = True
                                    btree.flag = True
                                else:
                                    btree.Risred = False
                            else:
                                btree.Lisred = False
                                btree.Risred = True
                            if btree.Lisred == True:
                                P[r] = 1
                                a = a + 1
                                listD[r].value = listD[r].value - 1
                                g_ma_graph[listD[r].id][listD[n -
                                                              1].id] = g_ma_graph[listD[n-1].id][listD[r].id] = 1
                                btree = btree.left
                            else:
                                btree = btree.right
                            r = r + 1
                        else:
                            # digui
                            temp = copy.deepcopy(listD[0:n-1])
                            temp.sort(key=lambda x: x.value, reverse=True)
                            re_part(temp, n-1)
                            downward = False
                            r = r - 1
                            if btree == btree.parent.left:
                                P[r] = 0
                                a = a - 1
                                listD[r].value = listD[r].value + 1
                                g_ma_graph[listD[r].id][listD[n -
                                                              1].id] = g_ma_graph[listD[n-1].id][listD[r].id] = 0
                            btree = btree.parent
                    else:
                        if btree.flag == True:
                            downward = True
                            btree.flag = False
                            btree = btree.right
                            r = r + 1
                        else:
                            if r == 0:
                                return
                            r = r - 1
                            if btree == btree.parent.left:
                                P[r] = 0
                                a = a - 1
                                listD[r].value = listD[r].value + 1
                                g_ma_graph[listD[r].id][listD[n -
                                                              1].id] = g_ma_graph[listD[n-1].id][listD[r].id] = 0
                            btree = btree.parent
# 生成给定度序列的所有图的主程序

    if program_init(D) == False:
        print("not graphic")
        exit()
    n = len(D)

    GL = []

    if LRG(n, D) == False:
        print("not graphic")
        exit()
    listD = getnode(D)
    listD.sort(key=lambda x: x.value, reverse=True)
    g_ma_graph = get_ma_graph(n)
    re_part(listD, n)
    return map(Graph, GL)


sage.graphs.graph_generators.graphs.Graphs_WithDegreeSequence = Graphs_WithDegreeSequence


def Hermitian_Adjacency_Matrix(DG):
    n = DG.order()
    A = DG.adjacency_matrix()
    B = matrix(QQ[I], n, lambda i, j: 1*(min(A[i, j], A[j, i])
                                         == 1)+I*(A[i, j] > A[j, i])-I*(A[i, j] < A[j, i]))
    return B


DiGraph.Hermitian_Adjacency_Matrix = Hermitian_Adjacency_Matrix


# 转换图的序列成 latex 文件
def GraphList2TeX(GL, nc=2, subfig=True, name=False, preview=False, src=True):
    engine = _Latex_prefs._option["engine"]
    tempfilename = tmp_filename()
    png_file = tempfilename + ".png"
    tex_file = tempfilename + ".tex"
    if nc != 0:
        L = Reshape(GL, nc=nc)
    else:
        L = GL
    nc = len(L[0])
    longtabu_header = r'''% !TEX TS-program = xelatex
\documentclass{article}
\usepackage{tikz}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{longtable,tabu}
'''
    longtabu_tabhead = r'''
\begin{document}
\begin{longtabu} to\textwidth {*%s{X[c]}}
''' % (nc)
    longtabu_tabtail = r'''
\end{longtabu}
\end{document}
'''
    longtabu_Body = ""
    for i in range(len(L)):
        li = list(enumerate(L[i]))
        Line = "&".join(map(lambda c: r'''
%s
%s
''' % ("%"*30+"Fig:(%s,%s)" % (i, c[0]+1)+"%"*30,
            c[1].Latex() if isinstance(c[1], (Hypergraph, Graph)) else latex(c[1])
       ), li))
        longtabu_Body = longtabu_Body+"\n"+"%"*20 + \
            "Line %s" % (i+1)+"%"*20+"\n"+Line+"\\\\\n"
        if name:
            longtabu_Body = longtabu_Body + \
                "\t&\t".join(map(lambda G: G.name(), L[i]))+"\\\\\n"

    header = r'''% !TEX TS-program = xelatex
\documentclass{article}
\usepackage{tikz}
\usepackage{subfigure}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{array}
\newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}
'''
    tabhead = r'''
\begin{document}
\begin{figure}[htbp]
\begin{tabular}{%s}
''' % (r"C{%s\textwidth}" % (round(0.9/nc, 2))*nc)
    tabtail = r'''
\end{tabular}
\caption{Several options}
\end{figure}
\end{document}
'''
    Body = ""
    for i in range(len(L)):
        li = list(enumerate(L[i]))
        Line = "&".join(map(lambda c: r'''
\subfigure [Subfigure %s] {
\resizebox{%s\textwidth}{!}{
%s
%s
}
}''' % ("(%s,%s)" % (i+1, c[0]+1), round(0.9/len(li), 2), "%"*30+"Fig:(%s,%s)" % (i, c[0]+1)+"%"*30,
            c[1].Latex() if isinstance(c[1], (Hypergraph, Graph)) else latex(c[1])), li))
        Body = Body+"\n"+"%"*20+"Line %s" % (i+1)+"%"*20+"\n"+Line+"\\\\\n"

    if subfig:
        latex_code = header+tabhead+Body+tabtail
    else:
        latex_code = longtabu_header+longtabu_tabhead+longtabu_Body+longtabu_tabtail

    open(tex_file, 'w').write(latex_code)
    if src:
        os.system('open "%s"' % (tex_file))

    if preview:
        open(tex_file, 'w').write(latex_code)
        e = _run_latex_(tex_file, png=True)
        if e.find("Error") == -1:
            return Image(filename=png_file)
        else:
            raise RuntimeError("编译错误，请检查TeX代码")


graphs_list.GraphList2TeX = GraphList2TeX


def GL2TeX(GL, nc=3, ratio_width=0.9, line_height=1.5, rotate=90):
    """
    Convert SageMath graphs (GL) into LaTeX code and images.

    Args:
    GL (list of SageMath objects): A list of SageMath graphs or mathematical objects to convert.
    nc (int): Number of columns in the LaTeX table (default is 3).
    ratio_width (float): Width ratio for the table (default is 0.9).
    line_height (float): Line height for LaTeX images (default is 1.5).
    rotate (int): Rotation angle for images (default is 90 degrees).

    Returns:
    None

    This function takes a list of SageMath objects (e.g., graphs) and generates LaTeX code and images for rendering them.
    The LaTeX code and images are stored in temporary files.

    Note:
    - You can customize the appearance of the LaTeX table by adjusting the LaTeX options in the 'opts' dictionary.
    - The generated LaTeX code is compiled into a temporary '.tex' file and opened using a LaTeX engine.

    Example:
    GL2TeX([Graph([(0, 1)]), Hypergraph([("A", "B")])])
    """
    opts = {"vertex_labels": False, "vertex_color": "black",
            "vertex_size": 10, "figsize": [2, 2]}
    name = "Gr"
    longtabu_header = r'''
    % !TEX TS-program = xelatex
    \documentclass[preview]{standalone}
    \usepackage{adjustbox,float}
    \usepackage{amsmath}
    \usepackage{amssymb}
    \usepackage{graphicx}
    \graphicspath{{fig/},{figures/},{img/}}
    \usepackage{longtable}[=v4.13]%需指定版本, 否则报错! 原因参加下面网页
    %https://tex.stackexchange.com/questions/600724/dimension-too-large-after-recent-longtable-update
    \usepackage{tabu}%
    '''
    longtabu_tabhead = r'''
    \begin{document}
    \begin{center}
    \vspace*{-8pt}
    \tabulinestyle{on 2pt gray}
    \begin{longtabu} to %s\textwidth {*%s{|X[c,m]}|}
    ''' % (RDF(ratio_width), nc)+r"\tabucline -"
    longtabu_tabtail = r'''
    \end{longtabu}
    \vspace*{-3mm}
    \end{center}
    \end{document}
    '''

    import random
    import string
    img_folder = "./temp/img/"
    if nc != 0:
        L = Reshape(GL, nc=nc)
    else:
        L = GL
    nc = len(L[0])
    from pathlib import Path
    Path(img_folder).mkdir(parents=True, exist_ok=True)
    Body = ""
    longtabu_Body = ""
    r = 0
    for row in L:
        r = r+1
        Row = []
        c = 0
        for G in row:
            c = c+1
            cell = ""
            if isinstance(G, (Hypergraph, Graph)):
                G = G.Normalize_Pos()
                P = G.plot(**opts)
                img = name+"%s-%s" % (r, c)
                codeimg = r'''
    \begin{adjustbox}{height=%s cm,center}
    %s
    \end{adjustbox}
    ''' % (RDF(line_height), r"\includegraphics[angle=%s]{%s}" % (rotate, img))
                Row.append(codeimg)
                P.save(img_folder+img+".pdf")
            else:
                if "_latex_" in dir(G):
                    Row.append("$%s$   " % (latex(G)))
                else:
                    Row.append(str(G)+"   ")

        Line = "\t\t"+("\t&\t".join(Row))[:-1]+r"\\\tabucline -"+"\n"
        longtabu_Body = longtabu_Body+Line

    latex_code = longtabu_header+longtabu_tabhead+longtabu_Body+longtabu_tabtail

    engine = _Latex_prefs._option["engine"]
    tempfilename = "./temp/" + \
        ''.join(random.choice(string.ascii_lowercase) for x in range(6))
    tex_file = tempfilename + ".tex"
    open(tex_file, 'w').write(latex_code)
    os.system('open "%s"' % (tex_file))



def get_centers(pos, lists):
    """
    Calculate the center positions of given sets of vertices.

    Parameters:
    pos (dictionary): The positions of vertices.
    lists (list of lists): The given sets of vertices.

    Returns:
    list: The center positions of each set of vertices.
    """
    centers = []
    for lst in lists:
        s = {v: pos[v] for v in lst}
        local_center = Mean(s.values())
        centers.append(local_center)
    return centers

def get_shifts(centers, index):
    """
    Calculate the shifts of positions given centers.

    Parameters:
    centers (list): The center positions of vertices.
    index (int): 0 for shifts on x-axis, 1 for shifts on y-axis.

    Returns:
    list: The shifts of each set of vertices.
    """
    cs = [cv[index] for cv in centers]
    length = max(cs) - min(cs)
    shifts = [length * (i - 0.5) for i in [0, 1 / (len(cs) - 1), .., 1]]
    return shifts

def get_centers(pos, lists):
    """
    Calculate the center positions of given sets of vertices.

    Parameters:
    pos (dictionary): The positions of vertices.
    lists (list of lists): The given sets of vertices.

    Returns:
    list: The center positions of each set of vertices.
    """
    centers = []
    for lst in lists:
        s = {v: pos[v] for v in lst}
        local_center = Mean(s.values())
        centers.append(local_center)
    return centers

def get_shifts(centers, index):
    """
    Calculate the shifts of positions given centers.

    Parameters:
    centers (list): The center positions of vertices.
    index (int): 0 for shifts on x-axis, 1 for shifts on y-axis.

    Returns:
    list: The shifts of each set of vertices.
    """
    if centers==[]:
        return []
    cs = [cv[index] for cv in centers]
    length = max(cs) - min(cs)
    shifts = [length * (i - 0.5) for i in [0, 1 / (len(cs) - 1), .., 1]]
    return shifts

def align_HV(Q, H=[], V=[],scale=[1,1]):
    """
    Align given sets of vertices in a hypergraph along horizontal and vertical directions.

    Parameters:
    Q (Hypergraph): The hypergraph whose vertices are to be aligned.
    H (list of lists, optional): The sets of vertices to be aligned horizontally. Default is an empty list.
    V (list of lists, optional): The sets of vertices to be aligned vertically. Default is an empty list.

    Returns:
        None

    Functionality:
    This function first calculates the center positions of the given vertices. Then, it calculates the shifts of 
    these centers' positions. The vertices are then moved based on these shifts.
    """
    G = copy(Q)
    pos = G.Pos()

    centers_H = get_centers(pos, H)
    centers_V = get_centers(pos, V)
    shifts_x = get_shifts(centers_V, 0)
    shifts_y = get_shifts(centers_H, 1)

    for i in range(len(V)):
        for v in V[i]:
            pos[v][0] = shifts_x[i]

    for i in range(len(H)):
        for v in H[i]:
            pos[v][1] = shifts_y[i]
    for v in G.vertices():
        G._pos[v]=G._pos[v]*diagonal_matrix(scale)
    return G

# HyperG.align_HV = align_HV

Graph.align_HV = align_HV




def CenterPos(G, c=None):
    H = copy(G)
    ndigits = 3
    c = None
    V = H.vertices()
    if isinstance(H, Graph):
        pos = H.Pos()
    else:
        pos = H._pos
    if c == None:
        center = Mean([pos[k] for k in pos])
    else:
        center = pos[c]
    pos = {k: pos[k]-center for k in pos}
    pos = {i: vector(map(lambda p: round(p, ndigits), pos[i])) for i in pos}
    return H


Graph.CenterPos = CenterPos


def rich_repr(self, display_manager, **kwds):
    prefs = display_manager.preferences
    if not Plot_Option.has_key('max_vertex_num_plot'):
        Plot_Option['max_vertex_num_plot'] = 50
    is_small = (0 < self.order() < Plot_Option)
    can_plot = (prefs.supplemental_plot != 'never')
    plot_graph = can_plot and (
        prefs.supplemental_plot == 'always' or is_small)
    # Under certain circumstances we display the plot as graphics
    if plot_graph:
        plot_kwds = dict(kwds)
        plot_kwds.setdefault('title', repr(self))
        output = self.plot(**plot_kwds)._rich_repr_(display_manager)
    if output is not None:
        return output
    # create text for non-graphical output
    if can_plot:
        text = '{0} (use the .plot() method to plot)'.format(repr(self))
    else:
        text = repr(self)
    # latex() produces huge tikz environment, override
    tp = display_manager.types
    if (prefs.text == 'latex' and tp.OutputLatex in display_manager.supported_output()):
        return tp.OutputLatex(r'\text{{{0}}}'.format(text))
    return tp.OutputPlainText(text)
    GenericGraph._rich_repr_ = rich_repr

#
# def findCycles(H):
#     B,C=H.blocks_and_cut_vertices()
#     maxlenC=max(map(len,B))
#     L=[]
#     for lenc in [3..maxlenC]:
#         L=L+list(H.subgraph_search_iterator(graphs.CycleGraph(lenc)))
#     Ecs=[]
#     for c in L:
#         lenc=len(c)
#         Ecs.append(Set([Set([c[i],c[(i+1)%lenc]]) for i in range(lenc)]))
#     Cycles=[Graph(map(tuple,s)) for s in Set(Ecs)]
#     return Cycles
# Graph.cosycles=findCycles
