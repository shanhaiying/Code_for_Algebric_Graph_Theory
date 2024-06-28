def magicsquare_normal_odd(n):
    r"""
    Generates nth odd normal cfor n greater than 1 using de la Loubere's method.

    EXAMPLES:
        sage: magicsquare_normal_odd(1)
        [8 1 6]
        [3 5 7]
        [4 9 2]
        sage: magicsquare_normal_odd(2)
        [17 24  1  8 15]
        [23  5  7 14 16]
        [ 4  6 13 20 22]
        [10 12 19 21  3]
        [11 18 25  2  9]

    AUTHOR: Timothy Clemans
    """
    square_length = 2*n+1
    square = [0 for i in range(square_length ^ 2)]

    def position(row, column):
        def row_id(h):
            return (h-1)*square_length
        return row_id(row)+column-1
    for i in range(1, square_length ^ 2+1):
        if i == 1:
            current_position = (1, n+1)
            last_position = current_position
            square[position(*current_position)] = i
        elif last_position == (1, square_length):
            current_position = (2, square_length)
            last_position = current_position
            square[position(*current_position)] = i
        elif last_position[0] == 1:
            current_position = (square_length, last_position[1]+1)
            last_position = current_position
            square[position(*current_position)] = i
        elif last_position[1] == square_length:
            current_position = (last_position[0]-1, 1)
            last_position = current_position
            square[position(*current_position)] = i
        elif square[position(*(last_position[0]-1, last_position[1]+1))] > 0:
            current_position = (last_position[0]+1, last_position[1])
            last_position = current_position
            square[position(*current_position)] = i
        else:
            current_position = (last_position[0]-1, last_position[1]+1)
            last_position = current_position
            square[position(*current_position)] = i
    return matrix(square_length, square)


def Siamese_magic_square(n):
    return [[(i+j+(n+1)/2) % n*n+(i+2*j+1) % n+1 for j in range(n)]
            for i in range(n)]








def my_print(self, *args):
    return "%sum_{"+latex(args[1])+"="+latex(args[2])+"}^{"+latex(args[3])+"} "+latex(args[0])
mysum = function('_sum', nargs=4, print_latex_func=my_print)
class MyIdentity:
    def __init__(self, item, k=var("k"),l=0, u=var("n")):
        self.item=item
        self.k=k
        self.l=l
        self.u=u
        self.left=sum(self.item,self.k,self.l,self.u)
        self.right=mysum(self.item,self.k,self.l,self.u)
    def __repr__(self):
        return "" 
    def __mul__(self, right):
        self.item=(right*self.item).simplify()
        self.left=sum(self.item,self.k,self.l,self.u)
        if self.left!=0:
            self.left=self.left.factor()
        return self
    def action(self,fun, both=False,**args):
        self.item=(fun(self.item)).simplify()
        if both:
            self.left=fun(self.left)
        else:
            self.left=sum(self.item,self.k,self.l,self.u)        
        if self.left!=0:
            self.left=self.left.factor()
        return self

    def subs(self,both=False,**args):
        self.item=self.item.subs(**args)
        if both:
            self.left=self.left.subs(**args)
        else:
            self.left=sum(self.item,self.k,self.l,self.u).subs(**args)
        if self.left!=0:
            self.left=self.left.factor()        
        return self
    def factor(self):
        self.left=self.left.factor()
        return self        
    def simplify(self):
        self.left=sum(self.item,self.k,self.l,self.u)
        self.item=self.item.simplify()
        return self 
    def show1(self):  
        return show(latex((self.left)==mysum(self.item,self.k,self.l,self.u)))
    def actionlist(self,al):
        re=[copy(self)]
        for ac in al:
            re=re+[copy(self.action(ac))]
        return re
        
    def show(self):  
        return show(latex(self.left==mysum(self.item,self.k,self.l,self.u)))

def divMatrix(n):
    return matrix(n,lambda k,i:((k+1)%(i+1)==0).real)

def getGridPath(x):   
    y=[1-i for i in x]
    u,v=[Word(i).partial_sums(0).to_integer_list() for i in [x,y]]
    opt=dict(plotjoined=True,ticks=[1, 1],marker=".",markerfacecolor="blue",color="red",markersize=12,figsize=[6,6], gridlines=True,zorder=-1,aspect_ratio=1,alpha=1,thickness=3)    
    G=list_plot(zip(u,v),frame=false, **opt)
    return G
    

def gridPath(seq,co="blue",withPoints=False):
    P=[(0,0)]
    for ite in seq:
        it=list(P[-1])
        if ite==1:
            it[0]=it[0]+1
        else:
            it[1]=it[1]+1
        P=P+[tuple(it)]
    pic=line(P,marker=".",aspect_ratio=1, frame=False,color=co)
    if withPoints==True:
        return (pic,P)
    return pic

def lattice_path(n,C,co="blue",withPoints=False):
    seq=[0]*n
    for i in C:seq[i]=1
    P=[(0,0)]
    for ite in seq:
        it=list(P[-1])
        if ite==1:
            it[0]=it[0]+1
        else:
            it[1]=it[1]+1
        P=P+[tuple(it)]
    pic=line(P,marker=".",aspect_ratio=1, frame=False,color=co)
    if withPoints==True:
        return (pic,P)
    return pic
    
def rectangle(minx,miny,maxx,maxy,co="black",f=True):
    if f==False:
        return line([(minx,miny),(maxx,miny),(maxx,maxy),(minx,maxy),(minx,miny)],color=co)
    else:
        return polygon([(minx,miny),(maxx,miny),(maxx,maxy),(minx,maxy),(minx,miny)],color=co)


def chess_board(M,cls,row_label=None,col_label=None):
    grids=Graphics()
    m,n=M.dimensions()
    for i in range(m):
        for j in range(n):
            l=line([(i,j),(i+1,j),(i+1,j+1),(i,j+1),(i,j)],color="black")
            rec=polygon([(i,j),(i+1,j),(i+1,j+1),(i,j+1)],color=cls(i,j))
            grids=grids+l+rec
    if row_label!=None and col_label!=None:
        for i in range(n):
            grids=grids+text(col_label[i], (i+0.5,-0.3))
        for j in range(n):
            grids=grids+text(row_label[i], (-0.3,i+0.5))
    return grids

def rook_polynomial(M):
    B=BipartiteGraph(M)
    f=B.matching_polynomial()
    f=(f*x^(-B.order())).expand()
    f=f.subs(x=x^(-1/2))
    return f.subs(x=-x)


def Partition_Number(n,k):
    if n==k or (n>0 and k==1): 
        return 1
    if k>n or (n>0 and k<1):
        return 0
    else:
        return add(Partition_Number(n-k,i) for i in (1..k))
#定义集合的包含关系矩阵
def zeta_mat(S):
    PS=S.subsets()
    n=PS.cardinality()
    M=matrix(n)
    for i in (0..n-1):
        for j in (0..n-1):
            M[i,j]=1*(PS[i].issubset(PS[j]))
    return M       

def BD_incidence_matrix(S,B):
    v=len(S)
    b = len(B)
    MS = MatrixSpace(ZZ, v, b)
    A = MS(0)
    for i in range(v):
        for j in range(b):
                A[i, j] = 1*(S[i] in B[j])
    return A

def divMatrix(n):
    return matrix(n,lambda k,i:((k+1)%(i+1)==0).real)

def mobuisMatrix(n):
    return matrix(n,lambda k,i:moebius((k+1)/(i+1)) if (k+1)%(i+1)==0 else 0)

def MyPowerSeries(f,R):
    x=R.gen()
    return eval(str(f))


#求解递推关系
def solve_rec(relation):
    maxima('load(solve_rec)')
    var('n')
    cmd='solve_rec('+relation+')'
    sol=maxima(cmd)
    return latex((sol.rhs().sage().simplify_full()))

    


#轮换类型为T的置换的计数(柯西Cauchy公式)
def num_cycletype(t):
    t=list(t)
    n=add(t)
    return factorial(n)/ prod([factorial(t.count(i))*(i^t.count(i)) for i in (1..n)])


#定义置换群的轮换指数多项式
def cycle_index_polynomial(G):    
    X=[var("x%s"%i) for i in G.domain()]
    G=map(Permutation,G)
    return add([prod([X[i-1]^list(g.cycle_type()).count(i) for i in (1..max(g.cycle_type()))]) for g in G])



#返回教材给定形式的轮换多项式
def getCycleIndex1(SG5):
    P = SG5.cycle_index()
    n=P.degree()
    R=PolynomialRing(QQ, 'x', n)
    mc=P.monomial_coefficients()
    return add([mc[d]*R.monomial(*(d+[0]*(n-len(d)))) for d in mc])

def Equivs(S,k):
    """
    从等价划分S中取出所有互不等价的k个代表元
    S=[[0, 5], [1, 4], [2], [3]]
    Equivs(S,2)
    输出：
    [[0, 5], [0, 1], [0, 2], [0, 3], [1, 4], [1, 2], [1, 3], [2, 3]]    
    """
    bk=map(len,S)
    ms=flatten([[i]*bk[i] for i in range(len(bk))])
    def lst2dic(it):
        theS=Set(it)
        K={i:it.count(i) for i in theS}
        return flatten([S[i][:K[i]] for i in K])
    L=Combinations(ms,k)
    return map(lst2dic,L)

    

#返回教材给定形式的轮换多项式
def getCycleIndex(SG5):
    P = SG5.cycle_index()
    n=P.degree()
    mc=P.monomial_coefficients()
    X=[var("x%d"%i) for i in [1..n]]
    return add([mc[d]*mul([X[i-1] for i in d]) for d in mc])
    
#定义算子的幂函数：将函数f与自己复合n次
power_operator=lambda f,n:reduce(compose, [f]*n, lambda args: args)
#定义一阶差分算子
diff_op=lambda f:lambda x:f(x+1)-f(x)
#通过一阶差分算子自身复合运算定义高阶差分；下面给出的是0至5阶的差分算子列表
diffop_list=[power_operator(diff_op,i) for i in  range(6)]  


#生成图G的匹配数m(G,k)的生成函数

def matching_generating_function(G):
    n=G.order()
    f=G.matching_polynomial()
    g=x^n*f.subs(x=1/x)
    g=g.power_series(ZZ).polynomial()
    g=g.subs(x=i*x)
    g=g.subs(x=x^(1/2))
    g=g.series(x,n)
    return g 
