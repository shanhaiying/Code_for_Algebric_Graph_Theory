
###############################分块矩阵###############
#块型转化为足标分组
Bs2Bk=lambda Bs:[[add(Bs[:i])..add(Bs[:i])+Bs[i]-1] for i in range(len(Bs))]

#块型转化为分割线标号
Bk2BL=lambda BK:[b[-1]+1 for b in BK[:-1]]

#分割线标号转化为块型(各块中的行不一定要求是连续的行)
def BL2Bk(n_Bk):
    n,rL=n_Bk
    Rind=range(n)
    rL=sorted(rL)
    if rL[-1]>=n:
        return
    return [Rind[:rL[0]]]+[Rind[rL[i-1]:rL[i]] for i in [1..len(rL)-1]]+[Rind[rL[-1]:]]


#对矩阵按行型列型进行分块
def Subdivide(A,rBk,cBk=None,bksize=True):
    M=copy(A)
    if type(rBk[0])==list:
        bksize=False
    if cBk==None:
        cBk=rBk
    if bksize:
        M.subdivide(*[Bk2BL(Bs2Bk(t)) for t in [rBk,cBk]])
    else:
        if type(rBk[0])==list:
            rind,cind=list(map(flatten,[rBk,cBk]))
            nrBk,ncBk=list(map(len,[rBk,cBk]))
            M=M[rind,cind]
            M.subdivide(*[Bk2BL(Bs2Bk(list(map(len,t)))) for t in [rBk,cBk]])
        else:
            M.subdivide(rBk,cBk)
    return M



def subdivide(N,rBk,cBk=None,bksize=True):
    M=copy(N)
    if cBk==None:
        cBk=rBk
    if bksize:
        rsd=Composition(rBk).partial_sums()[:-1]
        csd=Composition(cBk).partial_sums()[:-1]
        M.subdivide(rsd,csd)
    else:
        M.subdivide(rsd,csd)
    return M



def Submatrix(A,U,V=None):
    """返回矩阵的子阵或余子阵

    Args:
        A (Matrix): 矩阵
        U (list,tuple or int): 行指标集
        V (_type_, optional): 列指标集 Defaults to None.

    Returns:
        _type_: 子阵或余子阵
    """    
    if V==None:
        V=U
    m,n=A.dimensions()
    C=[]
    R=[]
    if isinstance(U,int):
        U=(U,)
    if isinstance(V,int):
        V=(V,)
    if isinstance(U,tuple):
        R=[i for i in range(m) if i not in U]
    if isinstance(U,list):
        R=[i for i in range(m) if i in U]
    
    if isinstance(V,tuple):
        C=[i for i in range(n) if i not in V]
    if isinstance(V,list):
        C=[i for i in range(n) if i in V]
    return A[R,C]



#对分块矩阵进行置换
def perm_block(A,R,C=None):
    if not C:
        C=R
    r,c=len(R),len(C)
    return block_matrix(r,c,[A.subdivision(R[i],C[j])  for i in range(len(R)) for j in  range(len(C))])


def bk_rop1(A,i,j):
    rBk,cBk=map(BL2Bk,zip(A.dimensions(),A.subdivisions()))
    rbk,cbk=map(len,[rBk,cBk])
    E=elementary_matrix(rbk,row1=i,row2=j)
    perm=E*vector(range(rbk))
    A=block_matrix(rbk,cbk,[A.subdivision(perm[i],j)  for i in range(rbk) for j in  range(cbk)])
    return A


def bk_rop2(A,i,p):
    rBk,cBk=map(BL2Bk,zip(A.dimensions(),A.subdivisions()))
    rbk,cbk=map(len,[rBk,cBk])
    if p in RR:
        p=MatrixSpace(QQ,len(rBk[i]))(p)
    A[rBk[i]]=p*A[rBk[i]]
    return A

def bk_rop3(A,i,j,p):
    rBk,cBk=map(BL2Bk,zip(A.dimensions(),A.subdivisions()))
    rbk,cbk=map(len,[rBk,cBk])
    if p in RR:
        p=MatrixSpace(QQ,len(rBk[i]))(p)
    A[rBk[i]]=A[rBk[i]]+p*A[rBk[j]]
    return A




######################################################




def show_equtions(equsys):
    header=r"\left\{\begin{array}{l,l}"
    tail=r"\end{array}\right."
    body=header
    for i in range(len(equsys)-1):
        body=body+latex(equsys[i])+r"\\"
    body=body+latex(equsys[len(equsys)-1])+tail
    return Math(body.replace("=","&="))


def stdequs(equsys):
    newequsys=[]
    for equ in equsys:
        equ=equ.rhs()-equ
        equ=equ.lhs().subs(dict(zip(equ.variables(),len(equ.variables())*[0])))-equ
        newequsys.append(equ)
    return newequsys

def movelines(equsys, x0,y0):
    delta_ys=[solve(equ.subs(x=x0),y,solution_dict=True) for equ in equsys]
    dys=flatten(delta_ys)
    newequsys=[equsys[i].subs(y=y-y0+dys[i][y]) for i in range(len(equsys))]
    return newequsys


def plotequs(equsys):
    p = implicit_plot(equsys[0],(x,-3,3), (y,-3,3), color=hue(0/len(equsys)))+ text('equ1:  $'+latex(equsys[0])+'$', (6,2.5),color=hue(0/len(equsys)))
    for n in range(1,len(equsys)):
        p += implicit_plot(equsys[n],(x,-3,3), (y,-3,3), color=hue(n*1.0/len(equsys)),frame=False,axes=True)+ text('equ'+str(n+1)+':  $'+latex(equsys[n])+'$', (6,2.5-0.4*n),color=hue(n*1.0/len(equsys)))
    return p


def plot3dequs(equsys):
    p = implicit_plot3d(equsys[0],(x,-3,3), (y,-3,3), (z,-3,3), color=hue(0/len(equsys)))+ text3d('equ1:  $'+latex(equsys[0])+'$', (6,2.5,4),color=hue(0/len(equsys)))
    for n in range(1,len(equsys)):
        p += implicit_plot3d(equsys[n],(x,-3,3), (y,-3,3), (z,-3,3), color=hue(n*1.0/len(equsys)),frame=False,axes=True)+ text3d('equ'+str(n+1)+':  $'+latex(equsys[n])+'$', (6,2.5-0.4*n,4),color=hue(n*1.0/len(equsys)))
    return p

def GenerateMatrix(equsys, vars):
    A=matrix([[equ.lhs().coefficient(v) for v in vars] for equ in equsys])
    b=matrix([[equ.rhs()] for equ in equsys])
    return (A,b)

def symbolMatrix(m,n,s):
    return matrix(m,n,[[var(s+str(i)+str(j),latex_name="{"+s+"_{"+str(i)+","+str(j)+"}}") for j in range(1,n+1)] for i in range(1,m+1)])


#定义初等变换(elementary operation)的Latex输出函数
def el_rop1(M,c,repl=False):
    pretty_print(LatexExpr(latex(M))+LatexExpr(r"\Aequiv{r_{%s} \leftrightarrow r_{%s}}"%tuple(map(operator.add,c,[1,1])))+LatexExpr(latex(M.with_swapped_rows(*c))))
    if repl:
        M.swap_rows(*c)

def el_rop2(M,c,repl=False):
    pretty_print(LatexExpr(latex(M))+LatexExpr(r"\Aequiv{{%s}\times r_{%s}}"%tuple(map(operator.add,c,[1,1])[::-1]))+LatexExpr(latex(M.with_rescaled_row(*c))))
    if repl:
        M.rescale_row(*c)

def el_rop3(M,c,repl=False):
    pretty_print(LatexExpr(latex(M))+LatexExpr(r"\Aequiv{" +'r_{{{0}}}'.format(c[0]+1)+(r"\,+\," if c[2] >=0 else r"\," )+(r'{1}\, \times \, r_{{{0}}}}}'.format(c[1]+1,c[2])))+LatexExpr(latex(M.with_added_multiple_of_row(*c))))
    if repl:
        M.add_multiple_of_row(*c)

def el_cop1(M,c,repl=False):
    pretty_print(LatexExpr(latex(M))+LatexExpr(r"\Aequiv{c_{%s} \leftrightarrow c_{%s}}"%tuple(map(operator.add,c,[1,1])))+LatexExpr(latex(M.with_swapped_columns(*c))))
    if repl:
        M.swap_columns(*c)

def el_cop2(M,c,repl=False):
    pretty_print(LatexExpr(latex(M))+LatexExpr(r"\Aequiv{{%s}\times c_{%s}}"%tuple(map(operator.add,c,[1,1])[::-1]))+LatexExpr(latex(M.with_rescaled_col(*c))))
    if repl:
        M.rescale_col(*c)

def el_cop3(M,c,repl=False):
    pretty_print(LatexExpr(latex(M))+LatexExpr(r"\Aequiv{" +'c_{{{0}}}'.format(c[0]+1)+(r"\,+\," if c[2] >=0 else r"\," )+(r'{1}\, \times \, c_{{{0}}}}}'.format(c[1]+1,c[2])))+LatexExpr(latex(M.with_added_multiple_of_column(*c))))
    if repl:
        M.add_multiple_of_column(*c)

def lst2ltx(L,Det=False):
    if Det:
        latex.matrix_delimiters("|", "|")
    else:
        latex.matrix_delimiters("(", ")")
    nL=[it if isinstance(it,str) else latex(it) for it in L]
    latex.matrix_delimiters("(", ")")
    itms=map(LatexExpr,nL)
    return LatexExpr(add(itms)[2:])

def sum_hold(L,Det=True):
    if Det:
        latex.matrix_delimiters("|", "|")
    else:
        latex.matrix_delimiters("(", ")")
    res=LatexExpr("+".join(map(latex,L)))
    latex.matrix_delimiters("(", ")")
    return res

#对行列式进行拆行变换
def split_row(A,i,Ls,disp=True):
    Ls=map(vector,Ls)
    Ls.append(A[i]-add(Ls))
    ML=[]
    for r in Ls:
        C=copy(A)
        C[i]=r
        ML.append(C)
    display(lst2ltx([A,"=",sum_hold(ML)],Det=True))
    if not disp:
        return [A]+ML

def split_vector(v):
    VL=[]
    R=v.base_ring()
    n=v.degree()
    for i in range(n):
        if v[i]!=0:
            tmpv=vector(R, n)
            tmpv[i]=v[i]
            VL.append(tmpv)
    return VL

def laplacian_expansion(A,i):
    Ls=split_vector(A[i])[:-1]
    rL=split_row(A,i,Ls,disp=False)
    return rL


def laplacian_expansion(A,i,minor=True,ltx=True):
    if not minor:
        Ls=split_vector(A[i])[:-1]
        rL=split_row(A,i,Ls,disp=False)
        return rL
    else:
        m,n=A.dimensions()
        AML=[alg_Mimor(A,i,j) for j in range(m)]

        if ltx:
            return LatexExpr(r"+".join(AML).replace("+-","-"))
        else:
            return AML




def alg_Mimor(A,i,j,ltx=True):
    R,C=map(range,A.dimensions())
    i=R.pop(i)
    j=C.pop(j)
    a,m=(-1)^(i+j)*A[i,j],A[R,C]
    ltxexp=lst2ltx([a,r"\cdot",m],Det=True)
    if ltx:
        return ltxexp
    else:
        return (-1)^(i+j)*A[i,j],A[R,C]

def toeplitz_matrix(R,c,r):
    return matrix(R, len(c), len(r)+1, lambda i,j: c[i-j] if i>=j else r[j-i-1])


#将具体矩阵转化为一般形式矩阵的Latex表示
def mat_generlize(A,r,c,D=None,ltx=True,Det=False):
    C=copy(A)
    m,n=C.dimensions()
    var("xtx",latex_name=r"\cdots")
    var("ytx",latex_name=r"\vdots")
    var("ztx",latex_name=r"\ddots")
    nA=C.T
    rr=n+r if r<0 else r
    cc=n+c if c<0 else c
    nA[cc]=[xtx]*m
    nA=nA.T
    nA[rr]=[ytx]*n
    nA[rr,cc]=ztx
    latex.matrix_column_alignment('c')
    if Det:
        latex.matrix_delimiters("|", "|")
        res=latex(nA)
        latex.matrix_delimiters("(", ")")
    else:
        latex.matrix_delimiters("(", ")")
        res=latex(nA)
    if D !=None:
        res=multiple_replace(D,res)
    if ltx:
        return res
    else:
        return Math(res)

#对矩阵进行初等变换
def do_Actions(A,ops,repl=False):
    if not repl:
        M=copy(A)
    else:
        M=A
    t=type(A)
    if "matrix" not in str(t):
        return
    el_oprs=[t.swap_rows,t.rescale_row,t.add_multiple_of_row]
    el_opcs=[t.swap_columns,t.rescale_col,t.add_multiple_of_column]
    elops=dict(list(zip(["r%d"%i for i in [1,2,3]],el_oprs))+list(zip(["c%d"%i for i in [1,2,3]],el_opcs)))
    for it in ops:
        elops[it[0]](M,*it[1:])
    return M


def Op_links(ops):
    opexps=[
    lambda op:LatexExpr(r"{\Large{r}}_{%s} \leftrightarrow {\Large{r}}_{%s}"%tuple(map(operator.add,op[1:],[1,1]))),
    lambda op:LatexExpr(r"{%s}\times {\Large{r}}_{%s}"%tuple(map(operator.add,op[1:],[1,0])[::-1])),
    lambda op:LatexExpr(r"{\Large{r}}_{{{%d}}}"%(op[1]+1)+(r"\,+\," if op[3] >=0 else r"\," )+r"%s\, \times \, {\Large{r}}_{%s}"%(op[3],op[2]+1)),
    lambda op:LatexExpr(r"{\Large{c}}_{%s} \leftrightarrow {\Large{c}}_{%s}"%tuple(map(operator.add,op[1:],[1,1]))),
    lambda op:LatexExpr(r"{%s}\times {\Large{c}}_{%s}"%tuple(map(operator.add,op[1:],[1,0])[::-1])),
    lambda op:LatexExpr(r"{\Large{c}}_{{{%d}}}"%(op[1]+1)+(r"\,+\," if op[3] >=0 else r"\," )+r"{0}\, \times %d {\Large{c}}_{{{1}}}".format(op[3],op[2]+1))]
    opexp=dict(zip(["r%d"%i for i in [1,2,3]]+["c%d"%i for i in [1,2,3]],opexps))
    exps=[opexp[op[0]](op) for op in ops]
    def links(exps):
        if len(exps)>1:
            return r"%stackrel{"+LatexExpr(r"\begin{subarray}{c}"+r"\\  ".join(map(latex,exps[:ceil(len(exps)/2)])) + r"\end{subarray}")+"}"+r"{\widetilde{" +LatexExpr(r"\begin{subarray}{c}"+r"\\  ".join(map(latex,exps[ceil(len(exps)/2):])) + r"\end{subarray}")+"}}"
        else:
            return r"{%stackrel{%s}{\widetilde{\hphantom{%s} } }}"%tuple(map(latex,exps*2))
    return links(exps)

def do_show_ops(A,ops):
    return latex(A)+Op_links(ops)+latex(do_Actions(A,ops,repl=True))


#将符号矩阵转化为多项式环上的矩阵
def Mat2PR(C,R):
    return matrix(C,PolynomialRing(R,C.arguments()))


#带误差限制的取整数函数:
round_with_eps=lambda x,eps=1e-10: x if abs(round(x)-x)>eps else round(x)


def matrix_round(M,ndigits=8,eps=1e-10):
    to_int=lambda x:int(x) if x.is_integer() else x
    return M.apply_map(lambda x:x if abs(round(x,ndigits)-x)>eps else to_int(round(x,ndigits)))

def map_round(M,ndigits=8,eps=1e-10):
    to_int=lambda x:int(x) if x.is_integer() else x
    my_round=lambda x:x if abs(round(x,ndigits)-x)>eps else to_int(round(x,ndigits))
    if "matrix" in str(type(M)):
        return M.apply_map(my_round)
    else:
        return map(my_round,M)


# #将列表显示为每行nc个项的图表
# def reshape(L,nc):
#     nL=L+["</pre>"]*(bool(len(L)%nc<>0)*nc-(len(L)%nc))
#     return [nL[i:i+nc] for i in range(0,len(nL),nc)]



#计算方阵的特征矩阵

def charmat(A,var=x):
    R.<x> = QQ[];
    A.is_square()
    n=A.dimensions()[0]
    return x*identity_matrix(R,n)-matrix(R,A)

def Jordan_block(t):
    J=block_diagonal_matrix([companion_matrix(R(t[0]))]*t[1])
    for j in range(R(t[0]).degree()*t[1]-1):
        J[j+1,j]=1
    return J
def Jordan_form(L):
    return block_diagonal_matrix(flatten([[Jordan_block(t[0].factor_list()[0])]*t[1] for t in L]))


#利用若当标准型计算矩阵函数
def matrix_function_Jordan(A, f):
    A=A.change_ring(QQbar)
    # returns jordan matrix J and invertible matrix P such that A = P*J*~P
    [J, P] = A.jordan_form(transformation=True)

    fJ = zero_matrix(SR, J.ncols())
    num_Jordan_blocks = 1+len(J.subdivisions()[0])
    fJ.subdivide(J.subdivisions())

    for k in range(num_Jordan_blocks):

        # get Jordan block Jk
        Jk = J.subdivision(k, k)

        # dimension of Jordan block Jk
        mk = Jk.ncols();

        fJk = zero_matrix(SR, mk, mk)

        # compute the first row of f(Jk)
        vk = [f.derivative(x, i)(Jk[i][i])/factorial(i) for i in range(mk)]

        # insert vk into each row (above the main diagonal)
        for i in range(mk):
            row_Jk_i = vector(SR, zero_vector(SR, i).list() + vk[0:mk-i])
            fJk.set_row(i, row_Jk_i)

        fJ.set_block(k, k, fJk)

    fA = P*fJ*~P

    return fA

def solve_linear_system(A, B):
    """
    Solve the linear system by smith form of coefficent matrix.
    We can use the function to solve systems of linear diophantine equations
    Example::
        A=matrix(3,5,lambda i,j:(i+1)^j)
        var("n d1")
        B=vector([n-1,4*(n-1)-d1,9*(n-1)+(n-1)^2-d1^2])
        res=solve_linear_system(A,B)

    """
    D, P, Q = A.smith_form()
    r, c = A.dimensions()
    Y = vector([var("y%i" % i) for i in range(c)])
    dif = D*Y-P*B
    su = solve(dif.list(), Y, solution_dict=True)
    Y1 = Y.subs(su[0])
    X = Q*Y1
    return X.column()


def Mat_Subs(M, D, ring="SR"):
    return eval("matrix(%s,%s,%s)" % (ring, M.dimensions()[0], multiple_replace(D, str(M.list()))))


def Compound_matrix(A,k):
    """计算矩阵A的k级乘性复合矩阵

    Args:
        A ([type]): [description]
        k ([type]): [description]

    Returns:
        [type]: A的k级乘性复合矩阵
    """    
    nr,nc=A.dimensions()
    Cr=Combinations(nr,k)
    Cc=Combinations(nc,k)
    return matrix(len(Cr),len(Cc),lambda i,j:A[Cr[i],Cc[j]].det())
    
def Additive_Compound_matrix(A,k):
    """计算矩阵A的k级加性复合矩阵

    Args:
        A ([type]): [description]
        k ([type]): [description]

    Returns:
        [type]: A的k级加性复合矩阵
    """    
    R=A.base_ring()
    nr=A.dimensions()[0]
    Q=R[x]
    B=(identity_matrix(nr)+x*A).change_ring(Q)
    CB=Compound_matrix(B,k)
    return matrix(CB.dimensions()[0],lambda i,j:CB[i,j][1])
    
    