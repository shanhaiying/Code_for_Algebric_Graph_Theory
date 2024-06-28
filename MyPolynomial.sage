

#多项式相关计算

#计算f的变号数
def num_chg(coefs):
    return add([mul(coefs[i:i+2])<0 for i in range(len(coefs)-1)])

def num_chg_of_polynomail(f):
    coefs=f.coefficients()
    return num_chg(coefs)
    
def strum_seq(f):
    strum=[f,f.diff(x)]
    while 1:
        q,r=strum[-2].quo_rem(strum[-1])
        if r==0: break
        strum.append(-r)
    return strum
    
def num_realroots(f,a,b):
    nca=list(map(lambda x:sgn(x(a)),strum_seq(f)))
    ncb=list(map(lambda x:sgn(x(b)),strum_seq(f)))
    return num_chg(nca)-num_chg(ncb)



def find_all_roots(f, a, b, eps=0.0000000001):
    roots = []
    intervals_to_check = [(a,b)]
    while intervals_to_check:
        start, end = intervals_to_check.pop()
        try:
            root = find_root(f, start, end)
        except RuntimeError:
            continue
        if root in roots:
            continue
        if abs(f(x=root)) < 1:
            roots.append(root)
        intervals_to_check.extend([(start, root-eps), (root+eps, end)])
    roots.sort()
    return roots