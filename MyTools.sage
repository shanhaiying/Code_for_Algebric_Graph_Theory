
# def reshape(L,dims,fil=0):
#     u"""
#     当数据不足时, 填充fil的值
#     """
#     if len(L)==0:
#         l=[]
#         for i in dims[::-1]:
#             l=[l*i]
#         return l[0]
#     l=flatten(L)
#     l=l+[fil]*(mul(dims)-len(l))
#     for i in dims[::-1]:
#         l=map(list,zip(*[iter(l)]*i))
#     return l[0]
#
#
def var_list(name,size):
    return [var("%s%i"%(name,i)) for i in range(size)]


def shortLex(s1, s2):
    """Shortlex less or equal comparator"""
    if len(s1) > len(s2):
        return False
    if len(s1) < len(s2):
        return True
    return s1 <= s2

def quickSort(nums):  # 这种写法的平均空间复杂度为 O(nlogn)
    if len(nums) <= 1:
        return nums
    pivot = nums[0]  # 基准值
    left = [nums[i] for i in range(1, len(nums)) if shortLex(nums[i], pivot)] 
    right = [nums[i] for i in range(1, len(nums)) if not shortLex(nums[i],pivot)]
    return quickSort(left) + [pivot] + quickSort(right)


def Majorization_Comp(p1,p2): #比较两序列的优超序， 若不等长，返回lNon, 若等长不可比返回Non， 若p1大于等于p2返回True， 若p1小于p2返回False
    if len(p1)!=len(p2):
        return "lNon"
    p1=sorted(p1,reverse=True)
    p2=sorted(p2,reverse=True)
    res=[]
    for k in range(len(p1)):
        res.append(add(p1[:k+1])>=add(p2[:k+1]))
    Res=list(Set(res))
    if len(Res)==2:
        return "Non"
    else:
        return Res[0]

        



def tab2tex(T):
    import pyperclip
    s=table._latex_(T)
    pyperclip.copy(s)
    show("数据已拷贝至剪切板！")
sage.misc.table.table.tab2tex=tab2tex


def subd_list(L,bk):
    """
    将列表划分为指定大小的子列表组。

    参数：
    L (list)：要划分的列表。
    bk (list)：指定的子列表大小列表。

    返回值：
    LB (list)：划分后的子列表组，每个子列表的大小由bk指定。

    """
    L=list(L)
    LB=[]
    for t in bk:
        bl=[]
        for _ in range(t):
            bl.append(L.pop(0))
        LB.append(bl)
    return LB



def Expr_add(L):
    res="";
    for i in map(LatexExpr,L):
        res+=i;
    return res

poly2sr=lambda f:sum([b*x^a for (a,b) in enumerate(f)])


import re
def multiple_replace(dict, text):
  # Create a regular expression  from the dictionary keys
  regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))
  # For each match, look-up corresponding value in dictionary
  return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)
#对表达式中的普遍性
#Gensubs=lambda exp,su:eval(multiple_replace(su,str(sage_input(exp))))


def Gensubs(exp,su):
    """对sagemath中的数学表达式或数学对象进行字符替换

    Args:
        exp ([type]): 数学表达式或数学对象
        su ([dict]): 字典

    Returns:
        [type]: 替换后的数学对象
    """    
    su= dict(zip(map(str,su.keys()),map(str,su.values())))
    if "matrix" in str(type(exp)):
        A=matrix(eval(multiple_replace(su, str([t for t in exp]))))
        return A.change_ring(PolynomialRing(ZZ,A.arguments()))
    else:
        return eval(multiple_replace(su,str(sage_input(exp))))




#带误差限制的取整数函数:
round_with_eps=lambda x,eps=1e-10: x if abs(round(x)-x)>eps else round(x)


def matrix_round(M,ndigits=8,eps=1e-10):
    to_int=lambda x:int(x) if x.is_integer() else x
    return M.apply_map(lambda x:x if abs(round(x,ndigits)-x)>eps else to_int(round(x,ndigits)))

def map_round(M,ndigits=8,eps=1e-10):
    to_int=lambda x:int(x) if x.is_integer() else x
    my_round=lambda x:round(x,ndigits) if abs(round(x,ndigits)-x)>eps else to_int(round(x,ndigits))
    if "matrix" in str(type(M)):
        return M.apply_map(my_round)
    else:
        return map(my_round,M)

def random_string(l):
    import random
    import string
    s=string.lowercase+string.digits
    return ''.join(random.sample(s,l))






#将ocr的latex矩阵代码转化为Sage矩阵

#先将公式复制到剪切板中

def l2math():
    import pyperclip
    import re
    T = pyperclip.paste()
    k1 = re.findall(r"[-]*\d+", T)
    L= map(Integer, k1)
    nc = len(re.findall(r"{array}{([r|c|l]+)}", T)[0])
    nr=len(L)/nc
    return matrix(nr,L)


def L2Math():
    import pyperclip
    import re
    T = pyperclip.paste()
    FD=re.findall(r"{array}{([r|c|l]+)}(.*)\\end{array}", T)
    opt,dt=FD[0]
    nc=len(opt)
    DT=re.findall(r"{([^\}]*)}", dt)
    dT=[T.replace(" ","*") for T in DT]
    vD=[sage_eval(T, locals={'x':x}) for T in dT]
    M=matrix(len(vD)/nc,vD)
    return M



# 这个函数生成一个包含所有长度为 k 的子列表的列表，子列表以元组形式保存
def generate_sublists(L, k):  
    result = []  
    for i in range(len(L) - k + 1):  
        sublist = tuple(L[i:i + k])
        result.append(sublist)  
    return result

# 这个函数抽取出 psline 命令中的坐标
def extract_coordinates(psline_command):
    import re
    # 正则表达式模式来匹配坐标，坐标是小数或整数，可以是负数
    num_pattern=r'[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?'
    coordinates_pattern =  f'\({num_pattern},{num_pattern}\)'   
    # 找出所有符合 pattern 的坐标
    coordinates = re.findall(coordinates_pattern, psline_command)
    # 使用 sage_eval （一个 SageMath 函数）评估字符串中的坐标，并返回一个坐标列表
    return [sage_eval(cor) for cor in coordinates]

# 该函数从文本中生成图
def genGraph_from_code(text=None):
    import subprocess
    # 如果没有传入文本，那么就从剪贴板获取文本
    if text==None:
        text = subprocess.check_output(['pbpaste'], universal_newlines=True).strip()
    # 提取出所有以 "\psline" 开始的行
    psline_commands = [line for line in text.split("\n") if line.strip().startswith("\\psline") ]
    # 提取出所有以 "\psdots" 开始的行
    psdots_commands = [line for line in text.split("\n") if line.strip().startswith("\\psdots") ]

    # 从 psline 命令中提取路径
    paths=[]
    for psline_command in psline_commands:
        paths.append(extract_coordinates(psline_command))

    # 从 psdots 命令中提取那些 'dots'
    dots=[]
    for psdot_command in psdots_commands:
        dots.append(extract_coordinates(psdot_command))    

    # 对路径使用 generate_sublists 函数来得到边
    Es=[generate_sublists(p,2) for p in paths]
    # 将嵌套的列表拉平，得到所有边
    Edges=list(Set(flatten(Es,max_level=1)))
    # 将嵌套的列表拉平，并删除重复的元素，得到所有顶点
    Vs=Set(flatten(dots,max_level=1))
    # 过滤掉边界情况和无效的边
    Edges=[e for e in Edges if Set(e).cardinality()>1 and Set(e).issubset(Vs)]
    # 使用这些边创建图
    G=Graph(Edges)
    V=G.vertices()
    # 设置图中每个顶点的位置为其自身的值
    G.set_pos({v:v for v in V})
    # 返回重标记后的图
    return G.Relabel()

