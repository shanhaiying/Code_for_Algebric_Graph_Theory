import types






class HyperMatrix():
    def __init__(self, *args):
        self.name = 'HyperMatrix' 
        narg=len(args)
        if all([type(1)==type(i) for i in args]):
            self.dims=list(args)
            self.l=reshape([0]*mul(self.dims),self.dims)
        else:
            if isinstance(args[-1], (types.FunctionType,types.BuiltinMethodType, types.LambdaType, types.MethodType)):
                self.dims=list(args[-2]) 
                self.l=reshape([args[-1](i) for i in mrange(self.dims,tuple)],self.dims)
            if any(k in str(type(args[-1])) for k in ["modules","list"]):
                if len(args)==1:
                    self.l=args[-1]
                    self.dims=[len(self.l)]
                else:
                    if isinstance(args[-2],(list,tuple)):
                        self.dims=list(args[-2])
                        self.l=reshape(args[-1],self.dims)
            if isinstance(args[-1],type(1)):
                self.dims=list(args[-2])
                self.l=reshape([args[-1]]*mul(self.dims),self.dims)
            if isinstance(args[-1],tuple):
                self.dims=list(args[-1])
                self.l=reshape([0]*mul(self.dims),self.dims) 
            if "matrix" in str(type(args[-1])):
                self.dims=list(args[-1].dimensions())
                self.l=reshape(args[-1].list(),self.dims)

    def order(self):
        return len(self.dims)

    def get_dims(self):
        t=self.l
        dim=[]
        while isinstance(t,list):
            dim.append(len(t))
            t=t[0]
        return  map(Integer,dim)

    def get_element(self, key):
        t=self.l
        if isinstance(key,(list,tuple)):
            if len(key)!=self.order():
                raise TypeError("dimension error!")                 
            for i in key:                
                t=t[i]
            return t
        else:
            raise TypeError("dimension error!")  

    def set_element(self, key,value):
        t=self.l
        if isinstance(key,(list,tuple)):
            if len(key)!=self.order():
                raise TypeError("dimension error!")                 
            for i in key[:-1]:                
                t=t[i]
            t[key[-1]]=value
        else:
            raise TypeError("dimension error!")  
    def __setitem__(self, key,value):
        self.set_element(key,value)
        return 
        
    def __getitem__(self, key):
        t=self.l
        if isinstance(key,(list,tuple)) and len(key)==self.order() and slice not in map(type, key):
            return self.get_element(key)        
        if isinstance(key,(list,tuple)):
            key=list(key)
            if len(key)>self.order():
                raise TypeError("dimension error!")
            else:
                for i in range(len(key)):
                    if type(key[i])==type(1):
                        key[i]=slice(key[i],key[i]+1)                        
                key=key+[slice(None)]*(self.order()-len(key))
                subindex=map(lambda i,j:range(i)[j],self.dims,key)
                subdims=map(len,subindex)
                indx=map(tuple,cartesian_product(subindex))
                newindx=mrange(subdims,tuple)
                inddict=dict(zip(newindx,indx))
                newHyperMatrix=HyperMatrix(tuple(subdims))
                for k in inddict.keys():
                    newHyperMatrix.set_element(k,self[inddict[k]])       
                return newHyperMatrix
        else:
            return t[key]
       
    def __repr__(self):
        return self.l.__repr__()

                
    def __eq__(self,B):
        return self.l==B.l 

    from operator import add               
    def __add__(self,B):
        if self.dims!=B.dims:
            raise TypeError("dimension error!")
        C=HyperMatrix(tuple(self.dims))  
        for ind in mrange(C.dims,tuple): 
            C[ind]=self[ind]+B[ind]
        return C

    def __rmul__(self,A):
        if A in CC:
            C=HyperMatrix(tuple(self.dims))  
            for ind in mrange(C.dims,tuple): 
                C[ind]=self[ind]*A
            return C 
    def apply_map(self,fun):
        C=HyperMatrix(tuple(self.dims))  
        for ind in mrange(C.dims,tuple): 
            C[ind]=fun(self[ind])
        return C 
        
    def __sub__(self,B):
        if self.dims!=B.dims:
            raise TypeError("dimension error!")
        C=HyperMatrix(tuple(self.dims))  
        for ind in mrange(C.dims,tuple): 
            C[ind]=self[ind]-B[ind]
        return C
    
    def __mul__(self,A):
        if "matrix" in str(type(A)):
            B=HyperMatrix(A)
        if "vector" in str(type(A)) or "free_module_element" in str(type(A)) :
            B=list(A)
            dim=len(B)
            B=HyperMatrix((dim,1),list(A))
        if A in CC:
            C=HyperMatrix(tuple(self.dims))  
            for ind in mrange(C.dims,tuple): 
                C[ind]=self[ind]*A
            return C            
        else:
            B=A
        # if min(B.dims)==1:
        #     B=reduced_order(B)
        conind=self.dims[1:]+B.dims[0:1]
        if Set(conind).cardinality()!=1:
            raise TypeError("dimension error!") 
        order1=self.order()
        order2=B.order()
        cdims=self.dims[0:1]+B.dims[1:]*(order1-1)
        C=HyperMatrix(tuple(cdims))  
        for ind in mrange(C.dims):
            cind=[ind[0]]+reshape(ind[1:],[order1-1,order2-1])      
            C[tuple(ind)]=add([self[tuple([cind[0]]+a)]*mul([B[tuple([a[k]]+cind[k+1])] for k in range(order1-1)]) for a in mrange(self.dims[1:])])
        if min(C.dims)==1:
            C=C.reduced_order()
        return C

    def is_square_tensor(self):
        return Set(self.dims).cardinality()==1
    
    def is_nonnegative(self):
        Sgn=map(sign,flatten(self.l))
        return not -1 in Sgn
    
    def is_positive(self):
        L=flatten(self.l)
        Sgn=map(sign,flatten(self.l))
        return not (-1 in Sgn or 0 in Sgn)
    
    def row_sum(self):
        return [add(flatten(i)) for i in self]
    
    def degree_tensor(self):
        if not self.is_square_tensor():
            raise TypeError("dimension error!") 
        return diagonal_tensor(self.order(),self.row_sum())
    
    def reduced_order(self):
        l=self.dims
        highdims=Set([i for i in range(len(l)) if l[i]>1])
        ndims=[l[k] for k in highdims]
        A=HyperMatrix(tuple(ndims))
        for ind in mrange(l):
            nind=[ind[k] for k in highdims]
            A[nind]=self[ind]
        return A
    
    def spectral_radius(self, niter=2000):
        if not self.is_nonnegative() or not self.is_square_tensor():
            raise TypeError("A must be a nonnegative square tensor")        
        r=self.order()
        n=self.dims[0]
        X=random_vector(RR,self.dims[0],0.5,1)
        X=HM([n],X)
        I=identity_hypermatrix(self.dims)
        eigval_old=+oo
        go_on=True
        while go_on:    
            Y=self*X
            IX=I*X
            D=diagonal_matrix(IX)
            eigval=max(HM(D^-1)*Y)
            if abs(eigval-eigval_old)<1e-18:
                go_on=False        
            X1=Y+IX
            Y1=X1.apply_map(lambda i:operator.pow(i,1/3))
            X=1/add(Y1)*Y1
            eigval_old=eigval
        return eigval,X
        
        

def ones_hypermatrix(dims):
    return HyperMatrix(dims,1)

def zeros_hypermatrix(dims):
    return HyperMatrix(dims,0) 

def identity_hypermatrix(dims):
    return HyperMatrix(dims,lambda k:(Set(k).cardinality()==1)*1)

def random_tensor(dims,lb,up):
    l=[ZZ.random_element(lb,up) for _ in range(mul(dims))]
    return HyperMatrix(dims,l)

def diagonal_tensor(order,lst):
    dims=tuple([len(lst)]*order)
    A=HyperMatrix(dims)
    for i in range(dims[0]):
        A[[i]*order]=lst[i]
    return A


    
HM=HyperMatrix