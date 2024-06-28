# -*- coding: utf-8 -*-

#########################################################################################
#      Here is  an enhancement to the part of  graph theory within sage                 #
#                                                                                       #
#       Copyright (C) 2014~2023 HaiyingShan <shan_haiying@tongji.edu.cn>                #
#                             Last Modified: 2023-11-10                                  #
#########################################################################################
import numpy as np
from sage.structure.element import is_Matrix,is_Vector
class Tensor:
    def __init__(self, input_data=None, shape=None, symmetric=False, info=None):
        if  isinstance(input_data, (np.ndarray,list,Tensor)) or is_Matrix(input_data) or is_Vector(input_data):
            self.array = np.asarray(input_data)
            self.shape = self.array.shape
            if shape is not None:
                self.array=self.array.reshape(shape)
                self.shape = shape
        elif  isinstance(input_data, (dict)):
            keys = list(input_data.keys())
            values = list(input_data.values())
            self.shape = np.max(keys, axis=0) + 1
            self.array = np.zeros(self.shape)
            for key, value in zip(keys, values):
                if symmetric==True:
                    for idx in list(Permutations(key)):
                        self.array[idx] = value
                else:
                    self.array[key] = value
        elif isinstance(input_data, RingElement):
            self.array = np.full(shape, input_data)
        elif input_data==None and shape !=None:
            self.array = np.full(shape, 0)
        else:
            raise ValueError("Either input_array or shape must be provided.")           
        self.shape = self.array.shape
        self.rank = len(self.shape)
        self.ndim = self.array.ndim
        self.symmetric=symmetric
        self.dtype = self.array.dtype
        self.info = info


    def size(self):
        return np.size(self.array)

    def transpose(self,args):
        return Tensor(self.array.transpose(args))

    def reshape(self, shape):
        return self.array.reshape(shape)
    
    def is_cubical_tensor(self):
        return all(dim == self.shape[0] for dim in self.shape)  

    def to_dict(self,symmetric=None):
        """
        Converts the tensor to a dictionary representation.

        Returns:
            dict: A dictionary representing the tensor, where the keys are tuples
                  representing the tensor indices and the values are the corresponding
                  tensor elements.
        """    
        Idx=[idx for idx in mrange(self.shape) if self[*idx]!=0]
        if  self.symmetric==None:
            self.symmetric=self.is_symmetric()
        if self.symmetric==True and symmetric==True:
                Idx=Set([sorted(idx) for idx in Idx])
        return {tuple((idx)):self[*idx] for idx in Idx} 
        
            

    def mean(self, axis=None):
        return self.array.mean(axis=axis)

    def tensordot(self, other, axes=2):
        return np.matmul(self.array, other)
    
    def tensor_power(self, power):
        result = self.array.copy()
        for _ in range(power - 1):
            result = np.tensordot(result, self.array, axes=0)
        return Tensor(result) 
    
    def shao_mul(A, B,start=1):
        """
        Perform a tensor multiplication between two arrays.

        Parameters:
        A (ndarray): The first array to be multiplied.
        B (ndarray): The second array to be multiplied.

        Returns:
        ndarray: The result of the tensor multiplication.

        Examples:
        >>> A = np.array([[1, 2], [3, 4]])
        >>> B = np.array([[5, 6], [7, 8]])
        >>> shao_mul(A, B)
        array([[19, 22],
            [43, 50]])
        """
        res = A
        for i in range(start, len(A.shape)):
            res = np.tensordot(res, B, [start, 0])
        return Tensor(res)

    def shao_mul1(A, LB):
        """
        Perform a tensor multiplication operation using the `einsum` function.

        Parameters:
        A (ndarray): The input multi-dimensional array.
        LB (list of ndarrays): The list of arrays to be multiplied with `A`.

        Returns:
        ndarray: The result of the tensor multiplication.

        """
        k = len(A.shape)
        pms = [A, list(range(k))]
        outaxes = [0]
        lstax=k
        for i in range(1, k):
            lnb=len(LB[i-1].shape)
            baxes = [i] + [lstax + j for j in range(lnb - 1)]
            lstax=lstax+lnb-1
            pms += [LB[i-1], baxes]
            outaxes += baxes[1:]
        pms.append(outaxes)
        return np.einsum(*pms)    

    def symmetrized(self):
        """
        Returns a generally symmetrized tensor, calculated by taking
        the sum of the tensor and its transpose with respect to all
        possible permutations of indices
        """
        perms = list(Permutations(range(self.rank)))
        return Tensor(np.sum([np.transpose(self, ind) for ind in perms],axis=0) / len(perms))
    
     
    def is_symmetric(self,epsilon=1e-6):
        # if len(self.shape) < 2:
        #     return True
        # if self.symmetric!=None:
        #     return self.symmetric
        # else:
        self.symmetric= np.max(self.symmetrized()-self)<epsilon
        return self.symmetric
        
    def is_diagonal_tensor(self):
        return len([key for key in self.to_dict() if len(Set(key))>1])==0
    

    def identity_tensor(self,order=2,dim=3):
        return Tensor({tuple([i]*order):1 for i in range(dim)})



    def sub_tensor(self, axis: list[list[int]]):
        if len(axis) != len(self.shape):
            raise ValueError("The number of indices must match the number of dimensions of the tensor.")
        Idx = cartesian_product(axis)
        shape = [len(ax) for ax in axis]
        return Tensor(np.array([self[tuple(idx)] for idx in Idx]).reshape(shape))    

    def lin_comb(self, c, k=0):
        """
        Calculate the linear combination of subtensor along axis k of arrays.

        Parameters:
            c (np.ndarray): The array to be multiplied with self.
            k (int, optional): The index of the array to be multiplied with c. Defaults to 0.

        Returns:
            np.ndarray: The result of the linear combination.
        """
        rk = len(self.shape)
        return np.einsum(self,[0..rk-1],c,[k])

    def SpectralRadius(self,epsilon=1e-6, max_iterations=100):
        if not self.is_cubical_tensor():
            raise TypeError("The tensor must be cubical.")
        order=self.shape[0]
        m=len(self.shape)
        def F(A,X,m):
            TX=Tensor(X).tensor_power(m-1)
            Y=np.tensordot(A.array,TX,[[1..m-1],[0..m-2]])
            lambda_1=np.dot(X, Y)
            return lambda_1,np.power((Y*X)/lambda_1,1/m)

        X=np.array([1-np.random.random() for _ in range(order)])
        abs_err,step=1,0
        while abs_err>epsilon and step < max_iterations:
            lambda_1,X1=F(self,X,m)
            X,abs_err,step=X1,vector(X1-X).norm(),step+1
        return lambda_1,vector(X)
    
    def save(self, filename):
        np.save(filename, self.array)

  #  @property


    

    @classmethod
    def load(cls, filename):
        array = np.load(filename)
        return cls(array)

    def __str__(self):
        return str(self.array)

    def __repr__(self):
        return repr(self.array)

    def __len__(self):
        return len(self.array)

    def __iter__(self):
        return iter(self.array)

    def __getitem__(self, key):
        return self.array[key]

    def __setitem__(self, key, value):
        self.array[key] = value

    def __eq__(self, other):
        return np.array_equal(self.array, other)

    def __ne__(self, other):
        return not np.array_equal(self.array, other)

    def __gt__(self, other):
        return np.greater(self.array, other)

    def __lt__(self, other):
        return np.less(self.array, other)

    def __ge__(self, other):
        return np.greater_equal(self.array, other)

    def __le__(self, other):
        return np.less_equal(self.array, other)

    def __abs__(self):
        return np.absolute(self.array)

    def __neg__(self):
        return Tensor(np.negative(self.array))

    def __add__(self, other):
        return Tensor(np.add(self.array, other))

    def __sub__(self, other):
        return Tensor(np.subtract(self.array, other))
    
    def __mul__(self, other):
        if isinstance(other, RingElement):
            return Tensor(self.array * other)
        elif isinstance(other, np.ndarray):
            return Tensor(self.array * other)
        elif isinstance(other, Tensor):
            return Tensor(self.array * other.array)
        else:
            raise TypeError("Unsupported operand type(s) for mul.")
    def __rmul__(self,other):
        return self.__mul__(other)
    

    def __pow__(self, power):
        if power<0:
            if self.is_diagonal_tensor():
                Diagonal=np.array([self[*[i]*self.rank] for i in range(self.shape[0])])
                if not np.any(Diagonal==0):
                    return Tensor(np.diag(1/Diagonal))
            else:       
                raise ValueError("The tensor can't be inverted.")
        else:
            return Tensor(np.power(self.array, power))
        

    def __matmul__(self, other):
        return np.matmul(self.array, other)

    def __mod__(self, other):
        return np.mod(self.array, other)

    def __and__(self, other):
        return np.bitwise_and(self.array, other)

    def __or__(self, other):
        return np.bitwise_or(self.array, other)

    def __xor__(self, other):
        return np.bitwise_xor(self.array, other)

    def __iadd__(self, other):
        self.array += other
        return self

    def __isub__(self, other):
        self.array -= other
        return self

    def __bool__(self):
        return bool(self.array)
    
    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, self.__str__())

    def __call__(self, *args, **kwargs):
        return self.array(*args, **kwargs)
    

def DiagonalTensor(D, rank=2):
    """
    Create a diagonal tensor with the given diagonal elements.

    Parameters:
        D (list): A list of diagonal elements.
        rank (int, optional): The rank of the tensor (default is 2).

    Returns:
        Tensor: A diagonal tensor with the specified diagonal elements.

    Example:
        >>> D = [1, 2, 3]
        >>> tensor = DiagonalTensor(D, rank=2)
        >>> print(tensor)
        [1 0 0]
        [0 2 0]
        [0 0 3]
    """
    shape = [len(D)] * rank
    DT = Tensor(shape=shape)
    for i in range(len(D)):
        DT[*([i] * rank)] = D[i]
    return DT