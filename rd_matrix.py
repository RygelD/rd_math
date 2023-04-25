# RD Matrix
# A simple (bad) matrix project
# version 0.1

RDM_EXP_DEPTH = 40 # Depth for power series of exp(M), set to -1 for NotImplemented
RDM_LOG_DEPTH = 100
import math

class Matrix:
    def __init__(self,shape=(0,0),elements=[],*,matrix=None):
        self.isabadmatrix= True
        if matrix != None:
            self.shape = matrix.shape
            self.elements = matrix.elements
        else:
            if len(shape) == 2:
                if len(elements) == shape[0] * shape[1]:
                    self.elements = list(elements)
                    self.shape = tuple(shape)
                else:
                    raise ValueError("Invalid shape - wrong number of element cells")
            else:
                raise ValueError("Only 2 dimentions are allowed for shape")

    def __getitem__(self,index):
        """ A[x,y] (0-indexed) """
        return self.elements[index[1]+index[0]*self.shape[1]]

    def __setitem__(self,index,result):
        """ A[x,y] (0-indexed) """
        self.elements[index[1]+index[0]*self.shape[1]] = result
        
    def __repr__(self):
        return "Matrix(" + repr(self.shape) + "," + repr(self.elements) + ")"

    def row_en(self):
        """ Enlist rows of matrix as a nested list"""
        rows=[]
        for i in range(self.shape[0]):
            rows.append(self.elements[i*self.shape[1]:(i+1)*self.shape[1]])
        return rows
    def col_en(self):
        """ Enlist columns of matrix as a nested list"""
        return self.transpose().row_en()
    def __str__(self):
        """ Prints represntations of values inside """
        string=""
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                if j != 0:
                    string += '\t'
                string += repr(self[i,j])
            string+= '\n'
        return string[0:-1]

    def __eq__(self,matrix):
        """ True if matrix type, shape and elements are equal """
        if type(self) == type(matrix): #STRONG EQUALS
            return (self.shape == matrix.shape) and (self.elements == matrix.elements)
        else:
            return False
    def __ne__(self,matrix):
        return not self==matrix
    def __add__(self,matrix):
        return type(self)(matrix=self.apl_ech(matrix,lambda a,b: a+b)) 
    def __neg__(self):
        return type(self)(matrix=self.apl_sch(lambda a: -a)) 
    def __sub__(self,matrix):
        return type(self)(matrix=self.apl_ech(matrix,lambda a,b: a-b)) 
    def __pos__(self):
        return self.transpose()
    def transpose(self):
        elem = []
        for i in range(self.shape[1]):
            for j in range(self.shape[0]):
                elem.append(self[j,i])
        return type(self)(matrix=Matrix((self.shape[1],self.shape[0]),elem))
    def __mul__(self,oth):
        if is_Matrix(oth):
            return self @ oth
        else:
            return type(self)(matrix=self.apl_sch(lambda a: a*oth))
    def __rmul__(self,oth):
        return self*oth
    def __matmul__(self,oth):
        return self.apl_oip(oth,lambda a,b: a+b,lambda a,b: a*b)
    def __truediv__(self,oth):
        if is_Matrix(oth):
            return NotImplemented # No inverses yet...
        else:
            return self * (1/oth)
    def apl_oip(self,oth,clp,op):
        """APL operator inner product"""
        if self.shape[1] == oth.shape[0]:
            elem = []
            shape = (self.shape[0],oth.shape[1])
            x = self.row_en()
            y = oth.col_en()
            for i in x:
                for j in y:
                    subelem =[]
                    for k in range(len(i)):
                        subelem.append(op(i[k],j[k]))
                    itr = subelem[-1]
                    for k in range(len(subelem)-1):
                        itr = clp(subelem[~(1+k)],itr)
                    elem.append(itr)
            return Matrix(shape,elem)
        else:
            raise ValueError("Matrices are of incompatible shapes")
    def apl_ech(self,oth,op):
        """APL each operator of one dyadic functions and two arguments"""
        if self.shape == oth.shape:
            elem = []
            for i in range(self.shape[0]):
                for j in range(self.shape[1]):
                    elem.append(op(self[i,j],oth[i,j]))
            return Matrix(self.shape,elem)
        else:
            raise ValueError("Matrices are of different shapes")
    def apl_sch(self,op):
        """APL each operator of one monadic function and argument"""
        elem = []
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                elem.append(op(self[i,j]))
        return Matrix(self.shape,elem)
    
    def is_square(self):
        return self.shape[0] == self.shape[1]

    def __pow__(self,oth):
        if self.is_square():
            if is_Matrix(oth):
                return NotImplemented   
            else:
                if oth == 0:
                    return self.identity()
                elif oth == round(oth) and oth > 0:
                    return self * (self ** (oth - 1))
                else:
                    return NotImplemented
        else:
            return NotImplemented
    def __rpow__ (self,oth): # DEPENDS ON RDM_EXP_DEPTH
        if RDM_EXP_DEPTH != -1 and self.is_square() and not is_Matrix(oth):
            if oth == 1:
                return self.identity()
            else:
                return (self*math.log(oth)).exp()
        else:
            return NotImplemented
    def exp(self): # DEPENDS ON RDM_EXP_DEPTH
        if RDM_EXP_DEPTH != -1 and self.is_square():
            m = self*0
            for i in range(RDM_EXP_DEPTH+1):
                m += (self**i)*(1/math.factorial(i))
            return m
        else:
            return NotImplemented
    def log(self,base=math.e): # DEPENDS ON RDM_LOG_DEPTH, DOESN'T WORK WELL
        if RDM_LOG_DEPTH != -1 and self.is_square():
            m = self*0
            for i in range(1,RDM_LOG_DEPTH+2):
                m += ((-1)**(i+1))*((self-self.identity())**i)/i
            return m
        else:
            return NotImplemented
            
    def identity(self):
        if self.is_square():
            m = self * 0
            for i in range(m.shape[0]):
                m[i,i] = 1
            return m
        else:
            return NotImplemented
        
    #def zero_mat(self):
    #   return self*0


    
class ComplexM(Matrix):
    # Complex numbers as a matrix
    def __init__(self,real=0,comp=0,*,matrix=None):
        if matrix != None:
            if matrix.shape == (2,2):
                if matrix[0,0] == matrix[1,1] and matrix[0,1] == -matrix[1,0]:
                        super().__init__(matrix=matrix)
                        self.real=matrix[0,0]
                        self.comp=matrix[1,0]
                else:
                    raise ValueError("Matrix doesn't represent a complex number")
            else:
                raise ValueError("Matrix doesn't represent a complex number")
        else:
            super().__init__(shape=(2,2),elements=[real,-comp,comp,real])
            self.real=real
            self.comp=comp

    def __str__(self):
        if self.bicomp == 0:
            return str(self.real)
        elif self.real == 0:
            return str(self.bicomp) + "i"
        elif self.bicomp > 0:
            return str(self.real) + "+" + str(self.bicomp) + "i"
        else:
            return str(self.real) + "-" + str(self.bicomp) + "i"
    def __invert__(self):
        return ComplexM(self.real,-self.comp)
    def __repr__(self):
        return "ComplexM("+str(self.real)+","+str(self.comp)+")"
    def __abs__(self):
        return (real**2 + comp**2)*0.5

class BicomplexM(Matrix):
    # Bicomplex numbers as a matrix
    def __init__(self,real=0,bicomp=0,*,matrix=None):
        if matrix != None:
            if matrix.shape == (2,2):
                if matrix[0,0] == matrix[1,1] and matrix[0,1] == matrix[1,0]:
                        super().__init__(matrix=matrix)
                        self.real=matrix[0,0]
                        self.bicomp=matrix[1,0]
                else:
                    raise ValueError("Matrix doesn't represent a bicomplex number")
            else:
                raise ValueError("Matrix doesn't represent a bicomplex number")
        else:
            super().__init__(shape=(2,2),elements=[real,bicomp,bicomp,real])
            self.real=real
            self.bicomp=bicomp

    def __str__(self):
        if self.bicomp == 0:
            return str(self.real)
        elif self.real == 0:
            return str(self.bicomp) + "j"
        elif self.bicomp > 0:
            return str(self.real) + "+" + str(self.bicomp) + "j"
        else:
            return str(self.real) + "-" + str(self.bicomp) + "j"
    def __invert__(self):
        return BicomplexM(self.real,-self.bicomp)
    def __repr__(self):
        return "BicomplexM("+str(self.real)+","+str(self.bicomp)+")"

class DualM(Matrix):
    # Bicomplex numbers as a matrix
    def __init__(self,real=0,dual=0,dual_t=1,*,matrix=None):
        if matrix != None:
            if matrix.shape == (2,2):
                if matrix[0,0] == matrix[1,1] and (matrix[1,0] == 0 or matrix[0,1] == 0):
                        super().__init__(matrix=matrix)
                        self.real=matrix[0,0]
                        self.dual=matrix[1,0] or matrix[0,1]
                        self.dual_t = (matrix[1,0] == 0) - (matrix[0,1] == 0)
                else:
                    raise ValueError("Matrix doesn't represent a dual number")
            else:
                raise ValueError("Matrix doesn't represent a dual number")
        else:
            if dual_t == 1:
                super().__init__(shape=(2,2),elements=[real,dual,0,real])
            else:
                super().__init__(shape=(2,2),elements=[real,0,dual,real])
            self.real=real
            self.dual=dual
            if dual == 0:
                self.dual_t = 0
            self.dual_t=dual_t
    def __str__(self):
        if self.dual == 0:
            return str(self.real)
        elif self.real == 0:
            return str(self.dual) + ["!","q","Q"][self.dual_t]
        elif self.dual > 0:
            return str(self.real) + "+" + str(self.dual) + ["!","q","Q"][self.dual_t]
        else:
            return str(self.real) + "-" + str(self.dual) + ["!","q","Q"][self.dual_t]
    def __invert__(self):
        return DualM(self.real,-self.dual,self.dual_t)
    def __repr__(self):
        return "DualM("+str(self.real)+","+str(self.dual)+["","",",-1"][self.dual_t]+")"

def is_Matrix(x):
    try:
        a = x.isabadmatrix
        return a
    except:
        return False
