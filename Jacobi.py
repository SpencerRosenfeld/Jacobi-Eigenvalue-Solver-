from math import *

class Matrix :
    def __init__(self,array):
        self.array = array
    def dimensions(self):
        return [len(self.array),len(self.array[0])]
    def __add__(self,Other) :
        m = self.dimensions()[0]
        n = self.dimensions()[1]
        Sum = list()
        for i in range(m):
            row = list() 
            for j in range(n) :
                row.append(self.array[i][j] + Other.array[i][j])
            Sum.append(row)
        return Matrix(Sum) 
    def __mul__(self,Other) : 
        m = self.dimensions()[0]
        n = Other.dimensions()[1]
        Product = list() 
        for i in range(m) :
            row = list() 
            for j in range(n) :
                sum = 0
                for s in range(self.dimensions()[1]) :
                    sum += self.array[i][s] * Other.array[s][j]
                row.append(sum)
            Product.append(row)
        return Matrix(Product)
    def transpose(self) :
        T = list()
        m = self.dimensions()[1]
        n = self.dimensions()[0]
        for i in range(m) :
            row = list() 
            for j in range(n) :
                row.append( self.array[j][i] )
            T.append(row)           
        return Matrix(T)
    def __repr__(self) :
        m = self.dimensions()[0]
        n = self.dimensions()[1]
        string = ""
        for i in range(m):
                for j in range(n) :
                    string += str(self.array[i][j]) + " "
                string += "\n"
        return string
    def col(self, i) :
        L = []
        for row in range(self.dimensions()[0]) :
            L.append([self.array[row][i]])
        return Matrix(L)
    def row( self, i) :
        return self.transpose().col(i).transpose()
        
def ematrix(x,i,j,n) :
    L = list() 
    for a in range(n):
        row = list() 
        for b in range(n):
            row.append(0)
        L.append(row)
    L[i][j]=x
    return Matrix(L)

def ident(n) : 
    L = list() 
    for a in range(n):
        row = list() 
        for b in range(n):
            row.append(0)
        L.append(row)
    for i in range(n) :
        L[i][i] = 1
    return Matrix(L)

def G(i,j,theta,n) :
    return ematrix(sin(theta),i,j,n) + ematrix(-sin(theta),j,i,n) + ematrix(cos(theta)-1,i,i,n) + ematrix(cos(theta)-1,j,j,n) + ident(n)

def K(M,i,j,n):
    if abs(M.array[j][j] - M.array[i][i]) > 1E-8 : 
        return G(i,j,0.5*atan(2*(M.array[i][j])/(M.array[j][j]-M.array[i][i])),n).transpose()*M*G(i,j,0.5*atan(2*(M.array[i][j])/(M.array[j][j]-M.array[i][i])),n)
    else :
        return G(i,j,0.5*pi/2,n).transpose()*M*G(i,j,0.5*pi/2,n)
        
def Jacobi(M) :
    n = M.dimensions()[0]
    V = ident(n)
    max = abs(M.array[0][1])
    imax = 0
    jmax = 1
    iteration_count = 0 
    for i in range(n) :
        for j in range(i) :
            if abs(M.array[i][j]) > max :
                max = abs(M.array[i][j] )
                imax = i
                jmax = j
    Result = Matrix(M.array) 
    while max > 1E-10 and iteration_count < 20 :
        iteration_count +=1 
        if( abs(Result.array[jmax][jmax] - Result.array[imax][imax]) > 1E-8 ) :
            V = V * G(imax,jmax,0.5*atan(2*(Result.array[imax][jmax])/(Result.array[jmax][jmax]-Result.array[imax][imax])),n)
        else :
            V = V * G(imax,jmax,pi/4,n)
        Result = K(Result,imax,jmax,n)
        max = abs(Result.array[0][1])
        imax = 0
        jmax = 1
        for i in range(n) :
            for j in range(i) :
                if abs(Result.array[i][j]) > max :
                    max = abs(Result.array[i][j] )
                    imax = i
                    jmax = j
    return Result, V 
        
