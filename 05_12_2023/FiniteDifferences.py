# -*- coding: utf-8 -*-


import numpy as np # Matrix Definitions and operations.
import scipy.sparse as sparse # Sparse Matrix Definitions and operations


class FiniteDifferences:
    

    def __init__(self,Grid):
        
        #Schemes for first-order (D)     i-1   i  i+1
        self.ForwardStencilD = np.array([-1.0,1.0,0.0])/(Grid.dx)
        self.CentralStencilD = np.array([-1.0,0.0,1.0])/(2.0*Grid.dx)
        self.BackwardStencilD= np.array([0.0,-1.0,1.0])/(Grid.dx)
        
        
        # second order (DD)              i-2  i-1  i  i+1  i+2
        #self.ForwardStencilDD = np.array([0.0,0.0,1.0,-2.0,1.0])/(Grid.dx**2)
        #self.CentralStencilDD = np.array([0.0,1.0,-2.0,1.0,0.0])/(Grid.dx**2)
        #self.BackwardStencilDD= np.array([1.0,-2.0,1.0,0.0,0.0])/(Grid.dx**2)
        self.CentralStencilDD = np.array([1.0,-2.0,1.0])/(Grid.dx**2)

        
        # Sparse Matrix Operators
        self.Identity=sparse.identity(Grid.Nx, dtype='float', format="csr")
    

        # first order
        self.DDXCentral=sparse.diags(self.CentralStencilD, [-1, 0, 1], shape=(Grid.Nx, Grid.Nx), dtype='float', format="csr")
        self.DDXCentral[-1,-3:]=self.BackwardStencilD   #Define right boundary stencil
        self.DDXCentral[0,:3]=self.ForwardStencilD      #Define left boundary stencil

        self.DDXBackward=sparse.diags(self.BackwardStencilD, [-2, -1, 0], shape=(Grid.Nx, Grid.Nx), dtype='float', format="csr" )
        self.DDXBackward[0,:3]=self.ForwardStencilD     #Define boundary stencil

        self.DDXForward=sparse.diags(self.ForwardStencilD, [0, 1, 2], shape=(Grid.Nx, Grid.Nx), dtype='float', format="csr")
        self.DDXForward[-1,-3:]=self.BackwardStencilD   #Define boundary stencil
        
        # second order
        self.D2DX2=sparse.diags(self.CentralStencilDD,  [-1, 0, 1], shape=(Grid.Nx, Grid.Nx), dtype='float', format="csr")
        self.D2DX2[-1,-3:]=self.CentralStencilDD       #Define right boundary stencil
        self.D2DX2[0,:3]=self.CentralStencilDD          #Define left boundary stencil
        

    # Efficient Implementation for 1D csr type FD matrix
    #do not Change implementation below!
    def SetDirichletLeft(self,M):
        M.data[[0,1,2]]=[1.0, 0.0, 0.0]

    def SetDirichletRight(self,M):
        M.data[[-3,-2,-1]]=[0.0, 0.0, 1.0]
    
    def SetNeumannLeft(self,M):
        M.data[[0,1,2]]=self.BackwardStencilD
    
    def SetNeumannRight(self,M):
        M.data[[-3,-2,-1]]=self.ForwardStencilD

