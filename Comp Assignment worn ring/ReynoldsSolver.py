# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:59:11 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""


import numpy as np # Matrix Definitions and operations.
import scipy.sparse as sparse # Sparse Matrix Definitions and operations
import scipy.sparse.linalg as linalg # Sparse Matrix Linear Algebra
import matplotlib.pyplot as plt
import VisualLib as vis
from scipy.sparse.linalg import spsolve #voor snel stelsels oplossen met veel 0 elementen


class ReynoldsSolver:
    def __init__(self,Grid,Time,Ops,FluidModel,Discretization):
        self.MaxIter=[]
        self.TolP=[]
        self.UnderRelaxP=[]
        self.SetSolver()
        self.VisualFeedbackLevel=0
        
        self.Grid=Grid
        self.Time=Time
        self.Ops=Ops
        self.FluidModel=FluidModel
        self.Discretization=Discretization
        self.VisualFeedbackLevel=0
    
    def SetSolver(self,MaxIter: int=10000,TolP: float=1.0e-5 ,UnderRelaxP: float=0.001,TolT: float=1.0e-5 ,UnderRelaxT: float=0.001, VisualFeedbackLevel: int=0):
        self.MaxIter=MaxIter
        self.TolP=TolP
        self.UnderRelaxP=UnderRelaxP
        self.TolT=TolT
        self.UnderRelaxT=UnderRelaxT
        self.VisualFeedbackLevel=VisualFeedbackLevel

 
#################
##### TO DO #####
#################       
    def SolveReynolds(self,StateVector,time): # StateVector is both in and output
        #in de statevector bevindt zich gwn alles van properties
        #onderverdeeld in tijdsstappen
        #statevector[time] is dus "de GRID (aangemaakt in main.py) waar alles uit gehaald
        # kan worden voor moment time"
        #1. reset convergence Hier gewoon een druk en temperatuursvector beginnende met stap 1

        epsP=np.zeros(self.MaxIter+1)
        epsP[0]=1.0
        epsT=np.zeros(self.MaxIter+1)
        epsT[0]=1.0

        #!!! de grid is dus gwn het 2D veld van de ring en de cylinderwand
        #elk punt in de grid heeft op dat moment een density en visc etc die verandert 
        #in de tijd en elk tijdsmoment wordt er een nieuwe grid aangemaakt met aangepaste waarden
        #dingen zoals h blijjven dus constant door de tijd maar verschillen voor elk punt in de grid
        #daarom wordt er altijd scalar vermenigvuldiging gedaan ipv matrixvermenigvuldiging

        
        #de volgende variabelen worden niet geupdatet
        #2. Predefine variables outside loop for Computational Efficiency
        DensityFunc    =self.FluidModel.Density
        ViscosityFunc  =self.FluidModel.DynamicViscosity
        SpecHeatFunc   =self.FluidModel.SpecificHeatCapacity
        ConducFunc     =self.FluidModel.ThermalConductivity      
        
        #fluidfraction??


        #hiervoor eerst die functie finitediff doen
        DDX=self.Discretization.DDXCentral
        DDXBackward=self.Discretization.DDXBackward
        DDXForward=self.Discretization.DDXForward
        D2DX2=self.Discretization.D2DX2
        SetDirichletLeft=self.Discretization.SetDirichletLeft
        SetDirichletRight=self.Discretization.SetDirichletRight
        SetNeumannLeft=self.Discretization.SetNeumannLeft
        SetNeumannRight=self.Discretization.SetNeumannRight
        
        #define your own when desired
        k=0
        I = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        UPLUSdt = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        UMINdt = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        Estart = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        PHI = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        DPHIDX = sparse.identity(self.Grid.Nx, dtype='float', format="csr")

        while ((((epsP[k]>self.TolP)) or (epsT[k]>self.TolT)) and (k<self.MaxIter)):
        
     
            #0. Calc Properties
            #nu op moment dat in de loop wordt gegaan is enkel T en p en de time daarvan gekend dus nu gaan we voor elke time step de eigenschappen 
            #opnieuw bepalen

            #rho enal ifv p,T die in statevector staan
            Density = DensityFunc(StateVector[time])
            SpecHeat = SpecHeatFunc(StateVector[time])
            Viscosity = ViscosityFunc(StateVector[time])
            Conduc = ConducFunc(StateVector[time])
            h = StateVector[time].h #gewoon filmthickness eruithalen

            # #statevector aanpassen
            # StateVector[time].Viscosity = Viscosity
            # StateVector[time].Density = Density
            # StateVector[time].SpecHeat = SpecHeat
            # StateVector[time].Conduc = Conduc


        
            #1. LHS Pressure



            #Voor LHS is phi dphi dx nodig dus begin van
            #start voor finite differences (unity grid)
            

            phi = ((h**3)/(12*Viscosity))*Density 
            #phi dus volledige grid met die waardes afh van plek

            #A is a diagonal matrix such that we need to trim phi:
            #altijd de phi grid updaten beginnend met allml eentjes
            PHI.data = phi
            DPHIDX.data = DDX @ phi #afleiden


            A = PHI @ D2DX2 #geupdatet PHI grid
            B = DPHIDX @ DDX
            M = A + B


            #RHS:
            #Snelheid is te vinden in de Ops class en wordt
            #constant geupdatet
            U = self.Ops.SlidingVelocity[time] #is speed time dependend???
            #!
            AA = (U/2) * (DDX @ (Density*h))

            hprev = StateVector[time-1].h
            Densityprev    =self.FluidModel.Density(StateVector[time-1])
            BB = (Density * h - Densityprev * hprev) / self.Time.dt
            #Time is gewoon een python class



            #zoals in LHS: M@p = RHS = AA (Reynolds) + BB (squeeze term)
            RHS = AA + BB
            #we willen pressure dus moeten we oplossen maar eerst BC zetten
   
            #3. Set Boundary Conditions Pressure with dirichlet functions that are given
            SetDirichletLeft(M) 
            SetDirichletRight(M)
            #in notes staat pcart-patm en pcc(psi)-patm???
            RHS[0] = self.Ops.AtmosphericPressure
            RHS[-1] = self.Ops.CylinderPressure[time] 



           
            
            #4. Solve System for Pressure + Update oplossen zoals pagina17 (67)
            #M*x = RHS
            pstar = spsolve(M, RHS)
            #pk is gewoon de oude druk van atm
            deltap = np.maximum(pstar,0) - StateVector[time].Pressure


            #Update pressure
            #what is theta, zeker iets met relaxationtheta???
            StateVector[time].Pressure += self.UnderRelaxP * deltap




            #5. LHS Temperature M = I+D+E
            #I eenheids uit finite diff
            #Uavg uitdrukking op p5 van taak
            uavg = -h**2 / (12 * Viscosity) * (DDX @ StateVector[time].Pressure) + self.Ops.SlidingVelocity[time] / 2
            uplus = np.maximum(uavg,0)
            umin = np.minimum(uavg,0)
            #grid aanmaken deze moeten nog uit de loop allemaal
            #UPLUSdt = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
            #UMINdt = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
            UPLUSdt.data = uplus * self.Time.dt
            UMINdt.data = umin * self.Time.dt

            D = UPLUSdt @ DDXBackward + UMINdt @ DDXForward
            #ook nog grid voor E
            #Estart = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
            Estart.data = (self.Time.dt*Conduc)/(Density * SpecHeat)
            E = -Estart @ D2DX2


            MTemp = I + D + E

            
            #6. RHS Temperature
            #formule voor Q op p5 in slides oo
            Tprev = StateVector[time-1].Temperature
            Qavg = (h**2 )/(12* Viscosity) * (DDX @ StateVector[time].Pressure)**2 + Viscosity *( U**2 /h**2 )
            RHS = Tprev + (self.Time.dt * Qavg) / (Density * SpecHeat)



            #BC:
            #vraag: mag dit met Neumann en dirichlet??


            if U <= 0:
                
                # MTemp[0,0:1] = [-1/self.Grid.dx, 1/self.Grid.dx] 
                # MTemp[0,3:] = 0
                # MTemp[-1,-1] = 1 
                # MTemp[-1, 1:-2] = 0

                #Mag dit ook??
                SetNeumannLeft(MTemp)
                SetDirichletRight(MTemp)
                
                RHS[0] = 0.0
                RHS[-1] = self.Ops.OilTemperature
            else:
                # MTemp[0,0] = 1     
                # MTemp[2:, 0] = 0
                # MTemp[-1,-2:] = [-1/self.Grid.dx, 1/self.Grid.dx]
                # MTemp[-1,1:-3] = 0

                #zelfde vraag
                SetDirichletLeft(MTemp)
                SetNeumannRight(MTemp)
                RHS[0] = self.Ops.OilTemperature
                RHS[-1] = 0.0
            #7. Solve System for Temperature + Update
            Tstar = spsolve(MTemp, RHS)

            #deltap = np.maximum(pstar,0 moet 0 vervangen worden???) - StateVector[time].Pressure
            #StateVector[time].Pressure += self.UnderRelaxP * deltap
            #voor de rest analoog als bij pressure
            deltaT = np.minimum(np.maximum(Tstar,self.Ops.OilTemperature),2.0*self.Ops.OilTemperature) - StateVector[time].Temperature
            StateVector[time].Temperature += self.UnderRelaxT * deltaT
            #Waarom wordt er vanaf 300 graden 300 genomen???
            #StateVector[time].Temperature = np.minimum(StateVector[time].Temperature, 300+273.15)
            
            #8. Calculate other quantities
            k += 1

            epsP[k] = np.linalg.norm(deltap / StateVector[time].Pressure) / self.Grid.Nx
            epsT[k] = np.linalg.norm(deltaT / StateVector[time].Temperature) / self.Grid.Nx
           
            
            # 9. Provide a plot of the solution
            if (k % 500 == 0):
                CFL=np.max(uavg)*self.Time.dt/self.Grid.dx
                print("ReynoldsSolver:: CFL", np.round(CFL,2) ,"Residual [P,T] @Time:",round(self.Time.t[time]*1000,5),"ms & Iteration:",k,"-> [",np.round(epsP[k],6),",",np.round(epsT[k],6),"]")
                if self.VisualFeedbackLevel>2:
                    fig=vis.Report_PT(self.Grid,StateVector[time]) 
                    plt.close(fig)

                
            if (epsP[k]<=self.TolP) and (epsT[k]<=self.TolT):
                print("ReynoldsSolver:: Convergence [P,T] to the predefined tolerance @Time:",round(self.Time.t[time]*1000,5),"ms & Iteration:",k,"-> [",np.round(epsP[k],6),",",np.round(epsT[k],6),"]")
                
            if k>=self.MaxIter:
                print("ReynoldsSolver:: Residual [P,T] @Time:",round(self.Time.t[time]*1000,5),"ms & Iteration:",k,"-> [",np.round(epsP[k],6),",",np.round(epsT[k],6),"]")
                print("ReynoldsSolver:: Warning: Maximum Iterations without converging to the predefined tolerance]")


            
        #11. Calculate other quantities (e.g. Wall Shear Stress, Hydrodynamic Load, ViscousFriction)
            StateVector[time].HydrodynamicLoad = np.trapz(StateVector[time].Pressure, self.Grid.x)
            WallShearStress_h = Viscosity * U/h  + (DDX @ StateVector[time].Pressure) * h /2
            StateVector[time].WallShearStress = WallShearStress_h
            StateVector[time].ViscousFriction = np.trapz(StateVector[time].WallShearStress, x=self.Grid.x)