# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 13:18:44 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""

import numpy as np
import scipy.special as special # Sparse Matrix Definitions and operations
import scipy.integrate as integral # Sparse Matrix Definitions and operations

class TriboContact:
    
    def __init__(self,Engine):
    
        self.Engine=Engine
        
        """ Equivalent Young's modulus of Hertzian contact"""
        self.YoungsModulus=1.0/((1.0-Engine.Cylinder.Material.PoissonModulus**2.0)/Engine.Cylinder.Material.YoungsModulus + (1.0-Engine.CompressionRing.Material.PoissonModulus**2)/Engine.CompressionRing.Material.YoungsModulus);
        self.Domain=np.array([-Engine.CompressionRing.Thickness/2,Engine.CompressionRing.Thickness/2])
        
        """ Roughness parameters """
        self.Roughness=np.sqrt(Engine.Cylinder.Roughness**2.0 + Engine.CompressionRing.Roughness**2.0);
        self.Zeta=97.0e9;
        self.Kappa=1.56e-6;
        self.Tau0=2.0e6;
        self.f_b=0.3;
        self.RoughnessParameter=self.Zeta*self.Kappa*self.Roughness
        self.RoughnessSlope=self.Roughness/self.Kappa
        self.Lambda_c = 2.2239
        
        """Wear Coefficients"""
        self.WearCoefficient_Cylinder=2.5e-10;
        self.WearCoefficient_CompressionRing=1.25e-10;
    
        """Geometry of the cylinder and ring"""
        self.L = 2 * np.pi * Engine.Cylinder.Radius
        self.b = Engine.CompressionRing.Thickness
        self.delta = Engine.CompressionRing.CrownHeight

    def I2(self,l): 
        I2 = (0.5*(l**2+1)*special.erfc(l/np.sqrt(2.0)) - (l/np.sqrt(2.0*np.pi))*np.exp(-l**2.0/2.0))/np.sqrt(l)
        return I2
    
    def I52(self,l):
        I52 = ((1.0/(8.0*np.sqrt(np.pi)))*np.exp(-l**2.0/4.0)*(l**(3.0/2.0))*((2.0*l**2.0+3.0)*special.kv(3.0/4.0,l**2.0/4.0)-(2.0*l**2.0+5.0)*special.kv(1.0/4.0,l**2.0/4.0)))/np.sqrt(l)
        return I52


#################
##### TO DO #####
#################

    def AsperityContact(self,StateVector,time):

        
        Lambda=StateVector[time].Lambda

        if Lambda < self.Lambda_c:  # contact is made, all formulas come from assignement
            StateVector[time].AsperityArea= np.pi**2 * self.RoughnessParameter ** 2 * self.L * np.sqrt(self.Roughness * (self.b ** 2) * 0.25 / self.delta) * integral.quad(self.I2, Lambda, self.Lambda_c,limit=100)[0]
            StateVector[time].AsperityLoad= 16/15 * np.sqrt(2) * np.pi * (self.RoughnessParameter ** 2) * np.sqrt(self.Roughness / self.Kappa) * self.YoungsModulus * np.sqrt(self.Roughness * (self.b ** 2)/(4* self.delta)) * integral.quad(self.I52, Lambda, self.Lambda_c,limit=100)[0]
            StateVector[time].AsperityFriction= self.Tau0 * StateVector[time].AsperityArea / self.L + self.f_b * StateVector[time].AsperityLoad
            StateVector[time].AsperityContactPressure= StateVector[time].AsperityLoad/StateVector[time].AsperityArea
              
        
        else: # no contact is made so no loads
            StateVector[time].AsperityArea= 0
            StateVector[time].AsperityLoad= 0
            StateVector[time].AsperityFriction= 0
            StateVector[time].AsperityContactPressure= 0

        StateVector[time].HertzianContactPressure= (np.pi / 4) * np.sqrt((StateVector[time].AsperityLoad * self.YoungsModulus) / (np.pi * self.Engine.CompressionRing.Curvature))


#################
##### TO DO #####
#################       
    def Wear(self,Ops,Time,StateVector,time):
        
        # Calculate Wear Depth on the Piston Ring  
        # p_t = StateVector[time].HertzianContactPressure
        # s_t = Ops.SlidingDistance[time - 1]
        # s_tmin = Ops.SlidingDistance[time-2] if time != 1 else 0.0
        # StateVector[time].WearDepthRing= StateVector[time-1].WearDepthRing +  self.WearCoefficient_CompressionRing * p_t/ self.Engine.CompressionRing.Material.Hardness *(s_t-s_tmin) # accumulated wear depth on the ring
        
        #alternatief
        ind=np.arange(0,time,1)
        p_t=np.array([StateVector[k].HertzianContactPressure for k in ind])
        StateVector[time].WearDepthRing=np.trapz(self.WearCoefficient_CompressionRing * p_t/ self.Engine.CompressionRing.Material.Hardness,Ops.SlidingDistance[ind])



        # Calculate The Wear Depth on the Cylinder wall
        StateVector[time].WearLocationsCylinder=np.unique(np.round(Ops.PistonPosition, 8)) # array of unique Positions where the pistion passes by  


        index=int(np.where(np.round(Ops.PistonPosition[time],8) == StateVector[time].WearLocationsCylinder)[0][0])
        DW= self.WearCoefficient_Cylinder*StateVector[time].HertzianContactPressure/self.Engine.Cylinder.Material.Hardness *np.abs(Ops.PistonVelocity[time])*Time.dt 
        StateVector[time].WearDepthCylinder[index] += DW #incremental wear depth on the positions in the array above