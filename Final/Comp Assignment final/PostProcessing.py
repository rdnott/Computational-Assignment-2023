#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:49:17 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""



"""
Import libraries
"""
import numpy as np # Matrix Definitions and operations.
import scipy.integrate as integral # Sparse Matrix Definitions and operations
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon



from EngineParts import Engine  #import all classes from file
from TriboContact import TriboContact #import all classes from file
from Grid import Grid #import all classes from file
from Time import Time #import all classes from file
from Ops import Ops #import all classes from file
from FluidLibrary import Liquid,Gas #import all classes from file
from TwoPhaseModel import CavitationModel 
from SolutionState import State #import all classes from file
from FiniteDifferences import FiniteDifferences
from IOHDF5 import IOHDF5
import VisualLib as vis

# from IPython import get_ipython
# get_ipython().magic('reset -sf')


"""I/O Operator"""
IO=IOHDF5()

""" Input Parameters"""
EngineType='VW 2.0 R4 16v TDI CR 103kW'
OilTemperature=95.0 #C
EngineRPM=2400.0 #rpm
EngineAcceleration=0.0;

""" Define Engine Geometry"""
Engine=Engine(EngineType)

"""Define Dry Contact parameters"""
Contact=TriboContact(Engine)


"""1D Computational Grid"""
Nodes=256;
Grid=Grid(Contact,Nodes)

"""Temporal Discretization"""
TimeStep=5e-5 # Choose Temperal Resolution 
EndTime=4.0*np.pi/(EngineRPM*(2.0*np.pi/60.0))
Time=Time(EndTime,TimeStep)

"""Define Operational Conditions""" 
Ops=Ops(Time,Engine,EngineRPM,EngineAcceleration,OilTemperature)


"""Define Two-Phase Lubricant-Vapour flow"""
Oil=Liquid('SAE5W40')
Vapour=Gas('SAE5W40')
Mixture=CavitationModel(Oil,Vapour)


"""Define the State Vector = List of All States over time"""
StateVector=[]
for t in range(Time.nt):
    StateVector.append(State(Grid))

""" "Spatial "Discretization by Finite Differences """
Discretization=FiniteDifferences(Grid)



"""Read Data"""
time=0
for time in range(1,Time.nt):
    FileName='Comp Assignment final\Data\Time_'+str(round(Time.t[time]*1000,4))+'ms.h5' 

    Data=IO.ReadData(FileName)
    StateVector[time].h0=float(Data['State']['h0'])
    StateVector[time].Lambda=float(Data['State']['Lambda'])
    StateVector[time].HydrodynamicLoad=float(Data['State']['HydrodynamicLoad'])
    StateVector[time].ViscousFriction=float(Data['State']['ViscousFriction'])
    StateVector[time].AsperityLoad=float(Data['State']['AsperityLoad'])
    StateVector[time].AsperityFriction=float(Data['State']['AsperityFriction'])
    StateVector[time].AsperityContactArea=float(Data['State']['AsperityContactArea'])
    StateVector[time].AsperityContactPressure=float(Data['State']['AsperityContactPressure'])
    StateVector[time].HertzianContactPressure=float(Data['State']['HertzianContactPressure'])
    StateVector[time].COF=float(Data['State']['COF'])
    StateVector[time].WearDepthRing=float(Data['State']['WearDepthRing'])
    StateVector[time].Viscosity=Data['State']['Viscosity']
    StateVector[time].VapourVolumeFraction=Data['State']['VapourVolumeFraction']
    StateVector[time].Density=Data['State']['Density']
    StateVector[time].SpecHeat = Data['State']['SpecHeat']
    StateVector[time].Conduc = Data['State']['Conduc']

    StateVector[time].Hersey = abs(Ops.SlidingVelocity[time]) * StateVector[time].Viscosity / ( StateVector[time].HydrodynamicLoad + StateVector[time].AsperityLoad)

    StateVector[time].h= Data['State']['h']
    StateVector[time].Pressure=Data['State']['Pressure']
    StateVector[time].Temperature=Data['State']['Temperature']
    StateVector[time].WallShearStress=Data['State']['WallShearStress']
    StateVector[time].WearLocationsCylinder=Data['State']['WearLocationsCylinder']
    StateVector[time].WearDepthCylinder=Data['State']['WearDepthCylinder']
    
    time+=1
    
 
    

"""Post-Processing"""
#################
##### TO DO #####
#################   

interesting_points = np.array([1, np.argmax(Ops.SlidingVelocity),np.argmin(Ops.SlidingVelocity), 500,  999-np.argmax(Ops.CompressionRingLoad[::-1]),718, 999])
#vis.Report_Ops(Time, Ops, interesting_points)
plt.close()

vline = (interesting_points/1000*4*np.pi).tolist()

# for i,e in enumerate(interesting_points):
#     print('time=%s, angle=%s' % (Time.t[e]*1000, interesting_points[i-1]/1000*4*np.pi))

#i Dimensionless film thickness (Lambda) evolution a.f.o. crank angle (psi) 

Lambdas = np.zeros(Time.nt-1)

for t in range(Time.nt-1):
    Lambdas[t] = StateVector[t].Lambda
color = 'tab:blue'
plt.plot(Ops.CranckAngle[1:], Lambdas, 'b-')
plt.xlabel('Crank angle ($\psi$) [rad]')
plt.ylabel('Dimensionless film thickness ($\Lambda$) [-]')
plt.xlim([-.5, 13.066])
plt.ylim([0,40])
psi = np.arange(0, 4 * np.pi + np.pi/2, step=(np.pi/2))
plt.xticks(psi,['0','π/2', 'π', '3π/2', '2π','5π/2', '3π', '7π/2', '4π'])

dot_x = []
dot_y = []
label_offset = 0.3

plt.hlines(y=1, xmin=-0.5,xmax=13.066, color='k',linestyles='--')
plt.hlines(y=2.5, xmin=-0.5,xmax=13.066, color='k',linestyles='--')

for i,line in enumerate(vline):
    plt.plot(line, Lambdas[interesting_points-1][i],'o', color=color)
    dot_x.append(line)
    dot_y.append(Lambdas[interesting_points-1][i])

for i, (x, y) in enumerate(zip(dot_x, dot_y)):
    label = f' {i+1}'
    plt.text(x + label_offset, y + label_offset, label, color='black', fontsize=12, ha='right', va='bottom')
    


plt.legend()
plt.show()
plt.close()

#ii The Stribeck curve, displaying the coefficient of friction versus the Hersey number

COFs = np.zeros(Time.nt-1)
Herseys = np.zeros(Time.nt-1)

for t in range(Time.nt-1):
    Herseys[t] = abs(np.mean(StateVector[t].Hersey))  
    COFs[t] = abs(StateVector[t].COF)

plt.plot(Herseys, COFs, 'b-')
plt.xlabel('Hersey number [-]')
plt.ylabel('Coefficient of Friction [-]')
plt.show()
plt.close()

#iii Characteristic Pressure & Temperature fields at interesting and relevant locations





def Report_PT1(Grid,State): # initiatlization

    f1, ax1 = plt.subplots()
    color = 'tab:blue'
    ax1.set_xlabel('$x [mm]$ at location %s (%.2f ms)' %(np.where(interesting_points == time)[0][0]+1, time/1000*50))
    ax1.set_ylabel('$P [MPa]$',color=color)
    ax1.plot(Grid.x[1:-2]*1000,State.Pressure[1:-2]/1e6,'x-', linewidth=1,color=color)
    ax1.tick_params(axis='y')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('$T [^\circ C]$',color=color)  # we already handled the x-label with ax1
    ax2.plot(Grid.x[1:-2]*1000,State.Temperature[1:-2]-273.15,'x-', linewidth=1,color=color)
    ax2.tick_params(axis='y')
    f1.tight_layout() # otherwise the right y-label is slightly clipped
    plt.show()
    return f1


for time in interesting_points:
    fig=Report_PT1(Grid,StateVector[time])
    plt.close(fig)

#iv Vapour Volume Fraction, Viscosity, Density,... at relevant locations
for time in interesting_points:

    plt.plot(Grid.x*1000, StateVector[time].VapourVolumeFraction)
    plt.ylabel( 'Vapour Volume Fraction '+str(chr(945)) + ' [-]')
    plt.xlabel('x [mm]')
    plt.title('Location %s (%.2f ms)' %(np.where(interesting_points == time)[0][0]+1, time/1000*50))
    plt.show()
    plt.close()



    plt.plot(Grid.x*1000, StateVector[time].Density)
    plt.ylabel( 'Density '+str(chr(961)) + '  [kg/m³]')
    plt.xlabel('x [mm]' )
    plt.title('Location %s (%.2f ms)' %(np.where(interesting_points == time)[0][0]+1, time/1000*50))
    plt.show()
    plt.close()
    
    plt.plot(Grid.x*1000, StateVector[time].Viscosity)
    plt.ylabel( 'Viscosity '+str(chr(956)) + '  Pa s')
    plt.xlabel('x [mm]')
    plt.title('Location %s (%.2f ms)' %(np.where(interesting_points == time)[0][0]+1, time/1000*50))
    plt.show()
    plt.close()

    plt.plot(Grid.x*1000, StateVector[time].SpecHeat)
    plt.ylabel( 'Specific Heat Capacity c [J/(K*kg)]')
    plt.xlabel('x [mm]')
    plt.title('Location %s (%.2f ms)' %(np.where(interesting_points == time)[0][0]+1, time/1000*50))
    plt.show()
    plt.close()    

    plt.plot(Grid.x*1000, StateVector[time].Conduc)
    plt.ylabel( 'Thermal Conductivity '+str(chr(954)) + ' [W/(m*K)]')
    plt.xlabel('x [mm]')
    plt.title('Location %s (%.2f ms)' %(np.where(interesting_points == time)[0][0]+1, time/1000*50))
    plt.show()
    plt.close()       

#v Velocity Field

    j = 1

for time in interesting_points:

    visc_x = StateVector[time].Viscosity
    density = Mixture.Density(StateVector[time])
    p = StateVector[time].Pressure
    p_x = Discretization.DDXCentral @ p
    DDX = Discretization.DDXCentral
    p_y = 0.0
    x_grid = Grid.x
    z_n = Grid.Nx
    h = StateVector[time].h
    
    u1 = 0
    u2 = Ops.SlidingVelocity[time] 
    Nx = Grid.Nx
    u_x = np.zeros((z_n, Nx))
    u_z = np.zeros((z_n, Nx))

    z_grid = np.linspace(0.0, h[0], z_n)
    for x in range(Grid.Nx):
        u_x[:, x] = np.array([1/visc_x[x] * p_x[x] * 1/2 * (z**2 - h[x] * z) if z < h[x] else 0 for z in z_grid])  # Poiseuille
        u_x[:, x] += np.array([(u2 - u1) / h[x] * z + u1 if z < h[x] else 0 for z in z_grid])                      # Couette
        
    skip = 10
    scale = 1
    if time in np.array([1, 500, 999]):
        skip = 12
        scale = 40
    skip1 = (slice(None, None, skip))
    skip2 = (slice(None, None, skip), slice(None, None, skip))

    ## Make grid for vectors
    X, Z = np.meshgrid(x_grid[skip1], z_grid[skip1])

    X_l, Z_l = np.meshgrid((x_grid[:Nx//2])[skip1], (z_grid)[skip1])
    X_r, Z_r = np.meshgrid((x_grid[Nx//2+1:])[skip1], (z_grid)[skip1])


    ## Make vector plot
    pts = [[-0.75,0.02]]
    for i in range(len(x_grid)):
        pts.append([x_grid[i]*1000,StateVector[time].h[i]*1000])
    pts.append([0.75,0.02])


    plt.quiver(X*1000,Z*1000,u_x[skip2]*scale,u_z[skip2],pivot='tail',minlength=0,scale=350) #scale=350 for v!=0

    plt.plot(x_grid*1000, StateVector[time].h*1000)
    plt.title('Location %s (%.2f ms)' %(np.where(interesting_points == time)[0][0]+1, time/1000*50))
    plt.xlabel('x [mm]')
    plt.ylabel('z [mm]')
    plt.xlim([-.8,.8])
    plt.ylim([0.0, 0.0175])
    plt.tight_layout()
    plt.show()
    plt.close()
    j += 1

#vi Wear of Compression Ring and Wear at Cylinder liner after 1 combustion cycle.

WearDepthRing_values = np.zeros(Time.nt - 1)
WearDepthCylinder_values = np.zeros(len(interesting_points))

for time in range(Time.nt - 1):
    WearDepthRing_values[time] = StateVector[time].WearDepthRing

plt.plot(Ops.CranckAngle[1:], WearDepthRing_values, '-',markersize=3)
plt.xlabel('Crank angle [rad]')

plt.ylabel('Weardepth ring [mm]')
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','π/2', 'π', '3π/2', '2π','5π/2', '3π', '7π/2', '4π'])
plt.show()
plt.close()


# for time in interesting_timestamps:
time = 998 # We are only interested in wear after a full combustion cycle
plt.plot(StateVector[time].WearLocationsCylinder*1000 - 95.5, StateVector[time].WearDepthCylinder, '-',markersize=3)
plt.xlabel('Location on cylinder liner [mm]')
plt.ylabel('Wear depth [m]')
plt.show()
plt.close()

print('Maximim wear depth on cylinder sleeve = ' + str(max(StateVector[time].WearDepthCylinder)))


# lifetime compression ring

WearDepth_one_comb_cycle = StateVector[998].WearDepthRing # constant wear rate assumed
print('Wear depth compression ring: ' + str(StateVector[999].WearDepthRing)+ 'm')
reduction = 0.2 * Engine.CompressionRing.CrownHeight
nr_comb_cycles = reduction / WearDepth_one_comb_cycle
rot = nr_comb_cycles * 2 #  1 combustion cycle = 2 rotations
km = rot / 1200 # 120 km/h @ 2400 rpm --> 1 km/30s --> 1km = 1200 rot
print('aantal km= ', km)

# lifetime cylinder liner

max_one_comb_cycle = np.max(StateVector[998].WearDepthCylinder)
nr_comb_cycles_2 = 0.000002 / max_one_comb_cycle
rot2 = nr_comb_cycles_2 * 2
km2 = rot2 / 1200
print('aantal km cylinder liner = ', km2)

