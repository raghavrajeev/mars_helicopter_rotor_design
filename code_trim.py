# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt

#  Density, Radius, root radius, number of blades, chord and mass respectively
rho=0.02
R=0.5
Rrc=0.0525
b=3
c=0.0508
m=2
g=3.72

# RPM, angular velocity and vertical velocity
rpm=2000
omega = 2*math.pi*rpm/60

# Lift coefficient terms:
a=8

# # Drag coefficient terms:
# cd0=0.0113
# cd1=1.25

# Coefficient of drag (fuselage) and wetted area
Cd_f = 0.4
S=0.9

# Forward velocity
V_inf=30 #100 km/hr = 27.7 m/s

# alpha_tpp calculation
drag = 0.5 * rho * (Cd_f) * V_inf**2 * S
weight=m*g
alpha_tpp = math.atan(drag/weight)

# Thrust, thrust coefficient and advance ratio calculation
T = weight/math.cos(alpha_tpp)
C_T = T/(rho*math.pi*(R**2)*(omega*R)**2)
mu=(V_inf*math.cos(alpha_tpp))/(omega*R)


# Inflow model
lamb_ig = C_T/(mu*2)
lamb_g =lamb_ig+(V_inf*math.sin(alpha_tpp)/(omega*R))

def inflow(r, psi): # non-uniform main rotor inflow
    lamb_i = lamb_ig*(1+(((4*mu)/3)/(mu+1.2*(lamb_g)))*(r/R)*math.cos(psi))
    lamb_i = lamb_i + ((V_inf*math.sin(alpha_tpp))/(omega*R))
    return (lamb_i)

# Angles (Pitch and flap)
# collective pitch, lateral and longitudinal cyclic
theta_0 = 27
theta_tw = 12
theta_1c = 5.61
theta_1s = -12.43
# flap angle, longitudinal and lateral tilt
beta_0 = 24

# converting to radians
theta_0 = theta_0*math.pi/180
theta_tw = theta_tw*math.pi/180
theta_1s = theta_1s*math.pi/180
theta_1c = theta_1c*math.pi/180
beta_0 = beta_0*math.pi/180

# Function to calculate Pitch values. 
def pitch(r, psi):
    return (theta_0 + theta_tw*(r/R) + theta_1c*math.cos(psi) + theta_1s*math.sin(psi))

# Angle of attach calculation
def aoa(r, psi):
    U_P = inflow(r,psi)*omega*R+V_inf*math.cos(psi)*math.sin(beta_0)
    U_T = omega*r + V_inf*math.sin(psi)*math.cos(alpha_tpp)
    return (pitch(r, psi) - math.atan(U_P/U_T))

rpts=[]
for i in range(2000):
    rpts.append(Rrc + i*(R-Rrc)/2000)
#Dividing rotor blade into 2000 segmets

psipts=[]
for i in range(360):
    psipts.append(i*math.pi/180)
# Dividing Azimuthal angles into 360 components


#Total Thrust Calculation
def thrust_func():
    total_thrust=0
    thrust_psi = []
    for psi in psipts:
        thrust_r = []
        for r in rpts:
            U_P = inflow(r,psi)*omega*R*+V_inf*math.cos(psi)*math.sin(beta_0)
            U_T = omega*r + V_inf*math.sin(psi)*math.cos(alpha_tpp)
            lift1 = 0.5*rho*c*a*aoa(r, psi)*(U_T**2 + U_P**2)
            thrust_r.append(np.sign(U_T)*lift1)
        lift_si = ((R-Rrc)/len(rpts))*(thrust_r[0]+2*sum(thrust_r[1:-1])+thrust_r[-1])/2
        # Integration of sectional lift at diff radial locations to get lift at azimuthal angle 
        thrust_psi.append(lift_si)
    total_thrust = b*(1/len(psipts))*(thrust_psi[0]+2*sum(thrust_psi[1:-1])+thrust_psi[-1])/(2*math.pi)
    # Integration of thrust at diff azimuthal locations to get total thrust 
    return(total_thrust)

print(thrust_func())

#Pitching Moment Calculation
def pm_func():
    pm=0
    pm_psi = []
    for psi in psipts:
        pm_r = []
        for r in rpts:
            U_P = inflow(r,psi)*omega*R*+V_inf*math.cos(psi)*math.sin(beta_0)
            U_T = omega*r + V_inf*math.sin(psi)*math.cos(alpha_tpp)
            lift1 = 0.5*rho*c*a*aoa(r, psi)*(U_T**2 + U_P**2)
            pm_r.append(np.sign(U_T)*lift1*r*math.cos(psi))
        pm1_psi = ((R-Rrc)/len(rpts))*(pm_r[0]+2*sum(pm_r[1:-1])+pm_r[-1])/2
        # Integration of PM at diff radial locations to get PM at azimuthal angle 
        pm_psi.append(pm1_psi)
    pm = b*(1/len(psipts))*(pm_psi[0]+2*sum(pm_psi[1:-1])+pm_psi[-1])/(2*math.pi)
    # Integration of PM at diff azimuthal locations to get total PM 
    
    return(pm)

print(pm_func())

#Rolling Moment Calculation
def rm_func():
    rm_psi = []
    for psi in psipts:
        rm_r = []
        for r in rpts:
            U_P = inflow(r,psi)*omega*R*+V_inf*math.cos(psi)*math.sin(beta_0)
            U_T = omega*r + V_inf*math.sin(psi)*math.cos(alpha_tpp)
            lift1 = 0.5*rho*c*a*aoa(r, psi)*(U_T**2 + U_P**2)
            rm_r.append(np.sign(U_T)*lift1*r*math.sin(psi))
        rm1_psi = ((R-Rrc)/len(rpts))*(rm_r[0]+2*sum(rm_r[1:-1])+rm_r[-1])/2
        # Integration of PM at diff radial locations to get PM at azimuthal angle 
        rm_psi.append(rm1_psi)
    rm = b*(1/len(psipts))*(rm_psi[0]+2*sum(rm_psi[1:-1])+rm_psi[-1])/(2*math.pi)
    # Integration of PM at diff azimuthal locations to get total PM 
    
    return(rm)

print(rm_func())