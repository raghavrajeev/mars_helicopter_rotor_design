# -*- coding: utf-8 -*-

import math
import numpy as np

# Density, Radius, root radius, number of blades, chord and mass respectively
rho=0.02
R=0.15
Rrc=0.0225
b=3
c=0.008
g=3.72

# RPM, angular velocity and vertical velocity
rpm=2000
omega = 2*math.pi*rpm/60

# Lift coefficient terms:
a=8

# Forward velocity
V_inf=30 #100 km/hr = 27.7 m/s

# alpha_tpp calculation
alpha_tpp = 0

# Thrust, thrust coefficient and advance ratio calculation
T = 0.282/0.5
C_T = T/(rho*math.pi*(R**2)*(omega*R)**2)
mu=(V_inf*math.cos(alpha_tpp))/(omega*R)


# Inflow model
lamb_ig = C_T/(mu*2)
lamb_g =lamb_ig+(V_inf*math.sin(alpha_tpp)/(omega*R))

def inflow(r, psi): # non-uniform main rotor inflow
    lamb_i = lamb_ig*(1+(((4*mu)/3)/(mu+1.2*(lamb_g)))*(r/R)*math.cos(psi))
    lamb_i = lamb_i + ((V_inf*math.sin(alpha_tpp))/(omega*R))
    return (lamb_i)

# Angles (Pitch)
theta_0 = 25

# converting to radians
theta_0 = theta_0*math.pi/180

def aoa(r, psi):
    U_P = inflow(r,psi)*omega*R
    U_T = omega*r + V_inf*math.sin(psi)*math.cos(alpha_tpp)
    return (theta_0 - math.atan(U_P/U_T))

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
            U_P = inflow(r,psi)*omega*R
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