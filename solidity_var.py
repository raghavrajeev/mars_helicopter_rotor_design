# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt

# Density, Radius, root radius, number of blades, root chord and    taper slope respectively
rho=0.02
R=0.5
Rrc=0.2
b=3
c0=0.5
c1=0

# Velocity  
rpm=100
omega = 2*math.pi*rpm/60
V=0

# Pitch terms
E=2*math.pi/180
F=0

solidity=0.02

# Lift coefficient terms:     
a=5.75

# Drag coefficient terms:
cd0=0.0113
cd1=1.25

# Chord, solidity, pitch and induced angle respectively
def c(r):
    return (c0+c1*(r-Rrc))


def theta_f(r):
    return((E+F*(r-Rrc)))

def phi_f(r):
    return math.atan(V/(omega*r))    
    
    
# Blade element momentum theory


# tip loss factor
def lamb_f(r):
    # find lambda from initial guess F; find f from lambda; calculate F from f; compare both F
    Fp=0.5                                   
    F1 = 1                                 
    while (abs(Fp-F1) > 0.00001):   
        lda1 = (a*solidity/(16*Fp))*((((1 - (8*Fp*V)/(a*R*omega*solidity))**2 + 
                    32*Fp*theta_f(r)*r/(a*R*solidity))**0.5) - (1 - (8*Fp*V)/(a*R*omega*solidity)))
        if lda1 == 0:          
            break       
            
        F1=Fp
        f = (b/2)*(1-r/R)/lda1
        Fp = (2/math.pi)*math.acos(math.e**(-f))
        
    # converged                      
    return lda1    
    
# coefficient of lift and drag respectively
def clpr(r):
     return (a*(theta_f(r)-math.atan(lamb_f(r)*R/r)))
    
def cdpr(r):
     return (cd0+cd1*(theta_f(r)-math.atan(lamb_f(r)*R/r))**2)
    
# array of radius
r_arr=[]
for i in range(10000):
    r_arr.append(Rrc + i*(R-Rrc)/10000)


# Lambda, CL and CD at each section
ldap=[]
clp=[]
cdp=[]

for r in r_arr:
    ldap.append(lamb_f(r))
    clp.append(clpr(r))
    cdp.append(cdpr(r))
    
# empty array for thrust and torque values from BEMT
tp=[]
pp=[]

sol_array=[]
for i in range(7):
    solidity=0.02+0.01*i
    sol_array.append(solidity)
    # Thrust and Torque sectional values
    dTp=[]
    dPp=[]
    
    # bemt sectional thrust and torque
    j=0
    for r in r_arr:
        dt = 0.5*rho*b*((omega*r)**2 + (omega*R*lamb_f(r))**2)*c(r)*(clpr(r)*math.cos((math.atan(lamb_f(r)))) 
                                                                - cdpr(r)*math.sin((math.atan(lamb_f(r)))))
        dTp.append(dt)
        dp = 0.5*rho*omega*b*r*((omega*r)**2 + (omega*R*lamb_f(r))**2)*c(r)*(cdpr(r)*math.cos((math.atan(lamb_f(r)))) 
                                                                     + clpr(r)*math.sin((math.atan(lamb_f(r)))))
        dPp.append(dp)
        j = j+1
    
    # total thrust
    tp1 = ((R-Rrc)/len(r_arr))*(dTp[0]+2*sum(dTp[1:-1])+dTp[-1])/2
    #thrust coeff
    ct1 = (tp1)/(rho*math.pi*(omega**2)*(R**4))
    tp.append(ct1)
    
    # total torque        
    pp1 = ((R-Rrc)/len(r_arr))*(dPp[0]+2*sum(dPp[1:-1])+dPp[-1])/2
    # torque coeff
    cp1 = (pp1)/(rho*math.pi*(omega**2)*(R**5))
    pp.append(cp1)
    
    

# plotting
fig,ax1=plt.subplots(dpi=300)
ax1.plot(sol_array,tp,'-', label = 'BEMT')
# ax1.plot(c1_array,Ct_array, '-o', label = 'BET')
# ax1.legend()
ax1.set(xlabel='Solidity',ylabel='CT')

fig,ax2=plt.subplots(dpi=300)
ax2.plot(sol_array,pp, '-', label = 'BEMT')
# ax2.plot(pitch,Cq_array, '-o', label = 'BET')
# ax2.legend()
ax2.set(xlabel='Solidity',ylabel='CP')

    