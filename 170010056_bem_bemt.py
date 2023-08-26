# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt

# Density, Radius, root radius, number of blades, root chord and    taper slope respectively
rho=0.02
R=0.762
Rrc=0.125
b=3
c0=0.0508
c1=0

# Velocity  
rpm=100
omega = 2*math.pi*rpm/60
V=0

# Pitch terms
E=0
F=0

# Lift coefficient terms:     
a=5.75

# Drag coefficient terms:
cd0=0.0113
cd1=1.25

# Chord, solidity, pitch and induced angle respectively
def c(r):
    return (c0+c1*(r-Rrc))

def solidity_f(r):
    return (((b*c0)/(math.pi*R))+((b*c1)/(math.pi*R))*(r-Rrc))

def theta_f(r):
    return((E+F*(r-Rrc)))

def phi_f(r):
    return math.atan(V/(omega*r))


# Blade element theory

# coefficient of lift and drag
def C_L(r):
    return a*(theta_f(r)-phi_f(r))

def C_D(r):
    return cd0+cd1*(theta_f(r)-phi_f(r))**2

# Thrust and Torque coefficients
def dT(r):
    ph=phi_f(r)
    Ut=omega*r
    Up=V
    return (b*rho*0.5*(Ut**2 + Up**2)*c(r)*(C_L(r)*math.cos(ph) - C_D(r)*math.sin(ph)))

def dQ(r):
    ph=phi_f(r)
    Ut=omega*r
    Up=V
    return (b*rho*0.5*r*(Ut**2 + Up**2)*c(r)*(C_D(r)*math.cos(ph) + C_L(r)*math.sin(ph)))    
    
    
    
# Blade element momentum theory


# tip loss factor
def lamb_f(r):
    # find lambda from initial guess F; find f from lambda; calculate F from f; compare both F
    Fp=0.5                                   
    F1 = 1                                 
    while (abs(Fp-F1) > 0.00001):        
        lda1 = (a*solidity_f(r)/(16*Fp))*((((1 - (8*Fp*V)/(a*R*omega*solidity_f(r)))**2 + 
                    32*Fp*theta_f(r)*r/(a*R*solidity_f(r)))**0.5) - (1 - (8*Fp*V)/(a*R*omega*solidity_f(r))))
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
    
    
    
    
# Calculating and storing values from BET

# empty array for thrust and torque values from BET
Ct_array=[]
Cq_array=[]

for i in range(7):
    E=((2)*math.pi/180)*i
    
    # array of sectional thrust and torque coefficient
    dT_array=[dT(i) for i in r_arr] 
    dQ_array=[dQ(i) for i in r_arr]

    dr=(R-Rrc)/10000 #section length

    T=0
    Q=0

    # total thrust and torque
    for i in range(len(r_arr)-1):
        Ti=(dT_array[i]+dT_array[i+1])/2
        T+=Ti*dr
        Qi=(dQ_array[i]+dQ_array[i+1])/2
        Q+=Qi*dr

    # coefficient of thrust and torque
    Ct=T/(rho*math.pi*(R**2)*(omega*R)**2)
    Cq=Q/(rho*math.pi*(R**3)*(omega*R)**2)
    
    Ct_array.append(Ct)
    Cq_array.append(Cq)
    
    


# calculating and storing values from BEMT


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
qp=[]


for i in range(7):
    E=((2)*math.pi/180)*i
    
    # Thrust and Torque sectional values
    dTp=[]
    dQp=[]
    
    # bemt sectional thrust and torque
    j=0
    for r in r_arr:
        dt = 0.5*rho*b*((omega*r)**2 + (omega*R*lamb_f(r))**2)*c(r)*(clpr(r)*math.cos((math.atan(lamb_f(r)))) 
                                                                - cdpr(r)*math.sin((math.atan(lamb_f(r)))))
        dTp.append(dt)
        dq = 0.5*rho*b*r*((omega*r)**2 + (omega*R*lamb_f(r))**2)*c(r)*(cdpr(r)*math.cos((math.atan(lamb_f(r)))) 
                                                                     + clpr(r)*math.sin((math.atan(lamb_f(r)))))
        dQp.append(dq)
        j = j+1
    
    # total thrust
    tp1 = ((R-Rrc)/len(r_arr))*(dTp[0]+2*sum(dTp[1:-1])+dTp[-1])/2
    #thrust coeff
    ct1 = (tp1)/(rho*math.pi*(omega**2)*(R**4))
    tp.append(ct1)
    
    # total torque        
    qp1 = ((R-Rrc)/len(r_arr))*(dQp[0]+2*sum(dQp[1:-1])+dQp[-1])/2
    # torque coeff
    cq1 = (qp1)/(rho*math.pi*(omega**2)*(R**5))
    qp.append(cq1)
    
    
    
# experimental values from paper   
pitch=[0,2,4,6,8,10,12]
ctpaper=[0,0.00051,0.00149,0.00274,0.004165,0.005625,0.00685]
cqpaper=[0.0000925,0.000103,0.00015,0.000237,0.000367,0.000528,0.000678]

# plotting
fig,ax1=plt.subplots(dpi=300)
ax1.plot(pitch,tp,'-o', label = 'BEMT')
ax1.plot(pitch,ctpaper,'-o', label = 'Experimental')
ax1.plot(pitch,Ct_array, '-o', label = 'BET')
ax1.legend()
ax1.set(xlabel='Pitch (deg)',ylabel='CT')

fig,ax2=plt.subplots(dpi=300)
ax2.plot(pitch,qp, '-o', label = 'BEMT')
ax2.plot(pitch,cqpaper, '-o', label = 'Experimental')
ax2.plot(pitch,Cq_array, '-o', label = 'BET')
ax2.legend()
ax2.set(xlabel='Pitch ( deg)',ylabel='CQ')

    