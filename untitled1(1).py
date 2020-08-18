# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 21:25:37 2020

@author: PPC
"""

import numpy as np
import matplotlib.pyplot as plt
T = 120
S =np.zeros([T])
E =np.zeros([T])
A =np.zeros([T])
I1 =np.zeros([T])
I2 = np.zeros([T])
I3 = np.zeros([T])
R = np.zeros([T])
Sq = np.zeros([T])
Eq = np.zeros([T])
Iq = np.zeros([T])
S[0] = 1*10**5
E[0] = 109
A[0] = 10
I1[0] = 10
I2[0] = 2
I3[0] = 0
R[0] = 0
Sq[0] = 0
Iq[0] = 21
N0 = 1368000
I = I1 +I2 + I3

B = 9.08*10**-5
p = 0.50
fai = 0.50
ruo = 0.50
a = 1/7.2
n = 0.9
p1 = 0.50
p2 = 0.35
p3 = 0.15
y = yA = yq = 0.06
theta1 = 0.20
theta2 = 0.80
theta3 = 1
w = 0.1

for t in range(1,T):
    dS = -(1-p)*B*fai*(ruo*E[t-1]+I[t-1]+A[t-1])*S[t-1]-p*B*(1-fai)*(ruo*E[t-1]+I[t-1]+A[t-1])*S[t-1]-p*B*fai*(ruo*E[t-1]+I[t-1]+A[t-1])*S[t-1]+w*Sq[t-1]
    S[t] = S[t-1]+dS
    dE = (1-p)*B*fai*(ruo*E[t-1]+I[t-1]+ruo*A[t-1])*S[t-1]-a*E[t-1]
    E[t] = E[t-1]+dE
    dA = a*(1-n)*E[t-1]-yA*A[t-1]
    A[t] = A[t-1]+dA
    dI1 = a*n*p1*E[t-1]-y*I1[t-1]
    I1[t] = I1[t-1]+dI1
    dI2 = a*n*p2*E[t-1]-y*I2[t-1]
    I2[t] = I2[t-1]+dI2
    dI3 = a*n*p3*E[t-1]-y*I3[t-1]
    I3[t] = I3[t-1]+dI3
    dR = yA*A[t-1]+(1-theta1)*y*I1[t-1]+(1-theta2)*y*I2[t-1]+(1-theta3)*y*I3[t-1]
    R[t] = R[t-1]+dR
    dSq = p*B*fai*(ruo*E[t-1]+I[t-1]+A[t-1])*S[t-1]-w*Sq[t-1]
    Sq[t] = Sq[t-1]+dSq
    dEq = p*B*(1-fai)*(ruo*E[t-1]+I[t-1]+A[t-1])*S[t-1]-a*Eq[t-1]
    Eq[t] = Eq[t-1]+dEq
    dIq = a*E[t-1]+theta1*y*I1[t-1]+theta2*y*I2[t-1]+theta3*y*I3[t-1]-yq*Iq[t-1]
    Iq[t] = Iq[t-1]+dIq
    I = I1 +I2 + I3
    
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(S, c='b', lw=2, label='S')
ax.plot(E, c='orange', lw=2, label='E')
ax.plot(I, c='r', lw=2, label='I')
ax.plot(R, c='g', lw=2, label='R')
ax.set_xlabel('Day',fontsize=20)
ax.set_ylabel('Infective Ratio', fontsize=20)
ax.grid(1)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend();
