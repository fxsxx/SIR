import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint as ode

def beta_H(t):
    beta_H = phi*np.sin(t*np.pi)**8
    return beta_H

def beta_M(t):
    beta_M = phi*np.sin(t*np.pi)**8
    return beta_M

def model(y,t,param):
    S_H,I_H1,I_H2,R_H,S_M,I_M = y
    mu_H,sigma_H,mu_M = param
    N_H = S_H + I_H1 + I_H2 + R_H
    N_M = S_M + I_M
    
    #Human Population
    dS_H = mu_H*N_H - beta_H(t)*I_M*S_H/N_M - mu_H*S_H
    dI_H1 = 0.9*beta_H(t)*I_M*S_H/N_M - sigma_H*I_H1 - mu_H*I_H1
    dI_H2 = 0.1*beta_H(t)*I_M*S_H/N_M - 0.1*sigma_H*I_H2 - mu_H*I_H2
    dR_H = sigma_H*I_H1 + 0.1*sigma_H*I_H2 - mu_H*R_H
    
    #Mosquito Population
    dS_M = mu_M*N_H - beta_M(t)*(I_H1 + 0.1*I_H2)*S_M/N_H - mu_M*S_M
    dI_M = beta_M(t)*(I_H1 + 0.1*I_H2)*S_M/N_H - mu_M*I_M
    
    return dS_H,dI_H1,dI_H2,dR_H,dS_M,dI_M

basepop = 10000 #sets baseline population size
relpop = 10 #determines relative population size of mosquitos to humans
phi = 36.5 #36.5 ~ 1 bite/day * 10% transmission probability

#Human parameters
mu_H = 1.0/75
sigma_H = 26 #equates to two week infection length

#Mosquito parameters
mu_M = 1.0/75

#Initial state
yinit = [0.01*basepop-10,5,5,0.99*basepop,(relpop*basepop)*0.9,(relpop*basepop)*0.1]

param = [mu_H,sigma_H,mu_M]
t=np.linspace(0,10,200)

y=ode(model,yinit,t,(param,))

fig, ax = plt.subplots(3,1,sharex=True)
ax[0].plot(t,y[:,0],'g')
ax[0].set_ylabel('susceptible')
ax[1].plot(t,y[:,1],'r',label='acute')
ax[1].plot(t,y[:,2],'r--',label='chronic')
ax[1].set_ylabel('infected')
ax[1].legend()
ax[2].plot(t,y[:,3],'b')
ax[2].set_ylabel('recovered')
ax[2].set_xlabel('time')