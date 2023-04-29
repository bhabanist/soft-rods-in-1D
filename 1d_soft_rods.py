#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 17:46:17 2022

@author: bhabanist
"""
import time
# here i am putting all my libraries that i need
import numpy as np
import matplotlib.pyplot as plt  # for plotting
# here i am defining all my parameters
t1 = time.time()

density = 0.2
num_particles = 50
box_length= int(num_particles/density)
    



Tmax = 5
dt = 0.1
sig =1
t = np.arange(0,Tmax,dt)
# lets initialize the position.
# I will consider one point of rod which is its centre of mass.
# lets make a array of position of 100 particle
x = np.zeros(num_particles)
x = np.linspace(0,box_length,num_particles)

# Lets initialize the velocity remeber that we have to use velocity from 
# Maxwell Boltzmann distribution so we will use Box_muller method
v = np.zeros(num_particles)    # Making an array of 
v1 = np.random.uniform(size=(num_particles))
v2 = np.random.uniform(size=(num_particles))
#Lets use Box Muller algorithm
v = 0.1*np.sqrt(-2*np.log(v1))*np.cos(2*np.pi*v2)
#print(x)
# Initialization of position  and velocity is done.

# Now we have to make an array of force. As how each particle is affected by other particle
# %%
def force(x):
    force= np.zeros(num_particles)
    for i in range(num_particles):
        for j in range(num_particles):
            r = x[i]-x[j]
            
            mr =abs(x[i]-x[j])
            
            if mr< sig and x[i] != x[j]:
                force[i] = (0.5/sig)*(1.0- (mr/sig))*(r/mr)
                
                #if i ==0 or i == num_particles-1:
                  #  force[i] = -12*(sig/r**13)
                #else:
                 #   force[i] = (0.5/sig)*(1.0- (mr/sig))*(r/mr)
                    
                    
                
            
            else:
                force[i]= 0
    return force
F = force(x)
def PE(x):
    pot= np.zeros(num_particles)
    for i in range(num_particles):
        for j in range(num_particles):
            
            mr =abs(x[i]-x[j])
            
            if mr< sig and x[i] != x[j]:
                pot[i] = 0.5*(1.0- (mr/sig))**2
                #if i ==0 or i == num_particles-1:
                    #pot[i] = sig/r**12 
                #else:
                   # pot[i] = 0.5*(1.0- (mr/sig))**2
            else:
                pot[i]= 0
    return pot
pe = PE(x)
# intialization of force is done lets check verlet algorithm
def KE(v):
    K_E = np.zeros(num_particles)
    K_E = 0.5*v**2
    return K_E
ke = KE(v)



def TE(pe,ke):
    T_E = np.zeros(num_particles)
    for i in range(num_particles):
        T_E[i] =pe[i] + ke[i]
    return T_E
te = TE(pe,ke)
t_final = int(Tmax/dt)
energy_time = np.zeros(t_final)
velocity_time = np.zeros(t_final)
KE_time = np.zeros(t_final)
PE_time = np.zeros(t_final)


z = np.zeros(num_particles)
z = np.linspace(0,box_length,num_particles)
def vel_verlet(x,v,dt,Tmax):
    t = 0
    
    x_new = x
    for t in range(0,t_final):
        
        
        
        for i in range(num_particles):
            
            Fold= force(x)
            x[i] += dt*v[i] +0.5*(Fold[i])*dt**2
            if x[i]<0:
                x[i] =abs(x[i])
                
                vel = v[i]
                Fold= force(x)
                x[i] += dt*vel +0.5*(Fold[i])*dt**2
                
            elif x[i] > box_length:
                x[i] =box_length-abs(x[i]-box_length)
                
                vel = v[i]
                Fold= force(x)
                x[i] += dt*vel +0.5*(Fold[i])*dt**2
            Fnew = force(x)
            v[i] += (dt/2) *(Fold[i] + Fnew[i])
            
            
                
                
                
                
            ke[i] = 0.5*v[i]**2
            PE_new = PE(x)
            pe[i] = PE_new[i]
            TE_new = TE(pe,ke)
            te[i] = TE_new[i]
          
        KE_const = np.sum(ke) 
        KE_time[t] =KE_const 
        PE_const = np.sum(pe)
        PE_time[t] =PE_const
        energy_const = np.sum(te)  
        energy_time[t] =energy_const
        velocity_const = np.sum(v)  
        velocity_time[t] =velocity_const

        x_new=np.append(x_new,x)
    
        
    return x,v,energy_time,velocity_time,KE_time,PE_time,x_new

position,velocity,t_KE,t_PE,t_energy,t_velocity,x_new =vel_verlet(x, v, dt, Tmax)

print(time.time()- t1)
# %%
'''
Result =[t_KE,t_PE,t_energy,t_velocity]
fig,ax = plt.subplots(2,2,figsize =(12,12),tight_layout = True)
ax[0,0].set_ylabel('Kinetic energy',size = 10)
ax[0,1].set_ylabel('Potential energy',size = 10)
ax[1,0].set_ylabel('total energy',size = 10)
ax[1,1].set_ylabel('Momentum',size = 10)
ax[0,0].set_xlabel('time',size = 10)
ax[0,1].set_xlabel('time',size = 10)
ax[1,1].set_xlabel('time',size = 10)
ax[1,0].set_xlabel('time',size = 10)  
ax[0,0].plot(Result[0],color = 'blue')
ax[0,1].plot(Result[1] ,color = 'green')
ax[1,0].plot(Result[2],color = 'red')
ax[1,1].plot(Result[3],color = 'magenta')

plt.plot(t_KE,color = "black",label = 'Kinetic Energy')
plt.plot(t_PE,color = 'blue',label = "potential Energy")
plt.plot(t_energy,color = "green",label = 'Total Energy')
plt.plot(t_velocity,color = 'magenta',label = "momentum")
'''

plt.hist(x_new,bins = 60)
# print(x_new)
# print(position)        
               
                





































