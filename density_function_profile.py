#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 19:27:52 2022

@author: bhabanist
"""
import time
import functionrs as rs
import numpy as np
import matplotlib.pyplot as plt
t1 = time.time()
box_length = 500
density =[0.2,0.3,0.5,0.7]
pos_list =[]
for i in density:
    num_particles = int(box_length*i)
    



    Tmax = 5
    dt = 0.001
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

    F = rs.force(x,num_particles,sig)
    pe = rs.PE(x,num_particles,sig)
    ke = rs.KE(v,num_particles)
    te = rs.TE(pe,ke,num_particles)
    t_final = int(Tmax/dt)
    energy_time = np.zeros(t_final)
    velocity_time = np.zeros(t_final)
    position,velocity,t_energy,t_velocity =rs.vel_verlet(x,v,dt,Tmax,t_final,num_particles,box_length,ke,pe,te,energy_time,velocity_time,sig)
    pos_list.append(position)

'''
fig,ax = plt.subplots(1,2)      
ax[0].plot(t_energy,color = "red",label = 'Total Energy')
ax[0].plot(t_velocity,color = 'blue',label = "momentum")

ax[0].legend()

'''
print(time.time()- t1)
# %%
fig,ax = plt.subplots(2,2,figsize =(9,9))
ax[0,0].set_title("density = 0.2")
ax[0,1].set_title("density = 0.3")
ax[1,0].set_title("density = 0.5")
ax[1,1].set_title("density = 0.7")
ax[0,0].set_ylabel('density profile',size = 10)
ax[0,1].set_ylabel('density profile',size = 10)
ax[1,0].set_ylabel('density profile',size = 10)
ax[1,1].set_ylabel('density profile',size = 10)
ax[0,0].set_xlabel('position',size = 10)
ax[0,1].set_xlabel('position',size = 10)
ax[1,1].set_xlabel('position',size = 10)
ax[1,0].set_xlabel('position',size = 10)    
ax[0,0].hist(pos_list[0],bins = 100,density = False,color = 'blue')
ax[0,1].hist(pos_list[1],bins = 100,density = False,color = 'green')
ax[1,0].hist(pos_list[2],bins = 100,density = False,color = 'red')
ax[1,1].hist(pos_list[3],bins = 100,density = False,color = 'magenta')




