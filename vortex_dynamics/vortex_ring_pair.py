#!/usr/bin/python

#Author: Michael Brothers
#Date: 10-22-2013
#Title: Vortex Ring Pair Interaction
#
#Note: must be run with four ranks

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import numpy.ma as ma
from mpl_toolkits.mplot3d import Axes3D
from mpi4py import MPI

comm = MPI.COMM_WORLD

rank = comm.Get_rank()

FLOOR = -30
CEILING = 30

class AnimatedScatter(object):
    def __init__(self, 
                 numpoints= 128,
                 Strength_One = 1.0,
                 Strength_Two = 1.0,
                 Radius_One = 10.0,
                 Radius_Two = 10.0,
                 X_Coord_One = -2.50,
                 X_Coord_Two = 2.50,
                 Timestep = 1.05,
                 Maxframes = 800,
		 Framerate = 20.0):

        self.numframes = 0
        self.maxframes = Maxframes
	self.framerate = Framerate
        self.timestep = Timestep
        self.numpoints = numpoints
        self.Radius_One = np.float(Radius_One)
        self.Radius_Two = np.float(Radius_Two)
        self.Strength_One = Strength_One
        self.Strength_Two = Strength_Two
        self.X_Coord_One = X_Coord_One
        self.X_Coord_Two = X_Coord_Two

        self.stream = self.data_stream()

        self.angle_sweep = np.linspace(0.0, (1.0 - 1.0/self.numpoints)*2.0*np.pi , self.numpoints)
        self.both_vortex_segment_centroids = np.zeros((2*self.numpoints,3))
        self.both_vortex_cart_velocities = np.zeros((2*self.numpoints, 3))
        self.set_Vortex_Segment_Centroids()
        self.recvbuff = np.zeros((8*self.numpoints, 3))
	self.scalarbuff = np.float(0.0)
	
	#rank 0 will write to disk, while rank 1 will show realtime output
	if(rank < 2):	
       		self.fig = plt.figure()
       		self.fig.canvas.mpl_connect('draw_event',self.forceUpdate)
       		self.ax = self.fig.add_subplot(111,projection = '3d')
       		self.ani = animation.FuncAnimation(self.fig, self.update, interval=1000.0/self.framerate, 
                                       init_func=self.setup_plot, frames= self.maxframes)
	#force the non-drawing ranks to iterate anyway
	if(rank > 1):
		while self.numframes < self.maxframes:
			data = next(self.stream)

    def set_Vortex_Segment_Centroids(self):

        for point in range(self.numpoints):
            self.both_vortex_segment_centroids[point, 0] = self.X_Coord_One
            self.both_vortex_segment_centroids[point, 1] = self.Radius_One*np.sin(self.angle_sweep[point])
            self.both_vortex_segment_centroids[point, 2] = self.Radius_One*np.cos(self.angle_sweep[point])
            self.both_vortex_segment_centroids[point + self.numpoints, 0] = self.X_Coord_Two
            self.both_vortex_segment_centroids[point + self.numpoints, 1] = self.Radius_Two*np.sin(self.angle_sweep[point])
            self.both_vortex_segment_centroids[point + self.numpoints, 2] = self.Radius_Two*np.cos(self.angle_sweep[point])


    def forceUpdate(self, event):
        self.scat.changed()

    def setup_plot(self):
        X = next(self.stream)
        color_one = ['r' for point in range(self.numpoints)]
        color_two = ['b' for point in range(self.numpoints)]
        colors = color_one + color_two

        self.scat = self.ax.scatter(X[:,0], X[:,1], X[:,2] , c=colors)

        self.ax.set_xlim3d(FLOOR, CEILING)
        self.ax.set_ylim3d(FLOOR, CEILING)
        self.ax.set_zlim3d(FLOOR, CEILING)

        return self.scat,

    def data_stream(self):
        data = self.both_vortex_segment_centroids
        
        while True and self.numframes < self.maxframes:
            self.numframes += 1

            if(rank == 0):
                self.compute_U_One_On_One_Vectorized()
            if(rank == 1):    
                self.compute_U_One_On_Two_Vectorized()
            if(rank == 2):
                self.compute_U_Two_On_One_Vectorized()
            if(rank == 3):
                self.compute_U_Two_On_Two_Vectorized()

            comm.Allgatherv([self.both_vortex_cart_velocities, MPI.DOUBLE], [self.recvbuff, MPI.DOUBLE])
            self.compute_Updated_Positions()

	    #print "rank " + str(rank) + " " + str(self.numframes)
            yield data
 
    def compute_U_One_On_One_Vectorized(self):
        m = self.angle_sweep
        s = self.angle_sweep
        dU = np.zeros((self.numpoints, self.numpoints -1, 3))
        dl = 0.0
        H_cubed = np.zeros(self.numpoints - 1)
        omega_direction = np.zeros((self.numpoints -1 , 3))

        x = np.zeros((self.numpoints -1, 3))
        x_prime = np.zeros((self.numpoints -1, 3))
        
        for self_point in range(self.numpoints):
            x = self.both_vortex_segment_centroids[self_point]
            x_prime[:self_point] = self.both_vortex_segment_centroids[:self_point]
            x_prime[self_point:] = self.both_vortex_segment_centroids[self_point +1:self.numpoints]

            omega_direction[:self_point,0] = 0.0 #x-dir magnitude
            omega_direction[:self_point,1] = np.cos(s[:self_point]) #y-dir magnitude
            omega_direction[:self_point,2] = -np.sin(s[:self_point]) #z-dir magnitude
            omega_direction[self_point:,0] = 0.0 #x-dir magnitude
            omega_direction[self_point:,1] = np.cos(s[self_point+1:]) #y-dir magnitude
            omega_direction[self_point:,2] = -np.sin(s[self_point+1:]) #z-dir magnitude

            dl = self.Radius_One* 2.0*np.pi * np.power(self.numpoints, -1.0)
           
            delta_x = x - x_prime
            H_cubed = np.power(np.sum(np.power(delta_x, 2.0), axis = 1), 3.0/2.0)

            dU[self_point] = self.Strength_One*np.power(4.0*np.pi, -1.0) * dl * np.power(H_cubed, -1.0)[:, np.newaxis] * np.cross(omega_direction, delta_x)

            #self.vortex_one_cart_velocities[self_point] = np.sum(dU[self_point], axis=0)
            self.both_vortex_cart_velocities[self_point] = np.sum(dU[self_point], axis=0)

        return 0

    def compute_U_One_On_Two_Vectorized(self):
        m = self.angle_sweep
        s = self.angle_sweep
        dU = np.zeros((self.numpoints, self.numpoints -1, 3))
        dl = 0.0
        H_cubed = np.zeros(self.numpoints - 1)
        omega_direction = np.zeros((self.numpoints -1 , 3))

        x = np.zeros((self.numpoints -1, 3))
        x_prime = np.zeros((self.numpoints -1, 3))
        
        for self_point in range(self.numpoints):
            x = self.both_vortex_segment_centroids[self_point + self.numpoints]
            x_prime[:self_point] = self.both_vortex_segment_centroids[:self_point]
            x_prime[self_point:] = self.both_vortex_segment_centroids[self_point +1:self.numpoints]

            omega_direction[:self_point,0] = 0.0 #x-dir magnitude
            omega_direction[:self_point,1] = np.cos(s[:self_point]) #y-dir magnitude
            omega_direction[:self_point,2] = -np.sin(s[:self_point]) #z-dir magnitude
            omega_direction[self_point:,0] = 0.0 #x-dir magnitude
            omega_direction[self_point:,1] = np.cos(s[self_point+1:]) #y-dir magnitude
            omega_direction[self_point:,2] = -np.sin(s[self_point+1:]) #z-dir magnitude

            dl = self.Radius_One* 2.0*np.pi * np.power(self.numpoints, -1.0)
           
            delta_x = x - x_prime
            H_cubed = np.power(np.sum(np.power(delta_x, 2.0), axis = 1), 3.0/2.0)

            dU[self_point] = self.Strength_One*np.power(4.0*np.pi, -1.0) * dl * np.power(H_cubed, -1.0)[:, np.newaxis] * np.cross(omega_direction, delta_x)

            #self.vortex_one_cart_velocities[self_point] = np.sum(dU[self_point], axis=0)
            self.both_vortex_cart_velocities[self_point + self.numpoints] = np.sum(dU[self_point], axis=0)

        return 0

    def compute_U_Two_On_One_Vectorized(self):
        m = self.angle_sweep
        s = self.angle_sweep
        dU = np.zeros((self.numpoints, self.numpoints -1, 3))
        dl = 0.0
        H_cubed = np.zeros(self.numpoints - 1)
        omega_direction = np.zeros((self.numpoints -1 , 3))

        x = np.zeros((self.numpoints -1, 3))
        x_prime = np.zeros((self.numpoints -1, 3))
        
        for self_point in range(self.numpoints):
            x = self.both_vortex_segment_centroids[self_point]
            x_prime[:self_point] = self.both_vortex_segment_centroids[self.numpoints :self_point +self.numpoints]
            x_prime[self_point:] = self.both_vortex_segment_centroids[self_point +1 +self.numpoints:]

            omega_direction[:self_point,0] = 0.0 #x-dir magnitude
            omega_direction[:self_point,1] = np.cos(s[:self_point]) #y-dir magnitude
            omega_direction[:self_point,2] = -np.sin(s[:self_point]) #z-dir magnitude
            omega_direction[self_point:,0] = 0.0 #x-dir magnitude
            omega_direction[self_point:,1] = np.cos(s[self_point+1:]) #y-dir magnitude
            omega_direction[self_point:,2] = -np.sin(s[self_point+1:]) #z-dir magnitude

            dl = self.Radius_Two* 2.0*np.pi * np.power(self.numpoints, -1.0)
           
            delta_x = x - x_prime
            H_cubed = np.power(np.sum(np.power(delta_x, 2.0), axis = 1), 3.0/2.0)

            dU[self_point] = self.Strength_Two*np.power(4.0*np.pi, -1.0) * dl * np.power(H_cubed, -1.0)[:, np.newaxis] * np.cross(omega_direction, delta_x)

            #self.vortex_one_cart_velocities[self_point] = np.sum(dU[self_point], axis=0)
            self.both_vortex_cart_velocities[self_point] = np.sum(dU[self_point], axis=0)

        return 0

    def compute_U_Two_On_Two_Vectorized(self):
        m = self.angle_sweep
        s = self.angle_sweep
        dU = np.zeros((self.numpoints, self.numpoints -1, 3))
        dl = 0.0
        H_cubed = np.zeros(self.numpoints - 1)
        omega_direction = np.zeros((self.numpoints -1 , 3))

        x = np.zeros((self.numpoints -1, 3))
        x_prime = np.zeros((self.numpoints -1, 3))
        
        for self_point in range(self.numpoints):
            x = self.both_vortex_segment_centroids[self_point + self.numpoints]
            x_prime[:self_point] = self.both_vortex_segment_centroids[self.numpoints :self_point + self.numpoints]
            x_prime[self_point:] = self.both_vortex_segment_centroids[self_point +1 + self.numpoints: ]

            omega_direction[:self_point,0] = 0.0 #x-dir magnitude
            omega_direction[:self_point,1] = np.cos(s[:self_point]) #y-dir magnitude
            omega_direction[:self_point,2] = -np.sin(s[:self_point]) #z-dir magnitude
            omega_direction[self_point:,0] = 0.0 #x-dir magnitude
            omega_direction[self_point:,1] = np.cos(s[self_point+1:]) #y-dir magnitude
            omega_direction[self_point:,2] = -np.sin(s[self_point+1:]) #z-dir magnitude

            dl = self.Radius_Two * 2.0*np.pi * np.power(self.numpoints, -1.0)
           
            delta_x = x - x_prime
            H_cubed = np.power(np.sum(np.power(delta_x, 2.0), axis = 1), 3.0/2.0)

            dU[self_point] = self.Strength_Two*np.power(4.0*np.pi, -1.0) * dl * np.power(H_cubed, -1.0)[:, np.newaxis] * np.cross(omega_direction, delta_x)

            #self.vortex_one_cart_velocities[self_point] = np.sum(dU[self_point], axis=0)
            self.both_vortex_cart_velocities[self_point + self.numpoints] = np.sum(dU[self_point], axis=0)

        return 0

    def compute_Updated_Positions(self):
        #self.vortex_one_segment_centroids += self.vortex_one_cart_velocities 
        #self.vortex_two_segment_centroids += self.vortex_two_cart_velocities 
        velocity = self.recvbuff[:2*self.numpoints] + self.recvbuff[2*self.numpoints: 4*self.numpoints] + self.recvbuff[4*self.numpoints: 6*self.numpoints] + self.recvbuff[6*self.numpoints: 8*self.numpoints] 
        self.both_vortex_segment_centroids += self.timestep*velocity

        #The new radius equals the average distance to the center of the ring

        self.Radius_One = np.linalg.norm(np.sum(np.power(np.power(self.both_vortex_segment_centroids[:self.numpoints] -
                                      np.sum(self.both_vortex_segment_centroids[:self.numpoints],
                                          axis = 0), 2.0), .5), axis = 0)[1:])/self.numpoints  
        self.Radius_Two = np.linalg.norm(np.sum(np.power(np.power(self.both_vortex_segment_centroids[self.numpoints:] -
                                      np.sum(self.both_vortex_segment_centroids[self.numpoints:],
                                          axis = 0), 2.0), .5), axis = 0)[1:])/self.numpoints 
                                      
        return 0

    def update(self, i):
        data = next(self.stream)
       	self.scat._offsets3d = ( np.ma.ravel(data[:,0]) , np.ma.ravel(data[:,1]) , np.ma.ravel(data[:,2]) )
       	return self.scat,

    def show(self):
        plt.show()

if __name__ == '__main__':
    a = AnimatedScatter() 
    if(rank==0):
    	a.ani.save("Rank" +str(rank) +".avi", writer='mencoder')
    if(rank==1):
    	a.show()

