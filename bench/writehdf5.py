#!/usr/bin/python3
#########################################################################
#    SPF - Stochastic Phase Field
#    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#    See the README file in the top-level SPF directory.
#########################################################################
# File: writehdf5.py
# Purpose: This script writes an hdf5 file to use as input to SPF. Modify
#   accordingly to create the field values you desire.
# Requirements: numpy, h5py

import os.path
#import argparse
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.colors as colors
#from matplotlib.gridspec import GridSpec
#from matplotlib.patches import Circle
#from matplotlib.patches import Rectangle
#import statistics as stat
#from scipy.io import netcdf
import h5py

def main():

    debug = False
    #debug = True

    #parser = argparse.ArgumentParser()
    #parser.add_argument("inFilePath", help="data file in netcdf format")

    #args = parser.parse_args()
    #inFilePath = str( args.inFilePath )
    #inFileName = os.path.split( inFilePath )[1]
    #inFileRef = os.path.splitext( inFileName )[0]
    #inFileDir = os.path.split( inFilePath )[0]

    fileName = "initial_fields/initial_field_4x4x4_mean1.h5"

    #Nt = 2
    #Nx = 3
    #Ny = 3
    #Nz = 3
    Nx = 4
    Ny = 4
    Nz = 4
    #Nx = 30
    #Ny = 30
    #Nz = 30

    #data = np.zeros((Nt,Nx,Ny,Nz))
    data = np.zeros((Nx,Ny,Nz))
    #xx = np.zeros(Nx)
    #yy = np.zeros(Ny)
    #zz = np.zeros(Nz)

    #for t in range(Nt):
    #data[2,2,2] = 0.4
    mu, sigma = 1.0, 1.00
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                #if ( (i,j,k) == (0,0,0)):
                #    data[i,j,k] = 1.0
                #else:
                #    data[i,j,k] = 0.0
                data[i,j,k] = np.max((0.0, np.min((np.random.normal(mu, sigma), 2.0))))
                data[i,j,k] = np.round( data[i,j,k])
                #if ((i-50)**2 + (j-50)**2 + (k-50)**2 < (5)**2 ):
                #    data[i,j,k] = 100
                #else: data[i,j,k] = 0
                #data[i,j,k] = np.sin(i*2*np.pi / Nx )#*np.sin(j*2*np.pi/Ny - np.pi)*np.sin(k*2*np.pi/Nz - np.pi)
                #data[i,j,k] = np.sin(i*2*np.pi/Nx - np.pi)*np.sin(j*2*np.pi/Ny - np.pi)*np.sin(k*2*np.pi/Nz - np.pi)
                #if ( i == Nx/2 ) and ( j == Ny/2 ) and ( k == Nz/2 ) :
                #    data[i,j,k] = 1
                #if ( np.abs(i - Nx/2) < 2) and ( np.abs(j - Ny/2) < 2) and ( np.abs(k - Nz/2) < 2) :
                #    data[i,j,k] = 100
                #else :
                #    data[i,j,k] = 100
                #data[i,j,k] = 1.0/(1.0+np.abs(0.5*Nx - i) \
                #        + np.abs(0.5*Ny - j)   \
                #        + np.abs(0.5*Nz - k))
                #data[i,j,k] = 1.0/(1.0+((np.cos(np.abs(0.5*Nx - i)))**2 \
                #        + (np.cos(np.abs(0.5*Ny - j)))**2   \
                #        + (np.cos(np.abs(0.5*Nz - k)))**2 ))
                #data[i,j,k] = (i*1.0/Nx)
                #data[i,j,k] = (1.0/(Nx*Ny*Nz))*(k + Nz*(j + Ny*i))
                ##data[i,j,k] = k + Nz*(j + Ny*i)
                #print( data[i,j,k] )
                ##data[0,i,j,k] = k + Nz*(j + Ny*i)
                ##data[1,i,j,k] = 1 + k + Nz*(j + Ny*i)
                ##print( data[0,i,j,k] )

    outFile = h5py.File( fileName, 'w')#, mmap=False)
    #phidataset = outFile.create_dataset("phi", (Nt,Nx,Ny,Nz), dtype='f')
    phidataset = outFile.create_dataset("phi", (Nx,Ny,Nz), dtype='f')

    #for t in range(Nt):
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz): 
                #phidataset[t,i,j,k] = data[t,i,j,k]
                phidataset[i,j,k] = data[i,j,k]
       
    outFile.close()
    return


if __name__ == "__main__":
    main()
