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

    fileName = "initial_fields/gl_initial_field_4x4x4_mean50.h5"

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
    #data = np.zeros((Nx,Ny,Nz))
    dataConc = 0.5*np.ones((Nx,Ny,Nz))
    dataPhi = 0.5*np.ones((Nx,Ny,Nz))
    dataT = 9.8*np.ones((Nx,Ny,Nz))
    #dataConc = np.zeros((Nx,Ny,Nz))
    #dataPhi = np.zeros((Nx,Ny,Nz))
    #dataT = np.zeros((Nx,Ny,Nz))
    #xx = np.zeros(Nx)
    #yy = np.zeros(Ny)
    #zz = np.zeros(Nz)

    #for t in range(Nt):
    #data[2,2,2] = 0.4
    mu, sigma = 50, 1.25
    maxConc = 100
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                #if ( (i,j,k) == (0,0,0)):
                #    data[i,j,k] = 1.0
                #else:
                #    data[i,j,k] = 0.0

                dataConc[i,j,k] = np.random.normal(mu, sigma)
                dataConc[i,j,k] = np.round( dataConc[i,j,k])/maxConc 

                #if i == round(Nx/2): 
                #    dataT[i,j,k] = 11.0
                #else:
                #    dataT[i,j,k] = 9.3
                #    
                #if j == round(Nx/2):
                #    dataPhi[i,j,k] = 1.0
                #else:
                #    dataPhi[i,j,k] = 0.0

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
    concdataset = outFile.create_dataset("conc", (Nx,Ny,Nz), dtype='f')
    phidataset = outFile.create_dataset("phi", (Nx,Ny,Nz), dtype='f')
    Tdataset = outFile.create_dataset("T", (Nx,Ny,Nz), dtype='f')

    #for t in range(Nt):
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz): 
                #phidataset[t,i,j,k] = data[t,i,j,k]
                concdataset[i,j,k] = dataConc[i,j,k]
                phidataset[i,j,k] = dataPhi[i,j,k]
                Tdataset[i,j,k] = dataT[i,j,k]
       
    outFile.close()
    return


if __name__ == "__main__":
    main()
