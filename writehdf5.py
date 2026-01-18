#!/usr/bin/python3
#########################################################################
#    SPF - Stochastic Phase Field
#    Copyright (C) 2025 
#    Peng Geng <penggeng@g.ucla.edu>
#    Nicholas Huebner Julian <njulian@ucla.edu>
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

import numpy as np
import h5py

def main():
    
    "Initial fields file name."
    fileName = "initial_fields/initial_field_60x60x60_C0.5T300.h5"
    
    "Input any size of voxel system along the x, y, and z directions."
    Nx = 60
    Ny = 60
    Nz = 60
    
    "Please choose one voxel field by commenting out the corresponding commands below."
    
    "1. Homogenious system (Example: Cr-50 at.% W at 300 K)"
    
    #Concentration
    data = 0.5*np.ones((Nx,Ny,Nz))
    #Temperature
    dataT = 300*np.ones((Nx,Ny,Nz))
    
    "2. Gradient system along y direction"
    "(Example: 6 layers concentration gradient range from Cr-30 at.% W to Cr-80 at.% W)"
    "(         Each concentration layer has 10 at.% W increment)"
    "(         60 layers temperature gradient range from 300-360 K)"
    "(         Each temperature layer has 1 K increment)"
    
    #Concentration gradient
    #data = 0.3*np.ones((Nx,Ny,Nz))

    #num_parts = 6
    #values_per_part = 10
    #increment_per_part = 0.1

    #for part in range(num_parts):
    #    start_idx = part * values_per_part
    #    end_idx = (part+1) * values_per_part
    #    data[:, start_idx:end_idx, :] += part * increment_per_part
    
    
    #Temperature gradient
    #dataT = 300*np.ones((Nx,Ny,Nz))
    
    #num_parts = 60
    #values_per_part = 1
    #increment_per_part = 1

    #for part in range(num_parts):
    #    start_idx = part * values_per_part
    #    end_idx = (part+1) * values_per_part
    #    dataT[:, start_idx:end_idx, :] += part * increment_per_part
    
    "set up output"
    outFile = h5py.File( fileName, 'w')
    phidataset = outFile.create_dataset("phi", (Nx,Ny,Nz), dtype='f')
    Tdataset = outFile.create_dataset("T", (Nx,Ny,Nz), dtype='f')

    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz): 
                phidataset[i,j,k] = data[i,j,k]
                Tdataset[i,j,k] = dataT[i,j,k]
       
    outFile.close()
    return

if __name__ == "__main__":
    main()
