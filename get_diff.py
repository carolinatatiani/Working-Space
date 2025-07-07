# by Carolina Tatiani
# email: carolina.tatiani@unesp.br
# Last Modified: 07/07/2025 
# #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 14:10:06 2024

@author: rafael
"""

# by Carolina Tatiani
# email: carolina.tatiani@unesp.br
# Last Modified: 03/07/2025

import mdtraj as md
import numpy as np

def parse_top_file(top_file):
    """Parse the .top file to extract angle and dihedral indices."""
    angles_indices = []
    dihedrals_indices = []
    angles=[]
    dihedrals=[]
    
    with open(top_file, 'r') as file:
        inside_angles = False
        inside_dihedrals = False
        for line in file:
            line = line.strip()
            if line.startswith("[ angles ]"):
                inside_angles = True
                inside_dihedrals = False
                continue
            elif line.startswith("[ dihedrals ]"):
                inside_dihedrals = True
                inside_angles = False
                continue
            elif line.startswith('['):  # Other sections
                inside_angles = False
                inside_dihedrals = False
            
            # Parse angles
            if inside_angles and line and not line.startswith(';'):
                angvalue = line.split()
                angles_indices.append([int(angvalue[0])-1, int(angvalue[1])-1, int(angvalue[2])-1])
                           # Convert to 0-indexed
                angles.append([(int(angvalue[0])-1), (int(angvalue[1])-1), (int(angvalue[2])-1), angvalue[3], angvalue[4], angvalue[5]])
            # Parse dihedrals
            if inside_dihedrals and line and not line.startswith(';'):
                dihvalue = line.split()
                dihedrals.append([int(dihvalue[0])-1, int(dihvalue[1])-1, int(dihvalue[2])-1, int(dihvalue[3])-1,dihvalue[4],dihvalue[5],dihvalue[6]])  # Store the line for later use
                dihedrals_indices.append([int(dihvalue[0])-1, int(dihvalue[1])-1, int(dihvalue[2])-1, int(dihvalue[3])-1])  # Convert to 0-indexed
    np.save('angvalues.npy', angles)
    np.save('dihvalues.npy', dihedrals)
    return np.array(angles_indices), np.array(dihedrals_indices),np.array(angles), np.array(dihedrals)

def diff(angles1, angles2):
    diff = angles2 - angles1
    min_diff = (diff + 180) % 360 - 180
    return np.abs(min_diff)

# # File paths
# top = 'closed.AA.top'
# conf1 = md.load('closed.AA.pdb')
# conf2 = md.load('opened.AA.pdb')


# angles_indices, dihedrals_indices,angvalues,dihvalues = parse_top_file(top)


# ang1 = np.rad2deg(md.compute_angles(conf1, angles_indices)).T
# ang2 = np.rad2deg(md.compute_angles(conf2, angles_indices)).T
# angd =  diff(ang1, ang2)

# # #Verificar e printar se ang1 for 20% maior que ang2
# # indices_ang_percent = []
# # for i in range(len(ang1)):
# #    if np.any(ang1[i] > 1.5 * ang2[i]):
# #        indices_ang_percent.append(i)
# #        print(f"Angle {i}: ang1 = {ang1[i]}, ang2 = {ang2[i]} (ang1 is 20% greater)")
# # np.save('sel_ang_1based.npy', angles_indices[indices_ang_percent])


# # #### usando unique pq funct 1 ocupa duas linhas e funct 2 apenas uma
# dihedrals_indices= np.unique(dihedrals_indices, axis=0)
# dih1 = np.rad2deg(md.compute_dihedrals(conf1, dihedrals_indices)).T
# dih2 = np.rad2deg(md.compute_dihedrals(conf2, dihedrals_indices)).T
# dihd =  diff (dih1, dih2)

# # Verificar e printar se dih1 for 50% maior que dih2
# #indices_dih_percent = []
# #for i in range(len(dih1)):
# # #    if np.any(dih1[i] > 1.5 * dih2[i]):
# # #        indices_dih_percent.append(i)
# # #        print(f"Dihedral {i}: dih1 = {dih1[i]}, dih2 = {dih2[i]} (dih1 is 20% greater)")
# # #np.save('indices_dih_50_percent.npy', dihedrals_indices[indices_dih_percent])


# print(len(np.where(dihd>54.5)[0]))
# ind_dih=dihedrals_indices[np.where(dihd>54.5)[0]]
# np.save('sel_dih_1based.npy',ind_dih+1)
# # Adjust dihedrals difference to account for periodicity (-pi to pi)
# #dihedrals_diff = (dihedrals_diff + np.pi) % (2 * np.pi) - np.pi

