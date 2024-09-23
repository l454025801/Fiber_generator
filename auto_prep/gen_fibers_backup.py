#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 16:46:19 2022

@author: leon
"""

import numpy as np
import sys

# import os
# import sys
# import inspect

# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0, parentdir) 

def write_inp(density=9, vacancy=0.3, height=5, filler="filler_h", ligand=False, ratio=50, layers=16, ligand_pos=None):
    """
    This function writes an .inp file for packmol to generate a fiber PDB file
    Parameters
    ----------
    density : int
        How many PA the fiber has per layer.
    vacancy : float
        The radius of the vacancy in the center of the fiber (unit A).
    height : float
        The distance between each layer (unit A).
    filler : string, optional
        The file name of the filler molecule. The default is "filler_h".
    ligand : string, optional
        The file name of the ligand molecule. The default is False.
    ratio : float, optional
        The ratio of filler:ligand. The default is 50.
    layers : int
        number of layers in the fiber. 
    ligand_pos : list
        index of ligands (where to place the ligands in the fiber).
        
    Returns
    -------
    None.

    """
    print("SETTING")
    print("The number of peptide amphiphile per layer is %i" %(density))
    print("The radius of the vacancy in the fiber center is %.2f A" %(vacancy))
    print("The distance between each layer is %.2f A" %(height))
    if ligand == False:
        print("No ligand in the fiber. Pure fillers only")
    else:
        print("The ligand in the fiber is %s" %(ligand))
        print("The ratio is %i : 1 [filler:ligand]" %(ratio))
    print ("-----------------------------------------------")
    print ("Writing packmol input file")
    print ("-----------------------------------------------")
    filename = "build_fiber.inp"
    inp = open(filename, "w") 
    rad = 2*3.1415926/density
    inp_script =\
"tolerance 0.5\n\
output fiber.pdb\n\
filetype pdb\n\
seed -1\n\
\n"
    if ligand == False:
        for j in range(int(layers)):
            if j % 2 == 0:
                for i in range(density):
                    inp_script += \
    "structure %s.pdb\n\
        number 1\n\
        fixed %.2f 0. %.1f 0. 0. 0.\n\
        fixed 0. 0. 0. 0. 0. %.3f \n\
    end structure\n\
    \n" %(filler, vacancy, j*height, i*rad)
            else:
                for i in range(density):
                    inp_script += \
    "structure %s.pdb\n\
        number 1\n\
        fixed %.2f 0. %.1f 0. 0. 0.\n\
        fixed 0. 0. 0. 0. 0. %.3f \n\
    end structure\n\
    \n" %(filler, vacancy, j*height, i*rad+rad/2)                    
        inp.write(inp_script)
        
    else:
        num_ligand = round(density*layers/ratio)
        #calculate how many ligands we need for the fiber
        if ligand_pos == None:
            ligand_pos = np.random.randint(0, density*layers, num_ligand)
        else:
            if len(ligand_pos) != num_ligand:
                sys.exit("Check the number of ligand needed")
        print("The indice for the ligand are")
        print(ligand_pos)
        print("----------------------------------------------------------")
        #give the index of ligands

        #create a string for ligand positions
        ligand_inp_script = ""
        for j in range((layers)):
            if j % 2 == 0:
                for i in range(density):
                    if density*j+i in ligand_pos:
                        ligand_inp_script += \
    "structure %s.pdb\n\
        number 1\n\
        fixed %.2f 0. %.1f 0. 0. 0.\n\
        fixed 0. 0. 0. 0. 0. %.3f \n\
    end structure\n\
    \n" %(ligand, vacancy, j*height, i*rad)
                    else:
                        inp_script += \
    "structure %s.pdb\n\
        number 1\n\
        fixed %.2f 0. %.1f 0. 0. 0.\n\
        fixed 0. 0. 0. 0. 0. %.3f \n\
    end structure\n\
    \n" %(filler, vacancy, j*height, i*rad)
            else:
                for i in range(density):
                    if density*j+i in ligand_pos:
                        ligand_inp_script += \
    "structure %s.pdb\n\
        number 1\n\
        fixed %.2f 0. %.1f 0. 0. 0.\n\
        fixed 0. 0. 0. 0. 0. %.3f \n\
    end structure\n\
    \n" %(ligand, vacancy, j*height, i*rad+rad/2)
                    else:
                        inp_script += \
    "structure %s.pdb\n\
        number 1\n\
        fixed %.2f 0. %.1f 0. 0. 0.\n\
        fixed 0. 0. 0. 0. 0. %.3f \n\
    end structure\n\
    \n" %(filler, vacancy, j*height, i*rad+rad/2)                
        inp_script += ligand_inp_script
        inp.write(inp_script)

    inp.close()
    print ("finish writing .inp file")
    return 

# write_inp(9,0.3,5,ligand="Z33_GG_h", layers=16)