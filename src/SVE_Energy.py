
#####################################################################
# 
# This code is based on the methodology described in the paper.
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 07/25/2017
#
#####################################################################

# Import built-in libraries =======================================================================

# import numpy as np

# Import classes ==================================================================================
import Input_Class
import Discretization_Class
#import Initialization

# Code begins =====================================================================================
print(" -------- 1D Saint-Venant Equation based on the face reconstruction --------")
print(" ---------- Developers:: Babak Poursartip/Ben R. Hodges/Frank Liu ----------")
print()
print()
print(" Simulation starts ...")
print()

# Input data ======================================================================================
# -- Read the name of the input file from address file
Address    = open("Address.txt","r")
Temp       = Address.readline().split("\n")  # 1
File       = Address.readline().split("\n")  # 2, Input file name
Temp       = Address.readline().split("\n")  # 3
Temp       = Address.readline().split("\n")  # 4
Input_Dir  = Address.readline().split("\n")  # 5
Temp       = Address.readline().split("\n")  # 6
Temp       = Address.readline().split("\n")  # 7
Output_Dir = Address.readline().split("\n")  # 8


InputFileName = Input_Dir[0] + "\\" + File[0] # InputFileName: the name of the input file
Output_Dir = Output_Dir[0]

print(" The input file name is: %s" % InputFileName)
print(" The output directory is: %s" % Output_Dir)

# Empty memory
Address.close()
del Temp
del File
del Input_Dir

# Reading data from the input file 
Experiment = Input_Class.Input_Info(InputFileName)
Input_Data = Experiment.Input() # Put the input data in Input_Data

# Initialization ==================================================================================












#def main():
#    print('Hello, World!')
#
#
#if __name__ == '__main__':
#    main()





