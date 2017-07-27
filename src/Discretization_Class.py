
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 07/27/2017
#
# Comment: This class reads the simulation data from file
#
#####################################################################

class Discretization:
  
  def __init__(self, Input_Info):
    
    print(" Discretization of the domain ...")


  def discreization_func(self):

    print(" Loop over reaches to discretize the domain ...")
    N_Cells = 0
    for ii in range(No_reaches):

      CntrlVolume_Length = (Reach_Length[ii]/Reach_Disc[ii])  



      N_Cells += Reach_Disc[ii] # Have a check on this