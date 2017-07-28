
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 07/28/2017
#
# Comment: This class Discretizes the domain.
#
#####################################################################

class Discretization:
  
  def __init__(self):
    # -- Import libs/classes
    import Input_Class

    # Reading data from the input file 
    Experiment = Input_Class.Input_Info()
    print(" Discretization of the domain ...")


  def discreization_func(self):
    # -- Import libs/classes
    import numpy as np
    print(" Loop over reaches to discretize the domain ...")

    # Find total number of cells in the domain:
    N_Cells = 0
    for ii in range(No_reaches):
      N_Cells += Reach_Disc[ii] # Have a check on this <Modify>

    # Finding the highest point in the domain:
    Max_Height = 0
    for ii in range(No_reaches):  
      Max_Height += Experiment.Reach_Slope[ii] * Experiment.Reach_Length[ii]

    Length_Cell  = np.zeros( N_Cells ) # Stores the length of each cell
    Z_Cell       = np.zeros( N_Cells ) # Stores bottom elevation of each cell
    Manning_Cell = np.zeros( N_Cells ) # Stores the Manning's number of each cell

    Cell_Counter = 0
    Height       = Max_Height
    for ii in range(No_reaches):
      CntrlVolume_Length = round(Reach_Length[ii]/Reach_Disc[ii],5)  # Control volume length, rounded to 5 decimal points. The length of the final cell in reach would be adjusted to fix the discretization, if necessary. 
      Z_loss             =  CntrlVolume_Length * Experiment.Reach_Slope[ii]
      Height             = Height - 0.5  * Z_loss
      for jj in range( Experiment.Reach_Disc[ii] - 1 ):
        Length_Cell[Cell_Counter]  = CntrlVolume_Length




        Z_Cell[Cell_Counter]       = Height
        Height                     -=  Z_loss
        Manning_Cell[Cell_Counter] = Experiment.Reach_Manning[ii]
        Cell_Counter += 1
      # The last cell
        
    if N_Cells != Cell_Counter:
      sys.exit("FATAL ERROR: Mismatch between the number of cells! Check the Discretization_Class.")
      


# Check on the discretization: if ..._Cell = 0 => ERROR

      