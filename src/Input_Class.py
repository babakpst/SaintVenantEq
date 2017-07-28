

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


class Input_Info:

  # -- Class initialization 
  def __init__(self):    # Initialization - constructor
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

    self.InputFileName = Input_Dir[0] + "\\" + File[0] # InputFileName: the name of the input file
    self.Output_Dir = Output_Dir[0]

    print(" The input file name is: %s" % self.InputFileName)
    print(" The output directory is: %s" % self.Output_Dir)

    # Empty memory
    Address.close()
    del Temp
    del File
    del Input_Dir

    print(" Opening the input file ...")
    self.File_Input = open(self.InputFileName,"r")
    print()

  # -- Opens and Reads data from the input file
  #def Input(self):
    # Import built-in libraries =======================================================================
    import numpy as np

    Temp = self.File_Input.readline().split("\n")  # 1
    Temp = self.File_Input.readline().split("\n")  # 2
    Temp = self.File_Input.readline().split("\n")  # 3
    Temp = self.File_Input.readline().split("\n")  # 4

    Temp = self.File_Input.readline().split("\n")  # 5
    Temp = self.File_Input.readline().split("\n")  # 6
    Total_Time = float(Temp[0])  # Total simulation time
    print(" The total simulation time is: %f" % Total_Time)
    print()

    Temp = self.File_Input.readline().split("\n") 
    Temp = self.File_Input.readline().split("\n") 
    Temp = self.File_Input.readline().split("\n") 
    Time_Step = float(Temp[0])  # Time step
    print(" The time step is: %f" % Time_Step)
    print()

    Temp = self.File_Input.readline().split("\n")  
    Temp = self.File_Input.readline().split("\n")  
    Temp = self.File_Input.readline().split("\n")  
    Q_Up = float(Temp[0])  # A constant flow rate at the upstream
    print(" Flow rate at the upstream is: %f" % Q_Up)
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    h_dw = float(Temp[0])  # Downstream water depth
    print(" Downstream water depth is: %f" % h_dw)
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    No_reaches  = int(Temp[0])  # Total number of reaches
    print(" Total number of reach(es) is(are): %d" % No_reaches)
    print()

    # Define arrays: 
    Reach_Length  = np.zeros( No_reaches ) # Stores the length of each reach
    Reach_Disc    = np.zeros( No_reaches ) # Stores the no. of control volume in each reach
    Reach_Slope   = np.zeros( No_reaches ) # Stores the slope of each reach
    Reach_Manning = np.zeros( No_reaches ) # Stores the Manning's number for each reach
    Reach_Width   = np.zeros( No_reaches ) # Stores the width of each reach


    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # Length of each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Length[ii] = float(Temp[0])
      print(" The length of reach %d is: %f" % (ii+1, Reach_Length[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # Total number of control volumes in each reach/ For now we have a constant discretization in each reach.
      Temp       = self.File_Input.readline().split("\n")
      Reach_Disc[ii] = int(Temp[0])
      print(" No. of discretization of reach %d is: %f" % (ii+1, Reach_Disc[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for i in range(No_reaches): # Slope of each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Slope[ii] = float(Temp[0])
      print(" The slope of reach %d is: %f" % (ii+1, Reach_Slope[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # The Manning's number for each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Manning[ii] = float(Temp[0])
      print(" The Manning's no. for reach %d is: %f" % (ii+1, Reach_Manning[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # The width of each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Width[ii] = float(Temp[0])
      print(" The width of reach %d is: %f" % (ii+1, Reach_Width[ii]))
    print()

  # -- Class destructor
  def __del__(self):   
    print(" Closing the input file ... ")
    self.File_Input.close()
    print(" End Input_Class. ")


