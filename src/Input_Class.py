

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
  def __init__(self, InputFileName0):    # Initialization - constructor
    print(" Reading data from the input file ...")
    self.InputFileName = InputFileName0
    print(" Opening the input file ...")
    self.File_Input = open(self.InputFileName,"r")
    print()

  # -- Opens and Reads data from the input file
  def Input(self):

    # <Modify> Substitute these with LISTS with array.
    Reach_Length  = [] # Stores the length of each reach
    Reach_Disc    = [] # Stores the no. of control volume in each reach
    Reach_Slope   = [] # Stores the slope of each reach
    Reach_Manning = [] # Stores the Manning's number for each reach
    Reach_Width   = [] # Stores the width of each reach
  

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
    Time_Step = float(Temp[0])  # Time step
    print(" The time step is: %f" % Time_Step)
    print()

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

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # Length of each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Length.append(float(Temp[0]))
      print(" The length of reach %d is: %f" % (ii+1,Reach_Length[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # Total number of control volumes in each reach/ For now we have a constant discretization in each reach.
      Temp       = self.File_Input.readline().split("\n")
      Reach_Disc.append(int(Temp[0]))
      print(" No. of discretization of reach %d is: %f" % (ii+1,Reach_Disc[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for i in range(No_reaches): # Slope of each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Slope.append(float(Temp[0]))
      print(" The slope of reach %d is: %f" % (ii+1,Reach_Slope[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # The Manning's number for each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Manning.append(float(Temp[0]))
      print(" The Manning's no. for reach %d is: %f" % (ii+1,Reach_Manning[ii]))
    print()

    Temp       = self.File_Input.readline().split("\n")
    Temp       = self.File_Input.readline().split("\n")
    for ii in range(No_reaches): # The width of each reach
      Temp       = self.File_Input.readline().split("\n")
      Reach_Width.append(float(Temp[0]))
      print(" The width of reach %d is: %f" % (ii+1,Reach_Width[ii]))
    print()


    '''
    Reach_Length  = [] # Stores the length of each reach
    Reach_Disc    = [] # Stores the no. of control volume in each reach
    Reach_Slope   = [] # Stores the slope of each reach
    Reach_Manning = [] # Stores the Manning's number for each reach
    Reach_Width   = [] # Stores the width of each reach

    Q_Up,
    h_dw,    
    No_reaches, 
    '''
#    return Input_Data []


  # -- Class destructor
  def __del__(self):   
    print(" Closing the input file ... ")
    self.File_Input.close()
    print(" End Input_Class. ")


