
#####################################################################
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
    print()

  # -- Opens and Reads data from the input file
  def Input(self):
    print(" Opening the input file ...",self.InputFileName)
    File_Input = open(self.InputFileName,"r")

    # <Modify> Substitute these with LISTS with array.
    Reach_Length  = [] # Stores the length of each reach
    Reach_Disc    = [] # Stores the no. of control volume in each reach
    Reach_Slope   = [] # Stores the slope of each reach
    Reach_Manning = [] # Stores the Manning's number for each reach
    Reach_Width   = [] # Stores the width of each reach
  

    Temp = File_Input.readline().split("\n")  # 1
    Temp = File_Input.readline().split("\n")  # 2
    Temp = File_Input.readline().split("\n")  # 3
    Temp = File_Input.readline().split("\n")  # 4

    Temp = File_Input.readline().split("\n")  # 5
    Temp = File_Input.readline().split("\n")  # 6
    Q_Up = float(Temp[0])  # A constant flow rate at the upstream
    print(" Flow rate at the upstream is: %f" % Q_Up)

    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    h_dw = float(Temp[0])  # Downstream water depth
    print(" Downstream water depth is: %f" % h_dw)

    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    No_reaches  = int(Temp[0])  # Total number of reaches
    print(" Total number of reach(es) is(are): %d" % No_reaches)

    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    for i in range(No_reaches): # Length of each reach
      Temp       = File_Input.readline().split("\n")
      Reach_Length.append(float(Temp[0]))
      print(" Length of reach %d is: %f" % (i+1,Reach_Length[i-1]))

    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    for i in range(No_reaches): # Total number of control volumes in each reach
      Temp       = File_Input.readline().split("\n")
      Reach_Disc.append(int(Temp[0]))
      print(" No. of discretization of reach %d is: %f" % (i+1,Reach_Length[i-1]))

    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    for i in range(No_reaches): # Slope of each reach
      Temp       = File_Input.readline().split("\n")
      Reach_Slope.append(float(Temp[0]))
      print(" The slope of reach %d is: %f" % (i,Reach_Length[i-1]))

    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    for i in range(No_reaches): # The Manning's number for each reach
      Temp       = File_Input.readline().split("\n")
      Reach_Manning.append(float(Temp[0]))

    Temp       = File_Input.readline().split("\n")
    Temp       = File_Input.readline().split("\n")
    for i in range(No_reaches): # The width of each reach
      Temp       = File_Input.readline().split("\n")
      Reach_Width.append(float(Temp[0]))
    
    #return Input_Data


  # -- Class destructor
  def __del__(self):   
    Experiment.close()
    print("End Input_Class. ")


