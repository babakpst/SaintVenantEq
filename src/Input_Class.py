
###############################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 07/27/2017
#
# Comment: This class reads the simulation data from file
#
###############################################################################

class Input_Info:

    # -- Class initialization 
    def __init__(self):    # Initialization - constructor
        # Import built-in libraries ===========================================
        import numpy as np
        import os

        # Input data ==========================================================
        # -- Read the name of the input file from address file

        print(" ========== Input Class ==========")
        Address = open("Address.txt","r")
        Temp = Address.readline().rstrip("\n")  # 1
        File = Address.readline().rstrip("\n")  # 2, Input file name
        Temp = Address.readline().rstrip("\n")  # 3
        Temp = Address.readline().rstrip("\n")  # 4
        Input_Dir  = Address.readline().rstrip("\n")  # 5
        Temp = Address.readline().rstrip("\n")  # 6
        Temp = Address.readline().rstrip("\n")  # 7
        Output_Dir = Address.readline().rstrip("\n")  # 8

        InputFileName = os.path.join(Input_Dir, File) 
        self.Output_Dir = os.path.join(Output_Dir, os.path.splitext(File)[0])

        print(" {0} {1}".format("The input file name is:",InputFileName))
        print(" {0} {1}".format(" The output path is:", self.Output_Dir))

        # Empty memory
        Address.close()
        del Temp
        del File
        del Input_Dir

        print(" Opening the input file ...")
        self.File_Input = open(InputFileName,"r")
        print()

        # -- Opens and Reads data from the input file
        Temp = self.File_Input.readline().rstrip("\n")  # 1
        Temp = self.File_Input.readline().rstrip("\n")  # 2
        Temp = self.File_Input.readline().rstrip("\n")  # 3
        Temp = self.File_Input.readline().rstrip("\n")  # 4

        Temp = self.File_Input.readline().rstrip("\n")  # 5
        Temp = self.File_Input.readline().rstrip("\n")  # 6
        self.Total_Time = float(Temp)  # Total simulation time
        print("{:40} {:f}".format(" The total simulation time is:", self.Total_Time))
 
        Temp = self.File_Input.readline().rstrip("\n") 
        Temp = self.File_Input.readline().rstrip("\n") 
        Temp = self.File_Input.readline().rstrip("\n") 
        self.Time_Step = float(Temp)  # Time step
        print("{:40} {:f}".format(" The time step is:", self.Time_Step))

        Temp = self.File_Input.readline().rstrip("\n")  
        Temp = self.File_Input.readline().rstrip("\n")  
        Temp = self.File_Input.readline().rstrip("\n")  
        self.Q_Up = float(Temp)  # A constant flow rate at the upstream
        print("{:40} {:f}".format(" Flow rate at the upstream is:", self.Q_Up))

        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        self.h_dw = float(Temp)  # Downstream water depth
        print(" Downstream water depth is:           %f" % self.h_dw)
        

        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        self.V_in = float(Temp)  # Control volume at the upstream
        print(" Downstream water depth is:           %f" % self.V_in)

        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        self.V_ratio = float(Temp)  # Control volume at the upstream
        print(" Downstream water depth is:           %f" % self.V_in)


        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        self.No_reaches  = int(Temp)  # Total number of reaches
        print(" Total number of reach(es) is(are):   %d" % self.No_reaches)
        
        # Define arrays: 
        self.Reach_Length  = np.zeros( self.No_reaches, dtype=np.float ) # Stores the length of each reach
        self.Reach_Disc    = np.zeros( self.No_reaches, dtype=np.int   ) # Stores the no. of control volume in each reach
        self.Reach_Slope   = np.zeros( self.No_reaches, dtype=np.float ) # Stores the slope of each reach
        self.Reach_Manning = np.zeros( self.No_reaches, dtype=np.float ) # Stores the Manning's number for each reach
        self.Reach_Width   = np.zeros( self.No_reaches, dtype=np.float ) # Stores the width of each reach


        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        for ii in range(self.No_reaches): # Length of each reach
            Temp       = self.File_Input.readline().rstrip("\n")
            self.Reach_Length[ii] = float(Temp)
            print("      The length of reach %d is:         %f" % (ii+1, self.Reach_Length[ii]))
        

        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        for ii in range(self.No_reaches): # Total number of control volumes in each reach/ For now we have a constant discretization in each reach.
            Temp       = self.File_Input.readline().rstrip("\n")
            self.Reach_Disc[ii] = int(Temp)
            print("      No. of discretization of reach %d is:   %f" % (ii+1, self.Reach_Disc[ii]))

        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        for ii in range(self.No_reaches): # Slope of each reach
            Temp       = self.File_Input.readline().rstrip("\n")
            self.Reach_Slope[ii] = float(Temp)
            print("      The slope of reach %d is:               %f" % (ii+1, self.Reach_Slope[ii]))

        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        for ii in range(self.No_reaches): # The Manning's number for each reach
            Temp       = self.File_Input.readline().rstrip("\n")
            self.Reach_Manning[ii] = float(Temp)
            print("      The Manning's no. for reach %d is:      %f" % (ii+1, self.Reach_Manning[ii]))
        

        Temp       = self.File_Input.readline().rstrip("\n")
        Temp       = self.File_Input.readline().rstrip("\n")
        for ii in range(self.No_reaches): # The width of each reach
            Temp       = self.File_Input.readline().rstrip("\n")
            self.Reach_Width[ii] = float(Temp)
            print("      The width of reach %d is:               %f" % (ii+1, self.Reach_Width[ii]))
        
        print(" ========== Input Class Ends ==========")
        print()

    # -- Class destructor
    def __del__(self):   
        print(" Closing the input file ... ")
        self.File_Input.close()
        print(" End Input_Class. ")


