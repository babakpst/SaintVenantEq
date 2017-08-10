
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 07/31/2017
#
# Comment: This class Discretizes the domain.
#
#####################################################################

class Discretization:
  
  def __init__(self):
    # -- Import libs/classes
    import Input_Class

    # Reading data from the input file 
    self.Experiment = Input_Class.Input_Info()

    print(" ========== Discretization Class ==========")
    print(" Discretization of the domain ...")

#  def discreization_func(self):
    # -- Import libs/classes
    import numpy as np
    print(" Loop over reaches to discretize the domain ...")

    # Find total number of cells in the domain:
    self.N_Cells = 0

    Temp_No_reaches = self.Experiment.No_reaches

    # Total number of cells
    print(" Calculating the total number of the cells in the domain ... ")
    for ii in range(Temp_No_reaches):
      self.N_Cells += self.Experiment.Reach_Disc[ii] # Have a check on this <Modify>
    print("   Total number of cells:   %d" % self.N_Cells)

    # Finding the highest point in the domain:
    print(" Calculating the highest point in the domain ... ")
    Max_Height = 0
    for ii in range(Temp_No_reaches):  
      Max_Height += self.Experiment.Reach_Slope[ii] * self.Experiment.Reach_Length[ii]
    print("   Maximum height is:       %d" % Max_Height)

    self.Length_Cell  = np.zeros( self.N_Cells, dtype=np.float ) # Stores the length of each cell
    self.Z_Cell       = np.zeros( self.N_Cells, dtype=np.float ) # Stores bottom elevation of each cell
    self.Manning_Cell = np.zeros( self.N_Cells, dtype=np.float ) # Stores the Manning's number of each cell
    self.Width_Cell   = np.zeros( self.N_Cells, dtype=np.float ) # Stores the Manning's number of each cell
    self.X_Disc       = np.zeros( self.N_Cells, dtype=np.float ) # Stores the Manning's number of each cell

    print(" Basic calculations ...")
    Cell_Counter = 0
    for ii in range(Temp_No_reaches):
      CntrlVolume_Length = round(self.Experiment.Reach_Length[ii]/self.Experiment.Reach_Disc[ii],5)  # Control volume length, rounded to 5 decimal points. The length of the final cell in reach would be adjusted to fix the discretization, if necessary. 
      print("   Cell length in the reach %d is: %f" % (ii+1, CntrlVolume_Length))
      Height             = Max_Height

      Z_loss             =  round(CntrlVolume_Length * self.Experiment.Reach_Slope[ii],8)
      Height            += round(0.5  * Z_loss,8)
      Total_Length       = 0

      X_distance = 0.5 * CntrlVolume_Length
      for jj in range(ii):
        X_distance  += self.Experiment.Reach_Length[ii]

      for jj in range( self.Experiment.Reach_Disc[ii] - 1 ):
        self.Length_Cell[Cell_Counter]  = CntrlVolume_Length
        self.X_Disc[Cell_Counter]       = X_distance
        X_distance                     += CntrlVolume_Length
        Total_Length                   += CntrlVolume_Length
        Height                         -= Z_loss

        self.Z_Cell[Cell_Counter]       = Height
        self.Manning_Cell[Cell_Counter] = self.Experiment.Reach_Manning[ii]
        self.Width_Cell[Cell_Counter]   = self.Experiment.Reach_Width[ii]
        Cell_Counter += 1

      # The last cell: we need to separate the last cell in each reach to adjust the numerical error in the total length of the reach
      self.Length_Cell[Cell_Counter]  = self.Experiment.Reach_Length[ii] - Total_Length
      X_distance                      = X_distance - 0.5 * CntrlVolume_Length + 0.5 * self.Length_Cell[Cell_Counter]
      self.X_Disc[Cell_Counter]       =  X_distance  
      Height                         -= ( 0.5*Z_loss + 0.5 * self.Length_Cell[Cell_Counter] * self.Experiment.Reach_Slope[ii] )

      self.Z_Cell[Cell_Counter]       = Height
      self.Manning_Cell[Cell_Counter] = self.Experiment.Reach_Manning[ii]
      self.Width_Cell[Cell_Counter]   = self.Experiment.Reach_Width[ii]
      Cell_Counter += 1
      Max_Height   -= self.Experiment.Reach_Length[ii] * self.Experiment.Reach_Slope[ii]

    self.Q_Up       = self.Experiment.Q_Up
    self.V_in       = self.Experiment.V_in
    self.Total_Time = self.Experiment.Total_Time
    self.Time_Step  = self.Experiment.Time_Step
    self.h_dw       = self.Experiment.h_dw
  

    # <delete> This part is for debugging purposes:
    #print(" Debugging section ...")
    #print(" X Coordinate: ")
    #print(self.X_Disc)
    #print(" Cell Length: ")
    #print(self.Length_Cell)
    #print(" Z Coordinate: ")
    #print(self.Z_Cell)
    #print(" Manning: ")
    #print(self.Manning_Cell)
    #print(" Width: ")
    #print(self.Width_Cell)

    if self.N_Cells != Cell_Counter:
      sys.exit("FATAL ERROR: Mismatch between the number of cells! Check the Discretization_Class.")
      
    print(" Discretization ends successfully. ")
    print(" ========== Discretization Class ==========")
    print()


# Check on the discretization: if ..._Cell = 0 => ERROR
