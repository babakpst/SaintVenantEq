
###############################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 11/20/2017
#
# Comment: This class Discretizes the domain.
#
###############################################################################

class Discretization:

    def __init__(self,argv):
        import Input_Class
        self.Experiment = Input_Class.Input_Info(argv)
        self.command_line_args=self.Experiment.command_line_args

    def Discretize(self):
        import numpy as np
        import Visualization_Class as Visual

        Draw = Visual.Visualization()

        #Func = lambda x: 0.2 - 0.05 * ((x-7)-10)**2
        #DFunc = lambda x: -( - 0.05 * 2 * ((x-7)-10))

        # For now, we define the geometry function directly here, this should be moved to the input file.
        Func = lambda x: 0.2 - 0.05 * (x-10)**2
        DFunc = lambda x: -( - 0.05 * 2 * (x-10))

        # Reading data from the input file 
    
        self.Experiment.Read_Data()
        self.Experiment.Read_Input()
        self.Output_Dir = self.Experiment.Output_Dir

        print(" ========== Discretization Class ==========")
        print(" Discretization of the domain ...")

        print(" Loop over reaches to discretize the domain ...")

        # Find total number of cells in the domain:
        self.N_Cells = 0

        Temp_No_reaches = self.Experiment.No_reaches

        # Total number of cells
        print(" Calculating the total number of the cells in the domain ... ")
        for ii in range(Temp_No_reaches):
            self.N_Cells += self.Experiment.Reach_Disc[ii] # Have a check on this <Modify>
        print("{} {:d}".format("   Total number of cells:   ", self.N_Cells))

        self.Length_Cell  = np.zeros(self.N_Cells, dtype=np.float) # Stores the length of each cell
        self.HLength_Cell = np.zeros(self.N_Cells, dtype=np.float) # Stores the length of each cell
        self.Z_Cell       = np.zeros(self.N_Cells, dtype=np.float) # Stores bottom elevation of each cell
        self.S_Cell       = np.zeros(self.N_Cells, dtype=np.float) # Stores bottom elevation of each cell
        self.Z_Full       = np.zeros(self.N_Cells*2+1, dtype=np.float) # Stores bottom elevation of each cell
        self.Manning_Cell = np.zeros(self.N_Cells, dtype=np.float) # Stores the Manning's number of each cell
        self.Width_Cell   = np.zeros(self.N_Cells, dtype=np.float) # Stores the Manning's number of each cell
        self.X_Disc       = np.zeros(self.N_Cells, dtype=np.float) # Stores the Manning's number of each cell
        self.X_Full       = np.zeros(self.N_Cells*2+1, dtype=np.float) # Stores the Manning's number of each cell

        # Finding the highest point in the domain:
        print(" Calculating the highest point in the domain ... ")
        Max_Height = 0
        for ii in range(Temp_No_reaches):  
            Max_Height += self.Experiment.Reach_Slope[ii] * self.Experiment.Reach_Length[ii]
        print("{} {:f}".format(" Maximum height is:", Max_Height))

        print(" Basic calculations ...")
        Cell_Counter = 0
        for ii in range(Temp_No_reaches):
            if self.Experiment.Reach_Type[ii]==0:
                # Control volume length, rounded to 5 decimal points. The length of the final cell in reach would be 
                #    adjusted to fix the discretization, if necessary. 
                CntrlVolume_Length = round(self.Experiment.Reach_Length[ii]/self.Experiment.Reach_Disc[ii],10)  
                print("   Cell length in the reach %d is: %f" % (ii+1, CntrlVolume_Length))
                Height = Max_Height
                Z_loss = round(CntrlVolume_Length * self.Experiment.Reach_Slope[ii],10)
                Height += round(0.5 * Z_loss,10)
                Total_Length = 0

                X_distance = 0.5 * CntrlVolume_Length
                for jj in range(ii):
                    X_distance  += self.Experiment.Reach_Length[jj]

                for jj in range( self.Experiment.Reach_Disc[ii] - 1 ):
                    self.HLength_Cell[Cell_Counter]  = CntrlVolume_Length
                    self.Length_Cell[Cell_Counter]   = CntrlVolume_Length * (1.0+self.Experiment.Reach_Slope[ii])
                    self.X_Disc[Cell_Counter]       = X_distance
                    self.X_Full[Cell_Counter*2]     = X_distance - 0.5 * CntrlVolume_Length
                    self.X_Full[Cell_Counter*2+1]   = X_distance
                    X_distance                     += CntrlVolume_Length
                    Total_Length                   += CntrlVolume_Length
                    Height                         -= Z_loss

                    self.Z_Cell[Cell_Counter]       = Height
                    self.Z_Full[Cell_Counter*2]     = Height +0.5 * Z_loss
                    self.Z_Full[Cell_Counter*2+1]   = Height
                    self.Manning_Cell[Cell_Counter] = self.Experiment.Reach_Manning[ii]
                    self.Width_Cell[Cell_Counter]   = self.Experiment.Reach_Width[ii]
                    Cell_Counter += 1

                # The last cell: we need to separate the last cell in each reach to adjust the numerical error in the total length of the reach
                self.HLength_Cell[Cell_Counter] = self.Experiment.Reach_Length[ii] - Total_Length
                self.Length_Cell[Cell_Counter]  = (self.Experiment.Reach_Length[ii] - Total_Length) * (1.0+self.Experiment.Reach_Slope[ii])
                X_distance                      = X_distance - 0.5 * CntrlVolume_Length + 0.5 * self.HLength_Cell[Cell_Counter]
                self.X_Disc[Cell_Counter]       = X_distance  
                self.X_Full[Cell_Counter*2]     = X_distance - 0.5 * self.HLength_Cell[Cell_Counter]
                self.X_Full[Cell_Counter*2+1]   = X_distance  
                self.X_Full[Cell_Counter*2+2]   = X_distance + 0.5 * self.HLength_Cell[Cell_Counter] 
                Height                         -= ( 0.5*Z_loss + 0.5 * self.HLength_Cell[Cell_Counter] * self.Experiment.Reach_Slope[ii] )

                self.S_Cell[Cell_Counter]       = self.Experiment.Reach_Slope[ii]
                self.Z_Cell[Cell_Counter]       = Height 
                self.Z_Full[Cell_Counter*2]     = Height + 0.5 * self.HLength_Cell[Cell_Counter] * self.Experiment.Reach_Slope[ii]
                self.Z_Full[Cell_Counter*2+1]   = Height 
                self.Z_Full[Cell_Counter*2+2]   = Height - 0.5 * self.HLength_Cell[Cell_Counter] * self.Experiment.Reach_Slope[ii]
                self.Manning_Cell[Cell_Counter] = self.Experiment.Reach_Manning[ii]
                self.Width_Cell[Cell_Counter]   = self.Experiment.Reach_Width[ii]
                Cell_Counter += 1
                Max_Height   -= self.Experiment.Reach_Length[ii] * self.Experiment.Reach_Slope[ii]

            elif self.Experiment.Reach_Type[ii]==1:
                Height = Max_Height
                Projection_Length = round(self.Experiment.Reach_Length[ii]/self.Experiment.Reach_Disc[ii],10) 
                X_distance = 0.5 * Projection_Length

                for jj in range(ii):
                    print("DEBUG:: X-distance: ",X_distance) # <delete>
                    X_distance  += self.Experiment.Reach_Length[jj]

                for jj in range( self.Experiment.Reach_Disc[ii] ):

                    Z_loss = Func ( X_distance+ 0.5 * Projection_Length) - Func ( X_distance-0.5 * Projection_Length ) 
                    #print("DEBUG:: X-distance: ",X_distance,Z_loss) # <delete>

                    self.HLength_Cell[Cell_Counter] = Projection_Length
                    self.Length_Cell[Cell_Counter] = (Projection_Length**2 + Z_loss**2)**0.5
                    
                    self.X_Disc[Cell_Counter] = X_distance
                    self.X_Full[Cell_Counter*2] = X_distance - 0.5 * Projection_Length
                    self.X_Full[Cell_Counter*2+1] = X_distance

                    self.S_Cell[Cell_Counter]       = DFunc(X_distance)
                    self.Z_Cell[Cell_Counter]       = Height + Func(X_distance)
                    self.Z_Full[Cell_Counter*2]     = Height + Func(X_distance - 0.5 * Projection_Length)
                    self.Z_Full[Cell_Counter*2+1]   = Height + Func(X_distance)
                    self.Manning_Cell[Cell_Counter] = self.Experiment.Reach_Manning[ii]
                    self.Width_Cell[Cell_Counter]   = self.Experiment.Reach_Width[ii]
                    X_distance += Projection_Length
                    Total_Length += Projection_Length
                    Cell_Counter += 1

                Max_Height   -= self.Experiment.Reach_Length[ii] * self.Experiment.Reach_Slope[ii]

        self.Q_Up       = self.Experiment.Q_Up
        self.V_in       = self.Experiment.V_in
        self.V_ratio    = self.Experiment.V_ratio
        self.Total_Time = self.Experiment.Total_Time
        self.Time_Step  = self.Experiment.Time_Step
        self.h_dw       = self.Experiment.h_dw
      
        if self.N_Cells != Cell_Counter: # Check whether the discretization is correct.
            sys.exit("FATAL ERROR: Mismatch between the number of cells! Check the Discretization_Class.")

        del self.Experiment

        # Plot the discretized domain -- delete after debugging
        # Title = "Discretized domain at the cell level"
        # Vis =Draw.Plot_Domain(self.N_Cells, self.X_Disc, self.Z_Cell, Title)
        # Vis =Draw.Plot_Domain(2*self.N_Cells+1, self.X_Full, self.Z_Full, Title)

        print(" Discretization ends successfully. ")
        print(" ======== Discretization Class Ends ========")
        print()
