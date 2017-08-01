

#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 07/31/2017
#
# Comment: This class processes the input data and provides the initial information for the simulation. 
#
#####################################################################

class Initialization:



  def __init__(self):
    # -- Import libs/classes
    import numpy as np
    import Discretization_Class
    self.Disc = Discretization_Class.Discretization()

    print(" Initialization ...")
    self.Q  = np.zeros( self.Disc.N_Cells, dtype=np.float64 )
    self.V  = np.zeros( self.Disc.N_Cells, dtype=np.float64 )
    
    for ii in range( self.Disc.N_Cells ):
      self.Q[ii] = self.Disc.Q_Up
      self.V[ii] = self.Disc.V_in

  def Geometry_func(self, V, L): # Returns the cross-section area: A
    return V/L

  def Elevation_func(self, A, B , Z):  # Returns the free surface elevation: eta
    return A/B + Z