

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

  import numpy as np

  def __init__(self):
    # -- Import libs/classes
    import Discretization_Class
    Disc = Discretization_Class.Discretization()
    Var = Disc.discreization_func()


    print(" Initialization ...")
    self.Q  = np.zeros( Var.N_Cells, dtype=np.float )
    self.V  = np.zeros( Var.N_Cells, dtype=np.float )
    
    for ii in range(Var.N_Cells):
      self.Q[ii] = Var.Q_Up
      self.V[ii] = Var.V_in

  def Geometry_func(self, V, L): # Returns the cross-section area: A
    return V/L

  def Elevation_func(self, A, B , Z):  # Returns the free surface elevation: eta
    return A/B + Z