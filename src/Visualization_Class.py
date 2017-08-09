
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    08/03/2017
# Latest update: 08/08/2017
#
# Comment: This class visualizes the results using matplotlib
#
#####################################################################

class Visualization:
    
  def __init__(self):
    pass




  def Plot(self, N, X, Y, Z, T):
    import numpy as np
    import matplotlib.pyplot as plt

    print(" This is the visualization class")
    
    X_Arr = np.zeros(N, dtype = np.float64)
    Y_Arr = np.zeros(N, dtype = np.float64)
    Z_Arr = np.zeros(N, dtype = np.float64)

    Title     = T
    X_Arr[:] = X[:]
    Y_Arr[:] = Y[:]
    Z_Arr[:] = Z[:]

    plt.Figure(figsize=(35,15) )
    plt.plot(X_Arr, Y_Arr, label ="Water elevation" , color = "c")
    plt.plot(X_Arr, Z_Arr, label ="Bottom elevation" , color = "r")

    plt.title(Title, fontsize = 16)
    plt.xlabel("Distance (m)", fontsize=12)
    plt.ylabel("Elevation (m)", fontsize=12)

    plt.legend(loc=0)
    plt.show() # <modify> See why the execution stops when the the command gets here. 