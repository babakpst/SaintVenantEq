

#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    08/03/2017
# Latest update: 08/03/2017
#
# Comment: This class visualizes the results using matplotlib
#
#####################################################################

class Visualization:
    
  def __init__(self):
    pass

  def Plot(self, N, X, Y, Z, T):
    import matplotlib.pyplot as plt

    X_Arr = np.zeros(N, dtype = np.float64)
    Y_Arr = np.zeros(N, dtype = np.float64)
    Z_Arr = np.zeros(N, dtype = np.float64)

    Title     = T
    X_Plot[:] = X[:]
    Y_Plot[:] = Y[:]
    Z_Plot[:] = Z[:]

    plt.Figure(figsize=(35,15) )
    plt.plot(X_Plot, Y_Plot, label ="Water elevation" , color = "c")
    plt.plot(X_Plot, Z_Plot, label ="Bottom elevation" , color = "r")

    plt.title(Title, fontsize = 16)
    plt.xlabel("X label", fontsize=12)
    plt.ylabel("Y label", fontsize=12)

    plt.legend(loc=0)
    plt.show()    