
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

    def Plot_at_Cell(self, N, X, Z, Q, V, Eta, U, E, A, T):
        import numpy as np
        import matplotlib.pyplot as plt

        print(" This is the visualization class")
        
        Q_Arr   = np.zeros(N, dtype = np.float64)
        V_Arr   = np.zeros(N, dtype = np.float64)
        Eta_Arr = np.zeros(N, dtype = np.float64)
        U_Arr   = np.zeros(N, dtype = np.float64)
        E_Arr   = np.zeros(N, dtype = np.float64)
        X_Arr   = np.zeros(N, dtype = np.float64)
        Z_Arr   = np.zeros(N, dtype = np.float64)
        A_Arr   = np.zeros(N, dtype = np.float64)

        Title      = T
        Q_Arr[:]   = Q[:]
        V_Arr[:]   = V[:]
        Eta_Arr[:] = Eta[:]
        U_Arr[:]   = U[:]
        E_Arr[:]   = E[:]
        X_Arr[:]   = X[:]
        Z_Arr[:]   = Z[:]
        A_Arr[:]   = A[:]

        plt.figure(1)
        plt.subplot(321)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, Q_Arr, label ="Water flow" , color = "c", linewidth = 2.0)

        plt.title("Water flow"+Title, fontsize = 16)
        plt.xlabel("Distance (m)",          fontsize=12)
        plt.ylabel("Flow rate (m^3/s)",     fontsize=12)

        plt.subplot(322)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, V_Arr, label ="Control Volume" , color = "c", linewidth = 2.0)

        plt.title("Control Volume"+Title, fontsize = 16)
        plt.xlabel("Distance (m)",              fontsize=12)
        plt.ylabel("Contral Volume (m^3)",      fontsize=12)

        plt.subplot(323)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, Eta_Arr, label ="Water Elevation" ,  color = "c", linewidth = 2.0)
        plt.plot(X_Arr, Z_Arr,   label ="Bottom Elevation" , color = "r", linewidth = 2.0)

        plt.title("Water Elevation"+Title, fontsize = 16)
        plt.xlabel("Distance (m)", fontsize=12)
        plt.ylabel("Elevation (m)", fontsize=12)
        plt.legend(loc=0)

        plt.subplot(324)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, U_Arr, label ="Velocity" , color = "c", linewidth = 2.0)

        plt.title("Water Velocity"+Title, fontsize = 16)
        plt.xlabel("Distance (m)", fontsize=12)
        plt.ylabel("Velocity (m/s)", fontsize=12)

        plt.subplot(325)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, E_Arr, label ="Energy" , color = "c", linewidth = 2.0)

        plt.title("Energy"+Title,   fontsize = 16)
        plt.xlabel("Distance (m)",  fontsize=12)
        plt.ylabel("Energy (m/s)",  fontsize=12)        

        plt.subplot(326)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, A_Arr, label ="Area" , color = "c", linewidth = 2.0)

        plt.title("Area"+Title,   fontsize = 16)
        plt.xlabel("Distance (m)",  fontsize=12)
        plt.ylabel("Area (m^2)",  fontsize=12)

        #mng = plt.get_current_fig_manager()
        #mng.frame.Maximize(True) 
        
        #mng = plt.get_current_fig_manager()
        #mng.resize(*manager.window.maxsize())

        plt.show() # <modify> See why the execution stops when the the command gets here. 

    def Plot_Full(self, N, X_F, Z_F, Q, Q_F, Eta, Eta_F, U, U_F, E, E_F, A, A_F, T):
        import numpy as np
        import matplotlib.pyplot as plt
        
        Q_Arr   = np.zeros(N*2+1, dtype = np.float64)
        V_Arr   = np.zeros(N*2+1, dtype = np.float64)
        Eta_Arr = np.zeros(N*2+1, dtype = np.float64)
        U_Arr   = np.zeros(N*2+1, dtype = np.float64)
        E_Arr   = np.zeros(N*2+1, dtype = np.float64)
        X_Arr   = np.zeros(N*2+1, dtype = np.float64)
        Z_Arr   = np.zeros(N*2+1, dtype = np.float64)
        A_Arr   = np.zeros(N*2+1, dtype = np.float64)

        Title      = T
        X_Arr[:]   = X_F[:]
        Z_Arr[:]   = Z_F[:]

        for ii in range(N):
            Q_Arr[ii*2]    = Q_F[ii]
            Q_Arr[ii*2+1]  = Q[ii]

            Eta_Arr[ii*2]  = Eta_F[ii]
            Eta_Arr[ii*2+1]= Eta[ii]

            U_Arr[ii*2]    = U_F[ii]
            U_Arr[ii*2+1]  = U[ii]
            
            E_Arr[ii*2]    = E_F[ii]
            E_Arr[ii*2+1]  = E[ii]

            A_Arr[ii*2]    = A_F[ii]
            A_Arr[ii*2+1]  = A[ii]

        Q_Arr[N*2]    = Q_F[N]
        Eta_Arr[N*2]  = Eta_F[N]
        U_Arr[N*2]    = U_F[N]
        E_Arr[N*2]    = E_F[N]
        A_Arr[N*2]    = A_F[N]

        plt.figure(1)
        plt.subplot(321)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, Q_Arr, label ="Water flow (Q)" , color = "c", linewidth = 2.0)

        plt.title("Water flow (Q)"+Title, fontsize = 16)
        plt.xlabel("Distance (m)",          fontsize=12)
        plt.ylabel("Flow rate (m^3/s)",     fontsize=12)

        plt.subplot(323)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, Eta_Arr, label ="Water Elevation (Eta)" ,  color = "c", linewidth = 2.0)
        plt.plot(X_Arr, Z_Arr,   label ="Bottom Elevation" , color = "r", linewidth = 2.0)

        plt.title("Water Elevation (Eta)"+Title, fontsize = 16)
        plt.xlabel("Distance (m)", fontsize=12)
        plt.ylabel("Elevation (m)", fontsize=12)
        plt.legend(loc=0)

        plt.subplot(324)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, U_Arr, label ="Velocity (U)" , color = "c", linewidth = 2.0)

        plt.title("Water Velocity (U)"+Title, fontsize = 16)
        plt.xlabel("Distance (m)", fontsize=12)
        plt.ylabel("Velocity (m/s)", fontsize=12)

        plt.subplot(325)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, E_Arr, label ="Energy (E)" , color = "c", linewidth = 2.0)

        plt.title("Energy (E)"+Title,   fontsize = 16)
        plt.xlabel("Distance (m)",  fontsize=12)
        plt.ylabel("Energy (m/s)",  fontsize=12)        

        plt.subplot(326)
        plt.Figure(figsize=(50,15) )
        plt.plot(X_Arr, A_Arr, label ="Area (A)" , color = "c", linewidth = 2.0)

        plt.title("Area (A)"+Title,   fontsize = 16)
        plt.xlabel("Distance (m)",  fontsize=12)
        plt.ylabel("Area (m^2)",  fontsize=12)

        mng = plt.get_current_fig_manager()
        mng.frame.Maximize(True)
        plt.show() # <modify> See why the execution stops when the the command gets here. 

    def Plot(self, N, X, Z, T):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker

        print(" This is the visualization class")
        
        X_Arr = np.zeros(N, dtype = np.float64)
        Z_Arr = np.zeros(N, dtype = np.float64)

        Title = T
        X_Arr[:]   = X[:]
        Z_Arr[:]   = Z[:]


        fig, ax = plt.subplots()
        ax.plot(X_Arr, Z_Arr, label ="Water flow" , color = "c", linewidth = 2.0)
        
        y_labels = ax.get_yticks()
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%15.13f'))
        plt.show()


        #plt.figure(1)
        #plt.Figure(figsize=(50,15) )
        #plt.plot(X_Arr, Z_Arr, label ="Water flow" , color = "c", linewidth = 2.0)
        #fmt = ".0f%%"
        #xticks = mtick.FormatStrFormatter(fmt)
        #plt.title(Title, fontsize = 16)
        #plt.xlabel("Distance (m)",          fontsize=12)
        #plt.ylabel(Title,     fontsize=12)
        #plt.show() # <modify> See why the execution stops when the the command gets here.         