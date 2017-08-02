
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/31/2017
# Latest update: 08/01/2017
#
# Comment: This class conducts the numerical solution of the Saint-Venant equation
#
#####################################################################

class Solver:

  def __init__(self):
    # -- Import libs/classes
    import numpy as np
    import Initialization_Class

    Gravity = 9.81

    print(" Solving the DE ...")

  #def RK2(self):
    print(" Second-order Runge-Kutta method: ...")
    Ex = Initialization_Class.Initialization()
    N_Step  = int(Ex.Total_Time/Ex.Time_Step)
    N_Cells = Ex.N_Cells
    h_dw    = Ex.h_dw

    Q    = np.zeros( N_Cells, dtype=np.float64 )
    V    = np.zeros( N_Cells, dtype=np.float64 )
    L    = np.zeros( N_Cells, dtype=np.float64 )
    Z    = np.zeros( N_Cells, dtype=np.float64 )
    M    = np.zeros( N_Cells, dtype=np.float64 )
    B    = np.zeros( N_Cells, dtype=np.float64 )

    A    = np.zeros( N_Cells, dtype=np.float64 )
    U    = np.zeros( N_Cells, dtype=np.float64 )
    Eta  = np.zeros( N_Cells, dtype=np.float64 )
    E    = np.zeros( N_Cells, dtype=np.float64 )
    l_P  = np.zeros( N_Cells, dtype=np.float64 )
    R_h  = np.zeros( N_Cells, dtype=np.float64 )
    C    = np.zeros( N_Cells, dtype=np.float64 )
    F    = np.zeros( N_Cells, dtype=np.float64 )

    #
    A_F   = np.zeros( N_Cells+1, dtype=np.float64 )
    Eta_F = np.zeros( N_Cells+1, dtype=np.float64 )
    U_F   = np.zeros( N_Cells+1, dtype=np.float64 )
    E_F   = np.zeros( N_Cells+1, dtype=np.float64 )
    Q_F   = np.zeros( N_Cells+1, dtype=np.float64 )

    # Initialization 
    Q[:] = Ex.Q[:]
    V[:] = Ex.V[:]
    L[:] = Ex.L[:]
    Z[:] = Ex.Z[:]
    M[:] = Ex.M[:]
    B[:] = Ex.B[:]

    for nn in range(N_Step):
      print(" Time step: %d out of %d " % (nn, N_Step))

      for ii in range(N_Cells):
        A[ii]   = V[ii] / L[ii]
        U[ii]   = Q[ii] / A[ii]
        Eta[ii] = A[ii] / B[ii] + Z[ii]
        E[ii]   = ((U[ii])**2) /2 +  Gravity*Eta[ii]
        l_P[ii] = B[ii] + 2 * Eta[ii]
        R_h[ii] = A[ii] / l_P[ii]
        C[ii]   = ((Eta[ii])**2) / ((R_h[ii])**(4/3))
        F[ii]   = Gravity * C[ii] * V[ii] * ((U[ii])**2)

      # Face reconstruction
      if nn == 0:
        for ii in range(N_Cells+1):
          if ii==0: # Boundary condition at face 1/2
            A_F[ii]   = A[ii]
            Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
            U_F[ii]   = Ex.Q_Up / A_F[ii]
            E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]
          elif ii != 1 and ii != N_Cells: # middle cells
            A_F[ii]   = (L[ii+1]*A[ii]   + L[ii]*A[ii+1]   )/( L[ii+1] + L[ii] )
            Eta_F[ii] = (L[ii+1]*Eta[ii] + L[ii]*Eta[ii+1] )/( L[ii+1] + L[ii] )
            E_F[ii]   = (L[ii+1]*E[ii]   + L[ii]*E[ii+1]   )/( L[ii+1] + L[ii] )
            U_F[ii]   = (2*(E_F[ii] - Gravity * Eta_F[ii] ) )**(0.5)
            Q_F[ii]   = A_F[ii] * U_F[ii]
          elif ii == N_Cells: # Boundary condition at face N+1/2
            Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
            Z_N1      = Z[ii-1] - Delta_Z
            Eta_F[ii] = h_dw
            A_F[ii]   = (Eta_F[ii] - Z_N1 ) * B[ii]
            Q_F[ii]   = (1.0/M[ii-1]) * A[ii-1] * ((R_h[ii-1])**(2.0/3.0))  * ( ( (Eta[ii-1]-Eta_F[ii])/(0.5*L[ii-1])  )**(1.0/2.0) )
            U_F[ii]   = Q_F[ii] / A_F[ii]
            E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]
      elif nn != 0:     #   hereeeeeeeeeeeeeeeeeeeeeeeeee
        for ii in range(N_Cells+1):
          if ii==0: # Boundary condition at face 1/2
            A_F[ii]   = A[ii]
            Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
            U_F[ii]   = Ex.Q_Up / A_F[ii]
            E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]
          elif ii != 1 and ii != N_Cells: # middle cells
            A_F[ii]   = (L[ii+1]*A[ii]   + L[ii]*A[ii+1]   )/( L[ii+1] + L[ii] )
            Eta_F[ii] = (L[ii+1]*Eta[ii] + L[ii]*Eta[ii+1] )/( L[ii+1] + L[ii] )
            E_F[ii]   = (L[ii+1]*E[ii]   + L[ii]*E[ii+1]   )/( L[ii+1] + L[ii] )
            U_F[ii]   = (2*(E_F[ii] - Gravity * Eta_F[ii] ) )**(0.5)
            Q_F[ii]   = A_F[ii] * U_F[ii]
          elif ii == N_Cells: # Boundary condition at face N+1/2
            Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
            Z_N1      = Z[ii-1] - Delta_Z
            Eta_F[ii] = h_dw
            A_F[ii]   = (Eta_F[ii] - Z_N1 ) * B[ii]
            Q_F[ii]   = (1.0/M[ii-1]) * A[ii-1] * ((R_h[ii-1])**(2.0/3.0))  * ( ( (Eta[ii-1]-Eta_F[ii])/(0.5*L[ii-1])  )**(1.0/2.0) )
            U_F[ii]   = Q_F[ii] / A_F[ii]
            E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]






      for ii in range(N_Cells+1): # To find k1 in the Runge-Kutta method.
        pass

      for ii in range(N_Cells+1): # To find cell variables at the n+1/2 time
        pass

      # 
      for ii in range(N_Cells+1):
        if ii==0: # Boundary condition at face 1/2
          pass
        elif ii != 1 and ii != N_Cells: # middle cells
          pass
        elif ii == N_Cells: # Boundary condition at face N+1/2
          pass

      for ii in range(N_Cells+1): # To find k2 in the Runge-Kutta method.
        pass

      # Visualization
    # End loop on the time steps

