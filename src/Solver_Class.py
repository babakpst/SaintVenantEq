
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
        import sys

        import Initialization_Class
        import Visualization_Class


        Gravity = 9.81
        Draw = Visualization_Class.Visualization()

        #def RK2(self):
        Ex = Initialization_Class.Initialization()

        print(" ========== Solver Class ==========")
        print(" Solving the DE ...")
        print(" Second-order Runge-Kutta method: ...")

        print(" Allocating memory ...")
        N_Steps = int(Ex.Total_Time/Ex.Time_Step)
        N_Cells = Ex.N_Cells
        h_dw    = Ex.h_dw
        DT      = Ex.Time_Step

        Q    = np.zeros(N_Cells, dtype=np.float64 )
        V    = np.zeros(N_Cells, dtype=np.float64 )
        L    = np.zeros(N_Cells, dtype=np.float64 )
        Z    = np.zeros(N_Cells, dtype=np.float64 )
        Z_F  = np.zeros(N_Cells*2+1, dtype=np.float64 )
        M    = np.zeros(N_Cells, dtype=np.float64 )
        B    = np.zeros(N_Cells, dtype=np.float64 )

        A    = np.zeros(N_Cells, dtype=np.float64 )
        U    = np.zeros(N_Cells, dtype=np.float64 )
        Eta  = np.zeros(N_Cells, dtype=np.float64 )
        E    = np.zeros(N_Cells, dtype=np.float64 )
        l_P  = np.zeros(N_Cells, dtype=np.float64 )
        R_h  = np.zeros(N_Cells, dtype=np.float64 )
        C    = np.zeros(N_Cells, dtype=np.float64 )
        F    = np.zeros(N_Cells, dtype=np.float64 )
        X    = np.zeros(N_Cells, dtype=np.float64 )
        X_F  = np.zeros(N_Cells*2+1, dtype=np.float64 )

        V_0  = np.zeros(N_Cells, dtype=np.float64 )
        Q_0  = np.zeros(N_Cells, dtype=np.float64 )
        E_0  = np.zeros(N_Cells, dtype=np.float64 )
        F_0  = np.zeros(N_Cells, dtype=np.float64 )
        U_0  = np.zeros(N_Cells, dtype=np.float64 )
        E_F_0= np.zeros(N_Cells, dtype=np.float64 )

        V_1   = np.zeros(N_Cells, dtype=np.float64 )
        Q_1   = np.zeros(N_Cells, dtype=np.float64 )
        Eta_1 = np.zeros(N_Cells, dtype=np.float64 )
        E_1   = np.zeros(N_Cells, dtype=np.float64 )
        A_1   = np.zeros(N_Cells, dtype=np.float64 )
        l_P_1 = np.zeros(N_Cells, dtype=np.float64 )
        R_h_1 = np.zeros(N_Cells, dtype=np.float64 )
        F_1   = np.zeros(N_Cells, dtype=np.float64 )
        C_1   = np.zeros(N_Cells, dtype=np.float64 )
        U_1   = np.zeros(N_Cells, dtype=np.float64 )

        #
        A_F   = np.zeros(N_Cells+1, dtype=np.float64 )
        Eta_F = np.zeros(N_Cells+1, dtype=np.float64 )
        U_F   = np.zeros(N_Cells+1, dtype=np.float64 )
        E_F   = np.zeros(N_Cells+1, dtype=np.float64 )
        Q_F   = np.zeros(N_Cells+1, dtype=np.float64 )

        A_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        Eta_F_1 = np.zeros(N_Cells+1, dtype=np.float64 )
        U_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        E_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        Q_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )

        k_1V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array
        k_1Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array
        k_2V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array
        k_2Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array

        # Initialization 
        print(" Initialization ... ")
        Q[:] = Ex.Q[:]
        V[:] = Ex.V[:]
        L[:] = Ex.L[:]
        Z[:] = Ex.Z[:]
        Z_F[:] = Ex.Z_F[:]
        M[:] = Ex.M[:]
        B[:] = Ex.B[:]
        X[:] = Ex.X[:]
        X_F[:] = Ex.X_F[:]

        print(" Time marching ... ")
        for nn in range(N_Steps):
            print(" Time step: %d out of %d " % (nn, N_Steps))

            for ii in range(N_Cells):
                A[ii]   = V[ii] / L[ii]
                U[ii]   = Q[ii] / A[ii]
                Eta[ii] = A[ii] / B[ii] + Z[ii]
                E[ii]   = ((U[ii])**2) /2 +  Gravity*Eta[ii]
                l_P[ii] = B[ii] + 2 * Eta[ii]
                R_h[ii] = A[ii] / l_P[ii]
                C[ii]   = ((M[ii])**2) / ((R_h[ii])**(4.0/3.0))
                F[ii]   = Gravity * C[ii] * V[ii] * ((U[ii])**2)

            # <delete>
            RealTime = nn*DT
            TITLE = " at time: " + str(RealTime)
            Draw.Plot_at_Cell(N_Cells, X, Z, Q, V, Eta, U, E, A, TITLE)

            # Face reconstruction
            # Important comment: The size of the face arrays (..._F) are "N_Cells + 1". Face i+1/2 is indicated by index i. For example, face 1/2 is ..._F[0], face 1+1/2 is ..._F[1]
            #Rlts    = open("Results.txt","w")
            #Rlts.write("A %38.20f \n" % A[N_Cells-1])
            #Rlts.write("A %38.20f \n" % A[N_Cells-2])
            #Rlts.write("A %38.20f \n" % A[N_Cells-3])        
            #Rlts.write("A %38.20f \n" % A[0])
            #Rlts.close()

            #TITLE = " Area "
            #Draw.Plot(N_Cells, X, A, TITLE)

            #TITLE = " velocity "
            #Draw.Plot(N_Cells, X, U, TITLE)

            if nn == 0:
                for ii in range(N_Cells+1):
                    if ii==0: # Boundary condition at face 1/2
                        A_F[ii]   = A[ii] # This means A_(1/2) = A_1
                        Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        U_F[ii]   = Ex.Q_Up / A_F[ii]
                        E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]
                        Q_F[ii]   = A_F[ii] * U_F[ii]            
                    elif ii != 0 and ii != N_Cells: # middle cells - The subtraction is due to the fact that Python numbering is starting from 0
                        A_F[ii]   = (L[ii]*  A[ii-1] + L[ii-1]*  A[ii] )/( L[ii] + L[ii-1] )
                        Eta_F[ii] = (L[ii]*Eta[ii-1] + L[ii-1]*Eta[ii] )/( L[ii] + L[ii-1] )
                        E_F[ii]   = (L[ii]*  E[ii-1] + L[ii-1]*  E[ii] )/( L[ii] + L[ii-1] )
                        U_F[ii]   = (2*(E_F[ii] - Gravity * Eta_F[ii] ) )**(0.5)
                        Q_F[ii]   = A_F[ii] * U_F[ii]
                    elif ii == N_Cells: # Boundary condition at face N+1/2
                        Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
                        Z_N1      = Z[ii-1] - Delta_Z
                        Eta_F[ii] = h_dw
                        A_F[ii]   = (Eta_F[ii] - Z_N1 ) * B[ii-1]
                        #if Eta[ii-1]<Eta_F[ii]:
                        #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta[ii-1],Eta_F[ii]))    
                        #    sys.exit()
                        #Q_F[ii]   = (1.0/M[ii-1]) * A_F[ii] * ((R_h[ii-1])**(2.0/3.0))  * ( ( (Eta[ii-1]-Eta_F[ii])/(0.5*L[ii-1])  )**(1.0/2.0) ) # <modify>
                        Q_F[ii]   = Q[ii-1]
                        U_F[ii]   = Q_F[ii] / A_F[ii]
                        E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]
            elif nn != 0:   
                for ii in range(N_Cells+1):
                    if ii==0: # Boundary condition at face 1/2
                        A_F[ii]   = A[ii]
                        Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        U_F[ii]   = Ex.Q_Up / A_F[ii]
                        E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]
                        Q_F[ii]   = A_F[ii] * U_F[ii]            
                    elif ii != 0 and ii != N_Cells: # middle cells
                        A_F[ii]   = (L[ii]*  A[ii-1] + L[ii-1]*  A[ii] )/( L[ii] + L[ii-1] )
                        Eta_F[ii] = ( L[ii]*Eta[ii-1] + L[ii-1]*Eta[ii] ) / ( L[ii] + L[ii-1] )
                        E_F[ii]   = ( 1.0 / (V[ii-1]+V[ii]) ) * ( E_F_0[ii] *( V_0[ii-1] + V_0[ii] ) - E[ii-1]*V[ii-1] - E[ii]*V[ii] + E_0[ii-1]*V_0[ii-1] + E_0[ii]*V_0[ii] + DT * ( Q[ii-1] * E[ii-1] - Q[ii] * E[ii] - (1.0/2.0) * ( U[ii-1] * F[ii-1] + U[ii] * F[ii]  ) + Q_0[ii-1] * E_0[ii-1] - Q_0[ii] * E_0[ii] - (1.0/2.0) * ( U_0[ii-1] * F_0[ii-1] + U_0[ii] * F_0[ii]  )   ) )
                        if (E_F[ii] - Gravity * Eta_F[ii]) < 0:
                            print(" Fatal Error in Energy: %d, %d, %f, %f" % (ii, nn, E_F[ii], Gravity * Eta_F[ii] ))
                            sys.exit()
                        U_F[ii]   = (2*(E_F[ii] - Gravity * Eta_F[ii] ) )**(0.5)
                        Q_F[ii]   = A_F[ii] * U_F[ii]
                    elif ii == N_Cells: # Boundary condition at face N+1/2
                        Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
                        Z_N1      = Z[ii-1] - Delta_Z
                        Eta_F[ii] = h_dw
                        A_F[ii]   = (Eta_F[ii] - Z_N1 ) * B[ii-1]
                        #if Eta[ii-1]<Eta_F[ii]:
                        #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta[ii-1],Eta_F[ii]))    
                        #    sys.exit()
                        #Q_F[ii]   = (1.0/M[ii-1]) * A_F[ii] * ((R_h[ii-1])**(2.0/3.0))  * ( ( (Eta[ii-1]-Eta_F[ii])/(0.5*L[ii-1])  )**(1.0/2.0) )
                        Q_F[ii]   = Q[ii-1]
                        U_F[ii]   = Q_F[ii] / A_F[ii]
                        E_F[ii]   = ((U_F[ii])**2)/2 + Gravity * Eta_F[ii]

            # <delete>
            RealTime = nn*DT            
            TITLE = " at time: " + str(RealTime)
            Draw.Plot_Full(N_Cells, X_F, Z_F, Q, Q_F, Eta, Eta_F, U, U_F, E, E_F, A, A_F, TITLE)

            for ii in range(N_Cells): # To find k1 in the Runge-Kutta method and find the solution at n+1/2
                k_1V[ii]  = DT * ( Q_F[ii] - Q_F[ii+1] )  # <modify> remove
                k_1Q[ii]  = (DT / L[ii]) * ( Q_F[ii] * U_F[ii] - Q_F[ii+1] * U_F[ii+1] + Gravity * A_F[ii]* Eta_F[ii] - Gravity * A_F[ii+1]* Eta_F[ii+1] -F[ii] ) # <modify> remove
                # Solution at "n+1/2"
                V_1[ii]   = V[ii] + 0.5* k_1V[ii]
                Q_1[ii]   = Q[ii] + 0.5* k_1Q[ii]  # <modify> We really don't need to define k_1v and k_1q, just for the clarity of the code.

            for ii in range(N_Cells):  # These are the variables at {n+1}
                A_1[ii]   = V_1[ii] / L[ii]
                U_1[ii]   = Q_1[ii] / A_1[ii]
                Eta_1[ii] = A_1[ii] / B[ii] + Z[ii]
                E_1[ii]   = ((U_1[ii])**2) /2 +  Gravity*Eta_1[ii]
                l_P_1[ii] = B[ii] + 2 * Eta_1[ii]  # <modify>
                R_h_1[ii] = A_1[ii] / l_P_1[ii] # <modify>
                C_1[ii]   = ((M[ii])**2) / ((R_h_1[ii])**(4.0/3.0)) # <modify>
                F_1[ii]   = Gravity * C_1[ii] * V_1[ii] * ((U_1[ii])**2)

            for ii in range(N_Cells+1):
                if ii==0: # Boundary condition at face 1/2
                    A_F_1[ii]   = A_1[ii]
                    Eta_F_1[ii] = Eta_1[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                    U_F_1[ii]   = Ex.Q_Up / A_F_1[ii]
                    E_F_1[ii]   = ((U_F_1[ii])**2)/2 + Gravity * Eta_F_1[ii]
                    Q_F_1[ii]   = A_F_1[ii] * U_F_1[ii]                                
                elif ii != 1 and ii != N_Cells: # middle cells
                    A_F_1[ii]   = (L[ii]*  A_1[ii-1] + L[ii-1]*  A_1[ii] )/( L[ii] + L[ii-1] )
                    Eta_F_1[ii] = (L[ii]*Eta_1[ii-1] + L[ii-1]*Eta_1[ii] )/( L[ii] + L[ii-1] )
                    # Temp: E_{i+1/2}^{n+1/2}
                    E_F_1[ii]   = ( 1.0 / (V_1[ii-1]+V_1[ii]) ) * ( E_F_0[ii] *( V_0[ii-1] + V_0[ii] ) - E_1[ii-1]*V_1[ii-1] - E_1[ii]*V_1[ii] + E_0[ii-1]*V_0[ii-1] + E_0[ii]*V_0[ii] + DT * ( Q_1[ii-1] * E_1[ii-1] - Q_1[ii] * E_1[ii] - (1.0/2.0) * ( U_1[ii-1] * F_1[ii-1] + U_1[ii] * F_1[ii]  ) + Q_0[ii-1] * E_0[ii-1] - Q_0[ii] * E_0[ii] - (1.0/2.0) * ( U_0[ii-1] * F_0[ii-1] + U_0[ii] * F_0[ii]  )   ) )
                    if (E_F_1[ii] - Gravity * Eta_F_1[ii]) < 0:
                        print(" Fatal Error in energy at n+half: temp %d, %d, %f, %f" % (ii, nn, E_F_1[ii], Gravity * Eta_F_1[ii] ))
                        sys.exit()
                    U_F_1[ii]   = ( 2*(E_F_1[ii] - Gravity * Eta_F_1[ii] ) )**(0.5)
                    Q_F_1[ii]   = A_F_1[ii] * U_F_1[ii]          
                elif ii == N_Cells: # Boundary condition at face N+1/2
                    Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
                    Z_N1      = Z[ii-1] - Delta_Z
                    Eta_F_1[ii] = h_dw
                    A_F_1[ii]   = (Eta_F_1[ii] - Z_N1 ) * B[ii-1]
                    #if Eta_1[ii-1]<Eta_F_1[ii]:
                    #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta_1[ii-1],Eta_F[ii]))    
                    #    sys.exit()
                    #Q_F_1[ii]   = (1.0/M[ii-1]) * A_F_1[ii] * ((R_h_1[ii-1])**(2.0/3.0))  * ( ( (Eta_1[ii-1]-Eta_F_1[ii])/(0.5*L[ii-1])  )**(1.0/2.0) )  # <modify>
                    Q_F_1[ii]   = Q_1[ii-1]
                    U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                    E_F_1[ii]   = ((U_F_1[ii])**2)/2 + Gravity * Eta_F_1[ii]

            # <delete>
            RealTime = nn*DT            
            TITLE = " at time: " + str(RealTime)
            Draw.Plot_Full(N_Cells, X_F, Z_F, Q_1, Q_F_1, Eta_1, Eta_F_1, U_1, U_F_1, E_1, E_F_1, A_1, A_F_1, TITLE)

            for ii in range(N_Cells): # To find k2 in the Runge-Kutta method and find the solution at n + 1
                V_0[ii]  = V[ii]
                Q_0[ii]  = Q[ii]
                E_0[ii]  = E[ii]
                F_0[ii]  = F[ii]
                U_0[ii]  = U[ii]
                E_F_0[ii]= E_F[ii]

                k_2V[ii]  = DT * ( Q_F_1[ii] - Q_F_1[ii+1] )  # <modify> remove
                k_2Q[ii]  = (DT / L[ii]) * ( Q_F_1[ii] * U_F_1[ii] - Q_F_1[ii+1] * U_F_1[ii+1] + Gravity * A_F_1[ii]* Eta_F_1[ii] - Gravity * A_F_1[ii+1]* Eta_F_1[ii+1] -F_1[ii] ) # <modify> remove

                V[ii] = V[ii] + k_2V[ii]
                Q[ii] = Q[ii] + k_2Q[ii]


            if (nn%10000) == 0:
                RealTime = nn*DT
                TITLE = " at time: " + str(RealTime)
                Draw.Plot_at_Cell(N_Cells, X, Z, Q, V, Eta, U, E, A, TITLE)
        # End loop on the time steps

        print(" ========== Solver Class ==========")
        print()