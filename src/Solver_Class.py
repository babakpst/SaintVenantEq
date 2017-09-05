
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/31/2017
# Latest update: 08/15/2017
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
        IntType = Ex.InterpolationType

        Q    = np.zeros(N_Cells, dtype=np.float64 )
        V    = np.zeros(N_Cells, dtype=np.float64 )
        L    = np.zeros(N_Cells, dtype=np.float64 )
        Z    = np.zeros(N_Cells, dtype=np.float64 )
        Z_F  = np.zeros(N_Cells*2+1, dtype=np.float64 )
        M    = np.zeros(N_Cells, dtype=np.float64 )
        B    = np.zeros(N_Cells, dtype=np.float64 )
        Fr   = np.zeros(N_Cells, dtype=np.float64 )

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

        slowness = 50000
        h_upstream = V[0]/(B[0]*L[0])

        print(" Time marching ... ")
        for nn in range(N_Steps):
            print(" Time step: %d out of %d " % (nn, N_Steps))

            if nn < slowness:
              Q_Upstream = Ex.Q_Up * (nn/ float(slowness))
            else:
              Q_Upstream = Ex.Q_Up        
        
            for ii in range(N_Cells):
                A[ii]   = V[ii] / L[ii]
                U[ii]   = Q[ii] / A[ii]
                Eta[ii] = A[ii] / B[ii] + Z[ii]
                E[ii]   = ((U[ii])**2.0) /(float(2)) +  Gravity*Eta[ii]
                l_P[ii] = B[ii] + 2.0 * (Eta[ii]-Z[ii])
                R_h[ii] = A[ii] / l_P[ii]
                C[ii]   = ((M[ii])**2.0) / ((R_h[ii])**(4.0/3.0))
                F[ii]   = Gravity * C[ii] * V[ii] * ((U[ii])**2.0)
                Fr[ii]  = U[ii]/((Gravity * (Eta[ii] - Z[ii]) )**(0.5))
                if Fr[ii] >= 1.0:
                  print("Flow is not subcritical")
                  check = input(" Error: press ENTER to exit ")
                  sys.exit()


            # <delete>
            if (nn%10000) == 0:
                print(Q_Upstream)
                RealTime = round(nn*DT,5)
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
                        if IntType == 3 or IntType == 4:
                            Q_F[ii]   = Q_Upstream
                            #A_F[ii]   = A[ii] # This means A_(1/2) = A_1
                            #U_F[ii]   = Q_F[ii] / A_F[ii]
                            #Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))

                            Z_0 = Z[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                            #Eta_F[ii] = h_upstream+ Z_0
                            Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                            A_F[ii]   =  (Eta_F[ii] - Z_0)* B[ii]

                            U_F[ii]   = Q_F[ii] / A_F[ii]
                            E_F[ii]   = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]
                    elif ii != 0 and ii != N_Cells: # middle cells - The subtraction is due to the fact that Python numbering is starting from 0
                        if IntType == 3:
                            A_F[ii]   = (L[ii]*  A[ii-1] + L[ii-1]*  A[ii] )/( L[ii] + L[ii-1] )
                            Eta_F[ii] = (L[ii]*Eta[ii-1] + L[ii-1]*Eta[ii] )/( L[ii] + L[ii-1] )
                            #E_F[ii]   = (L[ii]*  E[ii-1] + L[ii-1]*  E[ii] )/( L[ii] + L[ii-1] ) 
                            E_F[ii]   = Eta_F[ii] * Gravity # <delete>
                            if ( E_F[ii] < Gravity * Eta_F[ii] ):
                                print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], Gravity * Eta_F[ii] ))
                                print(" %60.50f" % (E_F_1[ii] - Gravity * Eta_F_1[ii]))                            
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
                            if ( E_F[ii] < E_F[ii-1] ):
                                print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], E_F[ii-1] ) )
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
                            #U_F[ii]   = ( 2*(E_F[ii] - Gravity * Eta_F[ii] ) )**(0.500)
                            U_F[ii]   = 0.0 # <delete>
                            Q_F[ii]   = A_F[ii] * U_F[ii]
                        elif IntType == 4:
                            A_F[ii]   = (L[ii]*  A[ii-1] + L[ii-1]*  A[ii] )/( L[ii] + L[ii-1] )
                            U_F[ii]   = (L[ii]*  U[ii-1] + L[ii-1]*  U[ii] )/( L[ii] + L[ii-1] )    
                            Q_F[ii]   = A_F[ii] * U_F[ii]
                            E_F[ii]   = (L[ii]*  E[ii-1] + L[ii-1]*  E[ii] )/( L[ii] + L[ii-1] ) 
                            if ( E_F[ii] < E_F[ii-1] ):
                                print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], E_F[ii-1] ) )
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
                            Eta_F[ii] = (1.0/(2 * Gravity)) * ( 2 * E_F[ii] - (U_F[ii])**(2.0) )
                    elif ii == N_Cells: # Boundary condition at face N+1/2
                        if IntType == 3 or IntType == 4:
                            Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
                            Z_N1      = Z[ii-1] - Delta_Z
                            Eta_F[ii] = h_dw
                            A_F[ii]   = (Eta_F[ii] - Z_N1 ) * B[ii-1]
                            #if Eta[ii-1]<Eta_F[ii]:
                            #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta[ii-1],Eta_F[ii]))    
                            #    check = input(" Error: press ENTER to exit ")
                            #    sys.exit()
                            #Q_F[ii]   = (1.0/M[ii-1]) * A_F[ii] * ((R_h[ii-1])**(2.0/3.0))  * ( ( (Eta[ii-1]-Eta_F[ii])/(0.5*L[ii-1])  )**(1.0/2.0) ) # <modify>
                            Q_F[ii]   = Q[ii-1]
                            U_F[ii]   = Q_F[ii] / A_F[ii]
                            E_F[ii]   = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]
                            if ( E_F[ii] < E_F[ii-1] ):
                                print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], E_F[ii-1] ) )
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
            elif nn != 0:   
                for ii in range(N_Cells+1):
                    if ii==0: # Boundary condition at face 1/2
                        if IntType == 3 or IntType == 4:
                            Q_F[ii]   = Q_Upstream
                            #A_F[ii]   = A[ii]
                            #U_F[ii]   = Q_F[ii] / A_F[ii]
                            #Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))

                            Z_0 = Z[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                            #Eta_F[ii] = h_upstream+ Z_0
                            Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                            A_F[ii]   =  (Eta_F[ii] - Z_0)* B[ii]

                            U_F[ii]   = Q_F[ii] / A_F[ii]
                            E_F[ii]   = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]

                    elif ii != 0 and ii != N_Cells: # middle cells
                        if IntType == 3:
                            A_F[ii]   = ( L[ii]*  A[ii-1] + L[ii-1]*  A[ii] )/( L[ii] + L[ii-1] )
                            Eta_F[ii] = ( L[ii]*Eta[ii-1] + L[ii-1]*Eta[ii] ) / ( L[ii] + L[ii-1] )

                            #if Eta_F[ii] != A_F[ii] / B[ii]  + Z_F[2*ii]:
                            #    print(" Fatal Error: inconsistency btw A and eta: %d %d %f %f"% (nn, ii, Eta_F[ii], A_F[ii] / B[ii]  + Z_F[2*ii] ))
                            #    check = input(" Error: press ENTER to exit ")
                            #    sys.exit()
                            E_F[ii]   = ( 1.0 / (V[ii-1]+V[ii]) ) * ( E_F_0[ii] *( V_0[ii-1] + V_0[ii] ) - E[ii-1]*V[ii-1] - E[ii]*V[ii] + E_0[ii-1]*V_0[ii-1] + E_0[ii]*V_0[ii] + DT * ( Q[ii-1] * E[ii-1] - Q[ii] * E[ii] - (1.0/2.0) * ( U[ii-1] * F[ii-1] + U[ii] * F[ii]  ) + Q_0[ii-1] * E_0[ii-1] - Q_0[ii] * E_0[ii] - (1.0/2.0) * ( U_0[ii-1] * F_0[ii-1] + U_0[ii] * F_0[ii]  )   ) )
                            if E_F[ii]  < 0:
                                print(" Fatal Error: Negative Energy: %d, %d, %f " % (nn, ii, E_F[ii] ))
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
                            if (E_F[ii] - Gravity * Eta_F[ii]) < 0.0000001:
                                E_F[ii] = Gravity * Eta_F[ii]
                            if (E_F[ii] - Gravity * Eta_F[ii]) < 0.0:
                                print(" Fatal Error in Energy: %d, %d, %20.10f, %20.10f" % (nn, ii, E_F[ii], Gravity * Eta_F[ii] ))
                                print(" %30.20f" % (E_F_1[ii] - Gravity * Eta_F_1[ii]))                            
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
                            if ( E_F[ii] < E_F[ii-1] ):
                                print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], E_F[ii-1] ) )
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
                            U_F[ii]   = (2.0*(E_F[ii] - Gravity * Eta_F[ii] ) )**(0.5)
                            Q_F[ii]   = A_F[ii] * U_F[ii]
                        elif IntType == 4:
                            A_F[ii]   = ( L[ii]*  A[ii-1] + L[ii-1]*  A[ii] )/( L[ii] + L[ii-1] ) # interpolation
                            U_F[ii]   = ( L[ii]*  U[ii-1] + L[ii-1]*  U[ii] )/( L[ii] + L[ii-1] ) # interpolation
                            Q_F[ii]   = A_F[ii] * U_F[ii]                                         # solve
                            E_F[ii]   = ( 1.0 / (V[ii-1]+V[ii]) ) * ( E_F_0[ii] *( V_0[ii-1] + V_0[ii] ) - E[ii-1]*V[ii-1] - E[ii]*V[ii] + E_0[ii-1]*V_0[ii-1] + E_0[ii]*V_0[ii] + DT * ( Q[ii-1] * E[ii-1] - Q[ii] * E[ii] - (1.0/2.0) * ( U[ii-1] * F[ii-1] + U[ii] * F[ii]  ) + Q_0[ii-1] * E_0[ii-1] - Q_0[ii] * E_0[ii] - (1.0/2.0) * ( U_0[ii-1] * F_0[ii-1] + U_0[ii] * F_0[ii]  )   ) )
                            Eta_F[ii] = (1.0/(2 * Gravity)) * ( 2 * E_F[ii] - (U_F[ii])**(2.0) )  # solve
                            if abs(Eta_F[ii] - (A_F[ii] / B[ii]  + Z_F[2*ii])) > 0.001:
                                print(" Fatal Error: inconsistency btw A and eta: %d %d %30.20f %30.20f"% (nn, ii, Eta_F[ii], A_F[ii] / B[ii]  + Z_F[2*ii] ))
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
                    elif ii == N_Cells: # Boundary condition at face N+1/2
                        if IntType == 3 or IntType == 4:
                            Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
                            Z_N1      = Z[ii-1] - Delta_Z
                            Eta_F[ii] = h_dw
                            A_F[ii]   = (Eta_F[ii] - Z_N1 ) * B[ii-1]
                            #if Eta[ii-1]<Eta_F[ii]:
                            #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta[ii-1],Eta_F[ii]))    
                            #    check = input(" Error: press ENTER to exit ")
                            #    sys.exit()
                            #Q_F[ii]   = (1.0/M[ii-1]) * A_F[ii] * ((R_h[ii-1])**(2.0/3.0))  * ( ( (Eta[ii-1]-Eta_F[ii])/(0.5*L[ii-1])  )**(1.0/2.0) )
                            Q_F[ii]   = Q[ii-1]
                            U_F[ii]   = Q_F[ii] / A_F[ii]
                            E_F[ii]   = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]
                            if ( E_F[ii] < E_F[ii-1] ):
                                print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], E_F[ii-1] ) )
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()
            # <delete>
            if (nn%100000) == 0:
                RealTime = round(nn*DT,5)
                TITLE = "k-1 at time: " + str(RealTime)
                Draw.Plot_Full(2,N_Cells, X_F, Z_F, Q, Q_F, Eta, Eta_F, U, U_F, E, E_F, A, A_F, TITLE)

            for ii in range(N_Cells): # To find k1 in the Runge-Kutta method and find the solution at n+1/2
                k_1V[ii]  =  DT * ( Q_F[ii] - Q_F[ii+1] )  # <modify> remove
                k_1Q[ii]  = (DT / L[ii]) * ( Q_F[ii] * U_F[ii] - Q_F[ii+1] * U_F[ii+1] + Gravity * A_F[ii]* Eta_F[ii] - Gravity * A_F[ii+1]* Eta_F[ii+1] -F[ii] ) # <modify> remove
                # Solution at "n+1/2"
                V_1[ii]   = V[ii] + 0.5* k_1V[ii]
                Q_1[ii]   = Q[ii] + 0.5* k_1Q[ii]  # <modify> We really don't need to define k_1v and k_1q, just for the clarity of the code.

            #print(k_1Q)

            for ii in range(N_Cells):  # These are the variables at {n+1}
                A_1[ii]   = V_1[ii] / L[ii]
                U_1[ii]   = Q_1[ii] / A_1[ii]
                Eta_1[ii] = A_1[ii] / B[ii] + Z[ii]
                E_1[ii]   = ((U_1[ii])**2.0) /2.0 +  Gravity*Eta_1[ii]
                l_P_1[ii] = B[ii] + 2.0 * (Eta_1[ii] - Z[ii])  # <modify>
                R_h_1[ii] = A_1[ii] / l_P_1[ii] # <modify>
                C_1[ii]   = ((M[ii])**2.0) / ((R_h_1[ii])**(4.0/3.0)) # <modify>
                F_1[ii]   = Gravity * C_1[ii] * V_1[ii] * ((U_1[ii])**2.0)

            if nn < slowness:
              Q_Upstream = Ex.Q_Up * ((nn+0.5)/ float(slowness))
            else:
              Q_Upstream = Ex.Q_Up


            # <delete>
            if (nn%100000) == 0:
                RealTime = round(nn*DT,5)
                TITLE = " at half time: " + str(RealTime)
                Draw.Plot_at_Cell(N_Cells, X, Z, Q_1, V_1, Eta_1, U_1, E_1, A_1, TITLE)

            for ii in range(N_Cells+1):
                if ii==0: # Boundary condition at face 1/2
                    if IntType == 3 or IntType == 4:
                        Q_F_1[ii]   = Q_Upstream
                        #A_F_1[ii]   = A_1[ii]
                        #U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                        #Eta_F_1[ii] = Eta_1[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))

                        Z_0 = Z[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        #Eta_F[ii] = h_upstream+ Z_0
                        Eta_F_1[ii] = Eta_1[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        A_F_1[ii]   = (Eta_F_1[ii] - Z_0)* B[ii]

                        U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                        E_F_1[ii]   = ((U_F_1[ii])**2.0)/2.0 + Gravity * Eta_F_1[ii]

                elif ii != 0 and ii != N_Cells: # middle cells
                    if IntType == 3:
                        A_F_1[ii]   = (L[ii]*  A_1[ii-1] + L[ii-1]*  A_1[ii] )/( L[ii] + L[ii-1] )
                        Eta_F_1[ii] = (L[ii]*Eta_1[ii-1] + L[ii-1]*Eta_1[ii] )/( L[ii] + L[ii-1] )
                        #if Eta_F_1[ii] != A_F_1[ii] / B[ii]  + Z_F[2*ii]:
                        #    print(" Fatal Error in the half step: inconsistency btw A and eta: %d %d %30.25f %30.25f %f " % (nn, ii, Eta_F_1[ii], A_F_1[ii] / B[ii]  + Z_F[2*ii] , Z_F[2*ii]))
                        #    check = input(" Error: press ENTER to exit ")
                        #    sys.exit()
                        # Temp: E_{i+1/2}^{n+1/2}
                        E_F_1[ii]   = ( 1.0 / (V_1[ii-1]+V_1[ii]) ) * ( E_F[ii] *( V[ii-1] + V[ii] ) - E_1[ii-1]*V_1[ii-1] - E_1[ii]*V_1[ii] + E[ii-1]*V[ii-1] + E[ii]*V[ii] + DT * ( Q_1[ii-1] * E_1[ii-1] - Q_1[ii] * E_1[ii] - (1.0/2.0) * ( U_1[ii-1] * F_1[ii-1] + U_1[ii] * F_1[ii]  ) + Q[ii-1] * E[ii-1] - Q[ii] * E[ii] - (1.0/2.0) * ( U[ii-1] * F[ii-1] + U[ii] * F[ii]  )   ) )
                        if E_F_1[ii]  < 0:
                            print(" Fatal Error: Negative Energy: %d, %d, %f " % (nn, ii, E_F_1[ii] ))
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()
                        if (E_F_1[ii] - Gravity * Eta_F_1[ii]) < 0.00000001:
                            E_F_1[ii] = Gravity * Eta_F_1[ii]
                        if (E_F_1[ii] - Gravity * Eta_F_1[ii]) < 0.0:
                            print(" Fatal Error in energy at n+half: %d, %d, %20.10f, %20.10f" % (nn, ii, E_F_1[ii], Gravity * Eta_F_1[ii] ))
                            print(" %30.20f" % (E_F_1[ii] - Gravity * Eta_F_1[ii]))
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()
                        if ( E_F_1[ii] < E_F_1[ii-1] ):
                            print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F_1[ii], E_F_1[ii-1] ) )
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()
                        U_F_1[ii]   = ( 2.0*(E_F_1[ii] - Gravity * Eta_F_1[ii] ) )**(0.5)
                        Q_F_1[ii]   = A_F_1[ii] * U_F_1[ii]          
                    elif IntType == 4:
                        A_F_1[ii]   = ( L[ii]*  A_1[ii-1] + L[ii-1]*  A_1[ii] )/( L[ii] + L[ii-1] ) # interpolation
                        U_F_1[ii]   = ( L[ii]*  U_1[ii-1] + L[ii-1]*  U_1[ii] )/( L[ii] + L[ii-1] ) # interpolation
                        Q_F_1[ii]   = A_F_1[ii] * U_F_1[ii]                                         # solve
                        E_F_1[ii]   = ( 1.0 / (V_1[ii-1]+V_1[ii]) ) * ( E_F[ii] *( V[ii-1] + V[ii] ) - E_1[ii-1]*V_1[ii-1] - E_1[ii]*V_1[ii] + E[ii-1]*V[ii-1] + E[ii]*V[ii] + DT * ( Q_1[ii-1] * E_1[ii-1] - Q_1[ii] * E_1[ii] - (1.0/2.0) * ( U_1[ii-1] * F_1[ii-1] + U_1[ii] * F_1[ii]  ) + Q[ii-1] * E[ii-1] - Q[ii] * E[ii] - (1.0/2.0) * ( U[ii-1] * F[ii-1] + U[ii] * F[ii]  )   ) )
                        Eta_F_1[ii] = (1.0/(2 * Gravity)) * ( 2 * E_F_1[ii] - (U_F_1[ii])**(2.0) )  # solve
                        if abs(Eta_F_1[ii] - (A_F_1[ii] / B[ii]  + Z_F[2*ii])) > 0.001:
                            print(" Fatal Error halfway: inconsistency btw A and eta: %d %d %30.20f %30.20f"% (nn, ii, Eta_F_1[ii], A_F_1[ii] / B[ii]  + Z_F[2*ii] ))
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()
                        if ( E_F_1[ii] < E_F_1[ii-1] ):
                            print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F_1[ii], E_F_1[ii-1] ) )
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()
                elif ii == N_Cells: # Boundary condition at face N+1/2
                    if IntType == 3 or IntType == 4:
                        Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )  # <remove> after debugging
                        Z_N1      = Z[ii-1] - Delta_Z  # <remove> after debugging
                        Eta_F_1[ii] = h_dw
                        A_F_1[ii]   = (Eta_F_1[ii] - Z_N1 ) * B[ii-1]
                        #if Eta_1[ii-1]<Eta_F_1[ii]:
                        #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta_1[ii-1],Eta_F[ii]))    
                        #    check = input(" Error: press ENTER to exit ")
                        #    sys.exit()
                        #Q_F_1[ii]   = (1.0/M[ii-1]) * A_F_1[ii] * ((R_h_1[ii-1])**(2.0/3.0))  * ( ( (Eta_1[ii-1]-Eta_F_1[ii])/(0.5*L[ii-1])  )**(1.0/2.0) )  # <modify>
                        Q_F_1[ii]   = Q_1[ii-1]
                        U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                        E_F_1[ii]   = ((U_F_1[ii])**2.0)/2.0 + Gravity * Eta_F_1[ii]
                        if ( E_F_1[ii] < E_F_1[ii-1] ):
                            print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F_1[ii], E_F_1[ii-1] ) )
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()

            # <delete>
            if (nn%100000) == 0:
                RealTime = round(nn*DT,5)
                TITLE = "k-2 at time: " + str(RealTime)
                Draw.Plot_Full(3, N_Cells, X_F, Z_F, Q_1, Q_F_1, Eta_1, Eta_F_1, U_1, U_F_1, E_1, E_F_1, A_1, A_F_1, TITLE)


            for ii in range(N_Cells): # To find k2 in the Runge-Kutta method and find the solution at n + 1   <remove> this for loop after debugging and substitute this with 
                V_0[ii]  = V[ii]
                Q_0[ii]  = Q[ii]
                E_0[ii]  = E[ii]
                F_0[ii]  = F[ii]
                U_0[ii]  = U[ii]
                E_F_0[ii]= E_F[ii]

                k_2V[ii]  = DT * ( Q_F_1[ii] - Q_F_1[ii+1] )  # <modify> remove
                k_2Q[ii]  = (DT / L[ii]) * ( Q_F_1[ii] * U_F_1[ii] - Q_F_1[ii+1] * U_F_1[ii+1] + Gravity * A_F_1[ii]* Eta_F_1[ii] - Gravity * A_F_1[ii+1]* Eta_F_1[ii+1]  -F_1[ii] ) # <modify> remove

                #print(k_2Q)

                V[ii] += k_2V[ii]
                Q[ii] += k_2Q[ii]


        # End loop on the time steps

        print(" ========== Solver Class ==========")
        print()




