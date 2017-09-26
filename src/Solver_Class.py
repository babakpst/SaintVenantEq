
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/31/2017
# Latest update: 09/20/2017
#
# Comment: This class conducts the numerical solution of the Saint-Venant equation
#
#####################################################################

class Solver:

    def __init__(self):
        # -- Import libs/classes
        import numpy as np
        import sys

        import Initialization_Class as Initial
        import Visualization_Class  as Visual

        Gravity = 9.81
        Draw = Visual.Visualization()

        #def RK2(self):
        Ex = Initial.Initialization()

        print(" ========== Solver Class ==========")
        print(" Solving the DE ...")
        print(" Second-order Runge-Kutta method: ...")

        print(" Allocating memory ...")
        N_Steps = int(Ex.Total_Time/Ex.Time_Step)
        N_Cells = Ex.N_Cells
        h_dw    = Ex.h_dw
        DT      = Ex.Time_Step

        Epsilon_E = 0.001

        L    = np.zeros(N_Cells,     dtype=np.float64 )
        Z    = np.zeros(N_Cells,     dtype=np.float64 )
        Z_F  = np.zeros(N_Cells*2+1, dtype=np.float64 )
        M    = np.zeros(N_Cells,     dtype=np.float64 )
        B    = np.zeros(N_Cells,     dtype=np.float64 )
        Fr   = np.zeros(N_Cells,     dtype=np.float64 )
        X    = np.zeros(N_Cells,     dtype=np.float64 )
        X_F  = np.zeros(N_Cells*2+1, dtype=np.float64 )
        
        Q    = np.zeros(N_Cells, dtype=np.float64 )
        V    = np.zeros(N_Cells, dtype=np.float64 )
        A    = np.zeros(N_Cells, dtype=np.float64 )
        U    = np.zeros(N_Cells, dtype=np.float64 )
        Eta  = np.zeros(N_Cells, dtype=np.float64 )
        E    = np.zeros(N_Cells, dtype=np.float64 )
        l_P  = np.zeros(N_Cells, dtype=np.float64 )
        R_h  = np.zeros(N_Cells, dtype=np.float64 )
        C    = np.zeros(N_Cells, dtype=np.float64 )
        Gamma= np.zeros(N_Cells, dtype=np.float64 )  # <modify> This needs to be modified if we are dealing with the variable channels.

        # <modify> double check the defined variables here. Delete if not necessary.

        Q_1   = np.zeros(N_Cells, dtype=np.float64 )
        V_1   = np.zeros(N_Cells, dtype=np.float64 )
        A_1   = np.zeros(N_Cells, dtype=np.float64 )
        U_1   = np.zeros(N_Cells, dtype=np.float64 )
        Eta_1 = np.zeros(N_Cells, dtype=np.float64 )
        E_1   = np.zeros(N_Cells, dtype=np.float64 )
        l_P_1 = np.zeros(N_Cells, dtype=np.float64 )
        R_h_1 = np.zeros(N_Cells, dtype=np.float64 )
        C_1   = np.zeros(N_Cells, dtype=np.float64 )

        Eta_F = np.zeros(N_Cells+1, dtype=np.float64 )
        A_F   = np.zeros(N_Cells+1, dtype=np.float64 )
        U_F   = np.zeros(N_Cells+1, dtype=np.float64 )
        E_F   = np.zeros(N_Cells+1, dtype=np.float64 )
        Q_F   = np.zeros(N_Cells+1, dtype=np.float64 )
        F_q   = np.zeros(N_Cells*2, dtype=np.float64 )

        A_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        Eta_F_1 = np.zeros(N_Cells+1, dtype=np.float64 )
        U_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        E_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        Q_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        F_q_1   = np.zeros(N_Cells*2, dtype=np.float64 )

        k_1V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array
        k_1Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array
        k_2V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array
        k_2Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove this array

        # Initialization 
        print(" Initialization ... ")
        Q[:]   = Ex.Q[:]
        V[:]   = Ex.V[:]
        L[:]   = Ex.L[:]
        Z[:]   = Ex.Z[:]
        Z_F[:] = Ex.Z_F[:]
        M[:]   = Ex.M[:]
        B[:]   = Ex.B[:]
        X[:]   = Ex.X[:]
        X_F[:] = Ex.X_F[:]

        slowness = 0
        Plot1 = 2000
        Plot2 = 2000
        h_upstream = V[0]/(B[0]*L[0])


        Del_Coe = 2.0  # <delete> this coefficient in the entire code once you finalize the system


        print(" Time marching ... ")
        for nn in range(N_Steps):
            print(" Time step: %d out of %d " % (nn, N_Steps))

            if nn < slowness:
              Q_Upstream = Ex.Q_Up * (nn/ float(slowness))
            else:
              Q_Upstream = Ex.Q_Up        
        
            #########################################################
            Check_Energy = "False"
            while Check_Energy == "False":
                for ii in range(N_Cells):
                    A[ii]   = V[ii] / L[ii]
                    U[ii]   = Q[ii] / A[ii]
                    Eta[ii] = A[ii] / B[ii] + Z[ii]
                    E[ii]   = ((U[ii])**2.0) /(float(2)) +  Gravity*Eta[ii]
                    l_P[ii] = B[ii] + 2.0 * (Eta[ii]-Z[ii])
                    R_h[ii] = A[ii] / l_P[ii]
                    C[ii]   = ((M[ii])**2.0) / ((R_h[ii])**(4.0/3.0))
  
                    Fr[ii]  = U[ii]/((Gravity * (Eta[ii] - Z[ii]) )**(0.5))
                    if Fr[ii] >= 1.0:
                        print("Flow is not subcritical")
                        check = input(" Error: press ENTER to exit ")
                        sys.exit()

                # <delete>
                if (nn%Plot1) == 0:
                    print(Q_Upstream)
                    RealTime = round(nn*DT,5)
                    TITLE1 = Ex.Output_Dir + "/" +  "Time_" + str(RealTime)+"_s"
                    TITLE2 = "at time: " + str(RealTime) + " s"
                    Draw.Plot_at_Cell(N_Cells, X, Z, Q, V, Eta, U, E, A, TITLE1, TITLE2)

                for ii in range(N_Cells+1):
                    if ii==0: # Boundary condition at face 1/2
                        Q_F[ii]   = Q_Upstream
                        Z_0 = Z[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        A_F[ii]   =  (Eta_F[ii] - Z_0)* B[ii]
                        U_F[ii]   = Q_F[ii] / A_F[ii]
                        E_F[ii]   = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]

                        # <modify> We assume that the boundary faces do not need any modifications and are fixed.
                        #Eta_F_hat[ii] = Eta_F[ii]
                        #Q_F_hat_S[ii] = (Q_F[ii])**2.0
                        #A_F_hat[ii]   = A_F[ii]
                        #B_F_hat[ii]   = B[ii]

                    elif ii != 0 and ii != N_Cells: # middle cells - The subtraction is due to the fact that Python numbering is starting from 0
                        E_F[ii]        = (L[ii]*  E[ii-1] + L[ii-1]*  E[ii] )/( L[ii] + L[ii-1] )
                        Eta_F_hat      = (L[ii]*Eta[ii-1] + L[ii-1]*Eta[ii] )/( L[ii] + L[ii-1] )
                        Q_F_hat_S      = (L[ii]*  ((Q[ii-1])**2.0) + L[ii-1]*  ((Q[ii])**2.0) )/( L[ii] + L[ii-1] )
                        A_F_hat        = (L[ii]*  A[ii-1] + L[ii-1]*  A[ii] )/( L[ii] + L[ii-1] ) # <modify> Modify this equation for a variable area
                        B_F_hat        = B[ii] # <modify> Modify this equation for a variable area
                        Gamma_F_hat    = (L[ii]*Gamma[ii-1] + L[ii-1]*  Gamma[ii] )/( L[ii] + L[ii-1] ) # <modify> Modify this equation for a variable area
                        H_hat          = A_F_hat / B_F_hat
                        a              = 1 + (2.0 * Gamma_F_hat * (H_hat**2.0) )/(A_F_hat) - ( Gravity * Eta_F_hat/E_F[ii]  ) * ( 1 + 2.0*Gamma_F_hat*(H_hat**2.0)/A_F_hat ) - 2.0 * Gravity * H_hat/E_F[ii]
                        b              = 2.0 - Gravity * ( 2.0 * Eta_F_hat + H_hat ) / E_F[ii]
                        c              = 1 - (Gravity * Eta_F_hat/E_F[ii] ) - Q_F_hat_S/( 2.0 * E_F[ii] * (A_F_hat**2.0) )

                        if ((b+1)**(2.0) - 4.0 * a * c )< 0.0:
                            print(" Fatal error: negative" )
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()

                        Eta_Epsilon1   = ( A_F_hat / B_F_hat ) * ( -b-1 + (  ( (b+1)**2.0 - 4.0 * a * c )**0.5 ) ) / ( 2 * a ) 
                        Eta_Epsilon2   = ( A_F_hat / B_F_hat ) * ( -b-1 - (  ( (b+1)**2.0 - 4.0 * a * c )**0.5 ) ) / ( 2 * a ) 
                    
                        if abs(Eta_Epsilon1) < abs(Eta_Epsilon2):
                            Eta_Epsilon    =  Eta_Epsilon1
                        else:
                            Eta_Epsilon    =  Eta_Epsilon2

                        #print("  Eta_Epsilon  %3d: %30.20f" % (ii,Eta_Epsilon))  # <delete> delete after debugging

                        Eta_F[ii]  = Eta_F_hat + Eta_Epsilon
                        A_F[ii]    = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon

                        x = Eta_Epsilon/H_hat
                        Alfa_Epsilon  = 2 * E_F[ii] * (A_F_hat**2.0) * ( a * (x**2.0) + b * x + c )

                        #print(" Alfa_Epsilon %d:  %25.20f %25.20f %25.20f %25.20f %25.20f %25.20f %25.20f" % (ii, E_F[ii], A_F_hat, a, b, c, x, Alfa_Epsilon))  # <delete> delete after debugging

                        if ( Q_F_hat_S + Alfa_Epsilon ) < 0.0:
                            print("  Warning %5d %5d" % (nn,ii ))
                            Alfa_Epsilon = 0.0
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()                        

                        Q_F[ii] = ( Q_F_hat_S + Alfa_Epsilon )**(0.5)
                        Q_check = ( 2*( (A_F[ii])**(2.0) ) * ( E_F[ii]-Gravity*Eta_F[ii] )   ) **(0.5)

                        if abs(Q_check- Q_F[ii]) >0.01:
                            print(' Error: Q at the face is not consistent, %5d %5d %30.20f %30.20f' % (nn,ii,Q_check, Q_F[ii]))
                            #check = input(" Error: press ENTER to exit ")
                            #sys.exit()
                        U_F[ii]   = Q_F[ii] / A_F[ii]

                    elif ii == N_Cells: # Boundary condition at face N+1/2
                        Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
                        Z_N1      = Z[ii-1] - Delta_Z
                        Z_N1      = 0.0   # <modify>
                        Eta_F[ii] = h_dw + Z_N1
                        A_F[ii]   = h_dw * B[ii-1]
                        #if Eta[ii-1]<Eta_F[ii]:
                        #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta[ii-1],Eta_F[ii]))    
                        #    check = input(" Error: press ENTER to exit ")
                        #    sys.exit()
                        #Q_F[ii]   = (1.0/M[ii-1]) * A_F[ii] * ((R_h[ii-1])**(2.0/3.0))  * ( ( (Eta[ii-1]-Eta_F[ii])/(0.5*L[ii-1])  )**(1.0/2.0) ) # <modify>
                        Q_F[ii]   = Q[ii-1]
                        U_F[ii]   = Q_F[ii] / A_F[ii]
                        E_F[ii]   = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]
                        #if ( E_F[ii] < E_F[ii-1] ):
                        #    print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], E_F[ii-1] ) )
                        #    check = input(" Error: press ENTER to exit ")
                        #    sys.exit()

                    # Modifying the Energy at the face

                Check_Energy = "True"
                for ii in range(N_Cells):
                    if E[ii] > E[ii-1]:
                        Check_Energy = "False"
                        print(" Modifying the solution")
                        Delta_Energy = E[ii-1] - E[ii] 
                        Delta_Q = (Epsilon_E - Delta_Energy) / ( DT * ( ( Del_Coe*U_F[ii] -U[ii-1] ) * ( U[ii-1]/V[ii-1] ) +  ( Del_Coe*U_F[ii] -U[ii] ) * ( U[ii]/V[ii] ) + Gravity * (1.0/(B[ii-1]*L[ii-1]) + 1.0/(B[ii]*L[ii])  ) ) )

                        Q[ii]   = Q[ii]  - Del_Coe * DT/L[ii] * Delta_Q * U_F[ii]
                        Q[ii-1] = Q[ii-1]+ Del_Coe * DT/L[ii] * Delta_Q * U_F[ii]

                        V[ii]   = V[ii]  - DT * Delta_Q 
                        V[ii-1] = V[ii-1]+ DT * Delta_Q 


            # <delete>
            if (nn%Plot2) == 0:
                RealTime = round(nn*DT,5)
                TITLE1 = Ex.Output_Dir + "/" +  "Full__time_" + str(RealTime)+"_s"
                TITLE2 = "Full results at time: " + str(RealTime) + " s"
                Draw.Plot_Full(2,N_Cells, X_F, Z_F, Q, Q_F, Eta, Eta_F, U, U_F, E, E_F, A, A_F, TITLE1, TITLE2)

            for ii in range(N_Cells):
                F_q[ii*2  ] = Gravity * C[ii] * V[ii] * ( ( U_F[ii] + U[ii]     )**2.0) / 8.0
                F_q[ii*2+1] = Gravity * C[ii] * V[ii] * ( ( U[ii]   + U_F[ii+1] )**2.0) / 8.0

            for ii in range(N_Cells): # To find k1 in the Runge-Kutta method and find the solution at n+1/2
                k_1V[ii]  =  DT * ( Q_F[ii] - Q_F[ii+1] )  # <modify> remove
                k_1Q[ii]  = ( DT / L[ii] ) * ( Q_F[ii] * U_F[ii] - Q_F[ii+1] * U_F[ii+1] + Gravity * A_F[ii]* Eta_F[ii] - Gravity * A_F[ii+1]* Eta_F[ii+1] - F_q[2*ii  ] - F_q[2*ii+1] ) # <modify> remove
                # Solution at "n+1/2"
                V_1[ii]   = V[ii] + 0.5* k_1V[ii]
                Q_1[ii]   = Q[ii] + 0.5* k_1Q[ii]  # <modify> We really don't need to define k_1v and k_1q, just for the clarity of the code.

            if nn < slowness:
                Q_Upstream = Ex.Q_Up * ((nn+0.5)/ float(slowness))
            else:
                Q_Upstream = Ex.Q_Up

            Check_Energy = "False"
            while Check_Energy == "False":
                for ii in range(N_Cells):  # These are the variables at {n+1}
                    A_1[ii]   = V_1[ii] / L[ii]
                    U_1[ii]   = Q_1[ii] / A_1[ii]
                    Eta_1[ii] = A_1[ii] / B[ii] + Z[ii]
                    E_1[ii]   = ((U_1[ii])**2.0) /2.0 +  Gravity*Eta_1[ii]
                    l_P_1[ii] = B[ii] + 2.0 * (Eta_1[ii] - Z[ii])  # <modify>
                    R_h_1[ii] = A_1[ii] / l_P_1[ii] # <modify>
                    C_1[ii]   = ((M[ii])**2.0) / ((R_h_1[ii])**(4.0/3.0)) # <modify>

                    Fr[ii]  = U_1[ii]/((Gravity * (Eta_1[ii] - Z[ii]) )**(0.5))
                    if Fr[ii] >= 1.0:
                        print("Flow is not subcritical")
                        check = input(" Error: press ENTER to exit ")
                        sys.exit()

                # <delete>
                #if (nn%Plot2) == 0:
                #    RealTime = round(nn*DT,5)
                #    TITLE = " at half time: " + str(RealTime)
                #    Draw.Plot_at_Cell(N_Cells, X, Z, Q_1, V_1, Eta_1, U_1, E_1, A_1, TITLE)

                for ii in range(N_Cells+1):
                    if ii==0: # Boundary condition at face 1/2
                        Q_F_1[ii]   = Q_Upstream
                        Z_0 = Z[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        Eta_F_1[ii] = Eta_1[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                        A_F_1[ii]   = (Eta_F_1[ii] - Z_0)* B[ii]
                        U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                        E_F_1[ii]   = ((U_F_1[ii])**2.0)/(2.0) + Gravity * Eta_F_1[ii]

                    elif ii != 0 and ii != N_Cells: # middle cells
                        E_F_1[ii]        = (L[ii]*  E_1[ii-1] + L[ii-1]*  E_1[ii] )/( L[ii] + L[ii-1] )
                        Eta_F_hat      = (L[ii]*Eta_1[ii-1] + L[ii-1]*Eta_1[ii] )/( L[ii] + L[ii-1] )
                        Q_F_hat_S      = (L[ii]*  ((Q_1[ii-1])**2.0) + L[ii-1]*  ((Q_1[ii])**2.0) )/( L[ii] + L[ii-1] )
                        A_F_hat        = (L[ii]*  A_1[ii-1] + L[ii-1]*  A_1[ii] )/( L[ii] + L[ii-1] ) # <modify> Modify this equation for a variable area
                        B_F_hat        = B[ii] # <modify> Modify this equation for a variable area
                        Gamma_F_hat    = (L[ii]*Gamma[ii-1] + L[ii-1]*  Gamma[ii] )/( L[ii] + L[ii-1] ) # <modify> Modify this equation for a variable area
                        H_hat          = A_F_hat / B_F_hat
                        a              = 1 + (2.0 * Gamma_F_hat * (H_hat**2.0) )/(A_F_hat) - (  (Gravity * Eta_F_hat)/E_F_1[ii]  ) * (1 + (2.0*Gamma_F_hat*(H_hat**2.0))/A_F_hat ) - 2.0 * Gravity * H_hat/E_F_1[ii]
                        b              = 2.0 - Gravity * ( 2.0 * Eta_F_hat + H_hat ) / E_F_1[ii]
                        c              = 1 - (Gravity * Eta_F_hat/E_F_1[ii] ) - Q_F_hat_S/( 2.0 * E_F_1[ii] * (A_F_hat**2.0) )

                        if ((b+1)**(2.0) - 4.0 * a * c) < 0.0:
                            print(" Fatal error: negative" )
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()

                        Eta_Epsilon1   = ( A_F_hat / B_F_hat ) * ( -b-1 + ( ((b+1)**(2.0) - 4.0 * a * c )**(0.5) ))  / ( 2 * a ) 
                        Eta_Epsilon2   = ( A_F_hat / B_F_hat ) * ( -b-1 - ( ((b+1)**(2.0) - 4.0 * a * c )**(0.5) ))  / ( 2 * a ) 
                    
                        if abs(Eta_Epsilon1) < abs(Eta_Epsilon2):
                            Eta_Epsilon    =  Eta_Epsilon1
                        else:
                            Eta_Epsilon    =  Eta_Epsilon2
                        
                        Eta_F_1[ii]  = Eta_F_hat + Eta_Epsilon
                        A_F_1[ii]    = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon

                        x = Eta_Epsilon/H_hat
                        Alfa_Epsilon  = 2 * E_F_1[ii] * (A_F_hat**2.0) * ( a * (x**2.0) + b * x + c )

                        if ( Q_F_hat_S + Alfa_Epsilon ) < 0.0:
                            print("  Warning at the half step %5d %5d" % (nn,ii ))
                            Alfa_Epsilon = 0.0
                            check = input(" Error: press ENTER to exit ")
                            sys.exit()                        

                        Q_F_1[ii] = ( Q_F_hat_S+Alfa_Epsilon )**(0.5)
                        Q_check = ( 2*((A_F_1[ii])**(2.0)) * ( E_F_1[ii]-Gravity*Eta_F_1[ii] )   ) **(0.5)

                        if abs(Q_check- Q_F[ii]) > 0.01:
                            print(' Error: Q at the face is not consistent, %5d %5d %30.20f %30.20f' % (nn,ii,Q_check, Q_F[ii]))
                            #check = input(" Error: press ENTER to exit ")
                            #sys.exit()
                        U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]

                    elif ii == N_Cells: # Boundary condition at face N+1/2
                        Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )  # <remove> after debugging
                        Z_N1      = Z[ii-1] - Delta_Z  # <remove> after debugging
                        Z_N1      = 0.0 # <modify>
                        Eta_F_1[ii] = h_dw + Z_N1
                        A_F_1[ii]   = h_dw * B[ii-1]
                        #if Eta_1[ii-1]<Eta_F_1[ii]:
                        #    print(" Fatal error: Downstream boundary condition: %d , %d , %f , %f" % (ii,nn,Eta_1[ii-1],Eta_F[ii]))    
                        #    check = input(" Error: press ENTER to exit ")
                        #    sys.exit()
                        #Q_F_1[ii]   = (1.0/M[ii-1]) * A_F_1[ii] * ((R_h_1[ii-1])**(2.0/3.0))  * ( ( (Eta_1[ii-1]-Eta_F_1[ii])/(0.5*L[ii-1])  )**(1.0/2.0) )  # <modify>
                        Q_F_1[ii]   = Q_1[ii-1]
                        U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                        #E_F_1[ii]   = ((U_F_1[ii])**2.0)/(2.0*Gravity) + Eta_F_1[ii]
                        E_F_1[ii]   = ((U_F_1[ii])**2.0)/2.0 + Gravity * Eta_F_1[ii]
                        #if ( E_F_1[ii] < E_F_1[ii-1] ):
                        #    print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F_1[ii], E_F_1[ii-1] ) )
                        #    check = input(" Error: press ENTER to exit ")
                        #    sys.exit()

                Check_Energy = "True"
                for ii in range(N_Cells):
                    if E_1[ii] > E_1[ii-1]:
                        print(" Modifying the solution at the half time step")
                        Check_Energy = "False"
                        Delta_Energy = E_1[ii-1] - E_1[ii] 
                        Delta_Q = (Epsilon_E - Delta_Energy) / ( DT * ( ( Del_Coe*U_F_1[ii] - U_1[ii-1] ) * ( U_1[ii-1]/V_1[ii-1] ) +  ( Del_Coe*U_F_1[ii] -U_1[ii] ) * ( U_1[ii]/V_1[ii] ) + Gravity * (1.0/(B[ii-1]*L[ii-1]) + 1.0/(B[ii]*L[ii])  ) ) )

                        Q_1[ii]   = Q_1[ii]  - Del_Coe * DT/L[ii] * Delta_Q * U_F_1[ii]
                        Q_1[ii-1] = Q_1[ii-1]+ Del_Coe * DT/L[ii] * Delta_Q * U_F_1[ii]

                        V_1[ii]   = V_1[ii]  - DT * Delta_Q 
                        V_1[ii-1] = V_1[ii-1]+ DT * Delta_Q 

            # <delete>
            #if (nn%Plot2) == 0:
            #    RealTime = round(nn*DT,5)
            #    TITLE = "k-2 at time: " + str(RealTime)
            #    Draw.Plot_Full(3, N_Cells, X_F, Z_F, Q_1, Q_F_1, Eta_1, Eta_F_1, U_1, U_F_1, E_1, E_F_1, A_1, A_F_1, TITLE)
            for ii in range(N_Cells):
                F_q_1[ii*2  ] = Gravity * C_1[ii] * V_1[ii] * ( ( U_F_1[ii] + U_1[ii]     )**2.0) / 8.0
                F_q_1[ii*2+1] = Gravity * C_1[ii] * V_1[ii] * ( ( U_1[ii]   + U_F_1[ii+1] )**2.0) / 8.0

            for ii in range(N_Cells): # To find k2 in the Runge-Kutta method and find the solution at n + 1   <remove> this for loop after debugging and substitute this with 
                k_2V[ii]  = DT * ( Q_F_1[ii] - Q_F_1[ii+1] )  # <modify> remove
                k_2Q[ii]  = (DT / L[ii]) * ( Q_F_1[ii] * U_F_1[ii] - Q_F_1[ii+1] * U_F_1[ii+1] + Gravity * A_F_1[ii]* Eta_F_1[ii] - Gravity * A_F_1[ii+1]* Eta_F_1[ii+1]   - F_q_1[2*ii] - F_q_1[2*ii+1] ) # <modify> remove
 
                #print(k_2Q)

                V[ii] += k_2V[ii]
                Q[ii] += k_2Q[ii]


        # End loop on the time steps

        print(" ========== Solver Class ==========")
        print()




