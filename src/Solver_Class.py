
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/31/2017
# Latest update: 11/20/2017
#
# Comment: This class conducts the numerical solution of the Saint-Venant equation
#
#######################################################################################################################

class Solver:

    def __init__(self):
        pass

    def sgn_func(self, x): return 1.0 if x >= 0 else -1.0

    def solve(self):
        import numpy as np
        import sys

        import Initialization_Class as Initial
        import Visualization_Class  as Visual

        Gravity = 9.81
        Draw = Visual.Visualization()

        Ex = Initial.Initialization()

        print(" ========== Solver Class ==========")
        print(" Solving SVE ...")
        print(" Second-order Runge-Kutta method: ...")

        print(" Allocating memory ...")
        N_Steps = int(Ex.Total_Time/Ex.Time_Step)
        N_Cells = Ex.N_Cells
        h_dw    = Ex.h_dw
        DT      = Ex.Time_Step

        Epsilon_E = 1.0e-3

        L    = np.zeros(N_Cells,     dtype=np.float64 )
        HL   = np.zeros(N_Cells,     dtype=np.float64 )
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
        Gamma= np.zeros(N_Cells, dtype=np.float64 ) #<modify> for variable area

        # <modify> double check.

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

        Eta_F_1 = np.zeros(N_Cells+1, dtype=np.float64 )
        A_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        U_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        E_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        Q_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        F_q_1   = np.zeros(N_Cells*2, dtype=np.float64 )

        k_1V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_1Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_2V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_2Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove

        # Initialization 
        print(" Initialization ... ")
        Q[:]   = Ex.Q[:]
        V[:]   = Ex.V[:]
        L[:]   = Ex.L[:]
        HL[:]  = Ex.HL[:]
        Z[:]   = Ex.Z[:]
        Z_F[:] = Ex.Z_F[:]
        M[:]   = Ex.M[:]
        B[:]   = Ex.B[:]
        X[:]   = Ex.X[:]
        X_F[:] = Ex.X_F[:]

        #slowness = 0 # 
        Plot1 = 20  # Plots the results at the cell center, every "Plot1" steps.
        Plot2 = 100 # Plots the full results at the cell center, every "Plot2" steps.
        h_upstream = V[0]/(B[0]*HL[0])

        #TITLE = "This is it." # <delete> after debugging
        #Draw.Plot_Domain(N_Cells, X, Z, TITLE)

        #check = input(" Enter to continue ...")

        print(" Time marching ... ")
        for nn in range(N_Steps):
            print(" Time step: %d out of %d " % (nn, N_Steps))

            #if nn < slowness:
            #    Q_Upstream = Ex.Q_Up * (nn/ float(slowness))
            #else:
            #    Q_Upstream = Ex.Q_Up        
            Q_Upstream = Ex.Q_Up        
        
            # Calculate the parameters at the cell center for time step n
            for ii in range(N_Cells): # see the paper:
                A[ii]   = V[ii] / HL[ii]
                U[ii]   = Q[ii] / A[ii]
                Eta[ii] = A[ii] / B[ii] + Z[ii]
                E[ii]   = (U[ii]**2.0) /2.0 + Gravity*Eta[ii]
                l_P[ii] = B[ii] + 2.0*(Eta[ii]-Z[ii])
                R_h[ii] = A[ii]/l_P[ii]
                C[ii]   = ((M[ii])**2.0) / (R_h[ii]**(4.0/3.0))

                Fr[ii]  = U[ii]/((Gravity * (Eta[ii] - Z[ii]) )**0.5)
                if Fr[ii] >= 1.0:
                    print("Flow is not subcritical %d %d %f" % (nn, ii, Fr[ii]))
                    check = input(" Error: press ENTER to exit ")
                    sys.exit()

            # Plot the results at the cell center at every "Plot1" steps at the time step n
            if (nn%Plot1) == 0:
                RealTime = round(nn*DT,5)
                #TITLE1 = Ex.Output_Dir + "/" +  "Time_" + str(RealTime)+"_s__Repeat_" + str(repeat)
                #TITLE2 = "at time: " + str(RealTime) +"_s__Repeat_" + str(repeat)
                TITLE1 = Ex.Output_Dir + "/" +  "Time_" + str(RealTime)+"_s"
                TITLE2 = "at time: " + str(RealTime) +"_s"
                Draw.Plot_at_Cell(N_Cells, X, Z, Q, V, Eta, U, E, A, TITLE1, TITLE2)

            #check = input(" Enter to continue ...") # <delete>
            # This loop evaluates all the required values at the faces at time step n
            for ii in range(N_Cells+1): 
                if ii==0: # Boundary condition at face 1/2
                    Q_F[ii] = Q_Upstream
                    Z_0 = Z[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                    Eta_F[ii] = Eta[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                    A_F[ii] =  (Eta_F[ii] - Z_0)* B[ii]
                    U_F[ii] = Q_F[ii] / A_F[ii]
                    E_F[ii] = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]

                    # <modify> We assume that the boundary faces do not need any modifications and are fixed.

                elif ii != 0 and ii != N_Cells: # middle cells (non-boundary cells)
                    E_F[ii]   = (L[ii]*  E[ii-1] + L[ii-1]*  E[ii])/(L[ii] + L[ii-1])
                    Eta_F_hat = (L[ii]*Eta[ii-1] + L[ii-1]*Eta[ii])/(L[ii] + L[ii-1])
                    #print("{:} {:20.15f}".format(" Eta_Hat", Eta_F_hat))
                    #A_F_hat   = (L[ii]*  A[ii-1] + L[ii-1]*  A[ii])/(L[ii] + L[ii-1]) # <modify>for variable area
                    Z_F0 = Z[ii] + L[ii] * (( Z[ii-1] - Z[ii]) / (L[ii-1] + L[ii]))
                    B_F_hat = B[ii] # <modify> Modify this equation for a variable area
                    A_F_hat   = (Eta_F_hat - Z_F0) * B_F_hat
                    Q_F_hat_S = (self.sgn_func(Q[ii-1])*L[ii]*((Q[ii-1])**2.0)  \
                                +self.sgn_func(Q[ii])*L[ii-1]*((Q[ii])**2.0))/(L[ii]+L[ii-1])
                    
                    Gamma_F_hat = (L[ii]*Gamma[ii-1] + L[ii-1]*Gamma[ii])/(L[ii]+L[ii-1]) #<modify>for variable area
                    H_hat = A_F_hat/B_F_hat
                    a = 1.0 \
                        + (2.0 * Gamma_F_hat * (H_hat**2.0) )/A_F_hat \
                        - (Gravity*Eta_F_hat/E_F[ii])*(1.0 + 2.0*Gamma_F_hat*(H_hat**2.0)/A_F_hat) \
                        - 2.0*Gravity*H_hat/E_F[ii]
                    b = 2.0 - Gravity * ( 2.0 * Eta_F_hat + H_hat ) / E_F[ii]
                    c = 1.0 - (Gravity * Eta_F_hat/E_F[ii] ) - Q_F_hat_S/( 2.0 * E_F[ii] * (A_F_hat**2.0) )

                    if ((b+1.0)**2.0 - 4.0 * a * c ) >= 0.0:
                        # print(" Positive/Zero argument a b c: %f %f %f" % (a,b,c)) # <delete> after debugging

                        Eta_Epsilon1 = H_hat*(-b-1.0 + (((b+1.0)**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                        Eta_Epsilon2 = H_hat*(-b-1.0 - (((b+1.0)**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                    
                        if abs(Eta_Epsilon1) < abs(Eta_Epsilon2): # Picks the closest to zero
                            Eta_Epsilon    =  Eta_Epsilon1
                        else:
                            Eta_Epsilon    =  Eta_Epsilon2

                        Eta_F[ii] = Eta_F_hat + Eta_Epsilon
                        A_F[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon

                        x = Eta_Epsilon/H_hat
                        Alfa_Epsilon  = 2.0 * E_F[ii] * (A_F_hat**2.0) * ( a * (x**2.0) + b * x + c )

                        if Q_F_hat_S > 0.0:
                            print(" Q is positive ", (nn, ii))
                            if -Q_F_hat_S < Alfa_Epsilon:
                                print(" Modification on Alpha")
                                Alfa_Epsilon = -Q_F_hat_S

                        elif Q_F_hat_S == 0.0:
                            #print(" Q is zero ", (nn, ii)) # <delete> after debugging
                            if Alfa_Epsilon < 0.0:
                                print(" Modification on Alpha- originally negative")
                                Alfa_Epsilon = 0.0
                        elif Q_F_hat_S < 0.0:
                            print(" Q is negative")
                            #check = input(" Error: press ENTER to exit ")
                            #sys.exit()

                        Q_F[ii] = self.sgn_func(Q_F_hat_S)*(Q_F_hat_S + Alfa_Epsilon)**0.5
                        Q_check = self.sgn_func(Q_F_hat_S)*A_F[ii]* (2*(E_F[ii]-Gravity*Eta_F[ii]))**(0.5)
                        
                        if abs(Q_check- Q_F[ii]) >0.01:
                            print("{40}{:5d}{:5d}{:30.20f}{:30.20f}".\
                                format("Error: Q at the face is not consistent",nn,ii,Q_check, Q_F[ii]))
                            #check = input(" Error: inconsistency in cell values - Enter to continue ...")
                            sys.exit()

                    elif ( (b+1.0)**2.0 - 4.0 * a * c ) < 0.0:
                        print(" Warning: negative argument- a b c: %f %f %f" % (a,b,c))
                        if (Gravity*Eta_F_hat > 0.5*Q_F_hat_S/(A_F_hat**2.0)) and (b**2.0 - 4.0 *a*c >= 0.0):

                            Eta_Epsilon1 = H_hat*(-b + ((b**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                            Eta_Epsilon2 = H_hat*(-b - ((b**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                        
                            if abs(Eta_Epsilon1) < abs(Eta_Epsilon2):
                                Eta_Epsilon    =  Eta_Epsilon1
                            else:
                                Eta_Epsilon    =  Eta_Epsilon2

                            Eta_F[ii] = Eta_F_hat + Eta_Epsilon
                            A_F[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon

                            #x = Eta_Epsilon/H_hat
                            #Alfa_Epsilon  = 2.0 * E_F[ii] * (A_F_hat**2.0) * ( a * (x**2.0) + b * x + c )
                            Alfa_Epsilon  = 0.0

                            Q_F[ii] = self.sgn_func(Q_F_hat_S)*(Q_F_hat_S + Alfa_Epsilon)**0.5
                            Q_check = self.sgn_func(Q_F_hat_S)*A_F[ii]* (2*(E_F[ii]-Gravity*Eta_F[ii]))**(0.5)
                            
                            if abs(Q_check- Q_F[ii]) >0.01:
                                print("{40}{:5d}{:5d}{:30.20f}{:30.20f}".\
                                    format("Error: Q at the face is not consistent",nn,ii,Q_check, Q_F[ii]))
                                #check = input(" Error: inconsistency in cell values - Enter to continue ...")
                                #check = input(" Error: press ENTER to exit ")
                                #sys.exit()
                        else:
                            Eta_Epsilon = 0.0
                            Eta_F[ii] = Eta_F_hat + Eta_Epsilon
                            A_F[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon
                            Alfa_Epsilon  = 2.0 * (A_F[ii]**2.0) * (E_F[ii] - Gravity*Eta_F[ii]) - Q_F_hat_S
                            if Alfa_Epsilon < - Q_F_hat_S:
                                Alfa_Epsilon =- Q_F_hat_S
                                Eta_Epsilon =  ( E_F[ii] - Gravity * Eta_F_hat ) / Gravity
                                Eta_F[ii] = Eta_F_hat + Eta_Epsilon
                                A_F[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon
                            
                            Q_F[ii] = self.sgn_func(Q_F_hat_S)*(Q_F_hat_S + Alfa_Epsilon)**0.5
                            Q_check = self.sgn_func(Q_F_hat_S)*A_F[ii]* (2*(E_F[ii]-Gravity*Eta_F[ii]))**(0.5)
                            
                            if abs(Q_check- Q_F[ii]) >0.01:
                                print("{40}{:5d}{:5d}{:30.20f}{:30.20f}".\
                                    format("Error: Q at the face is not consistent",nn,ii,Q_check, Q_F[ii]))
                                check = input(" Error: press ENTER to exit ")
                                sys.exit()

                    #print("{:} {:20.15f}".format(" estimated face: ", Eta_F[ii]) )
                    U_F[ii] = Q_F[ii]/A_F[ii]

                elif ii == N_Cells: # Boundary condition at face N+1/2
                    #Delta_Z   = L[ii-1] * ( ( Z[ii-2]-Z[ii-1] )/( L[ii-2]+L[ii-1] )  )
                    #Z_N1      = Z[ii-1] - Delta_Z
                    Z_N1      = 0.0   # <modify>
                    Eta_F[ii] = h_dw + Z_N1
                    A_F[ii]   = h_dw * B[ii-1]
                    Q_F[ii]   = Q[ii-1]
                    U_F[ii]   = Q_F[ii] / A_F[ii]
                    E_F[ii]   = ((U_F[ii])**2.0)/2.0 + Gravity * Eta_F[ii]
                    #if ( E_F[ii] < E_F[ii-1] ):
                    #    print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F[ii], E_F[ii-1] ) )
                    #    check = input(" Error: press ENTER to exit ")
                    #    sys.exit()


            # Plot the full results at every "Plot2" steps at the time step n
            if (nn%Plot2) == 0:
                RealTime = round(nn*DT,5)
                #TITLE1 = Ex.Output_Dir + "/" +  "Full__time_" + str(RealTime)+"_s__Repeat_" + str(repeat)
                #TITLE2 = "Full results at time: " + str(RealTime) +"_s__Repeat_" + str(repeat)
                TITLE1 = Ex.Output_Dir + "/" +  "Full__time_" + str(RealTime)+"_s"
                TITLE2 = "Full results at time: " + str(RealTime) +"_s"
                Draw.Plot_Full(2, N_Cells, X, X_F, Z_F, V, Q, Q_F, Eta, Eta_F, U, U_F, E, E_F, A, A_F, TITLE1, TITLE2)

            #check = input(" Enter to continue ...") # <delete> after debug

            for ii in range(N_Cells):
                F_q[ii*2  ] = Gravity * C[ii] * V[ii] * ( ( U_F[ii] + U[ii]     )**2.0) / 8.0
                F_q[ii*2+1] = Gravity * C[ii] * V[ii] * ( ( U[ii]   + U_F[ii+1] )**2.0) / 8.0

            # find k1 in the Runge-Kutta method and find the solution at n+1/2
            for ii in range(N_Cells): 
                k_1V[ii] =  DT * ( Q_F[ii] - Q_F[ii+1] )  # <modify> remove
                k_1Q[ii] = (DT/L[ii])*(Q_F[ii]*U_F[ii] - Q_F[ii+1]*U_F[ii+1] \
                           + Gravity*A_F[ii]*Eta_F[ii] - Gravity*A_F[ii+1]*Eta_F[ii+1] - F_q[2*ii] - F_q[2*ii+1]) 
                # Solution at "n+1/2"
                V_1[ii]   = V[ii] + 0.5* k_1V[ii]
                Q_1[ii]   = Q[ii] + 0.5* k_1Q[ii]  # <modify> not necessary

            #if nn < slowness:
            #    Q_Upstream = Ex.Q_Up * ((nn+0.5)/ float(slowness))
            #else:
            #    Q_Upstream = Ex.Q_Up
            #Q_Upstream = Ex.Q_Up        

            #check = input(" Enter to continue ...") # <delete>
            # Calculate the parameters at the cell center for time step n+1/2            
            for ii in range(N_Cells):  # These are the variables at {n+1}
                A_1[ii]   = V_1[ii] / HL[ii]
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

            # Plot the results at the cell center at every "Plot2" steps at the time step n+1/2
            if (nn%Plot2) == 0:
                RealTime = round(nn*DT,5)
                #TITLE1 = Ex.Output_Dir + "/" +  "Time_half_" + str(RealTime)+"_s__Repeat_" + str(repeat)
                #TITLE2 = "at time_half: " + str(RealTime) +"_s__Repeat_" + str(repeat)
                TITLE1 = Ex.Output_Dir + "/" +  "Time_half_" + str(RealTime)+"_s"
                TITLE2 = "at time_half: " + str(RealTime) +"_s"
                Draw.Plot_at_Cell(N_Cells, X, Z, Q_1, V_1, Eta_1, U_1, E_1, A_1, TITLE1, TITLE2)

            #check = input(" Enter to continue ...") # <delete>
            # This loop evaluates all the required values at the faces at time step n+1/2
            for ii in range(N_Cells+1):
                if ii==0: # Boundary condition at face 1/2
                    Q_F_1[ii]   = Q_Upstream
                    Z_0 = Z[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                    Eta_F_1[ii] = Eta_1[ii] + L[ii] * (( Z[ii] - Z[ii+1]) / (L[ii] + L[ii+1]))
                    A_F_1[ii]   = (Eta_F_1[ii] - Z_0)* B[ii]
                    U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                    E_F_1[ii]   = ((U_F_1[ii])**2.0)/(2.0) + Gravity * Eta_F_1[ii]

                elif ii != 0 and ii != N_Cells: # middle cells
                    E_F_1[ii] = (L[ii]*  E_1[ii-1] + L[ii-1]*  E_1[ii] )/( L[ii] + L[ii-1] )
                    Eta_F_hat = (L[ii]*Eta_1[ii-1] + L[ii-1]*Eta_1[ii] )/( L[ii] + L[ii-1] )
                    #A_F_hat   = (L[ii]*A_1[ii-1] + L[ii-1]*A_1[ii])/(L[ii] + L[ii-1]) # <modify> for a variable area
                    Z_F0 = Z[ii] + L[ii] * (( Z[ii-1] - Z[ii]) / (L[ii-1] + L[ii]))
                    B_F_hat = B[ii] # <modify> Modify this equation for a variable area
                    A_F_hat   = (Eta_F_hat - Z_F0) * B_F_hat
                    Q_F_hat_S = (self.sgn_func(Q[ii-1])*L[ii]*((Q_1[ii-1])**2.0)  \
                                +self.sgn_func(Q[ii])*L[ii-1]*(Q_1[ii]**2.0))/(L[ii]+L[ii-1])
                    Gamma_F_hat = (L[ii]*Gamma[ii-1] + L[ii-1]*Gamma[ii])/(L[ii] + L[ii-1]) #<modify>for a variable area
                    H_hat    = A_F_hat / B_F_hat
                    a        = 1.0 + (2.0*Gamma_F_hat*(H_hat**2.0))/(A_F_hat) \
                               - ((Gravity*Eta_F_hat)/E_F_1[ii])*(1.0 + (2.0*Gamma_F_hat*(H_hat**2.0))/A_F_hat) \
                               - 2.0*Gravity*H_hat/E_F_1[ii]
                    b        = 2.0 - Gravity*(2.0*Eta_F_hat + H_hat)/E_F_1[ii]
                    c        = 1.0 - (Gravity*Eta_F_hat/E_F_1[ii]) - Q_F_hat_S/(2.0*E_F_1[ii]*(A_F_hat**2.0))

                    if ((b+1.0)**2.0 - 4.0 * a * c ) >= 0.0:
                        print(" Positive/Zero argument %f %f %f" % (a,b,c))

                        Eta_Epsilon1 = H_hat*(-b-1.0 + (((b+1.0)**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                        Eta_Epsilon2 = H_hat*(-b-1.0 - (((b+1.0)**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                    
                        if abs(Eta_Epsilon1) < abs(Eta_Epsilon2):
                            Eta_Epsilon    =  Eta_Epsilon1
                        else:
                            Eta_Epsilon    =  Eta_Epsilon2

                        Eta_F_1[ii] = Eta_F_hat + Eta_Epsilon
                        A_F_1[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon

                        x = Eta_Epsilon/H_hat
                        Alfa_Epsilon  = 2.0 * E_F_1[ii] * (A_F_hat**2.0) * ( a * (x**2.0) + b * x + c )

                        if Q_F_hat_S > 0.0:
                            print(" Q^2 is positive")
                            if -Q_F_hat_S < Alfa_Epsilon:
                                print(" Modification on Alpha")
                                Alfa_Epsilon = -Q_F_hat_S
                        elif Q_F_hat_S == 0.0:
                            #print(" Q is zero ", (nn, ii)) # <delete> after debugging
                            if Alfa_Epsilon < 0.0:
                                print(" Modification on Alpha- originally negative")
                                Alfa_Epsilon = 0.0
                        elif Q_F_hat_S < 0.0:
                            print(" Q^2 is negative")
                            #check = input(" Error: press ENTER to exit ")
                            #sys.exit()

                        Q_F_1[ii] = self.sgn_func(Q_F_hat_S)*(Q_F_hat_S + Alfa_Epsilon)**0.5
                        Q_check = self.sgn_func(Q_F_hat_S)*A_F_1[ii]* (2*(E_F_1[ii]-Gravity*Eta_F_1[ii]))**(0.5)
                        if abs(Q_check- Q_F_1[ii]) >0.01:
                            #print("check ", ii, Q_check, Q_F_1[ii])
                            print("{40}{:5d}{:5d}{:30.20f}{:30.20f}".format("Error: Q at the face is not consistent",nn,ii,Q_check, Q_F_1[ii]))
                            #check = input(" Error: inconsistency in cell values - Enter to continue ...")
                            #sys.exit()

                    elif ((b+1.0)**2.0 - 4.0 * a * c ) < 0.0:
                        print(" Warning: negative argument %f %f %f" % (a,b,c))
                        if (Gravity*Eta_F_hat > 0.5*Q_F_hat_S/(A_F_hat**2.0)) and (b**2.0 - 4.0 *a*c >= 0.0):

                            Eta_Epsilon1 = H_hat*(-b + ((b**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                            Eta_Epsilon2 = H_hat*(-b - ((b**2.0 - 4.0*a*c)**0.5))/(2.0*a) 
                        
                            if abs(Eta_Epsilon1) < abs(Eta_Epsilon2):
                                Eta_Epsilon    =  Eta_Epsilon1
                            else:
                                Eta_Epsilon    =  Eta_Epsilon2

                            Eta_F_1[ii] = Eta_F_hat + Eta_Epsilon
                            A_F_1[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon

                            #x = Eta_Epsilon/H_hat
                            #Alfa_Epsilon  = 2.0 * E_F[ii] * (A_F_hat**2.0) * ( a * (x**2.0) + b * x + c )
                            Alfa_Epsilon  = 0.0

                            Q_F_1[ii] = self.sgn_func(Q_F_hat_S)*(Q_F_hat_S + Alfa_Epsilon)**0.5
                            Q_check = self.sgn_func(Q_F_hat_S)*A_F_1[ii]* (2*(E_F_1[ii]-Gravity*Eta_F_1[ii]))**0.5
                            
                            if abs(Q_check- Q_F_1[ii]) >0.01:
                                print("{40}{:5d}{:5d}{:30.20f}{:30.20f}".\
                                    format("Error: Q at the face is not consistent",nn,ii,Q_check, Q_F_1[ii]))
                                #check = input(" Error: inconsistency in cell values - Enter to continue ...")
                                #check = input(" Error: press ENTER to exit ")
                                #sys.exit()
                        else:
                            Eta_Epsilon = 0.0
                            Eta_F_1[ii] = Eta_F_hat + Eta_Epsilon
                            A_F_1[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon
                            Alfa_Epsilon  = 2.0 * (A_F_1[ii]**2.0) * ( E_F_1[ii] - Gravity * Eta_F_1[ii] ) - Q_F_hat_S
                            if Alfa_Epsilon < - Q_F_hat_S:
                                Alfa_Epsilon =- Q_F_hat_S
                                Eta_Epsilon =  ( E_F_1[ii] - Gravity * Eta_F_hat ) / Gravity
                                Eta_F_1[ii] = Eta_F_hat + Eta_Epsilon
                                A_F_1[ii] = A_F_hat + (B_F_hat + Gamma_F_hat*Eta_Epsilon) * Eta_Epsilon
                            
                            Q_F_1[ii] = self.sgn_func(Q_F_hat_S)*(Q_F_hat_S + Alfa_Epsilon)**0.5
                            Q_check = self.sgn_func(Q_F_hat_S)*A_F_1[ii]* (2*(E_F_1[ii]-Gravity*Eta_F_1[ii]))**0.5
                            
                            if abs(Q_check- Q_F_1[ii]) >0.01:
                                print("{40}{:5d}{:5d}{:30.20f}{:30.20f}".\
                                    format("Error: Q at the face is not consistent",nn,ii,Q_check, Q_F_1[ii]))
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
                    #    print("Fatal error Downstream BC: %d, %d, %f, %f"%(ii,nn,Eta_1[ii-1],Eta_F[ii]))    
                    #    check = input(" Error: press ENTER to exit ")
                    #    sys.exit()
                    #Q_F_1[ii]= (1.0/M[ii-1])*A_F_1[ii]*((R_h_1[ii-1])**(2.0/3.0)) \
                    #           *(((Eta_1[ii-1]-Eta_F_1[ii])/(0.5*L[ii-1]))**(1.0/2.0))  # <modify>
                    Q_F_1[ii]   = Q_1[ii-1]
                    U_F_1[ii]   = Q_F_1[ii] / A_F_1[ii]
                    E_F_1[ii]   = (U_F_1[ii]**2.0)/2.0 + Gravity*Eta_F_1[ii]
                    #if ( E_F_1[ii] < E_F_1[ii-1] ):
                    #    print(" Fatal Error in Energy: %d, %d, %40.30f, %40.30f" % (nn, ii, E_F_1[ii], E_F_1[ii-1] ) )
                    #    check = input(" Error: press ENTER to exit ")
                    #    sys.exit()

            # Plot the full results at every "Plot2" steps at the time step n+1/2
            if (nn%Plot2) == 0:
                RealTime = round(nn*DT,5)
                #TITLE1 = Ex.Output_Dir + "/" +  "Full__time_" + str(RealTime)+"_s__Repeat_" + str(repeat)
                #TITLE2 = "Full results at time: " + str(RealTime) +"_s__Repeat_" + str(repeat)
                TITLE1 = Ex.Output_Dir + "/" +  "Full__time_half_" + str(RealTime)+"_s"
                TITLE2 = "Full results at half time: " + str(RealTime) +"_s"
                Draw.Plot_Full(3, N_Cells, X, X_F, Z_F, V_1, Q_1, Q_F_1, Eta_1, Eta_F_1, U_1, U_F_1, \
                                E_1, E_F_1, A_1, A_F_1, TITLE1, TITLE2)

            #check = input(" Enter to continue ...")

            for ii in range(N_Cells):
                F_q_1[ii*2  ] = Gravity * C_1[ii] * V_1[ii] * ( ( U_F_1[ii] + U_1[ii]     )**2.0) / 8.0
                F_q_1[ii*2+1] = Gravity * C_1[ii] * V_1[ii] * ( ( U_1[ii]   + U_F_1[ii+1] )**2.0) / 8.0

            # find k2 in the Runge-Kutta method and find the solution at n+1
            for ii in range(N_Cells):
                k_2V[ii]  = DT * ( Q_F_1[ii] - Q_F_1[ii+1] )  
                k_2Q[ii]  = (DT/L[ii])*(Q_F_1[ii]*U_F_1[ii] - Q_F_1[ii+1]*U_F_1[ii+1] \
                                        + Gravity*A_F_1[ii]*Eta_F_1[ii] - Gravity*A_F_1[ii+1]*Eta_F_1[ii+1] \
                                        - F_q_1[2*ii] - F_q_1[2*ii+1]) # <modify>remove
            V[ii] += k_2V[ii]
            Q[ii] += k_2Q[ii]

            # delete this section after debugging -- below
            #for ii in range(N_Cells):
            #    A[ii]   = V[ii] / L[ii]
            #    U[ii]   = Q[ii] / A[ii]
            #    Eta[ii] = A[ii] / B[ii] + Z[ii]
            #    E[ii]   = (U[ii]**2.0) /2.0 +  Gravity*Eta[ii]
            #
            #RealTime = round(nn*DT,5)
            ##TITLE1 = Ex.Output_Dir + "/" +  "Time_" + str(RealTime)+"_s__Repeat_" + str(repeat)
            ##TITLE2 = "at time: " + str(RealTime) +"_s__Repeat_" + str(repeat)
            #TITLE1 = Ex.Output_Dir + "/" +  "Before_Correction_Time_" + str(RealTime)+"_s"
            #TITLE2 = "Before_Correction_at time: " + str(RealTime) +"_s"
            #Draw.Plot_at_Cell(N_Cells, X, Z, Q, V, Eta, U, E, A, TITLE1, TITLE2)

            #check = input(" Press enter to continute... ")

            # delete this section after debugging -- above

            # Discussion on energy
            repeat = 0 
            Check_Energy = "False"
            while Check_Energy == "False":
                repeat += 1
                print("{:40} {:d} {:d}".format(" Repeat at this step: ",repeat,ii)) 
                #check = input(" Enter to continue ...")
                Check_Energy = "True"
                for ii in range(N_Cells):
                    A[ii]   = V[ii]/HL[ii]
                    U[ii]   = Q[ii]/A[ii]
                    Eta[ii] = A[ii]/B[ii] + Z[ii]
                    E[ii]   = (U[ii]**2.0) /2.0 + Gravity*Eta[ii]

                RealTime = round(nn*DT,5) # <delete> after debugging
                #TITLE1 = Ex.Output_Dir + "/" +  "Time_" + str(RealTime)+"_s__Repeat_" + str(repeat)
                #TITLE2 = "at time: " + str(RealTime) +"_s__Repeat_" + str(repeat)
                TITLE1 = Ex.Output_Dir + "/" +  "Before_Correction_Time_" + str(RealTime)+"_s" + str(repeat) 
                TITLE2 = "Before_Correction_at time: " + str(RealTime) +"_s" + str(repeat) 
                Draw.Plot_at_Cell(N_Cells, X, Z, Q, V, Eta, U, E, A, TITLE1, TITLE2)

                #check = input(" Press enter to continute - before correction ... ")

                for ii in range(1,N_Cells):
                    if E[ii] - E[ii-1]>0.000:
                        print("{:22} {:d} {:d} {:d} {:15.12f} {:15.12f} {:e}"\
                               .format(" Energy modification: ", nn, ii-1, ii, E[ii-1], E[ii], E[ii-1]-E[ii]))
                        Check_Energy = "False"
                        Delta_Energy = E[ii-1] - E[ii] 
                        Delta_Q = (Epsilon_E - Delta_Energy) / (DT*(-(U[ii-1]-U_F_1[ii])*(U[ii-1]/V[ii-1]) \
                            - (U[ii]-U_F_1[ii])*(U[ii]/V[ii]) + (Gravity/(B[ii-1]*L[ii-1])) + (Gravity/(B[ii]*L[ii]))))

                        print("{:26} {:d} {:d} {:15.12f} {:15.12f} {:15.12f}" \
                              .format(" Modifying the solution: ", nn, ii, Delta_Q, Q[ii-1], Q[ii]))

                        #check = input(" OKAY ")

                        V[ii] -= DT*Delta_Q 
                        V[ii-1] += DT*Delta_Q 

                        Q[ii] -=  DT/L[ii] * Delta_Q * U_F_1[ii]
                        Q[ii-1] += DT/L[ii] * Delta_Q * U_F_1[ii]

                        # Modifying Energy in each 
                        A[ii-1] = V[ii-1]/HL[ii-1]
                        U[ii-1] = Q[ii-1]/A[ii-1]
                        Eta[ii-1] = A[ii-1]/B[ii-1] + Z[ii-1]
                        E[ii-1] = (U[ii-1]**2.0) /2.0 +  Gravity*Eta[ii-1]

                        A[ii]   = V[ii]/HL[ii]
                        U[ii]   = Q[ii]/A[ii]
                        Eta[ii] = A[ii]/B[ii] + Z[ii]
                        E[ii]   = (U[ii]**2.0) /2.0 +  Gravity*Eta[ii]
                        #check = input(" Enter to continue ...")
                        break

        print(" ========== Solver Class ==========")
        print()




