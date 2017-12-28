
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/31/2017
# Latest update: 12/27/2017
#
# Comment: This class conducts the numerical solution of the Saint-Venant equation
#
#######################################################################################################################

class Solver:

    # definitions:
    # the ii=0 cell center is first cell inside the upstream boundary
    # the ii=N_Cells-1 cell center is first cell inside the downstream boundary 
    # the cells outside the downstream boundary (ii=N_Cells) do not exist
    # faces take on the index of the downstream cell.
    # the ii=0 face is the upstream boundary
    # the ii=N_Cells face is the downstream boundary
    # Because of pythong indexing
    #   range(0,N_Cells) gives you all cells (which are interior)
    #   range(1,N_Cells) gives you all interior faces
    #   range(0,N_Cells+1) gives you all faces
    #   Var_F[N_Cells] = furthest downstream face value
    #   Var_E(N_Cells-1) = furthest downstream cell
    
    def __init__(self): 
        pass
   
    def solve(self,argv):

        import numpy as np
        import sys
        import os

        import Initialization_Class as Initial
        import Visualization_Class  as Visual
        import Info_Class as Information
        import RK4_Class as Solver_RK4
        import RK2_Class as Solver_RK2
        import Interpolation_Class
        import RootFinder_Class
        import Forward_Euler_Class as Solver_ForwardEuler
        import Values_Class as Values
        
        global el, Face_Arrays, Face_hat, geo, fge
        global Gravity
        global N_Cells, nn
        global DT, DT_max, DT_min, DT_variable
        global CFLu_max, CFLb_max
        global Fr_max
        
        global face_correct
        global face_interp
        
        global np, sys
        
        global ifig
        
        global quit_after_plot
        
        global depth_min
        global volume_min
        global face_area_min
        
        global error_check
                
        quit_after_plot = False
        
        error_check = True
        
        ifig = 1

        # Define classes
        Draw = Visual.Visualization()
        Ex = Initial.Initialization(argv)
        RK4 = Solver_RK4.RK4()
        RK2 = Solver_RK2.RK2()
        FEuler = Solver_ForwardEuler.ForwardEuler()
        Values = Values.Values()


        # Setting
        # for automatically selecting time step use DT_variable = 1
        DT_variable = 0 
        #Gravity = Setting.Gravity
        
        DT_min = Ex.command_line_args['DT_min']; #print("DT_min",DT_min)
        DT_max = Ex.command_line_args['DT_max']; #print("DT_max",DT_max)

        Plot_at_Cell = Ex.command_line_args['Plot_at_Cell']; #print("Plot_at_Cell",Plot_at_Cell)
        Plot_at_Face = Ex.command_line_args['Plot_at_Face']; #print("Plot_at_Face",Plot_at_Face)

        CFLu_max = Ex.command_line_args['CFLu_max']; #print("CFLu_max",CFLu_max)
        CFLb_max = Ex.command_line_args['CFLb_max']; #print("CFLb_max",CFLb_max)
        
        Fr_max = Ex.command_line_args['Fr_max']; #print("Fr_max",Fr_max)
        
        depth_min = Ex.command_line_args['depth_min']; #print("depth_min",depth_min)
        
        time_advance = Ex.command_line_args['time_advance']; #print("time_advance",time_advance)
                
        Gravity = Ex.command_line_args['Gravity']; #print("Gravity",Gravity)

        face_correct = Ex.command_line_args['face_correct']; #print("face_correct",face_correct)
        face_interp = Ex.command_line_args['face_interp']; #print("face_interp",face_interp)
 
        plotstart = 0 # time step to start the plots
        printout = 1 # iteration for printouts
        


        # Open an info file
        Info_File_Dir = os.path.join(Ex.Output_Dir,("Info_"+Ex.Model) ) 
        Info_File = open(Info_File_Dir,"w")

        Info = Information.Info(Info_File,Ex.Model)
        Info_File.close()

        if face_correct == 1 and face_interp != 2:
            print('Inconsistent setup. Face correction requires 2nd order interpolation')
            sys.exit()
                
        print(" ========== Solver Class ==========")
        print(" Solving SVE ...")
        if time_advance == "rk2": print(" Second-order Runge-Kutta method: ...")
        elif time_advance == "rk4": print(" Fourth-order Runge-Kutta method: ...")
        elif time_advance == "forward_euler": print(" First-order forward-Euler method: ...")
        else:
            print(" Wrong choice of the solver. See the setting file.")
            sys.exist()

        print(" Allocating memory ...")
        N_Steps = int(Ex.Total_Time/Ex.Time_Step)
        N_Cells = Ex.N_Cells
        h_dw    = Ex.h_dw
        DT      = Ex.Time_Step

        # fixed geometry of elements for solution
        L    = np.zeros(N_Cells,     dtype=np.float64 ) # element length
        HL   = np.zeros(N_Cells,     dtype=np.float64 ) # Horzontal element length (projection length on the X-axis)
        Z    = np.zeros(N_Cells,     dtype=np.float64 ) # element bottom elevation
        M    = np.zeros(N_Cells,     dtype=np.float64 ) # element Manning's n
        B    = np.zeros(N_Cells,     dtype=np.float64 ) # element bottom breadth (square channel)
        X    = np.zeros(N_Cells,     dtype=np.float64 ) # element X distance
        
        geo = {'L':L, 'HL':HL,'Z':Z, 'M':M, 'B':B, 'X':X}
               
       # geometry for plotting (double density)
        Z_Fp  = np.zeros(N_Cells*2+1, dtype=np.float64) # face bottom elevation
        X_Fp  = np.zeros(N_Cells*2+1, dtype=np.float64) # element X distance on face

        # dynamic variables on elements      
        A    = np.zeros(N_Cells, dtype=np.float64) # Cross-sectional area
        C    = np.zeros(N_Cells, dtype=np.float64) # Temporary variable
        E    = np.zeros(N_Cells, dtype=np.float64) # Energy
        Eta  = np.zeros(N_Cells, dtype=np.float64) # Free surface elevation
        Fup  = np.zeros(N_Cells, dtype=np.float64) # Friction on upstream 1/2 of cell
        Fdn  = np.zeros(N_Cells, dtype=np.float64) # Friction on downstream 1/2 of cell  
        Fr   = np.zeros(N_Cells, dtype=np.float64) # Froude number
        Gamma= np.zeros(N_Cells, dtype=np.float64) # Gradient of Top Width with Z
        l_P  = np.zeros(N_Cells, dtype=np.float64) # wetted perimeter
        R_h  = np.zeros(N_Cells, dtype=np.float64) # hydraulic radius
        Q    = np.zeros(N_Cells, dtype=np.float64) # flow rate
        U    = np.zeros(N_Cells, dtype=np.float64) # velocity
        V    = np.zeros(N_Cells, dtype=np.float64) # volume

        T    = np.zeros(N_Cells, dtype=np.float64) # Top Width
        H    = np.zeros(N_Cells, dtype=np.float64) # hydraulic depth
        CFLb = np.zeros(N_Cells, dtype=np.float64) # barotropic CFL
        CFLu = np.zeros(N_Cells, dtype=np.float64) # advective CFL

        # dynamic variables for RK solution (Cell variables at the mid-step)
        A_1   = np.zeros(N_Cells, dtype=np.float64) # Cross-sectional area
        C_1   = np.zeros(N_Cells, dtype=np.float64) # Temporary variable
        E_1   = np.zeros(N_Cells, dtype=np.float64) # Energy
        Eta_1 = np.zeros(N_Cells, dtype=np.float64) # Free-surface elevation
        Fup_1 = np.zeros(N_Cells, dtype=np.float64) # Friction on upstream 1/2 of cell
        Fdn_1 = np.zeros(N_Cells, dtype=np.float64) # Friction on downstream 1/2 of cell  
        Fr_1  = np.zeros(N_Cells, dtype=np.float64) # Froude Number
        Gamma_1 = np.zeros(N_Cells, dtype=np.float64) # Gradient of Top Width with Z
        l_P_1 = np.zeros(N_Cells, dtype=np.float64) # wetted perimeter
        R_h_1 = np.zeros(N_Cells, dtype=np.float64) # hydraulic radius
        Q_1   = np.zeros(N_Cells, dtype=np.float64) # flow rate
        U_1   = np.zeros(N_Cells, dtype=np.float64) # Velocity
        V_1   = np.zeros(N_Cells, dtype=np.float64) # Volume

        #T_ed = np.zeros(N_Cells, dtype=np.float64 ) # <modify> when finalized
        #T_eu = np.zeros(N_Cells, dtype=np.float64 ) # <modify> when finalized

        # Define the "el" dictionary that contains all parameters related to each cell.
        el =  {'A':A, 'C':C, 'CFLb':CFLb, 'CFLu':CFLu, 'E':E, 'Eta':Eta, 
               'Fr':Fr, 'Fdn':Fdn, 'Fup':Fup, 'Gamma':Gamma, 
               'H':H, 'l_P':l_P, 'R_h':R_h, 'Q':Q, 'T':T, 'U':U, 
               'V':V} #, 'T_ed':T_ed, 'T_eu':T_eu } # Element Dictionary that contains all information for the cell center
                             
        # geometry on face
        Z_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # cross-sectional area
        B_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # breadth (rectangular channel)
 
        fge = {'Z_F':Z_F, 'B_F':B_F} # Face Geometries (Dictionary)
        
        # dynamic variables on face
        A_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # cross-sectional area
        E_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # energy
        Eta_F = np.zeros(N_Cells+1, dtype=np.float64 ) # free surface elevation
        Q_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # flow rate
        U_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # velocity
        
        # Define the "Face_Arrays" dictionary that contains all parameters related faces.
        Face_Arrays = {'A_F':A_F, 'E_F':E_F, 'Eta_F':Eta_F, 'Q_F':Q_F, 'U_F':U_F}
        
        # dynamic variables on face for RK solution  <delete> after debugging, Substitute with the above arrays.
        # (This is identical to the previous set of variables. We need this as a temporary set of variables for RK time marching.)
        A_F_1   = np.zeros(N_Cells+1, dtype=np.float64 ) # cross-sectional area
        E_F_1   = np.zeros(N_Cells+1, dtype=np.float64 ) # energy
        Eta_F_1 = np.zeros(N_Cells+1, dtype=np.float64 ) # free surface elevation
        Q_F_1   = np.zeros(N_Cells+1, dtype=np.float64 ) # flow rate
        U_F_1   = np.zeros(N_Cells+1, dtype=np.float64 ) # velocity

        # Define the "Face_Arrays" dictionary that contains all parameters related faces (Required by RK solution).
        Face_Arrays_1 = {'A_F':A_F_1, 'E_F':E_F_1, 'Eta_F':Eta_F_1, 'Q_F':Q_F_1, 'U_F':U_F_1} # Face arrays for RK

        # estimated variables on face
        A_F_hat     = np.zeros(N_Cells+1, dtype=np.float64 ) # Estimated face area at the face
        Eta_F_hat   = np.zeros(N_Cells+1, dtype=np.float64 ) # Estimated elevation at the face
        Q_F_hat_S   = np.zeros(N_Cells+1, dtype=np.float64 ) # Estimated flow rate at the face
        Gamma_F_hat = np.zeros(N_Cells+1, dtype=np.float64 ) # Estimated channel shape derivative at the face
        H_hat       = np.zeros(N_Cells+1, dtype=np.float64 ) # Estimated depth at the face
        
        # We define the "Face_hat" dictionary, that contains all the estimated variables
        Face_hat ={'Eta_F_hat':Eta_F_hat, 'A_F_hat':A_F_hat, 'Q_F_hat_S':Q_F_hat_S, 
                    'Gamma_F_hat':Gamma_F_hat, 'H_hat':H_hat} 
        
        # RK solution space <delete> after debugging
        k_1V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_1Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_2V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_2Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_3V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_3Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_4V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_4Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove        
       
        rke = {'k_1V':k_1V, 'k_1Q':k_1Q, 'k_2V':k_2V, 'k_2Q':k_2Q, 
               'k_3V':k_3V, 'k_3Q':k_3Q, 'k_4V':k_4V, 'k_4Q':k_4Q,

               'A':A_1, 'C':C_1, 'CFLb':CFLb, 'CFLu':CFLu, 'E':E_1, 'Eta':Eta_1, 
               'Fr':Fr_1, 'Fdn':Fdn_1, 'Fup':Fup_1, 'Gamma':Gamma_1,
               'H':H, 'l_P':l_P_1, 'R_h':R_h_1, 'Q':Q_1, 'T':T, 'U':U_1, 
               'V':V_1}

        volume_min = np.zeros(N_Cells, dtype=np.float64)
        face_area_min = np.zeros(N_Cells+1, dtype=np.float64)
        
        #=======================================================================
        # Initialization 
        print(" Initialization ... ") # For the cell center only
        el['Q'][:]   = Ex.Q[:]
        el['V'][:]   = Ex.V[:]
        geo['L'][:]  = Ex.L[:]
        geo['HL'][:] = Ex.HL[:]
        geo['Z'][:]  = Ex.Z[:]
        geo['M'][:]  = Ex.M[:]
        geo['B'][:]  = Ex.B[:]
        geo['X'][:]  = Ex.X[:]
        X_Fp[:] = Ex.X_F[:]
        Z_Fp[:] = Ex.Z_F[:]
        
        el['Q'][:] = Ex.Q[:]
        
        # for error checking rectangular channel   
        volume_min[:] = depth_min * geo['L'][:] * geo['B'][:]  # double <check> 

        # rectangular channel: Setting the minimums for the face
        face_area_min[1:N_Cells] = depth_min*0.5*(geo['B'][1:N_Cells] + geo['B'][0:N_Cells-1]) # <check>
        face_area_min[0] = face_area_min[1] # <check>
        face_area_min[N_Cells] = face_area_min[N_Cells-1] # <check>


        # Geometry setup: (Evaluate the Z_F and B_F)
        fge = Values.face_upstream_bottom_elevation(fge, geo, N_Cells)
        fge = Values.face_downstream_bottom_elevation(fge, geo, N_Cells)
        fge = Values.face_interior_bottom_elevation(fge, geo, N_Cells)
        fge = Values.face_breadth_rectangular_channel(fge, geo, N_Cells)

        print(" Time marching ... ")
        
        for nn in range(N_Steps):
            #if (nn%printout) == 0:
            print("========================== Time step: %d out of %d " % (nn, N_Steps))
                
            # store the inflow boundary for this step     
            Q_Upstream = Ex.Q_Up     # Hack

            # Cell center auxiliary values at time n
            el = Values.element_values(el, geo, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max )            
                      
            # face values at time n
            Face_Arrays = Values.face_values(el, Face_Arrays, geo, fge, Q_Upstream, h_dw, depth_min, face_area_min, Gravity, N_Cells )  

            # friction on center requires face values
            el['Fdn'] = Values.get_friction_half_element(el['Fdn'], el, Face_Arrays, +1, N_Cells, Gravity)
            el['Fup'] = Values.get_friction_half_element(el['Fup'], el, Face_Arrays, 0, N_Cells, Gravity )
              
            if ((nn%Plot_at_Face) == 0 and nn >= plotstart) \
                or (quit_after_plot):
                RealTime = round(nn*DT,5)
                TITLE1 = Ex.Output_Dir + "/" +  "Full__time_" + str(RealTime)+"_s"
                TITLE2 = "Full results at time: " + str(RealTime) +"_s"
                Draw.Plot_Full_Results(ifig, N_Cells, geo['X'], X_Fp, Z_Fp, el['V'], el['Q'], Face_Arrays['Q_F'],  \
                               el['Eta'], Face_Arrays['Eta_F'], el['U'], Face_Arrays['U_F'],    \
                               el['E'], Face_Arrays['E_F'], el['A'], Face_Arrays['A_F'],        \
                               TITLE1, TITLE2)
                if quit_after_plot:
                    print('exit after plot due to error')
                    sys.exit()
                    
                ifig=ifig+1
                                
            if time_advance == 'rk2':    
                el = RK2.RK2(rke, Face_Arrays_1, el, Face_Arrays, Face_hat, geo, fge, Q_Upstream, h_dw, N_Cells, Gravity, DT)  
            elif time_advance == 'rk4':    
                el = RK4.RK4(el, Face_Arrays, geo, fge, rke, Face_Arrays_1, Q_Upstream, h_dw, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max, depth_min, face_area_min )
            elif time_advance == 'forward_euler':
                el = FEuler.forward_euler( el, Face_Arrays, geo )
            else:
                print('error, unknown time advance of ',time_advance)
                sys.exit()
         

        del Info_File
        print(" ========== Solver Class ==========")
        print()
