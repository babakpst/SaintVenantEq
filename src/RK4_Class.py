
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    12/08/2017
# Latest update: 12/27/2017
#
# Comment: This class contains all functions related to Runge-Kutta solver
#
#######################################################################################################################

class RK4:


    def RK4(self, el, Face_Arrays, geo, fge, rke, Face_Arrays_1, Q_Upstream, h_dw, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max, depth_min, face_area_min ):
                 
        import Values_Class as Values
        Values = Values.Values()

        for ii in range(0,N_Cells):
            rke['k_1V'][ii] = Values.continuity_RHS_explicit( Face_Arrays['Q_F'], ii )
            rke['k_1Q'][ii] = Values.momentum_RHS_explicit( el, Face_Arrays, geo, ii, Gravity ) 
            
        rke['V'][:] = el['V'][:] + 0.5 * DT * rke['k_1V'][:]
        rke['Q'][:] = el['Q'][:] + 0.5 * DT * rke['k_1Q'][:]
                
        rke = Values.element_values( rke, geo, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max )
        Face_Arrays_1 = Values.face_values( rke, Face_Arrays_1, geo, fge, Q_Upstream, h_dw, depth_min, face_area_min, Gravity, N_Cells )
        
        for ii in range(0,N_Cells):
            rke['k_2V'][ii] = Values.continuity_RHS_explicit( Face_Arrays_1['Q_F'], ii )
            rke['k_2Q'][ii] = Values.momentum_RHS_explicit( rke, Face_Arrays_1, geo,  ii, Gravity ) 
            
        rke['V'][:] = el['V'][:] + 0.5 * DT * rke['k_2V'][:]
        rke['Q'][:] = el['Q'][:] + 0.5 * DT * rke['k_2Q'][:]
       
        rke = Values.element_values( rke, geo, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max )
        Face_Arrays_1 = Values.face_values( rke, Face_Arrays_1, geo, fge, Q_Upstream, h_dw, depth_min, face_area_min, Gravity, N_Cells )

        for ii in range(0,N_Cells):
            rke['k_3V'][ii] = Values.continuity_RHS_explicit( Face_Arrays_1['Q_F'], ii )
            rke['k_3Q'][ii] = Values.momentum_RHS_explicit( rke, Face_Arrays_1, geo, ii, Gravity ) 
            
        rke['V'][:] = el['V'][:] + DT * rke['k_3V'][:]
        rke['Q'][:] = el['Q'][:] + DT * rke['k_3Q'][:]
            
        rke = Values.element_values( rke, geo, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max )
        Face_Arrays_1 = Values.face_values(rke, Face_Arrays_1, geo, fge, Q_Upstream, h_dw, depth_min, face_area_min, Gravity, N_Cells )

        for ii in range(0,N_Cells):
            rke['k_4V'][ii] = Values.continuity_RHS_explicit(Face_Arrays_1['Q_F'], ii)
            rke['k_4Q'][ii] = Values.momentum_RHS_explicit(rke, Face_Arrays_1, geo, ii, Gravity) 
            
        # Update Volume and Flow rate at n+1
        el['V'][:] = el['V'][:]                                                \
          + (DT/6.0) * (    rke['k_1V'][:] + 2*rke['k_2V'][:]                  \
                      + 2.0*rke['k_3V'][:] +   rke['k_4V'][:])

        el['Q'][:] = el['Q'][:]                                                \
          + (DT/6.0) * (  rke['k_1Q'][:] + 2*rke['k_2Q'][:]                    \
                      + 2.0*rke['k_3Q'][:]   + rke['k_4Q'][:])

        return el


