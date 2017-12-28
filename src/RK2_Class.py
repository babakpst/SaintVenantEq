
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    12/08/2017
# Latest update: 12/10/2017
#
# Comment: This class contains all functions related to Runge-Kutta solver
#
#######################################################################################################################

class RK2:



    def RK2(self, rke, Face_Arrays_1, el, Face_Arrays, Face_hat, geo, fge, Q_Upstream, h_dw, N_Cells, Gravity, DT):  
        ''' RK2 time advance '''    

        import Values_Class as Values
        Values = Values.Values()
        
        print('BRH 20171123 error not tested')
        sys.exit()
        
        rke = self.RK2_step1( rke,  el, Face_Arrays, geo ) 
              
        print('in RK2 after first step')
        for ii in range(4,10):
            print(ii,'rke1V=',rke['k_1V'][ii],'rke1Q=',rke['k_1Q'][ii])
        print()
        
        rke = Values.element_values( rke, geo, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max )
        Face_Arrays_1 = self.face_values( rke, Face_Arrays_1, Face_hat, geo, fge, Q_Upstream, h_dw, depth_min, face_area_min, Gravity )
        
        print('in RK2 after face values')
        for ii in range(4,10):
            print(ii,'Qup',rke['Q'][ii-1],'Qf=',Face_Arrays_1['Q_F'][ii], 'Qdn=',rke['Q'][ii])
        print()

        for ii in range(4,10):
            print(ii,'Aup',rke['A'][ii-1],'Af=',Face_Arrays_1['A_F'][ii], 'Adn=',rke['A'][ii])
        print()
       
        # friction update
        rke['Fdn'] = Values.get_friction_half_element(rke['Fdn'], rke, Face_Arrays_1, +1, N_Cells, Gravity )
        rke['Fup'] = Values.get_friction_half_element(rke['Fup'], rke, Face_Arrays_1,  0, N_Cells, Gravity )    
        rke = self.RK2_step2( rke, Face_Arrays_1, el, geo )   
            
        for ii in range(0,N_Cells):
            el['V'][ii] = rke['V'][ii]
            el['Q'][ii] = rke['Q'][ii] 
        
        return el
    
   
    def RK2_step1(self, rke,  el, Face_Arrays, geo):
        ''' first step of an RK2 time adavance'''

        import Values_Class as Values
        Values = Values.Values()

        for ii in range(0,N_Cells): 
            rke['k_1V'][ii] = Values.continuity_RHS_explicit( Face_Arrays['Q_F'], ii) 
                        
            rke['k_1Q'][ii] = Values.momentum_RHS_explicit( el, Face_Arrays, ii, Gravity)
                
        # Solution at n+1/2
        rke['V'][:]   = el['V'][:] + 0.5 * DT * rke['k_1V'][:] 
        rke['Q'][:]   = el['Q'][:] + 0.5 * DT * rke['k_1Q'][:] / geo['L'][:]  

        for ii in range(0,N_Cells): 
            rke['V'][ii] = min( rke['V'][ii], volume_min[ii] )

        return rke

    def RK2_step2(self,rke, Face_Arrays_1, el, geo):
        ''' second step of an RK2 time advance'''

        import Values_Class as Values
        Values = Values.Values()

        for ii in range(0,N_Cells):
            rke['k_2V'][ii]  = Values.continuity_RHS_explicit( Face_Arrays_1['Q_F'], ii)
            
            rke['k_2Q'][ii]  = Values.momentum_RHS_explicit( rke, Face_Arrays_1, ii, Gravity)

        # Solution at n+1
        rke['V'][:]   = el['V'][:] + DT * rke['k_2V'][:]
        rke['Q'][:]   = el['Q'][:] + DT * rke['k_2Q'][:] / geo['L'][:] 

        for ii in range(0,N_Cells): 
            rke['V'][ii] = min( rke['V'][ii], volume_min[ii] )
                       
        return rke




