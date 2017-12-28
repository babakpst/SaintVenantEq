
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    12/08/2017
# Latest update: 12/27/2017
#
# Comment: This class contains all the interpolation functions
#
#######################################################################################################################

class Interpolation:


    def sgn_func(self, x): 
        ''' returns sign of argument '''
        return 1.0 if x >= 0 else -1.0
    

    def sgn_func_array(self, xx):
        ''' returns array with sign of a numpy array '''
        xsgn = np.ones(xx.size,dtype=int)
        aa   = np.nonzero(xx < 0)
        xsgn[aa] = -1        
        return xsgn
    

    def linear_interp_to_face(self, cc, LL, ii):
        ''' generic linear interpolation from centers to face'''
        return ( LL[ii] * cc[ii-1] + LL[ii-1] * cc[ii] )  \
                / ( LL[ii] + LL[ii-1] ) 


    def linear_interp_to_face_array(self, cc, LL):
        ''' generic linear interpolation from centers to face'''
        return (LL[1:N_Cells]*cc[0:N_Cells-1] + LL[0:N_Cells-1]*cc[1:N_Cells]) \
              / (LL[1:N_Cells] + LL[0:N_Cells-1])
                

    def linear_interp_to_face_onevalue(self, cc1, cc2,  L1, L2 ):
        ''' generic linear interpolation from centers to face'''
        return ( L1 * cc2 + L2 * cc1 ) / ( L1 + L2 ) 
  

    def linear_extrapolate_to_upstream_face(self, cc, L):
        ''' generic linear extrapolation to upstream boundary'''       
        return cc[0] + L[0] *  ( cc[0] - cc[1]) / ( L[0] + L[1]) 
        

    def time_weighted_interp(self, phi_m, phi_p, TT_i, TT_ip, ii): # For a single variable
        ''' generic time-weighted interpolation'''
        return (TT_ip*phi_m + TT_i*phi_p) / (TT_i + TT_ip) 


    def energy_value(self, U, Eta, Gravity):
        '''Compute the energy from velocity and free surface elevation'''
        return (0.5 * U**2.0) + Gravity*Eta


    def exponential_of_abs_value_with_sign(self, value, power):

        for ii in range(0,value.size-1):
            value[ii] = self.sgn_func(value[ii]) * (abs(value[ii]))**power
        
        return value
