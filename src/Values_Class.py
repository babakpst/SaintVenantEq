
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    12/27/2017
# Latest update: 12/27/2017
#
# Comment: This class contains all the value computations
#
#######################################################################################################################

class Values:



    def __init__(self): 
        pass


    def element_values(self, el, geo, N_Cells, Gravity, DT, error_check, Fr_max, CFLb_max, CFLu_max ):
        '''Compute auxiliary values on element center
           Assumes Q, V and geometry are known
           Friction (Fdn, Fup) are computed later as they require face values
        '''
        import numpy as np

        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()

        #for ii in range(0,N_Cells):  # <delete> after debugging
        #    print(ii, el['V'][ii])            
        #sys.exit()    
        
        global quit_after_plot

        C_ed = np.zeros(N_Cells, dtype=np.float64 ) # <remove> after debugging
        C_eu = np.zeros(N_Cells, dtype=np.float64 ) # <remove> after debugging

        # depends on input
        # note that we assume V > volume_min, which must be checked when
        # V is updated
        el['A'][:]   = el['V'][:] / geo['HL'][:]
       
        # depends onrectangular cell assumption and input
        el['Eta'][:] = geo['Z'][:] + el['A'][:] / geo['B'][:]
        el['T'][:]   = geo['B'][:]
        
        # depends on A and method to get T
        el['H'][:]   = el['A'][:] / el['T'][:]
     
        # depends onrectangular cell assumption and input
        el['l_P'][:] = geo['B'][:] + 2.0*el['H'][:]
        el['Gamma'][:] = 0.0  # <modify later>

        # depends on A, cell cross-section, and input
        el['U'][:]   = self.velocity_value(el['Q'][:] , el['A'][:])


        # Finding the time scale
        C_ed[:] = el['U'][:] + (Gravity * el['H'][:])**0.5
        C_eu[:] = el['U'][:] - (Gravity * el['H'][:])**0.5

        #el['T_ed'][:]   = +geo['L'][:] / (2.0*C_ed[:])    # <modify> when finalized
        #el['T_eu'][:]   = -geo['L'][:] / (2.0*C_eu[:])    # <modify> when finalized
        # <modify>   ################## What is T_i

        el['R_h'][:] = el['A'][:] / el['l_P'][:]

        # depends on Eta and U
        el['E'][:]   = Inter.energy_value( el['U'][:], el['Eta'][:], Gravity)

        # depends on R_h
        el['C'][:]   = ((geo['M'][:])**2.0) / (el['R_h'][:]**(4.0/3.0))

        # non-dimensional numbers
        el['Fr'][:]  = abs(el['U'][:]) / ( (Gravity * el['H'][:] )**0.5)               
        el['CFLu'][:] = abs(el['U'][:]) * DT / geo['L'][:]
        el['CFLb'][:] = (Gravity * el['H'][:]) * DT / geo['L'][:]

        if error_check:
            for ii in range(0,N_Cells): 
                if el['Fr'][ii] >= Fr_max:
                    print(ii, el['Q'][ii], el['A'][ii], el['U'][ii], el['H'][ii])
                    print("Flow is above max Fr %d %d %f" % (nn, ii, el['Fr'][ii]))
                    quit_after_plot = True  
                
                if el['CFLb'][ii] > CFLb_max:
                    print("CFLb_max exceeded %d %d %f" % (nn, ii, el['CFLb'][ii]))
                    #if DT_variable == 0:
                    quit_after_plot = True
    
                if el['CFLu'][ii] > CFLu_max:
                    print('Cell=',ii,'Q=',el['Q'][ii],'A=',el['A'][ii],'U=',el['U'][ii],'H=',el['H'][ii])                
                    print("CFLu_max exceeded %d %d %f" % (nn, ii, el['CFLu'][ii]))
                    #if DT_variable == 0:
                    quit_after_plot = True
           
        return el

    def face_values(self, el, Face_Arrays, geo, fge, Q_Upstream, h_dw, depth_min, face_area_min, Gravity, N_Cells):
        '''interpolates values to faces'''
               
        Face_Arrays = self.face_upstream_boundary_values(el, Face_Arrays, geo, fge, Q_Upstream, depth_min, face_area_min, Gravity )
        Face_Arrays = self.face_downstream_boundary_values(el, Face_Arrays, fge, h_dw, N_Cells, depth_min, face_area_min, Gravity )
                        
#        if face_correct == 0:
        Face_Arrays = self.face_by_interpolation( el, Face_Arrays, geo, fge, depth_min, N_Cells, Gravity, face_area_min) # face values on all faces except boundaries
#        elif face_correct ==1:
#            Face_Arrays   = self.face_interior_energy_values( el, Face_Arrays, geo )
#            Face_hat = self.face_interior_hat_values( el, fge, Face_hat)
#            Face_Arrays   = self.face_corrections(Face_Arrays, Face_hat, fge)
#        else:
#            print('Unknown value of face_correct = ',face_correct)   
#            sys.exit()
                
#        for ii in range(0,N_Cells+1):
#            print(ii, Face_hat['A_F_hat'][ii])
#        print('at 108')    
#        sys.exit()
                       
        return Face_Arrays  
    
    def face_by_interpolation(self, el, Face_Arrays, geo, fge, depth_min, N_Cells, Gravity, face_area_min):
        '''Interpolation for all face values'''

        import Interpolation_Class as Interpolation
        import RootFinder_Class

        Inter = Interpolation.Interpolation()
        RootFinder = RootFinder_Class.Root_Finder()

        for ii in range(1,N_Cells):  # All the internal faces
            # direction of fluxes on adjacent cells
            sup = Inter.sgn_func( el['Q'][ii-1])  # Upstream cell
            sdn = Inter.sgn_func( el['Q'][ii])    # Downsream cell
            
            if sup * sdn == 1: # not a convergent or divergent flow: flow in both cells in the same direction
                #print("flow in both cells in the same direction")

                if sup == +1: # flow in nominal downstream direction
                    upi = ii-1
                    dni = ii
                elif sup == -1: # flow in nominal upstream direction
                    upi = ii
                    dni = ii-1
                else:
                    print('error in sign of argument')
                    sys.exit()
                
                #if el['Fr'][ii-upi] >= 1.0: # double check: what does it mean?
                if el['Fr'][upi] >= 1.0:  # supercritical: 
                    # upstream energy change from friction
                    SpecificEnergy = el['E'][upi] - geo['Z'][upi] * Gravity
                    Face_Arrays['Q_F'][ii]  = el['Q'][upi]
                elif el['Fr'][upi] < 1.0: 
                    # subcritical
                    SpecificEnergyUp = el['E'][upi]  - geo['Z'][upi] * Gravity
                    SpecificEnergyDn = el['E'][dni]  - geo['Z'][dni] * Gravity
                    SpecificEnergy = Inter.linear_interp_to_face_onevalue(SpecificEnergyUp, SpecificEnergyDn, \
                                                                         geo['L'][upi], geo['L'][dni])

                    Face_Arrays['Q_F'][ii] = Inter.linear_interp_to_face(el['Q'], geo['L'], ii)
                   
                # transverse width at face (need storage?) HACK   
                faT = Inter.linear_interp_to_face(el['T'], geo['L'], ii) 
                
                # compute roots of cubic to find area (a tuple with the coefficient of the third-degree equation)
                pp = (1, -SpecificEnergy*faT/Gravity, 0,+0.5*faT*(Face_Arrays['Q_F'][ii]**2.0)/Gravity)

                # cubic with reduced form using a=1 and c=0
                rootout = RootFinder.cubic_without_a_c(pp[1],pp[3])
                rootout_modified = []

                for Inumber in rootout:
                    if Inumber.imag == 0.0 and Inumber.real>0.0:
                        rootout_modified.append(Inumber.real)

                if el['Fr'][upi] >= 1.0:
                    Face_Arrays['A_F'][ii] = min(rootout_modified)
                elif el['Fr'][upi] < 1.0:
                    Face_Arrays['A_F'][ii] = max(rootout_modified)

                Face_Arrays['A_F'][ii] = max(face_area_min[ii], Face_Arrays['A_F'][ii])

                Face_Arrays['U_F'][ii] = self.velocity_value(Face_Arrays['Q_F'][ii], Face_Arrays['A_F'][ii])

                Face_Arrays['Eta_F'][ii] = fge['Z_F'][ii]+(SpecificEnergy - 0.5*Face_Arrays['U_F'][ii]**2.0)/Gravity
                    
                Face_Arrays['E_F'][ii] = Inter.energy_value(Face_Arrays['U_F'][ii], Face_Arrays['Eta_F'][ii], Gravity)
                    

                if Face_Arrays['E_F'][ii] < fge['Z_F'][ii]:
                    print('problem in face interpolation, insufficient energy')
                    sys.exit()
                if el['Eta'][ii-upi] <= fge['Z_F'][ii]:
                    print('Possible hydraulic jump required - not yet coded')
                    sys.exit()
       
            elif sup * sdn == -1: # flow is converging or diverging on face
                print("flow is converging or diverging on face")

                Face_Arrays['Q_F'][ii] = el['Q'][ii-1] + el['Q'][ii]
                
                Face_Arrays['A_F'][ii] = Inter.linear_interp_to_face(el['A'], geo['L'], ii)
                Face_Arrays['A_F'][ii] = max( Face_Arrays['A_F'][ii], face_area_min(ii))  
                Face_Arrays['Eta_F'][ii] = Inter.linear_interp_to_face(el['Eta'], geo['L'], ii)
                Face_Arrays['Eta_F'][ii] = max(Face_Arrays['Eta_F'][ii], depth_min + fge['Z_F'][ii])  
                
                Face_Arrays['U_F'][ii] = self.velocity_value(Face_Arrays['Q_F'][ii], Face_Arrays['A_F'][ii])
                Face_Arrays['E_F'][ii] = Inter.energy_value(Face_Arrays['U_F'][ii], Face_Arrays['Eta_F'][ii], Gravity)
                
        return Face_Arrays

    def face_upstream_bottom_elevation(self, fge, geo, N_Cells):
        '''upstream face Z by extrapolation'''
        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()
        fge['Z_F'][0] = Inter.linear_extrapolate_to_upstream_face(geo['Z'], geo['L'])      
                        
        return fge
 
    def face_downstream_bottom_elevation(self, fge, geo, N_Cells):
        ''' temporary hard code of bottom elevation as 0''' 
        fge['Z_F'][N_Cells] = 0.0 #HACK
       
        return fge

    def face_interior_bottom_elevation(self, fge, geo, N_Cells):
        '''linear interpolation of bottom elevation between cell centers'''
        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()
        # note this only needs to be called once in a solution
        for ii in range(1,N_Cells):
            fge['Z_F'][ii] = Inter.linear_interp_to_face(geo['Z'], geo['L'], ii)
                    
        return fge

    def face_breadth_rectangular_channel(self, fge, geo, N_Cells):
        ''' simple linear interpolation and boundary extrapolation'''
        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()
        # handles both interior and boundaries
        for ii in range(1,N_Cells):
            fge['B_F'][ii] = Inter.linear_interp_to_face(geo['B'], geo['L'], ii)

        fge['B_F'][0] =  geo['B'][0]
        fge['B_F'][N_Cells] = geo['B'][N_Cells-1]
        
        return fge



    def face_upstream_boundary_values(self, el, Face_Arrays, geo, fge, Q_Upstream, depth_min, face_area_min, Gravity ):
        '''Compute face values at boundary given Q upstream element values 
            at cell ii=0
        '''
        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()

        Face_Arrays['Q_F'][0] = Q_Upstream  # use input boundary condition
        
        # extrapolate bottom gradient to get surface elevation on upstream face
        # note this CANNOT be replaced with linear_extrapolate_to_upstream_face
        Face_Arrays['Eta_F'][0] = el['Eta'][0] + geo['L'][0]   \
            *  ( geo['Z'][0] - geo['Z'][1])           \
             / ( geo['L'][1] + geo['L'][2]) 
        Face_Arrays['Eta_F'][0] = max( Face_Arrays['Eta_F'][0], depth_min + fge['Z_F'][0]) # <check>

        # area for rectangular channel only
        Face_Arrays['A_F'][0] = (Face_Arrays['Eta_F'][0] - fge['Z_F'][0]) * fge['B_F'][0]
        Face_Arrays['A_F'][0] = max(Face_Arrays['A_F'][0], face_area_min[0]) # <check>
        
        # face velocity and energy
        Face_Arrays['U_F'][0] = self.velocity_value( Face_Arrays['Q_F'][0] , Face_Arrays['A_F'][0])
        Face_Arrays['E_F'][0] = Inter.energy_value( Face_Arrays['U_F'][0] , Face_Arrays['Eta_F'][0], Gravity)
        
        return Face_Arrays
      
    def face_downstream_boundary_values(self, el, Face_Arrays, fge, h_dw, N_Cells, depth_min, face_area_min, Gravity ):
        '''compute face values at boundary given downstream height (h_dw)
           and interior cell N_Cells'''

        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()

        Face_Arrays['Eta_F'][N_Cells] = h_dw + fge['Z_F'][N_Cells]
        Face_Arrays['Eta_F'][N_Cells] = max( Face_Arrays['Eta_F'][N_Cells],    \
                                    depth_min + fge['Z_F'][N_Cells] )
        
        # area for rectangular channel only
        Face_Arrays['A_F'][N_Cells] = h_dw * fge['B_F'][N_Cells]
        Face_Arrays['A_F'][N_Cells] = max(Face_Arrays['A_F'][N_Cells], face_area_min[N_Cells])
        
        Face_Arrays['Q_F'][N_Cells]   = el['Q'][N_Cells-1] # extrapolate flow from interior
                
        Face_Arrays['U_F'][N_Cells]   = self.velocity_value( Face_Arrays['Q_F'][N_Cells] , Face_Arrays['A_F'][N_Cells])        
        Face_Arrays['E_F'][N_Cells]   = Inter.energy_value( Face_Arrays['U_F'][N_Cells] , Face_Arrays['Eta_F'][N_Cells], Gravity)
        
        return Face_Arrays
 
    def face_interior_energy_values(self, el, Face_Arrays, geo):

        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()

        for ii in range(1,N_Cells): 
            Face_Arrays['E_F'][ii] = Inter.linear_interp_to_face( el['E'], geo['L'], ii)
        return Face_Arrays

    def face_interior_hat_values(self, el, fge, Face_hat):
        ''' estimated values on interior faces'''

        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()

        # Q^2 with sign 
        Qsquared = self.exponential_of_abs_value_with_sign \
                    (el['Q'], 2.0)
        
        for ii in range(1,N_Cells): 
            Face_hat['Eta_F_hat'][ii] = Inter.linear_interp_to_face(el['Eta'], geo['L'], ii)
            Face_hat['A_F_hat'][ii] = Inter.linear_interp_to_face(el['A'], geo['L'], ii)                
            Face_hat['Q_F_hat_S'][ii] = Inter.linear_interp_to_face(Qsquared, geo['L'],ii)    
            Face_hat['Gamma_F_hat'][ii] = Inter.linear_interp_to_face(el['Gamma'], geo['L'], ii)    
            Face_hat['H_hat'][ii] = Inter.linear_interp_to_face(el['H'], geo['L'], ii)
          
        return Face_hat
    
    def face_corrections(self, Face_Arrays, Face_hat, fge):
        ''' Perform epsilon corrections to face interpolation'''
        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()

        for ii in range(1,N_Cells): 
            
            # selector for imaginary root conditions
            selector = np.zeros(3, dtype=np.float64)
            
            # solution to ax^2 + (b+1)x + c = 0
            t1 = 2.0 * Face_hat['Gamma_F_hat'][ii] * Face_hat['H_hat'][ii]**2.0  \
                        / Face_hat['A_F_hat'][ii]
                                        
            a = 1.0 + t1                                                 \
                - (Gravity * Face_hat['Eta_F_hat'][ii] / Face_Arrays['E_F'][ii])      \
                * (1.0 + t1)                                             \
                - 2.0  * Gravity * Face_hat['H_hat'][ii] / Face_Arrays['E_F'][ii]
                
            b = 2.0 - Gravity * ( 2.0 * Face_hat['Eta_F_hat'][ii]            \
                                 + Face_hat['H_hat'][ii] ) / Face_Arrays['E_F'][ii]
            
            c = 1.0 - (Gravity * Face_hat['Eta_F_hat'][ii] / Face_Arrays['E_F'][ii] ) \
                     - (Face_hat['Q_F_hat_S'][ii]**2.0)                      \
                     / ( 2.0 * Face_Arrays['E_F'][ii] * (Face_hat['A_F_hat'][ii]**2.0) )
            
            # sqrt term of basic correction
            selector[0] =  (b+1.0)**2.0 - 4.0 * a * c 
            # sqrt term of correction for Eta only
            selector[1] =   b**2.0 - 4.0 *a*c
            # PE > KE for correction of Eta only
            selector[2] =  Gravity * Face_hat['Eta_F_hat'][ii] \
                - 0.5 * (Face_hat['Q_F_hat_S'][ii] / Face_hat['A_F_hat'][ii])**2.0
                
            # correction term for PE    
            Eta_Epsilon = self.get_eta_epsilon(Face_hat['H_hat'][ii], a, b, c, selector)
            
            # corrections involving Eta
            Face_Arrays['Eta_F'][ii] = Face_hat['Eta_F_hat'][ii] + Eta_Epsilon  
            
            Face_Arrays['A_F'][ii] = Face_hat['A_F_hat'][ii]                           \
                + ( fge['B_F'][ii] + Face_hat['Gamma_F_hat'][ii] * Eta_Epsilon)\
                * Eta_Epsilon
            
            # correction term for KE
            Alfa_Epsilon = self.get_alpha_epsilon                          \
                ( Face_hat['H_hat'][ii], Face_Arrays['E_F'][ii], Face_hat['A_F_hat'][ii],   \
                  Face_Arrays['A_F'][ii], Face_Arrays['Eta_F'][ii], Face_hat['Q_F_hat_S'][ii],   \
                      a, b, c, Eta_Epsilon, selector )    
      
            # corrections involving alpha
            Face_Arrays['Q_F'][ii] = Inter.sgn_func(Face_hat['Q_F_hat_S'][ii])                 \
                * (abs(Face_hat['Q_F_hat_S'][ii] + Alfa_Epsilon))**0.5
                
            Face_Arrays['U_F'][ii] = self.velocity_value(Face_Arrays['Q_F'][ii], Face_Arrays['A_F'][ii])    
              
        return Face_Arrays
    

    def get_alpha_epsilon(self, Hhat, EF, AFhat, AF, EtaF, QFhatS, a, b, c, Eta_Epsilon, selector):
        ''' computes the alpha_epsilon for the selector cases'''
        rAlfa = 0.5 # fraction of Q^2 to adjust, must between 0 and 1.        
        if selector[0] >= 0:
            # standard approach
            x = Eta_Epsilon / Hhat
            Alfa_Epsilon  = 2.0 * EF * (AFhat**2.0) \
                * ( a * (x**2.0) + b * x + c )
        elif selector[1] > 0 and selector[2] >= 0:
            # imaginary selector(0), correcting Eta only
            Alfa_Epsilon  = 0.0
        else:
            # imaginary selector(0), correcting Q only
            Alfa_Epsilon  = 2.0 * (AF**2.0) * (EF - Gravity * EtaF) - QFhatS          

        # limit the correction to prevent sign change
        if  Alfa_Epsilon < -rAlfa * QFhatS**2:
            Alfa_Epsilon = -rAlfa * QFhatS**2
        else:
            Alfa_Epsilon = rAlfa * QFhatS**2
                    
        return Alfa_Epsilon

    def get_eta_epsilon(self, Hhat, a, b, c, selector):
        ''' finds the root with the smallest magnitude'''
        if selector[0] >= 0:
            # non-imaginary roots
            t1 = ((b+1.0)**2.0 - 4.0*a*c)**0.5
            e1 =  -b - 1.0 + t1  
            e2 =  -b - 1.0 - t1  
        
        elif selector[1] > 0 and selector[2] >= 0:
            # imaginary roots, correcting Eta only
            t1 = (b**2.0 - 4.0*a*c)**0.5
            e1 = -b + t1 
            e2 = -b - t1 
        else:
            # imaginary roots, correct Q only
            e1 = 0.0    
            e2 = 0.0
            
        e1 = Hhat * e1 / (2.0 * a)    
        e2 = Hhat * e2 / (2.0 * a) 
         
        if abs(e1) < abs(e2): # Picks the closest to zero
            Eta_Epsilon    =  e1
        else:
            Eta_Epsilon    =  e2
    
        return Eta_Epsilon
    
    def get_friction_half_element(self, Ff, el, Face_Arrays, inc, N_Cells, Gravity):
        ''' friction term on 1/2 of an element, inc = 0 
            for upstream, +1 for downstream'''

        import Interpolation_Class as Interpolation
        Inter = Interpolation.Interpolation()

        for ii in range(0,N_Cells):     
            Ff[ii] = Gravity * el['C'][ii] * el['V'][ii]                     \
                * Inter.sgn_func(Face_Arrays['U_F'][ii+inc] + el['U'][ii]) \
                * ( ( Face_Arrays['U_F'][ii+inc] + el['U'][ii]     )**2.0) / 8.0
                
        return Ff       
    
    def velocity_value(self, QQ, AA):        
        return QQ / AA

    def continuity_RHS_explicit(self, QF, ii ):
        return ( QF[ii] - QF[ii+1] ) 
    
    def momentum_RHS_explicit(self, el, Face_Arrays, geo, ii, Gravity):
        ''' Common RHS to explicit Euler and RK solutions''' 
        return (  Face_Arrays['Q_F'][ii]   * Face_Arrays['U_F'][ii]            \
                - Face_Arrays['Q_F'][ii+1] * Face_Arrays['U_F'][ii+1]          \
                + Gravity * el['A'][ii]                                        \
                   * ( Face_Arrays['Eta_F'][ii] - Face_Arrays['Eta_F'][ii+1] ) \
                - el['Fdn'][ii] - el['Fup'][ii]) / geo['L'][ii] 
#        return DT * ( Face_Arrays['Q_F'][ii]   * Face_Arrays['U_F'][ii]                         \
#                    - Face_Arrays['Q_F'][ii+1] * Face_Arrays['U_F'][ii+1]                       \
#                    + Gravity * (  Face_Arrays['A_F'][ii]   * Face_Arrays['Eta_F'][ii]          \
#                                 - Face_Arrays['A_F'][ii+1] * Face_Arrays['Eta_F'][ii+1] )      \
#                    - el['Fdn'][ii] - el['Fup'][ii]) 


