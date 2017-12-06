
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

    # definitions:
    # the ii=0 cell center is first cell inside the upstream boundary
    # the ii=N_Cells-1 cell center is first cell inside the downstream boundary 
    # the cells outside the downstream boundary (ii=N_Cells) do not exist
    # faces take on the index of the downstream cell.
    # the ii=0 face is the upstream boundary
    # the ii=N_Cell face is the downstream boundary
    # Because of pythong indexing
    #   range(0,N_Cells) gives you all cells (which are interior)
    #   range(1,N_Cells) gives you all interior faces
    #   range(0,N_Cells+1) gives you all faces
    #   Var_F[N_Cell] = furthest downstream face value
    #   Var_E(N_Cell-1) = furthest downstream cell
    
    def __init__(self): 
        pass

    def sgn_func(self, x): 
        ''' returns sign of argument '''
        return 1.0 if x >= 0 else -1.0
    
    def sgn_func_array(self, xx):
        ''' returns array with sign of a numpy array '''
        xsgn = np.ones(xx.size,dtype=int)
        aa   = np.nonzero(xx < 0)
        xsgn[aa] = -1        
        return xsgn
    
    def linear_interp_to_face(self, cc, L, ii):
        ''' generic linear interpolation from centers to face'''
        return ( L[ii] * cc[ii-1] + L[ii-1] * cc[ii] )  \
                / ( L[ii] + L[ii-1] ) 

    def linear_interp_to_face_array(self, cc, LL):
        ''' generic linear interpolation from centers to face'''
        #for ii in range(1:N_Cells)
        #    xx[ii] = ( L[ii] * cc[ii-1] + L[ii-1] * cc[ii] )  \
        #            / ( L[ii] + L[ii-1] )        
        return (  LL[1:N_Cells]   * cc[0:N_Cells-1]                          \
                + LL[0:N_Cells-1] * cc[1:N_Cells]   )                        \
              / ( LL[1:N_Cells] + LL[0:N_Cells-1] )
                
    def linear_interp_to_face_onevalue(self, cc1, cc2,  L1, L2 ):
        ''' generic linear interpolation from centers to face'''
        return ( L1 * cc2 + L2 * cc1 ) / ( L1 + L2 ) 
  
    def linear_extrapolate_to_upstream_face(self, cc, L):
        ''' generic linear extrapolation to upstream boundary'''       
        return cc[0] + L[0] *  ( cc[0] - cc[1]) / ( L[0] + L[1]) 
        
    def energy_value(self, U, Eta):
        '''Compute the energy from velocity and free surface elevation'''
        return ( 0.5 * U**2.0 ) + Gravity * Eta

    def cbrt(self, x):
        '''Computes the cube root of a number.
        https://stackoverflow.com/questions/16270941/whats-wrong-with-this-python-function-to-solve-cubic-equations
        Copyright (c) 2013 user2330618
        This Source Code Form is subject to the terms of the Mozilla Public
        License, v. 2.0. If a copy of the MPL was not distributed with this
        file, you can obtain one at http://www.mozilla.org/MPL/2.0/.
        '''      
        import cmath 
        
        if x.imag != 0.0:    # Finding the cubic root of a complex number
            return cmath.exp(cmath.log(x) / 3)
        else:
            if x < 0:        # Finding the cubic root of a real number
                d = (-x) ** (1 / 3)
                return -d
            elif x >= 0:
                return x ** (1 / 3)
            
    def cubic(self, a, b, c, d):
        '''Returns the real roots to cubic equations in expanded form.
        https://stackoverflow.com/questions/16270941/whats-wrong-with-this-python-function-to-solve-cubic-equations
        Copyright (c) 2013 user2330618
        This Source Code Form is subject to the terms of the Mozilla Public
        License, v. 2.0. If a copy of the MPL was not distributed with this
        file, you can obtain one at http://www.mozilla.org/MPL/2.0/.
        '''
        import cmath
        
        # cutoff for small imaginary parts to be neglected
        phase_tolerance = 1.0e-14
        round_small_values = True
        xx = np.array([0,0,0], dtype=np.complex128)

        # Define the discriminants
        D = (18 * a * b * c * d) - (4 * (b ** 3) * d) + ((b ** 2) * (c ** 2)) - \
        (4 * a * (c ** 3)) - (27 * (a ** 2) * d ** 2)
        D0 = (b ** 2) - (3 * a * c)
        # Test for some special cases
        if D == 0 and D0 == 0:
            xx[0] = -(b / (3 * a))
        elif D == 0 and D0 != 0:
            xx[0] = [((b * c) - (9 * a * d)) / (-2 * D0), ((b ** 3) - (4 * a * b * c) 
            + (9 * (a ** 2) * d)) / (-a * D0)]
        else:
            D1 = (2 * (b ** 3)) - (9 * a * b * c) + (27 * (a ** 2) * d)
            # More special cases
            if D != 0 and D0 == 0 and D1 < 0:
                C = self.cbrt((D1 - cmath.sqrt((D1 ** 2) - (4 * (D0 ** 3)))) / 2)
            else:
                C = self.cbrt((D1 + cmath.sqrt((D1 ** 2) - (4 * (D0 ** 3)))) / 2)
                u_2 = (-1 + (1j * cmath.sqrt(3))) / 2
                u_3 = (-1 - (1j * cmath.sqrt(3))) / 2

                xx[0] = (-(b + C + (D0 / C))) / (3 * a)
                xx[1] = (-(b + (u_2 * C) + (D0 / (u_2 * C)))) / (3 * a)
                xx[2] = (-(b + (u_3 * C) + (D0 / (u_3 * C)))) / (3 * a)
                
        if round_small_values:
            for ii in range(0,3):
                if     abs(cmath.phase(xx[ii])) == cmath.pi \
                    or abs(cmath.phase(xx[ii])) < phase_tolerance:                  
                    xx[ii] = xx[ii].real
                  
        if D > 0:
            return xx
        else:
            return xx[0]    
            
    def cubic_without_a_c(self, b,d):
        '''Returns the roots to cubic equation.
        for cubic ax^3+bx^2+cx+d = 0 where a=1 and c=0
        https://stackoverflow.com/questions/16270941/whats-wrong-with-this-python-function-to-solve-cubic-equations
        Copyright (c) 2013 user2330618
        This Source Code Form is subject to the terms of the Mozilla Public
        License, v. 2.0. If a copy of the MPL was not distributed with this
        file, you can obtain one at http://www.mozilla.org/MPL/2.0/.
        '''
        import cmath
        
        #print("Find roots...")

        # cutoff for small imaginary parts to be neglected brh20171126
        phase_tolerance = 1.0e-14
        round_small_values = True
        xx = np.array([0,0,0], dtype=np.complex128)

        # Define the discriminants
        D = - (4 * (b ** 3) * d) - (27  * d ** 2)
        D0 = (b ** 2)
        # Test for some special cases
        if D == 0 and D0 == 0:
            xx[0] = -(b / (3))
        elif D == 0 and D0 != 0:
            xx[0] = [( - (9 *  d)) / (-2 * D0), ((b ** 3) + (9 * d)) / (- D0)]
        else:
            D1 = (2 * (b ** 3))  + (27  * d)
            # More special cases
            if D != 0 and D0 == 0 and D1 < 0:
                C = self.cbrt((D1 - cmath.sqrt((D1 ** 2) - (4 * (D0 ** 3)))) / 2)
            else:
                C = self.cbrt((D1 + cmath.sqrt((D1 ** 2) - (4 * (D0 ** 3)))) / 2)
                u_2 = (-1 + (1j * cmath.sqrt(3))) / 2
                u_3 = (-1 - (1j * cmath.sqrt(3))) / 2

                xx[0] = (-(b + C + (D0 / C))) / (3 )
                xx[1] = (-(b + (u_2 * C) + (D0 / (u_2 * C)))) / (3)
                xx[2] = (-(b + (u_3 * C) + (D0 / (u_3 * C)))) / (3)
         
        #brh added 20171126    
        if round_small_values:
            for ii in range(0,3):
                if     abs(cmath.phase(xx[ii])) == cmath.pi \
                    or abs(cmath.phase(xx[ii])) < phase_tolerance:                  
                    xx[ii] = xx[ii].real
                  
        if D > 0:
            return xx
        else:
            return xx[0]    
   
    def element_values(self, el, geo):
        '''Compute auxiliary values on element center
           Assumes Q, V and geometry are known
           Friction (Fdn, Fup) are computed later as they require face values
        '''
        #for ii in range(0,N_Cells):
        #    print(ii, el['V'][ii])            
        #sys.exit()    
        
        global quit_after_plot


        # depends on input
        # note that we assume V > volume_min, which must be checked when
        # V is updated
        el['A'][:]   = el['V'][:] / geo['L'][:]
        #for ii in range(0,N_Cells): 
        #    el['A'][ii]   = el['V'][ii] / geo['L'][ii]
       
        # depends onrectangular cell assumption and input
        el['Eta'][:] = geo['Z'][:] + el['A'][:] / geo['B'][:]
        el['T'][:]   = geo['B'][:]
        #for ii in range(0,N_Cells): 
        #    el['Eta'][ii] = geo['Z'][ii] + el['A'][ii] / geo['B'][ii]
        #    el['T'][ii]   = geo['B'][ii]
        
        # depends on A and method to get T
        el['H'][:]   = el['A'][:] / el['T'][:]
        #for ii in range(0,N_Cells): 
        #    el['H'][ii]   = el['A'][ii] / el['T'][ii]
     
        # depends onrectangular cell assumption and input
        el['l_P'][:] = geo['B'][:] + 2.0*el['H'][:]
        el['Gamma'][:] = 0.0 
        #for ii in range(0,N_Cells): 
        #    el['l_P'][ii] = geo['B'][ii] + 2.0*el['H'][ii]
        #    el['Gamma'][ii] = 0.0 

        # depends on A, cell cross-section, and input
        el['U'][:]   = self.velocity_value( el['Q'][:] , el['A'][:])
        el['R_h'][:] = el['A'][:] / el['l_P'][:]
        #for ii in range(0,N_Cells): 
        #    el['U'][ii]   = self.velocity_value( el['Q'][ii] , el['A'][ii])
        #    el['R_h'][ii] = el['A'][ii] / el['l_P'][ii]

        # depends on Eta and U
        el['E'][:]   = self.energy_value( el['U'][:], el['Eta'][:] )
        #for ii in range(0,N_Cells): 
        #    el['E'][ii]   = self.energy_value( el['U'][ii], el['Eta'][ii] )

        # depends on R_h
        el['C'][:]   = ((geo['M'][:])**2.0) / (el['R_h'][:]**(4.0/3.0))
        #for ii in range(0,N_Cells): 
        #    el['C'][ii]   = ((geo['M'][ii])**2.0) / (el['R_h'][ii]**(4.0/3.0))

        # non-dimensional numbers
        el['Fr'][:]  = abs(el['U'][:]) / ( (Gravity * el['H'][:] )**0.5)               
        el['CFLu'][:] = abs(el['U'][:]) * DT / geo['L'][:]
        el['CFLb'][:] = (Gravity * el['H'][:]) * DT / geo['L'][:]
        #for ii in range(0,N_Cells): 
        #    el['Fr'][ii]  = abs(el['U'][ii]) / ( (Gravity * el['H'][ii] )**0.5)               
        #    el['CFLu'][ii] = abs(el['U'][ii]) * DT / geo['L'][ii]
        #    el['CFLb'][ii] = (Gravity * el['H'][ii]) * DT / geo['L'][ii]
        
        
        # compact version that loops over cells for debugging    
        #for ii in range(0,N_Cells):                        
            #el['A'][ii]   = el['V'][ii] / geo['L'][ii]                                  
            #el['Eta'][ii] = geo['Z'][ii] + el['A'][ii] / geo['B'][ii]
            #el['T'][ii]   = geo['B'][ii]
            #el['H'][ii]   = el['A'][ii] / el['T'][ii]            
            #el['l_P'][ii] = geo['B'][ii] + 2.0*el['H'][ii]
            #el['Gamma'][ii] = 0.0 
            #el['U'][ii]   = self.velocity_value( el['Q'][ii] , el['A'][ii])
            #el['R_h'][ii] = el['A'][ii] / el['l_P'][ii]
            #el['E'][ii]   = self.energy_value( el['U'][ii], el['Eta'][ii] )            
            #el['C'][ii]   = ((geo['M'][ii])**2.0) / (el['R_h'][ii]**(4.0/3.0))
            #el['Fr'][ii]  = abs(el['U'][ii]) / ( (Gravity * el['H'][ii] )**0.5)               
            #el['CFLu'][ii] = abs(el['U'][ii]) * DT / geo['L'][ii]
            #el['CFLb'][ii] = (Gravity * el['H'][ii]) * DT / geo['L'][ii]


        if error_check:
            for ii in range(0,N_Cells): 
                if el['Fr'][ii] >= Fr_max:
                    print(ii, el['Q'][ii], el['A'][ii], el['U'][ii], el['H'][ii])
                    print("Flow is above max Fr %d %d %f" % (nn, ii, el['Fr'][ii]))
                    #check = input(" Error: press ENTER to exit ")
                    quit_after_plot = True  
                
                if el['CFLb'][ii] > CFLb_max:
                    print("CFLb_max exceeded %d %d %f" % (nn, ii, el['CFLb'][ii]))
                    #check = input(" Error: press ENTER to exit ")
                    #if DT_variable == 0:
                    quit_after_plot = True
    
                if el['CFLu'][ii] > CFLu_max:
                    print('Cell=',ii,'Q=',el['Q'][ii],'A=',el['A'][ii],'U=',el['U'][ii],'H=',el['H'][ii])                
                    print("CFLu_max exceeded %d %d %f" % (nn, ii, el['CFLu'][ii]))
                    #check = input(" Error: press ENTER to exit ")
                    #if DT_variable == 0:
                    quit_after_plot = True
           
        return el
    
    def exponential_of_abs_value_with_sign(self, value, power):

        for ii in range(0,value.size-1):
            value[ii] = self.sgn_func(value[ii]) * (abs(value[ii]))**power
        
        return value

    def face_values(self, el, fa, geo, fge, Q_Upstream, h_dw):
        '''interpolates values to faces'''
               
        fa = self.face_upstream_boundary_values( el, fa, geo, fge, Q_Upstream)
        
        fa = self.face_downstream_boundary_values( el, fa, fge, h_dw )
                        
#        if face_correct == 0:
        fa = self.face_by_interpolation( el, fa, geo, fge)
#        elif face_correct ==1:
#            fa   = self.face_interior_energy_values( el, fa, geo )
#            fhat = self.face_interior_hat_values( el, fge, fhat)
#            fa   = self.face_corrections(fa, fhat, fge)
#        else:
#            print('Unknown value of face_correct = ',face_correct)   
#            sys.exit()
                
#        for ii in range(0,N_Cells+1):
#            print(ii, fhat['A_F_hat'][ii])
#        print('at 108')    
#        sys.exit()
                       
        return fa  
    
    def face_by_interpolation(self, el, fa, geo, fge):
        '''Interpolation for all face values'''

        import time

#        if face_interp == 1:                
        for ii in range(1,N_Cells):
            # direction of fluxes on adjacent cells
            sup = self.sgn_func( el['Q'][ii-1])  
            sdn = self.sgn_func( el['Q'][ii])
            
            if sup * sdn != -1: # not a convergent or divergent flow
                # flow in both cells in the same direction
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
                
                if el['Fr'][ii-upi] >= 1.0:
                    # supercritical: 
                    # upstream energy change from friction
                    SpecificEnergy = el['E'][upi] - geo['Z'][upi] * Gravity
                    fa['Q_F'][ii]  = el['Q'][upi]
                else:
                    # subcritical
                    SpecificEnergyUp = el['E'][upi]  - geo['Z'][upi] * Gravity
                    SpecificEnergyDn = el['E'][dni]  - geo['Z'][dni] * Gravity
                    SpecificEnergy = self.linear_interp_to_face_onevalue     \
                            (SpecificEnergyUp, SpecificEnergyDn,             \
                             geo['L'][upi], geo['L'][dni] )
                    fa['Q_F'][ii] = self.linear_interp_to_face            \
                            ( el['Q'], geo['L'], ii )    
                   
                # transverse width at face (need storage?) HACK   
                faT = self.linear_interp_to_face( el['T'], geo['L'], ii) 
                
                # compute roots of cubic to find area
                pp = ( 1,                                             \
                      -SpecificEnergy * faT / Gravity,                \
                       0,                                             \
                      +0.5 * faT * (fa['Q_F'][ii]**2.0) / Gravity )


                #start_time = time.time()
                #print(" before mask ...") # <delete>
                # python root finding (slow)
                #rootout = np.roots(pp) 
                # general cubic
                #rootout = self.cubic(pp[0], pp[1], pp[2] ,pp[3]) 
                # cubic with reduced form using a=1 and c=0
                rootout = self.cubic_without_a_c(pp[1],pp[3])
                #print(1,rootout)

                #elapsed_time = time.time() - start_time
                #print(" solve the third-order: ",elapsed_time)

                # remove any complex roots
                #print()
                #print(ii)
                #print("Roots before mask",rootout) 
                #rootout2 = rootout

                #start_time = time.time()
                #rootout = np.ma.masked_where(rootout.imag != 0.0, rootout, copy=False)
                #print("Roots after mask", rootout)
                #rootout = rootout.real
                #print("Roots after real", rootout)
                #print(2,rootout)
                # remove any negative roots
                #rootout = np.ma.masked_where(rootout <= 0, rootout, copy=False)
                #print("Roots make positive", rootout)
                #rootout = np.ma.compressed(rootout)
                #elapsed_time = time.time() - start_time
                #print(" mask                 : ",elapsed_time)
                #print("Roots after compressing", rootout)
                #check = input("check")
                #print(3,rootout)
                #print(min(rootout.data))
                #print(max(rootout.data))


                rootout_modified = []
                #print()
                #print("Roots before mask2",rootout) 
                #start_time = time.time()
                for Inumber in rootout:
                    if Inumber.imag == 0.0 and Inumber.real>0.0:
                        rootout_modified.append(Inumber.real)

                #elapsed_time = time.time() - start_time
                #print(" mask                 : ",elapsed_time)
                #print("New Roots", rootout_modified)
                #check = input("check")

                if el['Fr'][upi] >= 1.0:
                    fa['A_F'][ii] = min(rootout_modified)
                else:
                    fa['A_F'][ii] = max(rootout_modified)
                    
                fa['A_F'][ii] = max( face_area_min[ii], fa['A_F'][ii] ) 
                
                fa['U_F'][ii] = self.velocity_value(fa['Q_F'][ii], fa['A_F'][ii])
                
                fa['Eta_F'][ii] = fge['Z_F'][ii]                          \
                    +(SpecificEnergy - 0.5 * fa['U_F'][ii]**2.0) / Gravity
                    
                fa['E_F'][ii] = self.energy_value(fa['U_F'][ii], fa['Eta_F'][ii])
                    
                if fa['E_F'][ii] < fge['Z_F'][ii]:
                    print('problem in face interpolation, insufficient energy')
                    sys.exit()
                if el['Eta'][ii-upi] <= fge['Z_F'][ii]:
                    print('Possible hydraulic jump required - not yet coded')
                    sys.exit()
       
            else: 
                # flow is converging or diverging on face
                print("flow is converging or diverging on face")
                check  = input(" Hello")
                fa['Q_F'][ii] = el['Q'][ii-1] + el['Q'][ii]
                
                fa['A_F'][ii] = self.linear_interp_to_face                  \
                             ( el['A'], geo['L'], ii )                            
                fa['A_F'][ii] = max( fa['A_F'][ii], face_area_min(ii))  
                
                fa['Eta_F'][ii] = self.linear_interp_to_face                \
                             ( el['Eta'], geo['L'], ii )
                fa['Eta_F'][ii] = max( fa['Eta_F'][ii],                     \
                                       depth_min + fge['Z_F'][ii])  
                
                fa['U_F'][ii] = self.velocity_value(fa['Q_F'][ii], fa['A_F'][ii])
                fa['E_F'][ii] = self.energy_value(fa['U_F'][ii], fa['Eta_F'][ii])
                

#        print('=========================== in face interp')
#        for ii in range(955,966):
#            print(ii, 'A:', el['A'][ii-1], fa['A_F'][ii], el['A'][ii])     
#
#        print()                       
        return fa
           
    def face_interior_bottom_elevation(self, fge, geo):
        '''linear interpolation of bottom elevation between cell centers'''
        # note this only needs to be called once in a solution
        for ii in range(1,N_Cells):
            fge['Z_F'][ii] = self.linear_interp_to_face                       \
                            ( geo['Z'], geo['L'], ii )    
                    
        return fge

    def face_upstream_bottom_elevation(self, fge, geo):
        '''upstream face Z by extrapolation'''
        fge['Z_F'][0] = self.linear_extrapolate_to_upstream_face              \
                        (geo['Z'], geo['L'])      
                        
        return fge
 
    def face_downstream_bottom_elevation(self, fge, geo):
        ''' temporary hard code of bottom elevation as 0''' 
        fge['Z_F'][N_Cells] = 0.0 #HACK
       
        return fge

    def face_breadth_rectangular_channel(self, fge, geo):
        ''' simple linear interpolation and boundary extrapolation'''
        # handles both interior and boundaries
        for ii in range(1,N_Cells):
            fge['B_F'][ii] = self.linear_interp_to_face                       \
                            ( geo['B'], geo['L'], ii)
        fge['B_F'][0] =  geo['B'][0]
        fge['B_F'][N_Cells] = geo['B'][N_Cells-1]
        
        return fge
       
    def face_upstream_boundary_values(self, el, fa, geo, fge, Q_Upstream):
        '''Compute face values at boundary given Q upstream element values 
            at cell ii=0
        '''
        fa['Q_F'][0] = Q_Upstream  # use input boundary condition
        
        # extrapolate bottom gradient to get surface elevation on upstream face
        # note this CANNOT be replaced with linear_extrapolate_to_upstream_face
        fa['Eta_F'][0] = el['Eta'][0] + geo['L'][0]   \
            *  ( geo['Z'][0] - geo['Z'][1])           \
             / ( geo['L'][1] + geo['L'][2]) 
        fa['Eta_F'][0] = max( fa['Eta_F'][0], depth_min + fge['Z_F'][0]) 

        # area for rectangular channel only
        fa['A_F'][0] = (fa['Eta_F'][0] - fge['Z_F'][0]) * fge['B_F'][0]
        fa['A_F'][0] = max(fa['A_F'][0], face_area_min[0])
        
        # face velocity and energy
        fa['U_F'][0] = self.velocity_value( fa['Q_F'][0] , fa['A_F'][0])
        fa['E_F'][0] = self.energy_value( fa['U_F'][0] , fa['Eta_F'][0])
        
        return fa
      
    def face_downstream_boundary_values(self, el, fa, fge, h_dw):
        '''compute face values at boundary given downstream height (h_dw)
           and interior cell N_Cells'''
        fa['Eta_F'][N_Cells] = h_dw + fge['Z_F'][N_Cells]
        fa['Eta_F'][N_Cells] = max( fa['Eta_F'][N_Cells],                     \
                                    depth_min + fge['Z_F'][N_Cells] )
        
        # area for rectangular channel only
        fa['A_F'][N_Cells] = h_dw * fge['B_F'][N_Cells] 
        fa['A_F'][N_Cells] = max(fa['A_F'][N_Cells], face_area_min[N_Cells])
        
        fa['Q_F'][N_Cells]   = el['Q'][N_Cells-1] # extrapolate flow from interior
                
        fa['U_F'][N_Cells]   = self.velocity_value( fa['Q_F'][N_Cells] , fa['A_F'][N_Cells])        
        fa['E_F'][N_Cells]   = self.energy_value( fa['U_F'][N_Cells] , fa['Eta_F'][N_Cells])
        
        return fa
 
    def face_interior_energy_values(self, el, fa, geo):
        for ii in range(1,N_Cells): 
            fa['E_F'][ii] = self.linear_interp_to_face( el['E'], geo['L'], ii)
        return fa

    def face_interior_hat_values(self, el, fge, fhat):
        ''' estimated values on interior faces'''
        # Q^2 with sign 
        Qsquared = self.exponential_of_abs_value_with_sign \
                    (el['Q'], 2.0)
        
        for ii in range(1,N_Cells): 
            fhat['Eta_F_hat'][ii] = self.linear_interp_to_face     \
                                    ( el['Eta'], geo['L'], ii)
                                    
            fhat['A_F_hat'][ii] = self.linear_interp_to_face \
                                    ( el['A'], geo['L'], ii)                
            
            fhat['Q_F_hat_S'][ii] = self.linear_interp_to_face \
                                    ( Qsquared, geo['L'],ii)    
            fhat['Gamma_F_hat'][ii] = self.linear_interp_to_face \
                                    ( el['Gamma'], geo['L'], ii)    
            fhat['H_hat'][ii] = self.linear_interp_to_face \
                                    ( el['H'], geo['L'], ii)
          
        return fhat
    
    def face_corrections(self, fa, fhat, fge):
        ''' Perform epsilon corrections to face interpolation'''
        for ii in range(1,N_Cells): 
            
            # selector for imaginary root conditions
            selector = np.zeros(3, dtype=np.float64)
            
            # solution to ax^2 + (b+1)x + c = 0
            t1 = 2.0 * fhat['Gamma_F_hat'][ii] * fhat['H_hat'][ii]**2.0  \
                        / fhat['A_F_hat'][ii]
                                        
            a = 1.0 + t1                                                 \
                - (Gravity * fhat['Eta_F_hat'][ii] / fa['E_F'][ii])      \
                * (1.0 + t1)                                             \
                - 2.0  * Gravity * fhat['H_hat'][ii] / fa['E_F'][ii]
                
            b = 2.0 - Gravity * ( 2.0 * fhat['Eta_F_hat'][ii]            \
                                 + fhat['H_hat'][ii] ) / fa['E_F'][ii]
            
            c = 1.0 - (Gravity * fhat['Eta_F_hat'][ii] / fa['E_F'][ii] ) \
                     - (fhat['Q_F_hat_S'][ii]**2.0)                      \
                     / ( 2.0 * fa['E_F'][ii] * (fhat['A_F_hat'][ii]**2.0) )
            
            # sqrt term of basic correction
            selector[0] =  (b+1.0)**2.0 - 4.0 * a * c 
            # sqrt term of correction for Eta only
            selector[1] =   b**2.0 - 4.0 *a*c
            # PE > KE for correction of Eta only
            selector[2] =  Gravity * fhat['Eta_F_hat'][ii] \
                - 0.5 * (fhat['Q_F_hat_S'][ii] / fhat['A_F_hat'][ii])**2.0
                
            # correction term for PE    
            Eta_Epsilon = self.get_eta_epsilon \
                (fhat['H_hat'][ii], a, b, c, selector)
            
            # corrections involving Eta
            fa['Eta_F'][ii] = fhat['Eta_F_hat'][ii] + Eta_Epsilon  
            
            fa['A_F'][ii] = fhat['A_F_hat'][ii]                           \
                + ( fge['B_F'][ii] + fhat['Gamma_F_hat'][ii] * Eta_Epsilon)\
                * Eta_Epsilon
            
            # correction term for KE
            Alfa_Epsilon = self.get_alpha_epsilon                          \
                ( fhat['H_hat'][ii], fa['E_F'][ii], fhat['A_F_hat'][ii],   \
                  fa['A_F'][ii], fa['Eta_F'][ii], fhat['Q_F_hat_S'][ii],   \
                      a, b, c, Eta_Epsilon, selector )    
      
            # corrections involving alpha
            fa['Q_F'][ii] = self.sgn_func(fhat['Q_F_hat_S'][ii])                 \
                * (abs(fhat['Q_F_hat_S'][ii] + Alfa_Epsilon))**0.5
                
            fa['U_F'][ii] = self.velocity_value(fa['Q_F'][ii], fa['A_F'][ii])    
              
        return fa
    

    def get_alpha_epsilon(self, Hhat, EF, AFhat, AF, EtaF, QFhatS,
                          a, b, c, Eta_Epsilon, selector ):
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
    
    def get_friction_half_element(self, Ff, el, fa, inc):
        ''' friction term on 1/2 of an element, inc = 0 
            for upstream, +1 for downstream'''
        for ii in range(0,N_Cells):     
            Ff[ii] = Gravity * el['C'][ii] * el['V'][ii]                     \
                * self.sgn_func(fa['U_F'][ii+inc] + el['U'][ii]) \
                * ( ( fa['U_F'][ii+inc] + el['U'][ii]     )**2.0) / 8.0
                
        return Ff       
    
    def velocity_value(self, QQ, AA):        
        return QQ / AA

    def continuity_RHS_explicit(self, QF, ii ):
        return ( QF[ii] - QF[ii+1] ) 
    
    def momentum_RHS_explicit(self, el, fa, geo, ii):
        ''' Common RHS to explicit Euler and RK solutions''' 
        return (  fa['Q_F'][ii]   * fa['U_F'][ii]                         \
                - fa['Q_F'][ii+1] * fa['U_F'][ii+1]                       \
                + Gravity * el['A'][ii]                                   \
                   * ( fa['Eta_F'][ii] - fa['Eta_F'][ii+1] )              \
                - el['Fdn'][ii] - el['Fup'][ii]) / geo['L'][ii] 
#        return DT * ( fa['Q_F'][ii]   * fa['U_F'][ii]                         \
#                    - fa['Q_F'][ii+1] * fa['U_F'][ii+1]                       \
#                    + Gravity * (  fa['A_F'][ii]   * fa['Eta_F'][ii]          \
#                                 - fa['A_F'][ii+1] * fa['Eta_F'][ii+1] )      \
#                    - el['Fdn'][ii] - el['Fup'][ii]) 



    def RK4(self, el, fa, geo, fge, rke, rkf, Q_Upstream, h_dw):
                       
        for ii in range(0,N_Cells):
            rke['k_1V'][ii] = self.continuity_RHS_explicit( fa['Q_F'], ii )
            rke['k_1Q'][ii] = self.momentum_RHS_explicit( el, fa, geo, ii ) 
            
        rke['V'][:] = el['V'][:] + 0.5 * DT * rke['k_1V'][:]
        rke['Q'][:] = el['Q'][:] + 0.5 * DT * rke['k_1Q'][:]
        #for ii in range(0,N_Cells):
        #    rke['V'][ii] = el['V'][ii] + 0.5 * DT * rke['k_1V'][ii]
        #    rke['Q'][ii] = el['Q'][ii] + 0.5 * DT * rke['k_1Q'][ii]
    
                
        rke = self.element_values( rke, geo )
        rkf = self.face_values( rke, rkf, geo, fge, Q_Upstream, h_dw )
        
        for ii in range(0,N_Cells):
            rke['k_2V'][ii] = self.continuity_RHS_explicit( rkf['Q_F'], ii )
            rke['k_2Q'][ii] = self.momentum_RHS_explicit( rke, rkf, geo,  ii ) 
            
        rke['V'][:] = el['V'][:] + 0.5 * DT * rke['k_2V'][:]
        rke['Q'][:] = el['Q'][:] + 0.5 * DT * rke['k_2Q'][:]
        #rke['V'][:] = el['V'][:] + 0.5 * DT * rke['k_2V'][:]    WHY TWICE???
        #rke['Q'][:] = el['Q'][:] + 0.5 * DT * rke['k_2Q'][:]
        #for ii in range(0,N_Cells):
        #    rke['V'][ii] = el['V'][ii] + 0.5 * DT * rke['k_2V'][ii]
        #    rke['Q'][ii] = el['Q'][ii] + 0.5 * DT * rke['k_2Q'][ii]
        #    rke['V'][ii] = el['V'][ii] + 0.5 * DT * rke['k_2V'][ii]
        #    rke['Q'][ii] = el['Q'][ii] + 0.5 * DT * rke['k_2Q'][ii]
       
        rke = self.element_values( rke, geo )
        rkf = self.face_values( rke, rkf, geo, fge, Q_Upstream, h_dw )

        for ii in range(0,N_Cells):
            rke['k_3V'][ii] = self.continuity_RHS_explicit( rkf['Q_F'], ii )
            rke['k_3Q'][ii] = self.momentum_RHS_explicit( rke, rkf, geo, ii ) 
            
        rke['V'][:] = el['V'][:] + DT * rke['k_3V'][:]
        rke['Q'][:] = el['Q'][:] + DT * rke['k_3Q'][:]
        #for ii in range(0,N_Cells):
        #    rke['V'][ii] = el['V'][ii] + DT * rke['k_3V'][ii]
        #    rke['Q'][ii] = el['Q'][ii] + DT * rke['k_3Q'][ii]
            
        rke = self.element_values( rke, geo )
        rkf = self.face_values( rke, rkf, geo, fge, Q_Upstream, h_dw )

        for ii in range(0,N_Cells):
            rke['k_4V'][ii] = self.continuity_RHS_explicit( rkf['Q_F'], ii )
            rke['k_4Q'][ii] = self.momentum_RHS_explicit( rke, rkf, geo,  ii ) 
            
        el['V'][:] = el['V'][:]                                         \
          + (DT/6.0) * (    rke['k_1V'][:] + 2*rke['k_2V'][:]                \
                      + 2.0*rke['k_3V'][:] +   rke['k_4V'][:])
        #for ii in range(0,N_Cells):
        #    el['V'][ii] = el['V'][ii]                                         \
        #      + (DT/6) * (    rke['k_1V'][ii] + 2*rke['k_2V'][ii]                \
        #                  + 2*rke['k_3V'][ii] +   rke['k_4V'][ii])

        el['Q'][:] = el['Q'][:]                                         \
          + (DT/6.0) * (  rke['k_1Q'][:] + 2*rke['k_2Q'][:]                \
                      + 2.0*rke['k_3Q'][:]   + rke['k_4Q'][:])
        #for ii in range(0,N_Cells):
        #    el['Q'][ii] = el['Q'][ii]                                         \
        #      + (DT/6) * (    rke['k_1Q'][ii] + 2*rke['k_2Q'][ii]                \
        #                  + 2*rke['k_3Q'][ii]   + rke['k_4Q'][ii])

        return el
    
    def RK2(self, rke, rkf, el, fa, fhat, geo, fge, Q_Upstream, h_dw):  
        ''' RK2 time advance '''    
        
        print('BRH 20171123 error not tested')
        sys.exit()
        
        rke = self.RK2_step1( rke,  el, fa, geo ) 
              
        print('in RK2 after first step')
        for ii in range(4,10):
            print(ii,'rke1V=',rke['k_1V'][ii],'rke1Q=',rke['k_1Q'][ii])
        print()
        
        rke = self.element_values( rke, geo )
        rkf = self.face_values( rke, rkf, fhat, geo, fge, Q_Upstream, h_dw )
        
        print('in RK2 after face values')
        for ii in range(4,10):
            print(ii,'Qup',rke['Q'][ii-1],'Qf=',rkf['Q_F'][ii], 'Qdn=',rke['Q'][ii])
        print()

        for ii in range(4,10):
            print(ii,'Aup',rke['A'][ii-1],'Af=',rkf['A_F'][ii], 'Adn=',rke['A'][ii])
        print()
       
        # friction update
        rke['Fdn'] = self.get_friction_half_element(rke['Fdn'], rke, rkf, +1 )
        rke['Fup'] = self.get_friction_half_element(rke['Fup'], rke, rkf,  0 )    
        rke = self.RK2_step2( rke, rkf, el, geo )   
            
        for ii in range(0,N_Cells):
            el['V'][ii] = rke['V'][ii]
            el['Q'][ii] = rke['Q'][ii] 
        
        return el
    
   
    def RK2_step1(self, rke,  el, fa, geo):
        ''' first step of an RK2 time adavance'''
        for ii in range(0,N_Cells): 
            rke['k_1V'][ii] = self.continuity_RHS_explicit( fa['Q_F'], ii) 
                        
            rke['k_1Q'][ii] = self.momentum_RHS_explicit( el, fa, ii)
                
        # Solution at n+1/2
        rke['V'][:]   = el['V'][:] + 0.5 * DT * rke['k_1V'][:] 
        rke['Q'][:]   = el['Q'][:] + 0.5 * DT * rke['k_1Q'][:] / geo['L'][:]  

        for ii in range(0,N_Cells): 
            rke['V'][ii] = min( rke['V'][ii], volume_min[ii] )

        return rke

    def RK2_step2(self,rke, rkf, el, geo):
        ''' second step of an RK2 time advance'''
        for ii in range(0,N_Cells):
            rke['k_2V'][ii]  = self.continuity_RHS_explicit( rkf['Q_F'], ii)
            
            rke['k_2Q'][ii]  = self.momentum_RHS_explicit( rke, rkf, ii)

        # Solution at n+1
        rke['V'][:]   = el['V'][:] + DT * rke['k_2V'][:]
        rke['Q'][:]   = el['Q'][:] + DT * rke['k_2Q'][:] / geo['L'][:] 

        for ii in range(0,N_Cells): 
            rke['V'][ii] = min( rke['V'][ii], volume_min[ii] )
                       
        return rke
        
    def forward_euler(self, el, fa, geo):
        print('BRH 20171123 error not tested')
        sys.exit()
        for ii in range(0,N_Cells):
            el['V'][ii] = el['V'][ii] + DT*( fa['Q_F'][ii] - fa['Q_F'][ii+1] )      
                
            el['Q'][ii] = el['Q'][ii] + ( DT / geo['L'][ii] )         \
                * (  fa['Q_F'][ii]   * fa['U_F'][ii]                  \
                   - fa['Q_F'][ii+1] * fa['U_F'][ii+1]                \
                   + Gravity * ( fa['A_F'][ii]   * fa['Eta_F'][ii]    \
                                -fa['A_F'][ii+1] * fa['Eta_F'][ii+1] )\
                   - el['Fdn'][ii] - el['Fup'][ii] )                            
        return el
    
#=============================================================================
    def solve(self,argv):
        import numpy as np
        import sys

        import Initialization_Class as Initial
        import Visualization_Class  as Visual
        
        
        global el, fa, fhat, geo, fge
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

        Draw = Visual.Visualization()
        Ex = Initial.Initialization(argv)

        ########################################################################
        ########################################################################
        ########################################################################
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
 
        ########################################################################
        ########################################################################
        ########################################################################


        plotstart = 0 # time step to start the plots
        printout = 1 # iteration for printouts
        
        # control for correcting face = 0 for no correcting
        face_correct = 0
        
        # control face interpolation 
        # 1 = 1st order upwind
        # 2 = 2nd order linear
        face_interp = 1
        
        
        
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
        HL  = np.zeros(N_Cells,     dtype=np.float64 ) # element length
        Z    = np.zeros(N_Cells,     dtype=np.float64 ) # element bottom elevation
        M    = np.zeros(N_Cells,     dtype=np.float64 ) # element Manning's n
        B    = np.zeros(N_Cells,     dtype=np.float64 ) # element bottom breadth (square channel)
        X    = np.zeros(N_Cells,     dtype=np.float64 ) # element X distance
        
        geo = {'L':L, 'HL':HL,'Z':Z, 'M':M, 'B':B, 'X':X}
               
       # geometry for plotting (double density)
        Z_Fp  = np.zeros(N_Cells*2+1, dtype=np.float64 ) # face bottom elevation
        X_Fp  = np.zeros(N_Cells*2+1, dtype=np.float64 ) # element X distance on face

        # dynamic variables on elements      
        A    = np.zeros(N_Cells, dtype=np.float64 ) # cross-sectional area
        C    = np.zeros(N_Cells, dtype=np.float64 ) #
        CFLb = np.zeros(N_Cells, dtype=np.float64 ) # barotropic CFL
        CFLu = np.zeros(N_Cells, dtype=np.float64 ) # advective CFL
        E    = np.zeros(N_Cells, dtype=np.float64 ) # energy
        Eta  = np.zeros(N_Cells, dtype=np.float64 ) # free surface elevation
        Fup  = np.zeros(N_Cells, dtype=np.float64 ) # Friction on upstream 1/2 of cell
        Fdn  = np.zeros(N_Cells, dtype=np.float64 ) # Friction on downstream 1/2 of cell  
        Fr   = np.zeros(N_Cells, dtype=np.float64 ) # Froude number
        Gamma= np.zeros(N_Cells, dtype=np.float64 ) # Gradient of Top Width with Z
        H    = np.zeros(N_Cells, dtype=np.float64 ) # hydraulic depth
        l_P  = np.zeros(N_Cells, dtype=np.float64 ) # wetted perimeter
        R_h  = np.zeros(N_Cells, dtype=np.float64 ) # hydraulic radius
        Q    = np.zeros(N_Cells, dtype=np.float64 ) # flow rate
        T    = np.zeros(N_Cells, dtype=np.float64 ) # Top Width
        U    = np.zeros(N_Cells, dtype=np.float64 ) # velocity
        V    = np.zeros(N_Cells, dtype=np.float64 ) # volume

        el =  {'A':A, 'C':C, 'CFLb':CFLb, 'CFLu':CFLu, 'E':E, 'Eta':Eta, 
               'Fr':Fr, 'Fdn':Fdn, 'Fup':Fup, 'Gamma':Gamma, 
               'H':H, 'l_P':l_P, 'R_h':R_h, 'Q':Q, 'T':T, 'U':U, 
               'V':V}             
                             
        # geometry on face
        Z_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # cross-sectional area
        B_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # breadth (rectangular channel)
 
        fge = {'Z_F':Z_F, 'B_F':B_F} 
        
        # dynamic variables on face
        A_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # cross-sectional area
        E_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # energy
        Eta_F = np.zeros(N_Cells+1, dtype=np.float64 ) # free surface elevation
        Q_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # flow rate
        U_F   = np.zeros(N_Cells+1, dtype=np.float64 ) # velocity
        
        fa = {'A_F':A_F, 'E_F':E_F, 'Eta_F':Eta_F, 'Q_F':Q_F, 'U_F':U_F}
        
        # estimated variables on face
        A_F_hat     = np.zeros(N_Cells+1, dtype=np.float64 )
        Eta_F_hat   = np.zeros(N_Cells+1, dtype=np.float64 )  
        Q_F_hat_S   = np.zeros(N_Cells+1, dtype=np.float64 )
        Gamma_F_hat = np.zeros(N_Cells+1, dtype=np.float64 )
        H_hat       = np.zeros(N_Cells+1, dtype=np.float64 )
        
        fhat ={'Eta_F_hat':Eta_F_hat, 'A_F_hat':A_F_hat, 
               'Q_F_hat_S':Q_F_hat_S , 'Gamma_F_hat':Gamma_F_hat,
               'H_hat':H_hat}
        
        # dynamic variables on face for RK solution
        A_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        E_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        Eta_F_1 = np.zeros(N_Cells+1, dtype=np.float64 )
        Q_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )
        U_F_1   = np.zeros(N_Cells+1, dtype=np.float64 )

        rkf = {'A_F':A_F_1, 'E_F':E_F_1, 'Eta_F':Eta_F_1, 'Q_F':Q_F_1, 
               'U_F':U_F_1,  }
        
        # RK solution space
        k_1V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_1Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_2V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_2Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_3V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_3Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_4V  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove
        k_4Q  = np.zeros(N_Cells, dtype=np.float64 ) # <modify> remove        

        # dynamic variables for RK solution
        A_1   = np.zeros(N_Cells, dtype=np.float64 )
        C_1   = np.zeros(N_Cells, dtype=np.float64 )
        E_1   = np.zeros(N_Cells, dtype=np.float64 )
        Eta_1 = np.zeros(N_Cells, dtype=np.float64 )
        Fup_1 = np.zeros(N_Cells, dtype=np.float64 ) # Friction on upstream 1/2 of cell
        Fdn_1 = np.zeros(N_Cells, dtype=np.float64 ) # Friction on downstream 1/2 of cell  
        Fr_1  = np.zeros(N_Cells, dtype=np.float64 )
        Gamma_1 = np.zeros(N_Cells, dtype=np.float64 )
        l_P_1 = np.zeros(N_Cells, dtype=np.float64 )
        Q_1   = np.zeros(N_Cells, dtype=np.float64 )
        R_h_1 = np.zeros(N_Cells, dtype=np.float64 )
        U_1   = np.zeros(N_Cells, dtype=np.float64 )
        V_1   = np.zeros(N_Cells, dtype=np.float64 )
       
        rke = {'k_1V':k_1V, 'k_1Q':k_1Q, 'k_2V':k_2V, 'k_2Q':k_2Q,
               'k_3V':k_3V, 'k_3Q':k_3Q, 'k_4V':k_4V, 'k_4Q':k_4Q,
               'A':A_1, 'C':C_1, 'CFLb':CFLb, 'CFLu':CFLu, 'E':E_1, 'Eta':Eta_1, 
               'Fr':Fr_1, 'Fdn':Fdn_1, 'Fup':Fup_1, 'Gamma':Gamma_1,
               'H':H, 'l_P':l_P_1, 'R_h':R_h_1, 'Q':Q_1, 'T':T, 'U':U_1, 
               'V':V_1}

        volume_min = np.zeros(N_Cells, dtype=np.float64)
        face_area_min = np.zeros(N_Cells+1, dtype=np.float64)
        
#==========================================================================
        # Initialization 
        print(" Initialization ... ")
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
        
        #print('BRH TEST: Initializing Q everywhere with upstream BC==========')
        el['Q'][:] = Ex.Q[:]
        #el['Q'][:] = Ex.Q_Up 
        #for ii in range(0,N_Cells):
        #    el['Q'][ii] = Ex.Q_Up 

        #print('BRH TEST: RESETTING L TO UNIFORM =============================') 
        #geo['L'][:] = 0.0125  # <delete>     We already did the same thing in the initialization class.
        #tempH = 0.33 - geo['Z'][:]
        #el['V'][:] = geo['L'][:] * tempH * geo['B'][:]
        
        # for error checking rectangular channel
        volume_min[:] = depth_min * geo['L'][:] * geo['B'][:]
        # rectangular channel
        face_area_min[1:N_Cells] = depth_min *  0.5 *         \
            ( geo['B'][1:N_Cells] + geo['B'][0:N_Cells-1] )
        face_area_min[0]       = face_area_min[1]
        face_area_min[N_Cells] = face_area_min[N_Cells-1]

        # Geometry setup        
        fge = self.face_downstream_bottom_elevation( fge, geo )
        fge = self.face_upstream_bottom_elevation( fge, geo )
        fge = self.face_interior_bottom_elevation( fge, geo )
        fge = self.face_breadth_rectangular_channel( fge, geo )

        print(" Time marching ... ")
        
        for nn in range(N_Steps):
            #if (nn%printout) == 0:
            print("========================== Time step: %d out of %d " % (nn, N_Steps))
                
            # store the inflow boundary for this step     
            Q_Upstream = Ex.Q_Up     #HACK
                    
            # Cell center auxiliary values at time n
            el = self.element_values( el, geo )            
                      
            # face values at time n
            fa = self.face_values(el, fa, geo, fge, Q_Upstream, h_dw)  
                                    
            # friction on center requires face values
            el['Fdn'] = self.get_friction_half_element                       \
                (el['Fdn'], el, fa, +1)
            el['Fup'] = self.get_friction_half_element                       \
                (el['Fup'], el, fa, 0)
 
#            print('after start')
#            for ii in range(955,966):
#                #print(ii,'Vc=',el['V'][ii])
#                print(ii,'Vc=',el['V'][ii],'QFu=' , fa['Q_F'][ii]  , 'DeltaE=',el['E'][ii-1] - el['E'][ii]  )
#                #print(ii,'QFu='  ,fa['Q_F'][ii]  ,'Qc='  ,el['Q'][ii]  ,'QFd=' ,fa['Q_F'][ii+1] )
#                #print(ii,'EtaFu=',fa['Eta_F'][ii],'Etac=',el['Eta'][ii],'Etad=',fa['Eta_F'][ii+1] )
#                print()
              
            if ((nn%Plot_at_Face) == 0 and nn >= plotstart) \
                or (quit_after_plot):
                RealTime = round(nn*DT,5)
                TITLE1 = Ex.Output_Dir + "/" +  "Full__time_" + str(RealTime)+"_s"
                TITLE2 = "Full results at time: " + str(RealTime) +"_s"
                Draw.Plot_Full_Results(ifig, N_Cells, geo['X'], X_Fp, Z_Fp, el['V'], el['Q'], fa['Q_F'],  \
                               el['Eta'], fa['Eta_F'], el['U'], fa['U_F'],    \
                               el['E'], fa['E_F'], el['A'], fa['A_F'],        \
                               TITLE1, TITLE2)
                if quit_after_plot:
                    print('exit after plot due to error')
                    sys.exit()
                    
                ifig=ifig+1
                                
            if time_advance == 'rk2':    
                el = self.RK2(rke, rkf, el, fa, fhat, geo, fge, Q_Upstream, h_dw)  
            elif time_advance == 'rk4':    
                el = self.RK4(el, fa, geo, fge, rke, rkf, Q_Upstream, h_dw)
            elif time_advance == 'forward_euler':
                el = self.forward_euler( el, fa, geo )
            else:
                print('error, unknown time advance of ',time_advance)
                sys.exit()
         
        print(" ========== Solver Class ==========")
        print()
