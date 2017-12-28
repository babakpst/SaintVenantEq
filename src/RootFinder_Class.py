
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    12/08/2017
# Latest update: 12/27/2017
#
# Comment: This class conducts the numerical solution of the Saint-Venant equation
# 
#
#######################################################################################################################

class Root_Finder:

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
        import numpy as np
        
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
