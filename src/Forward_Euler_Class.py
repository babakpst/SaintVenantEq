

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

class ForwardEuler:


    # This function has not been debugged yet.
    def forward_euler(self, el, Face_Arrays, geo):
        print('BRH 20171123 error not tested')
        sys.exit()

        for ii in range(0,N_Cells):
            el['V'][ii] = el['V'][ii] + DT*( Face_Arrays['Q_F'][ii] - Face_Arrays['Q_F'][ii+1] )      
                
            el['Q'][ii] = el['Q'][ii] + ( DT / geo['L'][ii] )         \
                * (  Face_Arrays['Q_F'][ii]   * Face_Arrays['U_F'][ii]                  \
                   - Face_Arrays['Q_F'][ii+1] * Face_Arrays['U_F'][ii+1]                \
                   + Gravity * ( Face_Arrays['A_F'][ii]   * Face_Arrays['Eta_F'][ii]    \
                                -Face_Arrays['A_F'][ii+1] * Face_Arrays['Eta_F'][ii+1] )\
                   - el['Fdn'][ii] - el['Fup'][ii] )                            
        return el

