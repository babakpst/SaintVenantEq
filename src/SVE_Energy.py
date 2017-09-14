
#!/usr/bin/python

# Tell the linux machine to run the script with python

#####################################################################
# 
# This code is based on the methodology described in the paper.
# The code is based on Python 3.6, NumPy 1.12, and MatLibPlot.
#
# 
# 
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 08/08/2017
#
#####################################################################

def main():


    # Import built-in libraries =======================================================================
    # import numpy as np --delete

    # Import classes ==================================================================================
    import Solver_Class


    # Code begins =====================================================================================
    print()
    print('''
    
     -------- 1D Saint-Venant Equation based on the face reconstruction --------
     ---------- Developers:: Babak Poursartip/Ben R. Hodges/Frank Liu ----------
    
     ''')
    print()
    print()
    print(" Simulation starts ...")
    print()

    Results = Solver_Class.Solver()
    #Results2=Results.RK2()

    print(" Simulation was conducted successfully. ")
    print()

if __name__ == '__main__':
     main()





