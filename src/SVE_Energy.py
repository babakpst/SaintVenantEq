

#####################################################################
# 
# This code is based on the methodology described in the paper "Write the final title and the Journal name"
# The code is based on Python 3.6, NumPy 1.12, and MatLibPlot.
#
# 
# 
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 11/20/2017
#
#####################################################################

def main(argv):

    # Import built-in libraries =======================================================================
    # import numpy as np --delete

    # Import classes ==================================================================================
    import sys
    import Solver_Class
    import math  
    import os  
    import time  
    from datetime import datetime

    # Code begins =====================================================================================
    print()
    print("{:^80}".format("-------- 1D Saint-Venant Equation based on the face reconstruction --------"))
    print("{:^80}".format("---------- Developers:: Babak Poursartip/Ben R. Hodges/Frank Liu ----------"))
    print()
    print("{:^80}".format(" Simulation starts ..."))
    print()

    Results = Solver_Class.Solver()
    Results.solve()

    print("{:^80}".format("---------- Simulation was conducted successfully ----------"))
    print()

if __name__ == '__main__':
    import sys    
    main(sys.argv)





