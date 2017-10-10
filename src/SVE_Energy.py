
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

def main(arg):


    # Import built-in libraries =======================================================================
    # import numpy as np --delete

    # Import classes ==================================================================================
    import sys
    import Solver_Class
    import math  # <delete>
    import os  # <delete>  You can create/del dir using this module
    import time  # time.sleep(2)
    from datetime import datetime  # Use

    # Code begins =====================================================================================
    print()
    print("{:^80}".format("-------- 1D Saint-Venant Equation based on the face reconstruction --------"))
    print("{:^80}".format("---------- Developers:: Babak Poursartip/Ben R. Hodges/Frank Liu ----------"))
    print()
    print("{:^80}".format(" Simulation starts ..."))
    print()

    Results = Solver_Class.Solver()

    print("{:80}".format("---------- Simulation was conducted successfully ----------"))
    print()

if __name__ == '__main__':
    import sys    
    main(sys.argv)





