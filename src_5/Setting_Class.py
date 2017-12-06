

###############################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    11/27/2017
# Latest update: 11/27/2017
#
# Comment: This class defines the various setting required for this code.
#
###############################################################################

class Setting:

    def __init__(self):

        self.DT_min = 1.0e-6
        self.DT_max = 1

        self.Plot_at_Cell = 5 # Plots the results at the cell center, every "Plot_at_Cell" steps.
        self.Plot_at_Face = 10 # Plots the results at the cell center, every "Plot_at_Face" steps.

        self.CFLu_max = 1.0
        self.CFLb_max = 1.0

        self.Fr_max = 10000.0

        self.depth_min = 0.001

        self.time_advance = 'rk4' # Options are: rk4 rk2 forward_euler

        self.Gravity = 9.81