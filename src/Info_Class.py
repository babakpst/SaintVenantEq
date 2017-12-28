
#######################################################################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    12/6/2017
# Latest update: 12/20/2017
#
# Comment: 
#
#######################################################################################################################


class Info:



    def __init__(self,OutPutFile,name):
        
        import datetime
        print(" This is the information class ")

        OutPutFile.write("\n")
        OutPutFile.write("{0} {1} \n".format(" The model name is:",name))
        OutPutFile.write("{} {} \n".format(" Simulation starts at:",datetime.datetime.now() ))
        OutPutFile.write("\n")
        