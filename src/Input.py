
#####################################################################
# 
# This class reads the simulation data from file
# Start date:    07/18/2017
# Latest update: 07/25/2017
#
#####################################################################


class Input_Class:

  # -- Class initialization 
  def __init__(self, InputFileName):
    print("Reading data from the input file ...")
    self.InputFileName = InputFileName
    print()

  # -- Opens and Reads data from the input file
  def Input(self, ):
    File_Input = open(self.InputFileName,"r")
    print("The input file opened successfully.")

    



