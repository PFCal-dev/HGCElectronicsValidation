import ROOT
import sys

from utils import *


if len( sys . argv ) != 3:
    print('USAGE : %s <input file > <output file >'%( sys . argv [0]))
    sys . exit (1)

inFileName  = sys . argv [1]
outFileName = sys . argv [2]
print (" Reading from ", inFileName , "and writing to", outFileName)


# works
inFile = ROOT . TFile . Open ( inFileName ," READ ")
inFile.cd("ana;1")


# magic to get all directories
keyList = inFile.GetKeyNames("ana")
# print "\nKeys in file:", keyList



outFile = ROOT.TFile.Open(outFileName, 'recreate')
run_demo( outFile )
outFile.Close()

# get all directories present in ana

# map wafers of first sextant to wafers outsde of it

## meeds a tough class handling all this
# find a way of loading all histograms within one folder/WAFER into a class

# find way of adding histograms of two wafers

# write everytying out
