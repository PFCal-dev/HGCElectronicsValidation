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
print "\nKeys in file:", len(keyList)

def split_key(ll):
    return map(lambda s: s.split('_'), ll)


import pandas as pd
def keys_to_df(lol):
    uu = pd.DataFrame(lol, columns=['det','lay','U','V'])
    # drop the last 6 rows which are spurious text
    #print( uu.shape )
    #print( uu.tail(5))
    #print( uu.tail(5).index)
    # print( uu.drop(uu.tail(5).index,inplace=True).shape )
    # print( uu.drop(uu.tail(5).index,inplace=True).shape )
    return uu[:-6]
    # return  uu.drop(uu.tail(5).index,inplace=True)

# print "\nKeys in file:", split_key(keyList)
keys_df = keys_to_df( split_key(keyList) )

print "\ndf head:", keys_df.head()
print()
print "\ndf tail:", keys_df.tail()

# print "\n0 Keys in file:", keyList[0]
# print "\n110 Keys in file:", keyList[110]
# print "\n-1 Keys in file:", keyList[-1]



outFile = ROOT.TFile.Open(outFileName, 'recreate')
run_demo( outFile )
outFile.Close()

# get all directories present in ana

# map wafers of first sextant to wafers outsde of it

## meeds a tough class handling all this
# find a way of loading all histograms within one folder/WAFER into a class

# find way of adding histograms of two wafers

# write everytying out
