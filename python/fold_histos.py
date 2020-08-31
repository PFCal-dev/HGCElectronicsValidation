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
print("\nexample of one key:", keyList[0])
print("example of one key:", keyList[-1],'\n')


# print "\nKeys in file:", split_key(keyList)
keys_df = keys_to_df(keyList)

print "\ndf head:\n", keys_df.head()
print "\ndf tail:\n", keys_df.tail()
print('\n\n Keys in file:', len(keyList), keys_df.shape)

# print('\n\n types',keys_df.dtypes)

print('\n\n Making groups: \n\n')
keys_df_groups = keys_df.groupby( ['det','lay','U_rot','V_rot'] )

# see the groups
# print_groups(keys_df_groups)

# get all key sets defining one group
# print(keys_df_groups.groups)
print( keys_df_groups.get_group( (0, '11', 6, 0) ))
print( type(keys_df_groups.get_group( (0, '11', 6, 0) ) ) )

print( keys_df_groups.get_group( (1, '5', 8, 7)  ))
print( keys_df_groups.get_group( (1, '5', 8, 7)  )['U'])

# k,i = keys_df_groups.get_group(0)
# print(keys_df_groups.get_group(k))
# print( type(keys_df_groups.get_group(k)) )



outFile = ROOT.TFile.Open(outFileName, 'recreate')
run_demo( outFile )
outFile.Close()

# get all directories present in ana

# map wafers of first sextant to wafers outsde of it

## meeds a tough class handling all this
# find a way of loading all histograms within one folder/WAFER into a class

# find way of adding histograms of two wafers

# write everytying out
