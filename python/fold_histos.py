import ROOT
import sys

from utils import *
from histo_classes import *


# run a test:
# cd   /afs/cern.ch/user/f/franzoni/work/CMSSW_11_2_0_pre3_afterTalk2/src/UserCode/HGCElectronicsValidation/python; esr ; python fold_histos.py /eos/home-f/franzoni/www/CMS/HGCSample/SiOptim/2020-08-21-nofold/ttbar_D49_1120pre1_PU200_eolupdate_qua_20200723_thr0p5_thrbxm12p5_flfalse.root  rrrr.root


if len( sys . argv ) != 3:
    print('USAGE : %s <input file > <output file >'%( sys . argv [0]))
    sys . exit (1)

inFileName  = sys . argv [1]
outFileName = sys . argv [2]
print (" Reading from ", inFileName , "and writing to", outFileName)



# works
inFile = ROOT . TFile . Open ( inFileName ," READ")
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
print( keys_df_groups.get_group( (0, 11, 6, 0) ))
print( type(keys_df_groups.get_group( (0, 11, 6, 0) ) ) )

print( keys_df_groups.get_group( (1, 5, 8, 7)  ))
print()
print( keys_df_groups.get_group( (1, 5, 8, 7)  )['txt_key'])
print( type(keys_df_groups.get_group( (1,5, 8, 7)  )['txt_key']) )

#keys_df_groups_kk = keys_df_groups.keys
print( len(keys_df_groups.keys)  )

print()
print()

outFile = ROOT.TFile.Open(outFileName, 'recreate')
run_demo( outFile )



counter=0

# loop over groups of wavers; members of a group will be folded together
for name, group in keys_df_groups:
    if counter==3:
        break
    counter+=1

    first_key = group.loc[group['first']==True]['txt_key'].to_numpy()[0]
    # first_key = ( group.loc[group['first']==True] ).at[0,'txt_key']
#    print('------|||')
#    print( type(group) )
#    print()
#    print( type(group.loc[group['first']==True])  )
#    print( group.loc[group['first']==True] )
#    print( group.loc[group['first']==True]['txt_key'].to_numpy()[0] )
#    print('------|||')
#    print('name is:') # GF first key needs work
#    print(name) # GF first key needs work
#    print('key value is: %s'%(first_key)) # GF first key needs work
#    print(group)
#    print("\n")
#
    print()
    first_histos = Histos(first_key)
    first_histos.set_infile_name(inFileName)
    first_histos.get_histos()

    print('++ main nun histos: %d'%len(first_histos.histos))
#    print('++ main histo found called: %s with entries %d'%(first_histos.histos['adc'].GetName(), first_histos.histos['adc'].GetEntries()))
#     print('++ main histo found called: %s with entries %d'%(first_histos.histos['busycounts'].GetName(), first_histos.histos['busycounts'].GetEntries()))
    print('++ main nun histos: %d'%len(first_histos.histos))

    first_histos.check_histos()


outFile.Close()

# logical plan

# get all directories present in ana

# map wafers of first sextant to wafers outsde of it

## meeds a tough class handling all this
# find a way of loading all histograms within one folder/WAFER into a class

# find way of adding histograms of two wafers

# write everytying out
