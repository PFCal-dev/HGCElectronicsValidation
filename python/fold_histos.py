import ROOT
import sys

from utils import *
from histo_classes import *


# run a test:
# cd   /afs/cern.ch/user/f/franzoni/work/CMSSW_11_2_0_pre3_afterTalk2/src/UserCode/HGCElectronicsValidation/python; esr ; python fold_histos.py /eos/home-f/franzoni/www/CMS/HGCSample/SiOptim/2020-08-21-nofold/ttbar_D49_1120pre1_PU200_eolupdate_qua_20200723_thr0p5_thrbxm12p5_flfalse.root  rrrr.root


if len( sys . argv ) < 2:
    print('USAGE : %s <input file > [ <output file > ]'%( sys . argv [0]))
    sys . exit (1)

inFileName  = sys . argv [1]
if len (sys . argv) == 2:
    outFileName = None
else:
    outFileName = sys . argv [2]
print (" Reading from: %s\n writing to: %s"%(inFileName, outFileName))


# works
inFile = ROOT . TFile . Open ( inFileName ," READ")
inFile.cd("ana;1")
# magic to get all directories
keyList = inFile.GetKeyNames("ana")


# cast the text keys into a dataframe, with additional columns for later operations
keys_df = keys_to_df(keyList)

# print "\ndf head:\n", keys_df.head()
# print "\ndf tail:\n", keys_df.tail()
print('\n\n Keys in file:', len(keyList), keys_df.shape)

# a group is formed by all wafers which will be folded together,
# with one representative element which is the only one in the first sextant (EE) / thirdtant (HE)
print('\n\n Making groups which will be folded together:... ')
keys_df_groups = keys_df.groupby( ['det','lay','U_rot','V_rot'] )
num_groups = len(keys_df_groups)
print('  ... %d  groups formed\n\n'%num_groups)


counter=-1
# loop over groups of wavers; members of a group will be folded together
for name, group in keys_df_groups:

    counter+=1
    if counter%100 ==0:
        print('\n %d groups foled out of %d'%(counter,num_groups))

    if group.size == 0:
        print('group: %s has 0 size... strange; continuing'%name)
        continue

    # first is the wafer within the group (located in the first sextant/thirdtant ) over which the other will be Add()-ed
    first_key_df = group.loc[group['first']==True]['txt_key'].to_numpy()  # extract the actul element
    if first_key_df.size ==0:
        print('group: %s has 0 elements within 1st sextant/third-tant... strange; continuing'%name)
        continue

    first_key = first_key_df[0]
    first_histo = Histos(first_key, inFileName, outFileName)
    first_histo.set_infile_name(inFileName)
    first_histo.get_histos()

    # print('++ main nun histos: %d'%len(first_histo.histos))
    # for key in first_histo.histos:
    #    print('==> LOOP MAIN check_histos histo found called: %s with entries %d'%(first_histo.histos[key].GetName(), first_histo.histos[key].GetEntries()))
    # first_histo.check_histos()

    # loop on other members of the group and Add() their histogram
    # to the instances of the first (representative) wafer
    # loop over (the other) members of the group, which will be added to the first
    other_histos_df = group.loc[group['first']==False]['txt_key']
    if other_histos_df.size ==0:
        print('group: %s has 0 elements outside 1st sextant/third-tant... strange; continuing'%name)
        continue

    for other_key in other_histos_df:
        # print('other_key is: %s'%other_key)
        other_histo = Histos(other_key, inFileName, )
        other_histo.get_histos()
        first_histo.add_histo(other_histo)

    first_histo.write_histos()

# outFile.Close()

# logical plan

# get all directories present in ana

# map wafers of first sextant to wafers outsde of it

## meeds a tough class handling all this
# find a way of loading all histograms within one folder/WAFER into a class

# find way of adding hroups = len(keys_df_groups.keys)istograms of two wafers

# write everytying out
