import ROOT
from   ROOT import TFile
from   ROOT import gDirectory
import sys

# this include can go away once we'llre move run_demo
from histo_classes import *

# this below is a lottle demo, lots of crap / needs to be exported and organised
#      ana/0_lay1_2_0/0_lay1_2_0_busycounts
hisN='ana/0_lay1_-11_-6/0_lay1_-11_-6_busycounts'
hisNtwo='ana/0_lay1_-11_-5/0_lay1_-11_-5_busycounts'
def run_demo(outFile, 
             # inFileName = '/data/franzoni/ttbar_D49_1120pre1_PU200_eolupdate_qua_20200723_thr0p5_thrbxm12p5_flfalse.root'
             inFileName = '/eos/home-f/franzoni/www/CMS/HGCSample/SiOptim/2020-08-21-nofold/ttbar_D49_1120pre1_PU200_eolupdate_qua_20200723_thr0p5_thrbxm12p5_flfalse.root'
         ):
    with HistogramFile( inFileName ) as f:
        hist_1 = f.get_histogram(hisN)
        hist_2 = f.get_histogram(hisNtwo)
        # print(hist_1.GetName())
        # print(hist_1.GetEntries())
        # how to write a hitogram to ad different file than it's been red from
        outFile.cd()
        outFile.mkdir('ana/0_lay1_-11_-6/')
        outFile.cd('ana/0_lay1_-11_-6/')
        hist_1.Write()
        outFile.mkdir('ana/0_lay1_-11_-5/')
        outFile.cd('ana/0_lay1_-11_-5/')
        hist_2.Write()
        hist_1.Add(hist_2)
        hist_3=hist_1.Clone()
        hist_3.SetName('puppa')
        outFile.mkdir('ana/SUM/')
        outFile.cd('ana/SUM/')
        hist_3.Write()


#
def rot(waferU, waferV):
    ''' rotate a (U,V) 
    pair by 60 degrees, clock-wise
     rotation by 60 degrees (skype on April 27, 2020 - 14:48)
        U' = u -v
        V' = u
    '''

    u_old = waferU
    v_old = waferV

    waferU = u_old - v_old;
    waferV = u_old;
    
    return waferU,waferV

#
def rotat(subdet, waferU, waferV):
    return rotate(subdet, waferU, waferV)[0]
        
def rotate(subdet, waferU, waferV):
    ''' 
    chose wether to rotate 
    by one sextant (if EE, 60 deg)
    or by one thirdtant (if HE, 120 deg)
    '''
    assert subdet in [0,1]
    
    if   subdet==0:    # 60 degrees in EE
        waferU,waferV = rot(waferU,waferV)
    elif subdet==1: # 120 = x2 60 degrees in HE
        waferU,waferV = rot(waferU,waferV)
        waferU,waferV = rot(waferU,waferV)
    else:
        pass

    return waferU,waferV


#
def isFirstSextant(waferU, waferV):
    ''' true if wafer in first sextant'''
    # condition to be in the first sextant
    if (waferV >= 0 and waferU >waferV) :
        return True
    return False


#
def isFirstThirdtant(waferU, waferV):
    ''' true if wafer in first thirdtant'''

    if isFirstSextant(waferU, waferV):
        return True

    waferU, waferV  = rot(waferU, waferV)
    # if wafer is in 1st sextant after rotating by 60deg => it's in 1st thirdant
    if isFirstSextant(waferU, waferV):
        return True

    return False

#
def isFirstSexORThirdTant(subdet,U,V):
    if   subdet==0:
        return isFirstSextant(U,V)
    elif subdet==1:
        return isFirstThirdtant(U,V)
    else:
        pass

#
def remapUV(subdet, waferU, waferV):
    '''
    rotate an UV pair far enogh
    to map it into the first sextant/thirdtant (in EE/EH, resp.)
    '''
    assert subdet in [0,1]

    if      subdet==0:    # EE: rotate to first sextant
        while not isFirstSextant(waferU, waferV):
            waferU,waferV = rotate(subdet,waferU,waferV)
    
    elif subdet==1:     # //HE: rotate to first third-tant
        while not isFirstThirdtant(waferU, waferV):
            waferU,waferV = rotate(subdet,waferU,waferV)

    else:
        pass

    return waferU,waferV



def split_key(ll):
    '''
    split the root txt keys into indices
    '''
    return map(lambda s: s.split('_'), ll)


import pandas as pd
def keys_to_df_raw(lol):
    '''
    turn list of list into a df, each row is unique ID of a wafer
    remove the last 6 raws which are garbage
    '''
    uu = pd.DataFrame(lol, columns=['det','lay_txt','U','V'])
    return uu[:-6]

#
def make_txt_key(det,lay,U,V):
    ''' 0_lay1_-11_-6 '''
    return str(det)+'_lay'+str(lay)+'_'+str(U)+'_'+str(V)


def keys_to_df(keyList):
    vv = keys_to_df_raw( split_key(keyList) )

    # force the indices to be integers (natively they're object type)
    vv = vv.astype({'U': 'int32'})
    vv = vv.astype({'V': 'int32'})
    vv = vv.astype({'det': 'int32'})

    # create a numeric column for layer, which is natively a string like 'lay22'
    vv['lay'] = vv['lay_txt'].str.replace('lay','')
    vv = vv.astype({'lay': 'int32'})
    
    # add the remapped UV, which will be used for grouping
    vv['U_rot'] = vv.apply( lambda row: remapUV(row['det'],row['U'],row['V'])[0], axis=1)
    vv['V_rot'] = vv.apply( lambda row: remapUV(row['det'],row['U'],row['V'])[1], axis=1)

    vv['txt_key'] = vv.apply( lambda row: make_txt_key(row['det'],row['lay'],row['U'],row['V']), axis=1)

    # mark the wafers which are already inside the folding region
    vv['first'] = vv.apply( lambda row: isFirstSexORThirdTant(row['det'],row['U'],row['V']), axis=1)

    return vv

def print_groups(keys_df_groups):
    '''
    print all groups identified by previous groupby
    '''
    for key, item in keys_df_groups:
        print(keys_df_groups.get_group(key), "\n\n")
