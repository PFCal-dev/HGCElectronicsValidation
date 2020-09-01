import ROOT
from   ROOT import TFile
from   ROOT import gDirectory
import sys


def GetKeyNames( self, dir = "" ):
    '''some magic to get all directories within dir'''
    self.cd(dir)
    return [key.GetName() for key in gDirectory.GetListOfKeys()]
TFile.GetKeyNames = GetKeyNames



class Histos(object):
    '''
    holds hisotograms of a director and 
    carries out basic operations
    '''
    def __init__(self, key):
        self.name   = key
        self.histos = {}
        
    def set_infile_name(self,n):
        self.filename = n
        print('\n\nHitos: filename set to %s'%self.filename)

    def get_histos(self):
        print('Histo class named: %s'%self.name)
        histo_names = ['busycounts','adc']
        # some root magic to make sure that cloning persists the histos in the dictionary
        # https://root-forum.cern.ch/t/pyroot-typecast-histograms-stored-in-dict/24744/11
        ROOT.TH1.AddDirectory(0) 
        with HistogramFile( self.filename ) as f:
            for histo_name in  histo_names:
                nn = 'ana/' + self.name +'/' + self.name + '_' + histo_name
                print('** st')
                print('looking for histo named: %s'%nn)
                # print(nn)
                one_histo = f.get_histogram( nn )
                # one_histo = f.get_histogram( 'ana/0_lay1_2_0/0_lay1_2_0_busycounts' )
                print('histo found called: %s with entries %d'%(one_histo.GetName(), one_histo.GetEntries()))
                # put it in dict
                self.histos[histo_name] = one_histo.Clone()
                print('==> RIGHT AFTER histo found called: %s with entries %d and key %s'%(self.histos[histo_name].GetName(), self.histos[histo_name].GetEntries(), histo_name) )
                print('en **')

    def check_histos(self):
        print('check_histos showing the dic')
        print(self.histos)
        for key in self.histos:
            print('==> LOOP LOOP check_histos histo found called: %s with entries %d'%(self.histos[key].GetName(), self.histos[key].GetEntries()))
        print('==> OUT OUT check_histos histo found called: %s with entries %d'%(self.histos['adc'].GetName(), self.histos['adc'].GetEntries()))

class HistogramFile(object):
    '''
    basic interface to a root file
    '''
    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        '''http://effbot.org/zone/python-with-statement.htm'''
        self.file = ROOT.TFile.Open(self.filename, 'read')
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.file.Close()

    def get_histogram(self, name):
        """Return the histogram identified by name from the file.
        """
        # The TFile::Get() method returns a pointer to an object stored in a ROOT file.
        hist = self.file.Get(name)
        if hist:
            return hist
        else:
            raise RuntimeError('Unable to retrieve histogram named {0} from {1}'.format(name, self.filename))
