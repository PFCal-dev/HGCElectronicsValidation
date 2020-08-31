import ROOT
from   ROOT import TFile
from   ROOT import gDirectory
import sys


def GetKeyNames( self, dir = "" ):
    '''some magic to get all directories within dir'''
    self.cd(dir)
    return [key.GetName() for key in gDirectory.GetListOfKeys()]
TFile.GetKeyNames = GetKeyNames



class Hitos(object):
    '''
    holds hisotograms of a director and 
    carries out basic operations
    '''
    def __init__(self, key):
        self.name = key
        
    def set_infile_name(self,n):
        self.filename = n
        print('Hitos: filename set to %s'%self.filename)

    def get_histos(self):
        histos = {}
        histo_names = ['busycounts','adc']
        with HistogramFile( self.filename ) as f:
            for histo_name in  histo_names:
                nn = 'ana/' + self.name +'/' + self.name + '_' + histo_name
                print('looking for histo named: %s'%nn)
                print('**')
                print(self.name)
                print('-')
                print(nn)
                print('**')
                one_histo = f.get_histogram( nn )
                # one_histo = f.get_histogram( 'ana/0_lay1_2_0/0_lay1_2_0_busycounts' )
                print('histo found is: %s with entries %d'%(one_histo.GetName(), one_histo.GetEntries()))
                # put it in dict


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
