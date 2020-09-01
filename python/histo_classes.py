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
        self.name        = key
        self.histos      = {}
        self.histo_types = ['busycounts','adc']
        self.out_file    = None

    def set_infile_name(self,n):
        self.filename = n
        print('\n\nHitos: filename set to %s'%self.filename)


    def get_histos(self):
        print('Histo class named: %s - loading histos '%self.name)
        # some root magic to make sure that cloning persists the histos in the dictionary
        # https://root-forum.cern.ch/t/pyroot-typecast-histograms-stored-in-dict/24744/11
        ROOT.TH1.AddDirectory(0) 

        # fill a dictionary with histogram objects, then re-close the input file
        with HistogramFile( self.filename ) as f:
            for histo_type in  self.histo_types:
                nn = 'ana/' + self.name +'/' + self.name + '_' + histo_type
                one_histo = f.get_histogram( nn )
                # put it in dict
                self.histos[histo_type] = one_histo.Clone()

    def write_histos(self):
        print('Histo class named: %s - writing histos '%self.name)
        # TAKE CARE OF UI for filename
        xx = '/afs/cern.ch/user/f/franzoni/work/CMSSW_11_2_0_pre3_afterTalk2/src/UserCode/HGCElectronicsValidation/python/out.root'
        # fill a dictionary with histogram objects, then re-close the input file
        with HistogramFile( xx, rw='recreate' ) as f:
            dir = 'ana/' + self.name +'/'
            print('write_histos: about to loop')
            for histo_type in  self.histo_types:
                # nn = 'ana/' + self.name +'/' + self.name + '_' + histo_type
                f.write_histogram( dir, self.histos[histo_type] )




    def check_histos(self):
        print('++ check_histos for object %s showing the dic'%self.name)
        print(self.histos)
        for key in self.histos:
            print('==> check_histos: LOOP LOOP check_histos histo found called: %s with entries %d'%(self.histos[key].GetName(), self.histos[key].GetEntries()))
        #print('==> check_histos: OUT OUT check_histos histo found called: %s with entries %d'%(self.histos['adc'].GetName(), self.histos['adc'].GetEntries()))

    #def histo_types():
    #    return self.histo_types

    def add_histo(self, other):
        print('accessing other %s from self %s'%(other.name,self.name) )
        # for histo_type in  other.histo_types:
        #    print('==> add_histo: from other %s histo found called: %s with entries %d'%(other.name,other.histos[histo_type].GetName(), other.histos[histo_type].GetEntries()))

        for histo_type in self.histo_types:
            print('==> add_histo BEF: I am class %s, histo called %s with entries %d'%(self.name,self.histos[histo_type].GetName(), self.histos[histo_type].GetEntries()))
            self.histos[histo_type].Add(other.histos[histo_type])
            print('==> add_histo AFT: I am class %s, histo called %s with entries %d'%(self.name,self.histos[histo_type].GetName(), self.histos[histo_type].GetEntries()))



class HistogramFile(object):
    '''
    basic interface to a root file
    '''
    def __init__(self, filename, rw='read'):
        self.filename = filename
        self.rw       = rw

    def __enter__(self):
        '''http://effbot.org/zone/python-with-statement.htm'''
        self.file = ROOT.TFile.Open(self.filename, self.rw)
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

    def write_histogram(self, dir, histo):
        """write the histogram identified by name to the file.
        """
        self.file.cd()
        # self.file.mkdir('ana/0_lay1_-11_-6/')
        self.file.mkdir(dir)
        TFdir = self.file.GetDirectory(dir);
        # GFGF this needs work
        histo.SetDirectory(TFdir)
        histo.Write()
