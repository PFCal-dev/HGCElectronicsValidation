import ROOT
from   ROOT import TFile
from   ROOT import gDirectory
import sys



def GetKeyNames( self, dir = "" ):
    '''some magic to get all directories within dir'''
    self.cd(dir)
    return [key.GetName() for key in gDirectory.GetListOfKeys()]
TFile.GetKeyNames = GetKeyNames




########################################################
class Histos(object):
    '''
    holds hisotograms of a directory/wafer and 
    carries out basic operations: reading from inflie, adding, writing to outfile
    '''

    def __init__(self, key, in_file, out_file=None):
        self.name        = key
        self.histos      = {}
        self.histo_types = ['busycounts','adc']
        self.in_file     = in_file
        self.out_file    = out_file
        self.consolidate_filename()


    def set_histo_types(self,ll):
        assert type(ll)  == list
        self.histo_types = ll 


    def consolidate_filename(self):
        '''
        if filename not specified in constructor
        set it to ./ plus the input filename interposig the suffix 'folded'
        '''
        if self.out_file == None:
            tmp1 = self.in_file.split('/')[-1].split('.')
            tmp2 = tmp1[0:-1] + ['folded'] + tmp1[-1:]
            tmp3 = './' + '.'.join(tmp2)
            self.out_file = tmp3

    def set_infile_name(self,n):
        self.in_file = n
        # print('\n\nHitos: in_file set to %s'%self.in_file)


    def get_histos(self):
        # print('Histo class named: %s - loading histos '%self.name)
        # some root magic to make sure that cloning persists the histos in the dictionary
        # https://root-forum.cern.ch/t/pyroot-typecast-histograms-stored-in-dict/24744/11
        ROOT.TH1.AddDirectory(0) 

        # fill a dictionary with histogram objects, then re-close the input file
        with HistogramFile( self.in_file, rw='read') as f:
            for histo_type in  self.histo_types:
                nn = 'ana/' + self.name +'/' + self.name + '_' + histo_type
                one_histo = f.get_histogram( nn )
                # put it in dict
                self.histos[histo_type] = one_histo.Clone()


    def write_histos(self):
        # print('\n\nentering TOP')
        dir = 'ana/' + self.name +'/'
        # it would be better not to open and close the outout file 
        # for every sensor group... next project ;)
        # print('BEF looping TOP')
        with HistogramFile( self.out_file, rw='update' ) as f:
            for histo_type in  self.histo_types:
                # print('\t TOP loop with: %s'%histo_type)
                f.write_histogram( dir, self.histos[histo_type] )
                # print('\t TOP after write_histogram -  with: %s'%histo_type)

    
    def check_histos(self):
        '''
        was useful during development to see if root histograms properly dealt with in memory
        '''
        print('++ check_histos for object %s showing the dic'%self.name)
        print(self.histos)
        for key in self.histos:
            print('==> check_histos: LOOP LOOP check_histos histo found called: %s with entries %d'%(self.histos[key].GetName(), self.histos[key].GetEntries()))


    def add_histo(self, other):
        '''
        receives anotehr object of this same class
        and add its histograms to those of the current object
        type-by-type
        '''
        # print('accessing other %s from self %s'%(other.name,self.name) )

        for histo_type in self.histo_types:

            one_bins   = self.histos[histo_type].GetNbinsX()
            other_bins = other.histos[histo_type].GetNbinsX()
            if (one_bins != other_bins):
                continue # we don't worry about the partials which have different number of bins

            self.histos[histo_type].Add(other.histos[histo_type])







########################################################
class HistogramFile(object):
    '''
    basic interface to a root file:
    lets you both to read and to write to a file
    '''
    def __init__(self, filename, rw='read'):
        self.filename = filename
        self.rw       = rw


    def __enter__(self):
        '''http://effbot.org/zone/python-with-statement.htm'''
        self.file = ROOT.TFile.Open(self.filename, self.rw)
        return self


    def __exit__(self, exception_type, exception_value, traceback):
        # if (self.rw in ['recreate','RECREATE','write',] ):
            # self.file.Write()
        self.file.Close()


    def get_histogram(self, name):
        """
        Return the histogram identified by name from the file.
        """

        hist = self.file.Get(name)
        if hist:
            return hist
        else:
            raise RuntimeError('Unable to retrieve histogram named {0} from {1}'.format(name, self.filename))


    def write_histogram(self, dir, histo):
        """
        write the histogram identified by name to the file.
        """
        assert self.rw in ['update']#'recreate','RECREATE','write',]
        # print('entering INNER')
        self.file.cd()
        self.file.mkdir(dir)
        # print('AFT mkdir INNER')
        self.file.cd(dir)
        # print('AFT cd INNER')
        histo.Write()
        # self.file.Write()
        # print('AFT write INNER')
        
