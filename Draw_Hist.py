import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--filename', type=str, required=True)
parser.add_argument('--histname', type=str, required=True)
args = parser.parse_args()

class HistogramFile(object):
    def __init__(self, filename):
        self.filename = filename
    def __enter__(self):
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

with HistogramFile(args.filename) as f:
    canvas = ROOT.TCanvas('canvas', '', 500, 500)
    hist = f.get_histogram(args.histname)
    hist.Draw()
    ROOT.gPad.Update()
    raw_input('Please press enter to continue.')
