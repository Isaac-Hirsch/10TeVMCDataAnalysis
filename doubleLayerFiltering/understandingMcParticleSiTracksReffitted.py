#Imports are meant to be run on snowmass21 server
from array import array
from pyLCIO import IOIMPL, EVENT, UTIL
import ROOT
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath
import glob
from optparse import OptionParser

#Code needed to download data to a root file
Bfield = 3.56  # T
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_hits_SiTracksNOBIB.root',
                  type=str, default='ntup_hits_SiTracksNOBIB.root')
(options, args) = parser.parse_args()

tree = TTree("tracks_tree", "tracks_tree")

fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_0_50/muonGun_pT_0_50_reco_2?00.slcio")

for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    for event in reader:
        MCPSiRe = event.getCollection('MCParticle_SiTracks_Refitted')
        for i in MCPSiRe:
            print(dir(i))
        
    