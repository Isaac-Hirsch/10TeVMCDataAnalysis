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
parser.add_option('-o', '--outFile', help='--outFile ntup_hits.root',
                  type=str, default='ntup_hits.root')
(options, args) = parser.parse_args()

tree = TTree("tracks_tree", "tracks_tree")

#Initialize the branches for the info to be stored in:
x_pos = array('d', [0])
y_pos = array('d', [0])
z_pos = array('d', [0])
time = array ('d', [0])
barOrEnd = array('d', [0]) #0 if in the barrel, 1 if in the endcap
location = array('d',[0]) #0 in vertex, 1 in inner layer, 2 in outer layer
system = array ('d', [0])
layer = array ('d', [0])
side = array('d', [0])

#Initialzing branches for each layer
tree.Branch("x",  x_pos,  'var/D')
tree.Branch("y",  y_pos,  'var/D')
tree.Branch("z", z_pos, 'var/D')
tree.Branch("t", time, 'var/D')
tree.Branch("barOrEnd", barOrEnd, "var/D")
tree.Branch("location", location, "var/D")
tree.Branch("module", system, 'var/D')
tree.Branch("layer", layer, 'var/D')
tree.Branch("side", side, 'var/D')

#Comment out one of the two fnames definitions to run the other
#No BIB input files
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl_v1/muonGun/reco/*.slcio")

#BIB input files
fnames = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/recoBIB/muonGun_reco_1[123]0.slcio")

# Loop over files
i = 0 #keep track of which event we are on for readout purposes
for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events
    for event in reader:
        print(i)
        if i %100 ==1:
            print("Reading event ",i)
        # setting decoder
        hitsCollection = event.getCollection("SiTracks")
        encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        #Get hits within the collection
        for track in event.getCollection("SiTracks"):
            #Writing pointers to all the data for each hit
            if i %100 ==1:
                print("Tracks")
                dir(track)
                for hits in track:
                    print("Hits")
                    dir(hits)
            #Filling the data from the pointers into the tree
            tree.Fill()


output_file = TFile(options.outFile, 'RECREATE')
tree.Write()