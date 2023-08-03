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

#Initialize the branches for the info to be stored in:
x_pos1 = array('d', [0])
y_pos1 = array('d', [0])
z_pos1 = array('d', [0])
time1 = array ('d', [0])
theta1 = array('d', [0])
phi1 = array('d', [0])
barOrEnd1 = array('d', [0]) #0 if in the barrel, 1 if in the endcap
location1 = array('d',[0]) #0 in vertex, 1 in inner layer, 2 in outer layer
system1 = array ('d', [0])
layer1 = array ('d', [0])
side1 = array('d', [0])

x_pos2 = array('d', [0])
y_pos2 = array('d', [0])
z_pos2 = array('d', [0])
time2 = array ('d', [0])
theta2 = array('d', [0])
phi2 = array('d', [0])
barOrEnd2 = array('d', [0]) #0 if in the barrel, 1 if in the endcap
location2 = array('d',[0]) #0 in vertex, 1 in inner layer, 2 in outer layer
system2 = array ('d', [0])
layer2 = array ('d', [0])
side2 = array('d', [0])

#Initialzing branches for each layer
tree.Branch("x1",  x_pos1,  'var/D')
tree.Branch("y1",  y_pos1,  'var/D')
tree.Branch("z1", z_pos1, 'var/D')
tree.Branch("t1", time1, 'var/D')
tree.Branch("theta1", theta1, 'var/D')
tree.Branch("phi1", phi1, 'var/D')
tree.Branch("module1", system1, 'var/D')
tree.Branch("layer1", layer1, 'var/D')
tree.Branch("side1", side1, 'var/D')

tree.Branch("x2",  x_pos2,  'var/D')
tree.Branch("y2",  y_pos2,  'var/D')
tree.Branch("z2", z_pos2, 'var/D')
tree.Branch("t2", time2, 'var/D')
tree.Branch("theta2", theta2, 'var/D')
tree.Branch("phi2", phi2, 'var/D')
tree.Branch("module2", system2, 'var/D')
tree.Branch("layer2", layer2, 'var/D')
tree.Branch("side2", side2, 'var/D')

#Comment out one of the two fnames definitions to run the other
#No BIB input files
fnames = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/recoBIB/*100.slcio")

#BIB input files
#fnames = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/recoBIB/muonGun_reco_1[123]0.slcio")

# Loop over files
i = 0 #keep track of which event we are on for readout purposes
for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events
    for event in reader:
        # setting decoder
        hitsCollection = event.getCollection("SiTracks")
        encoding = event.getCollection("VBTrackerHits").getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        #Get tracks within the collection
        for track in hitsCollection:
            #Writing pointers to all the data for the first hit in each track
            hits= [i for i in track.getTrackerHits()]
            x_pos1[0]=hits[0].getPositionVec().X()
            y_pos1[0]=hits[0].getPositionVec().Y()
            z_pos1[0]=hits[0].getPositionVec().Z()
            time1[0]=hits[0].getTime()
            theta1[0]=hits[0].getPositionVec().Theta()
            phi1[0]=hits[0].getPositionVec().Phi()

            #Decoder
            cellID1 = int(hits[0].getCellID0())
            decoder.setValue(cellID1)
            layer1[0] = decoder['layer'].value()
            system1[0] = decoder["system"].value()
            side1[0] = decoder["side"].value()

            #Second hit
            x_pos2[0]=hits[1].getPositionVec().X()
            y_pos2[0]=hits[1].getPositionVec().Y()
            z_pos2[0]=hits[1].getPositionVec().Z()
            time2[0]=hits[1].getTime()
            theta2[0]=hits[0].getPositionVec().Theta()
            phi2[0]=hits[0].getPositionVec().Phi()

            #Decoder
            cellID2 = int(hits[1].getCellID0())
            decoder.setValue(cellID2)
            layer2[0] = decoder['layer'].value()
            system2[0] = decoder["system"].value()
            side2[0] = decoder["side"].value()

            #Filling the data from the pointers into the tree
            tree.Fill()


output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
