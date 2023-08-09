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
parser.add_option('-o', '--outFile', help='--outFile ntup_hits_AllTracksBIB.root',
                  type=str, default='ntup_hits_AllTracksBIB.root')
(options, args) = parser.parse_args()

tree = TTree("tracks_tree", "tracks_tree")

#Initialize the branches for the info to be stored in:
x_pos1 = array('d', [0])
y_pos1 = array('d', [0])
z_pos1 = array('d', [0])
time1 = array ('d', [0])
theta1 = array('d', [0])
phi1 = array('d', [0])
system1 = array ('d', [0])
layer1 = array ('d', [0])
side1 = array('d', [0])

x_pos2 = array('d', [0])
y_pos2 = array('d', [0])
z_pos2 = array('d', [0])
time2 = array ('d', [0])
theta2 = array('d', [0])
phi2 = array('d', [0])
system2 = array ('d', [0])
layer2 = array ('d', [0])
side2 = array('d', [0])

x_pos3 = array('d', [0])
y_pos3 = array('d', [0])
z_pos3 = array('d', [0])
time3 = array ('d', [0])
theta3 = array('d', [0])
phi3 = array('d', [0])
system3 = array ('d', [0])
layer3 = array ('d', [0])
side3 = array('d', [0])

x_pos4 = array('d', [0])
y_pos4 = array('d', [0])
z_pos4 = array('d', [0])
time4 = array ('d', [0])
theta4 = array('d', [0])
phi4 = array('d', [0])
system4 = array ('d', [0])
layer4 = array ('d', [0])
side4 = array('d', [0])

#Initialzing branches for each layer
tree.Branch("x1",  x_pos1,  'var/D')
tree.Branch("y1",  y_pos1,  'var/D')
tree.Branch("z1", z_pos1, 'var/D')
tree.Branch("t1", time1, 'var/D')
tree.Branch("theta1", theta1, 'var/D')
tree.Branch("phi1", phi1, 'var/D')
tree.Branch("system1", system1, 'var/D')
tree.Branch("layer1", layer1, 'var/D')
tree.Branch("side1", side1, 'var/D')

tree.Branch("x2",  x_pos2,  'var/D')
tree.Branch("y2",  y_pos2,  'var/D')
tree.Branch("z2", z_pos2, 'var/D')
tree.Branch("t2", time2, 'var/D')
tree.Branch("theta2", theta2, 'var/D')
tree.Branch("phi2", phi2, 'var/D')
tree.Branch("system2", system2, 'var/D')
tree.Branch("layer2", layer2, 'var/D')
tree.Branch("side2", side2, 'var/D')

tree.Branch("x3",  x_pos3,  'var/D')
tree.Branch("y3",  y_pos3,  'var/D')
tree.Branch("z3", z_pos3, 'var/D')
tree.Branch("t3", time3, 'var/D')
tree.Branch("theta3", theta3, 'var/D')
tree.Branch("phi3", phi3, 'var/D')
tree.Branch("system3", system3, 'var/D')
tree.Branch("layer3", layer3, 'var/D')
tree.Branch("side3", side3, 'var/D')

tree.Branch("x4",  x_pos4,  'var/D')
tree.Branch("y4",  y_pos4,  'var/D')
tree.Branch("z4", z_pos4, 'var/D')
tree.Branch("t4", time4, 'var/D')
tree.Branch("theta4", theta4, 'var/D')
tree.Branch("phi4", phi4, 'var/D')
tree.Branch("system4", system4, 'var/D')
tree.Branch("layer4", layer4, 'var/D')
tree.Branch("side4", side4, 'var/D')

#Comment out one of the two fnames definitions to run the other
#No BIB input files
#fnames = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/reco/*.slcio")

#BIB input files
fnames = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/recoBIB/muonGun_reco_1[123]0.slcio")

# Loop over files
i = 0 #keep track of which event we are on for readout purposes
for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events
    for event in reader:
        # setting decoder
        tracksCollection = event.getCollection("AllTracks")
        encodCol = event.getCollection("IBTrackerHits")
        encoding = encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        #Get tracks within the collection
        for track in tracksCollection:
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
            theta2[0]=hits[1].getPositionVec().Theta()
            phi2[0]=hits[1].getPositionVec().Phi()

            #Decoder
            cellID2 = int(hits[1].getCellID0())
            decoder.setValue(cellID2)
            layer2[0] = decoder['layer'].value()
            system2[0] = decoder["system"].value()
            side2[0] = decoder["side"].value()

            #Third hit
            if len(hits) >= 3:
                x_pos3[0]=hits[2].getPositionVec().X()
                y_pos3[0]=hits[2].getPositionVec().Y()
                z_pos3[0]=hits[2].getPositionVec().Z()
                time3[0]=hits[2].getTime()
                theta3[0]=hits[2].getPositionVec().Theta()
                phi3[0]=hits[2].getPositionVec().Phi()

                #Decoder
                cellID3 = int(hits[2].getCellID0())
                decoder.setValue(cellID3)
                layer3[0] = decoder['layer'].value()
                system3[0] = decoder["system"].value()
                side3[0] = decoder["side"].value()

            else:
                x_pos3[0]=0
                y_pos3[0]=0
                z_pos3[0]=0
                time3[0]=0
                theta3[0]=0
                phi3[0]=0

                #Decoder
                layer3[0] = 0
                system3[0] = 0
                side3[0] = 0

            

            #Fourth hit
            if len(hits) >= 4:
                x_pos4[0]=hits[3].getPositionVec().X()
                y_pos4[0]=hits[3].getPositionVec().Y()
                z_pos4[0]=hits[3].getPositionVec().Z()
                time4[0]=hits[3].getTime()
                theta4[0]=hits[3].getPositionVec().Theta()
                phi4[0]=hits[3].getPositionVec().Phi()

                #Decoder
                cellID4 = int(hits[3].getCellID0())
                decoder.setValue(cellID4)
                layer4[0] = decoder['layer'].value()
                system4[0] = decoder["system"].value()
                side4[0] = decoder["side"].value()
            
            else:
                x_pos4[0]=0
                y_pos4[0]=0
                z_pos4[0]=0
                time4[0]=0
                theta4[0]=0
                phi4[0]=0

                #Decoder
                layer4[0] = 0
                system4[0] = 0
                side4[0] = 0

            #Filling the data from the pointers into the tree
            tree.Fill()


output_file = TFile(options.outFile, 'RECREATE')
tree.Write()
