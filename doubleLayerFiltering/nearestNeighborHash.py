import numpy as np
import glob
from optparse import OptionParser
import json
from pyLCIO import IOIMPL, EVENT, UTIL

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile hitsPerLayer',
                  type=str, default='BIBEndcapAnalysis')
(options, args) = parser.parse_args()

#Gather all the files you want to run over
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/photonGun_pT_0_50/photonGun_pT_0_50_reco_2??0.slcio")

#Sorting function for sorting the second doublet hits
nPseudoRap=300
nPhi=300
sorting=[]
for i in range(6):
    sorting.append([])
    for j in range(nPseudoRap):
        sorting[i].append([])

        for k in range(nPhi):
            sorting[i][j].append([])

#Collections we will be looking at
collections=["VBTrackerHits", "VETrackerHits"]

for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events
    for event in reader:
        
        #Looking at the only doublet layer in the vertex barrel
        tracksCollection = event.getCollection("VBTrackerHits")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encoding=tracksCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        firstLayerHit=[]

        for hit in tracksCollection:
            #Decoder
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()

            if layer==1:
                PseudoRapidity=hit.getPositionVec().PseudoRapidity()
                phi=hit.getPositionVec().Phi()
                sorting[0][((2.2+PseudoRapidity)*nPseudoRap)//4.4][((np.pi+phi)*nPhi)//(2*np.pi)].append((PseudoRapidity,phi))

            elif layer==0:
                firstLayerHit.append((hit.getPositionVec().PseudoRapidity(),hit.getPositionVec().Phi()))
        
        for (psuedoRap,phi) in firstLayerHit:
            #Impletement finding nearest neighbor using delta R
            minUnknown=min([(2.2+PseudoRapidity)//4.4,(np.pi+phi)//(2*np.pi), 1-(2.2+PseudoRapidity)//4.4,1-(np.pi+phi)//(2*np.pi)])
            #need to keep track of which of the min cases it is so that you know which box to check next

        #Looking at the four doublet layer in the vertex endcaps
        tracksCollection = event.getCollection("VETrackerHits")
        firstLayerHit=[]

        for hit in event.getCollection(tracksCollection):
            #Decoder
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()
            if (layer==1) | (layer==3) | (layer==5) | (layer==7):
                side = decoder['side'].value()
                PseudoRapidity=hit.getPositionVec().PseudoRapidity()
                phi=hit.getPositionVec().Phi()
                sorting[layer//2+4*(side==1)][((2.2+PseudoRapidity)*nPseudoRap)//4.4][((np.pi+phi)*nPhi)//(2*np.pi)].append((PseudoRapidity,phi))

            else:
                firstLayerHit.append((hit.getPositionVec().PseudoRapidity(),hit.getPositionVec().Phi()))

            for hit in firstLayerHit:
                #Impletement finding nearest neighbor using delta R
                NotImplemented