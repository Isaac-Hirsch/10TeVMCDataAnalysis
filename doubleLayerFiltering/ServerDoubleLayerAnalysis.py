from array import array
from pyLCIO import IOIMPL, EVENT, UTIL
import ROOT
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath
import glob
from optparse import OptionParser
import numpy as np
import json

Bfield=3.56  # T

#Boiler plate to export files from server to be analyzed locally
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile doubleLayerHits',
                  type=str, default='doubleLayerHits')
(options, args) = parser.parse_args()

#Gathering all muonGun files without BIB
noBIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/reco/*.slcio")

#Gathering the 29 files from 10-290 in the muon gun samples
BIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/photonGun/recoBIB/photonGun_reco_?0.slcio")

#initializing lists of lists to store the results in. For the outer list, 0-3 will be doublets in the -z endcaps, 4-7 are doublets in the barrel, and 8-11 are doublets in the +z endcaps
noBIBDeltaTheta=[]
noBIBDeltaPhi=[]

#Initialize an array for each doublet:
for i in range(12):
    noBIBDeltaTheta.append([])
    noBIBDeltaPhi.append([])

#Calculating the delta phi and delta theta for all noBIB samples for all SiTracks and vertex doublets
for file in noBIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)


    for event in reader:
        #We will be using SiTracks collection to find which hits coorespond to the same track.
        #It would be signficantly better if we could find a way to pair hits in the doublet layers prior to making tracks
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        #Get tracks within the collection
        for track in tracksCollection:
            
            #Creating an array to remember if the tracks has a hit on each layer of the vertex (since the vertex is the silicon detector with doublet layers)
            #thetaZeroHit will track the theta of hits on the first layer of the doublets with 0-3 being in the -z endcap, 4-7 being in the +z endcap, and 8-11 being in the barrel
            #thetaOnetHit will do the same with the second layer of the doublets
            zeroHit=np.zeros(12, dtype=bool)
            oneHit=np.zeros(12, dtype=bool)
            thetaZeroHit=np.zeros(12)
            thetaOneHit=np.zeros(12)

            #phiZeroHit and phiOnetHit will do the same with phi
            phiZeroHit=np.zeros(12)
            phiOneHit=np.zeros(12)

            #Gathering the hits in each track
            for hit in track.getTrackerHits():
                #Get information on which detector the hit came from
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()

                #Check whether the hit was in the vertex barrel (system==1) or in the vertex endcaps (system==2)
                if ((system==1) | (system==2)):
                    #Layer counts which layer of the detector you are in, with 0 being the smallest radius barrel layer or smallest z endcap
                    layer = decoder['layer'].value()
                    #side is 0 if the hit was int he -z endcap, 1 if in the barrel, and 2 if in the +z endcap
                    side = decoder["side"].value() + 1

                    #retrieve the position vector for each hit
                    position=hit.getPositionVec()
    
                    # setting the thetas and phi of the first layer of the doublets
                    if (layer % 2) ==0:
                        #Using a hashing function to map each layer to a spot in the array
                        id=4*side+layer//2
                        zeroHit[id]=True
                        thetaZeroHit[id]=position.Theta()
                        phiZeroHit[id]=position.Phi()
                    
                    #doing it for the second layer
                    if (layer % 2) ==1:
                        #Using a hashing function to map each layer to a spot in the array
                        id=4*side+layer//2
                        oneHit[id]=True
                        thetaOneHit[id]=position.Theta()
                        phiOneHit[id]=position.Phi()
            doubleID=zeroHit & oneHit
            for i in range(12):
                if doubleID[i]:
                    noBIBDeltaTheta[i].append(thetaOneHit[i]-thetaZeroHit[i])
                    noBIBDeltaPhi[i].append(phiOneHit[i]-phiZeroHit[i])
reader.close()

#initializing lists of lists to store the results in. For the outer list, 0-3 will be doublets in the -z endcaps, 4-7 are doublets in the barrel, and 8-11 are doublets in the +z endcaps
BIBDeltaTheta=[]
BIBDeltaPhi=[]
#Initialize a list for each doublet:
for i in range(12):
    BIBDeltaTheta.append([])
    BIBDeltaPhi.append([])

#Repeating the previous study fro BIB files     
for file in BIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)


    for event in reader:
        #We will be using SiTracks collection to find which hits coorespond to the same track.
        #It would be signficantly better if we could find a way to pair hits in the doublet layers prior to making tracks
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        #Get tracks within the collection
        for track in tracksCollection:
            
            #Creating an array to remember if the tracks has a hit on each layer of the vertex (since the vertex is the silicon detector with doublet layers)
            #thetaZeroHit will track the theta of hits on the first layer of the doublets with 0-3 being in the -z endcap, 4-7 being in the +z endcap, and 8-11 being in the barrel
            #thetaOnetHit will do the same with the second layer of the doublets
            zeroHit=np.zeros(12, dtype=bool)
            oneHit=np.zeros(12, dtype=bool)
            thetaZeroHit=np.zeros(12)
            thetaOneHit=np.zeros(12)

            #phiZeroHit and phiOnetHit will do the same with phi
            phiZeroHit=np.zeros(12)
            phiOneHit=np.zeros(12)

            #Gathering the hits in each track
            for hit in track.getTrackerHits():
                #Get information on which detector the hit came from
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()

                #Check wether the hit was in the vertex barrel (system==1) or in the vertex endcaps (system==2)
                if ((system==1) | (system==2)):
                    #Layer counts which layer of the detector you are in, with 0 being the smallest radius barrel layer or smallest z endcap
                    layer = decoder['layer'].value()
                    #side is 0 if the hit was int he -z endcap, 1 if in the barrel, and 2 if in the +z endcap
                    side = decoder["side"].value() + 1

                    #retrieve the position vector for each hit
                    position=hit.getPositionVec()
    
                    # setting the thetas and phi of the first layer of the doublets
                    if (layer % 2) ==0:
                        #Using a hashing function to map each layer to a spot in the array
                        id=4*side+layer//2
                        zeroHit[id]=True
                        thetaZeroHit[id]=position.Theta()
                        phiZeroHit[id]=position.Phi()
                    
                    #doing it for the second layer
                    if (layer % 2) ==1:
                        #Using a hashing function to map each layer to a spot in the array
                        id=4*side+layer//2
                        oneHit[id]=True
                        thetaOneHit[id]=position.Theta()
                        phiOneHit[id]=position.Phi()
            doubleID=zeroHit & oneHit
            for i in range(12):
                if doubleID[i]:
                    BIBDeltaTheta[i].append(thetaOneHit[i]-thetaZeroHit[i])
                    BIBDeltaPhi[i].append(phiOneHit[i]-phiZeroHit[i])
reader.close()

output={
    "BIB/negZTheta" : BIBDeltaTheta[:4],
    "BIB/posZTheta" : BIBDeltaTheta[4:8],
    "BIB/barTheta" : BIBDeltaTheta[8:],
    "BIB/negZPhi" : BIBDeltaPhi[:4],
    "BIB/posZPhi" : BIBDeltaPhi[4:8],
    "BIB/barPhi" : BIBDeltaPhi[8:],
    "noBIB/negZTheta" : noBIBDeltaTheta[:4],
    "noBIB/posZTheta" : noBIBDeltaTheta[4:8],
    "noBIB/barTheta" : noBIBDeltaTheta[8:],
    "noBIB/negZPhi" : noBIBDeltaPhi[:4],
    "noBIB/posZPhi" : noBIBDeltaPhi[4:8],
    "noBIB/barPhi" : noBIBDeltaPhi[8:]
}

output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
            json.dump(output, fp)