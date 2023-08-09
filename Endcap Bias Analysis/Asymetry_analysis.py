#Imports are meant to be run on snowmass21 server
from array import array
from pyLCIO import IOIMPL, EVENT, UTIL
import ROOT
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath
import glob
from optparse import OptionParser
import numpy as np
import json

#Boiler plate to export files from server to be analyzed locally
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile hitsPerLayer',
                  type=str, default='hitsPerLayer')
(options, args) = parser.parse_args()


#Gathering all muonGun files without BIB
muonNoBIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/reco/*.slcio")

#Gathering all pionGun files without BIB
photonNoBIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/photonGun/reco/*.slcio")

#Gathering all pionGun files without BIB
pionNoBIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/pionGun/reco/*.slcio")

#Gathering all muonGun files with BIB
muonBIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/recoBIB/?*.slcio")

#Gathering all muonGun files with BIB
photonBIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/photonGun/recoBIB/1?*.slcio")

#Gathering all muonGun files with BIB
pionBIBFiles=glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/pionGun/recoBIB/2?*.slcio")

muonNoBIBSiTrack=[np.zeros(8),np.zeros(8),np.zeros(8),np.zeros(7),np.zeros(3),np.zeros(7),np.zeros(4),np.zeros(3),np.zeros(4)]

for file in muonNoBIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        for track in tracksCollection:
            
            for hit in track.getTrackerHits():
                #Gathering data to calssify hits
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()
                side = decoder["side"].value() +1
                layer = decoder['layer'].value()

                #Counting the hit
                muonNoBIBSiTrack[3*(system-1)//2+side][layer]+=1
reader.close()

muonBIBSiTrack=[np.zeros(8),np.zeros(8),np.zeros(8),np.zeros(7),np.zeros(3),np.zeros(7),np.zeros(4),np.zeros(3),np.zeros(4)]

for file in muonBIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        for track in tracksCollection:
            
            for hit in track.getTrackerHits():
                #Gathering data to calssify hits
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()
                side = decoder["side"].value() +1
                layer = decoder['layer'].value()

                #Counting the hit
                muonBIBSiTrack[3*(system-1)//2+side][layer]+=1
reader.close()

photonNoBIBSiTrack=[np.zeros(8),np.zeros(8),np.zeros(8),np.zeros(7),np.zeros(3),np.zeros(7),np.zeros(4),np.zeros(3),np.zeros(4)]

for file in photonNoBIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        for track in tracksCollection:
            
            for hit in track.getTrackerHits():
                #Gathering data to calssify hits
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()
                side = decoder["side"].value() +1
                layer = decoder['layer'].value()

                #Counting the hit
                photonNoBIBSiTrack[3*(system-1)//2+side][layer]+=1
reader.close()

photonBIBSiTrack=[np.zeros(8),np.zeros(8),np.zeros(8),np.zeros(7),np.zeros(3),np.zeros(7),np.zeros(4),np.zeros(3),np.zeros(4)]

for file in photonBIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        for track in tracksCollection:
            
            for hit in track.getTrackerHits():
                #Gathering data to calssify hits
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()
                side = decoder["side"].value() +1
                layer = decoder['layer'].value()

                #Counting the hit
                photonBIBSiTrack[3*(system-1)//2+side][layer]+=1
reader.close()

pionNoBIBSiTrack=[np.zeros(8),np.zeros(8),np.zeros(8),np.zeros(7),np.zeros(3),np.zeros(7),np.zeros(4),np.zeros(3),np.zeros(4)]

for file in pionNoBIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        for track in tracksCollection:
            
            for hit in track.getTrackerHits():
                #Gathering data to calssify hits
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()
                side = decoder["side"].value() +1
                layer = decoder['layer'].value()

                #Counting the hit
                pionNoBIBSiTrack[3*(system-1)//2+side][layer]+=1
reader.close()

pionBIBSiTrack=[np.zeros(8),np.zeros(8),np.zeros(8),np.zeros(7),np.zeros(3),np.zeros(7),np.zeros(4),np.zeros(3),np.zeros(4)]

for file in pionBIBFiles:
    reader=IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)

    for event in reader:
        tracksCollection = event.getCollection("SiTracks")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encodCol=event.getCollection("IBTrackerHits")
        encoding=encodCol.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        for track in tracksCollection:
            
            for hit in track.getTrackerHits():
                #Gathering data to calssify hits
                cellID=int(hit.getCellID0())
                decoder.setValue(cellID)
                system=decoder["system"].value()
                side = decoder["side"].value() +1
                layer = decoder['layer'].value()

                #Counting the hit
                pionBIBSiTrack[3*(system-1)//2+side][layer]+=1
reader.close()

output={
    "photonBIB" : photonBIBSiTrack,
    "photonNoBIB" : photonNoBIBSiTrack,
    "muonBIB" : muonBIBSiTrack,
    "muonNoBIB" : muonNoBIBSiTrack,
    "pionBIB" : pionBIBSiTrack,
    "pionNoBIB" : pionNoBIBSiTrack
}

output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
            json.dump(output, fp)