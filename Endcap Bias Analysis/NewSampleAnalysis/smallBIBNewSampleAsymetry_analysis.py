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
                  type=str, default='BIBEndcapAnalysis')
(options, args) = parser.parse_args()


#Gathering all muonGun files with BIB
muon250BIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/muonGun_pT_250_1000/muonGun_pT_250_1000_reco_[12]?00.slcio")

#Gathering all muonGun files with BIB
muon0BIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/muonGun_pT_0_50/muonGun_pT_0_50_reco_[78]?00.slcio")

#Gathering all muonGun files with BIB
muon1000BIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/muonGun_pT_1000_5000/muonGun_pT_1000_5000_reco_[56]?00.slcio")

photonBIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/photonGun_pT_0_50/photonGun_pT_0_50_reco_[34]?00.slcio")

muon250BIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

for file in muon250BIBFiles:    
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
                muon250BIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

muon0BIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

for file in muon0BIBFiles:
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
                muon0BIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

muon1000BIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

for file in muon1000BIBFiles:
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
                muon1000BIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

photonBIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

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
                photonBIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

output={
    "muon0BIB" : muon0BIBSiTrack,
    "muon250BIB" : muon250BIBSiTrack,
    "muon1000BIB" : muon1000BIBSiTrack,
    "photonBIB" : photonBIBSiTrack
}

output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
    json.dump(output, fp)