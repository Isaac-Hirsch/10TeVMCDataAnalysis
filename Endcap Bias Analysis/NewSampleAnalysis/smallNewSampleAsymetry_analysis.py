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
parser.add_option('-o', '--outFile', help='--outFile NoBIBHitsPerLayer',
                  type=str, default='NoBIBHitsPerLayer')
(options, args) = parser.parse_args()


#Gathering all muonGun files without BIB
muon250NoBIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_250_1000/muonGun_pT_250_1000_reco_*.slcio")

muon250NoBIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

for file in muon250NoBIBFiles:
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
                muon250NoBIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

del muon250NoBIBFiles

#Gathering all pionGun files without BIB
muon0NoBIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_0_50/muonGun_pT_0_50_reco_*.slcio")

muon0NoBIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

for file in muon0NoBIBFiles:
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
                muon0NoBIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

del muon0NoBIBFiles

#Gathering all pionGun files without BIB
muon1000NoBIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_1000_5000/muonGun_pT_1000_5000_reco_*.slcio")

muon1000NoBIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

for file in muon1000NoBIBFiles:
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
                muon1000NoBIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

del muon1000NoBIBFiles

photonNoBIBFiles=glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/photonGun_pT_0_50/photonGun_pT_0_50_reco_*.slcio")

photonNoBIBSiTrack=[np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(8,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(7,dtype=int).tolist(),np.zeros(4,dtype=int).tolist(),np.zeros(3,dtype=int).tolist(),np.zeros(4,dtype=int).tolist()]

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
                photonNoBIBSiTrack[3*((system-1)//2)+side][layer]+=1
reader.close()

del photonNoBIBFiles

output={
    "muon0NoBIB" : muon0NoBIBSiTrack,
    "muon250NoBIB" : muon250NoBIBSiTrack,
    "muon1000NoBIB" : muon1000NoBIBSiTrack,
    'photonNoBIB' :  photonNoBIBSiTrack
}

output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
            json.dump(output, fp)
