#Imports are meant to be run on snowmass21 server
from array import array
from pyLCIO import IOIMPL, EVENT, UTIL
import ROOT
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath
import glob
from optparse import OptionParser
import numpy as np
from matplotlib import pyplot as plt
import json

#Code needed to export dats
Bfield = 5  # T
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile occupancyRate',
                  type=str, default='occupancyRate')
(options, args) = parser.parse_args()

#BIB input files
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/photonGun_pT_0_50/photonGun_pT_0_50_reco_2[12]?0.slcio")

#List of collections
collections=[
    "VBTrackerHits",
    "VETrackerHits",
    "IBTrackerHits",
    "IETrackerHits",
    "OBTrackerHits",
    "OETrackerHits"
]

#Keep track of number of events so occupancy cans be normalized to being per event
numEvents=0

#layers for design v0
layNum=[8,16,3,14,3,8]
layers=[]
for i in layNum:
    sys=[]
    for j in range(i):
          sys.append(0)
    layers.append(sys)

#Loop for through all
for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events
    for event in reader:
        #Count up the number of events and provide progress report
        numEvents+=1
        if (numEvents %100) ==0:
                print("Reading event ",numEvents)

        #Loop over collections
        for icol, collection in enumerate(collections):
            # setting decoder
            hitsCollection = event.getCollection(collection)
            encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)

            #Get hits within the collection
            for hit in event.getCollection(collection):
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                layer = decoder['layer'].value()
                system = decoder["system"].value()
                side = decoder["side"].value()
                layers[system-1][layer*(1+((side==1)&(system==2))*8+((side==1)&(system==4))*7)+((side==1)&(system==6))*4]+=1
reader.close()

#Normalizing to being per event
for num, i in enumerate(layNum):
    for j in range(i):
        layers[num][j]=layers[num][j]/numEvents

#Outputting file to a json
output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
            json.dump(layers, fp)