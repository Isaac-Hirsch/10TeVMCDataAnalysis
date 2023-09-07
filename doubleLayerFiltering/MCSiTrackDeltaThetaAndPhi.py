import numpy as np
import glob
from optparse import OptionParser
import json
from pyLCIO import IOIMPL, EVENT, UTIL

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile nearestPairBIB',
                  type=str, default='nearestPairBIB')
(options, args) = parser.parse_args()

fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_0_50/muonGun_pT_0_50_reco_1000.slcio")

#Loop over every file
for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events in the file
    for event in reader:
        tracksCollection = event.getCollection("MCParticle_SiTracks_Refitted")
        print("Collection:")
        print(dir(tracksCollection))
        for track in tracksCollection:
            print("Track:")
            print(dir(track))
            for hit in track:
                print("Hit")
                print(dir(hit))