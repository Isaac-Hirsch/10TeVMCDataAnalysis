import numpy as np
import glob
from optparse import OptionParser
import json
from pyLCIO import IOIMPL, EVENT, UTIL

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile MCSi0NoBIB',
                  type=str, default='MCSi0NoBIB')
(options, args) = parser.parse_args()

fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_0_50/muonGun_pT_0_50_reco_*.slcio")

deltaPhi=[]
deltaTheta=[]
totalHits=[]

#Loop over every file
for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events in the file
    for event in reader:
        particlesCollection = event.getCollection("MCParticle_SiTracks_Refitted")

        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encoding=particlesCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)
        for particle in particlesCollection:
            track=particle.getTo()

            phi=[]
            theta=[]
            numHits=[]
            for i in range(9):
                phi.append(0)
                theta.append(0)
                numHits.append(0)
            

            for hit in track.getTrackerHits():
                #Decoder
                cellID = int(hit.getCellID0())
                if (decoder['system'].value()<3): #Requiring this be in the vertex layers
                    #Getting info on detector of hit
                    decoder.setValue(cellID)
                    layer = decoder['layer'].value()
                    side = decoder['side'].value()

                    #Index 0 is the barrel doublet, 1-4 are the -z endcaps from innermost to outer most and 5-6 are the +z endcaps also starting with innermost
                    index=(1+layer//2)*(side!=0)+4*(side==1)
                    
                    #Getting hit info
                    position=hit.getgetPositionVec()
                    phi[index]=((-1)**(layer%2))*position.Phi()
                    theta[index]=((-1)**(layer%2))*position.Theta()
                    numHits[index]+=1
            #Appending the particle data to the list
            deltaPhi.append(phi)
            deltaTheta.append(theta)
            totalHits.append(numHits)

#Wrapping data into a dictionary that will be exported as a json
output={
    "deltaPhi" : deltaPhi,
    "deltaTheta" : deltaTheta,
    "totalHits" : totalHits
}

output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
    json.dump(output, fp)