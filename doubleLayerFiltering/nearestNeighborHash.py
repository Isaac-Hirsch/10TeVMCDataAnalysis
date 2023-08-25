import numpy as np
import glob
from optparse import OptionParser
import json
from pyLCIO import IOIMPL, EVENT, UTIL

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile nearestPair',
                  type=str, default='nearestPair')
(options, args) = parser.parse_args()

#Gather all the files you want to run over
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/photonGun_pT_0_50/photonGun_pT_0_50_reco_2[012]00.slcio")

#Setting number of bins for the sorting function
nPseudoRap=200
nPhi=200
#Creating output data storage
deltaPseudo=[]
deltaPhi=[]
deltaR=[]
nBox=[]
#Filling them with the correct number and dimension of empty lists
for i in range(9):
    deltaPseudo.append([])
    deltaPhi.append([])
    deltaR.append([])
    nBox.append(0)

for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    #Loop over events
    for event in reader:

        #Creating bins to sort hits from the second layer into.
        #The Pseudonapidity spans -2.2 to 2.2 and is split into nPsuedoRap evenly sized bins
        #The Phi spans -pi to pi and is split into nPhi evenly sized bins
        sorting=[]
        for i in range(nPseudoRap):
            sorting.append([])
            for k in range(nPhi):
                sorting[i].append([])

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
            
            #This is the second layer of the only vertex barrel doublet
            if layer==1:
                pseudoRapidity=hit.getPositionVec().PseudoRapidity()
                phi=hit.getPositionVec().Phi()
                #Hashing the hit into the correct bins
                print("Pseudo:" + str(int(((2.2+pseudoRapidity)*nPseudoRap)/4.4)))
                print("Phi:" + str(int(((np.pi+phi)*nPhi)/(2*np.pi))))
                sorting[int(((2.2+pseudoRapidity)*nPseudoRap)/4.4)][int(((np.pi+phi)*nPhi)/(2*np.pi))].append((pseudoRapidity,phi))


            elif layer==0:
                firstLayerHit.append((hit.getPositionVec().PseudoRapidity(),hit.getPositionVec().Phi()))
        
        #Run through every hit in the first layer of the doublet and find its nearest neighbor in terms of delta R on the second layer
        for (pseudoRap,phi) in firstLayerHit:
            #Navive search to be fixed later. Does not work well if the nearest neighbor is not in the same box the hit is in.
            box=sorting[int(((2.2+pseudoRap)*nPseudoRap)/4.4)][int(((np.pi+phi)*nPhi)/(2*np.pi))]
            #Checking if there is anything in the box, if not its skipped and added to nBox which is a measure of how good of an approximation searching 1 box is
            if len(box) !=0:
                minPseudo=box[0][0]
                minPhi=box[0][1]
                minRad=np.sqrt(minPseudo**2+minPhi**2)
                for (boxPseduo, boxPhi) in box[1:]:
                    boxRad=np.sqrt(boxPseduo**2+boxPhi**2)
                    if minRad > boxRad:
                        minRad=boxRad
                        minPhi=boxPhi
                        minPseudo=boxPseduo
                deltaPseudo[0].append(minPseudo)
                deltaPhi[0].append(minPhi)
                deltaR[0].append(minRad)
                
            else:
                nBox[0]+=1
                
        #Reseting sorting function
        sorting=[]
        for j in range(8):
            sorting.append([])
            for i in range(nPseudoRap):
                sorting[j].append([])
                for k in range(nPhi):
                    sorting[j][i].append([])
        
        #Looking at the four doublet layer in the vertex endcaps
        tracksCollection = event.getCollection("VETrackerHits")
        firstLayerHit=[]

        for hit in event.getCollection(tracksCollection):
            #Decoder
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()
            #Identifying any hit that is the second layer of any of the four doublets
            if (layer==1) | (layer==3) | (layer==5) | (layer==7):
                side = decoder['side'].value()
                pseudoRapidity=hit.getPositionVec().PseudoRapidity()
                phi=hit.getPositionVec().Phi()
                #layer/2+4*(side==1) uniquely hashes each outer doublet endcap into a value of 0-7
                sorting[int(layer/2+4*(side==1))][int(((2.2+pseudoRapidity)*nPseudoRap)/4.4)][int(((np.pi+phi)*nPhi)/(2*np.pi))].append((pseudoRapidity,phi))

            #All other hits are in the first layer of a doublet
            else:
                #layer/2+4*(side==1) uniquely hashes each inner doublet endcap into a value of 0-7
                firstLayerHit.append((hit.getPositionVec().PseudoRapidity(),hit.getPositionVec().Phi(),int(layer/2+4*(side==1))))

            for (psuedoRap,phi, pixel) in firstLayerHit:
                #Navive search to be fixed later
                #Pixel represents the endcap hash we should be looking at
                box=sorting[pixel][int(((2.2+pseudoRapidity)*nPseudoRap)/4.4)][int(((np.pi+phi)*nPhi)/(2*np.pi))]
                if len(box) !=0:
                    minPseudo=box[0][0]
                    minPhi=box[0][1]
                    minRad=np.sqrt(minPseudo**2+minPhi**2)
                    for (boxPseduo, boxPhi) in box[1:]:
                        boxRad=np.sqrt(boxPseduo**2+boxPhi**2)
                        if minRad > boxRad:
                            minRad=boxRad
                            minPhi=boxPhi
                            minPseudo=boxPseduo
                    #Add one to hash in delta{} because they need to account for the barrel doublet in the first spot
                    deltaPseudo[1+pixel].append(minPseudo)
                    deltaPhi[1+pixel].append(minPhi)
                    deltaR[1+pixel].append(minRad)
                
                else:
                    nBox[pixel+1]+=1

#Wrapping data into a dictionary that will be exported as a json
output={
    "deltaPesudorapdity" : deltaPseudo,
    "deltaPhi" : deltaPhi,
    "deltaR" : deltaR,
    "emptyBoxes" : nBox
}

output_json = options.outFile+".json"
with open(output_json, 'w') as fp:
    json.dump(output, fp)