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

#Gather all the files you want to run over.
#Comment out all but 1 fnames
#BIB
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/recoBIB/photonGun_pT_0_50/photonGun_pT_0_50_reco_2[123]00.slcio")
#0-50 pt muons
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_0_50/muonGun_pT_0_50_reco_*.slcio")
#250-1000 pt muons
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_250_1000/muonGun_pT_250_1000_reco_*.slcio")

#Setting number of bins for the sorting function
nPseudoRap=50
nPhi=50
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
        #The Pseudonapidity spans -2.4 to 2.4 and is split into nPsuedoRap evenly sized bins
        #The Phi spans -pi to pi and is split into nPhi evenly sized bins
        sorting=[]
        for i in range(nPseudoRap):
            sorting.append([])
            for k in range(nPhi):
                sorting[i].append([])

        #Looking at the only doublet layer in the vertex barrel
        hitsCollection = event.getCollection("VBTrackerHits")
        #creating a decoder that will be used layer to trace a hit back to its system and layer
        encoding=hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder=UTIL.BitField64(encoding)

        firstLayerHit=[]

        for hit in hitsCollection:
            #Decoder
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()
            
            #This is the second layer of the only vertex barrel doublet
            if layer==1:
                pseudoRapidity=hit.getPositionVec().PseudoRapidity()
                phi=hit.getPositionVec().Phi()
                #Hashing the hit into the correct bins
                sorting[int(((2.4+pseudoRapidity)*nPseudoRap)/4.8)][int(((np.pi+phi)*nPhi)/(2*np.pi))].append((pseudoRapidity,phi))


            elif layer==0:
                firstLayerHit.append((hit.getPositionVec().PseudoRapidity(),hit.getPositionVec().Phi()))
        
        #Run through every hit in the first layer of the doublet and find its nearest neighbor in terms of delta R on the second layer
        for (pseudoRap,phi) in firstLayerHit:
            #Navive search to be fixed later. Does not work well if the nearest neighbor is not in the same box the hit is in.
            box=sorting[int(((2.4+pseudoRap)*nPseudoRap)/4.8)][int(((np.pi+phi)*nPhi)/(2*np.pi))]
            #Checking if there is anything in the box, if not its skipped and added to nBox which is a measure of how good of an approximation searching 1 box is
            if len(box) !=0:
                minPseudo=np.abs(box[0][0]-pseudoRap)
                minPhi=np.abs(box[0][1]-phi)
                minRad=np.sqrt(minPseudo**2+minPhi**2)
                for (boxPseduo, boxPhi) in box[1:]:
                    deltaPseudoPos=np.abs(boxPseduo-pseudoRap)
                    deltaPhiPos=np.abs(boxPhi-phi)
                    deltaRadPos=np.sqrt(deltaPseudoPos**2+deltaPhiPos**2)
                    if minRad > deltaRadPos:
                        minRad=deltaRadPos
                        minPhi=deltaPhiPos
                        minPseudo=deltaPseudoPos
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
        hitsCollection = event.getCollection("VETrackerHits")
        firstLayerHit=[]

        for hit in hitsCollection:
            #Decoder
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            layer = decoder['layer'].value()
            side = decoder['side'].value()
            #Identifying any hit that is the second layer of any of the four doublets
            if (layer==1) | (layer==3) | (layer==5) | (layer==7):
                pseudoRapidity=hit.getPositionVec().PseudoRapidity()
                phi=hit.getPositionVec().Phi()
                #layer/2+4*(side==1) uniquely hashes each outer doublet endcap into a value of 0-7
                sorting[int(layer/2+4*(side==1))][int(((2.4+pseudoRapidity)*nPseudoRap)/4.8)][int(((np.pi+phi)*nPhi)/(2*np.pi))].append((pseudoRapidity,phi))

            #All other hits are in the first layer of a doublet
            else:
                #layer/2+4*(side==1) uniquely hashes each inner doublet endcap into a value of 0-7
                print("Input: "+str(int(layer/2+4*(side==1))))
                firstLayerHit.append((hit.getPositionVec().PseudoRapidity(),hit.getPositionVec().Phi(),int(layer/2+4*(side==1))))

        for (psuedoRap,phi, pixel) in firstLayerHit:
            #Navive search to be fixed later
            #Pixel represents the endcap hash we should be looking at
            box=sorting[pixel][int(((2.4+pseudoRapidity)*nPseudoRap)/4.8)][int(((np.pi+phi)*nPhi)/(2*np.pi))]
            if len(box) !=0:
                minPseudo=np.abs(box[0][0]-pseudoRap)
                minPhi=np.abs(box[0][1]-phi)
                minRad=np.sqrt(minPseudo**2+minPhi**2)
                for (boxPseduo, boxPhi) in box[1:]:
                    deltaPseudoPos=np.abs(boxPseduo-pseudoRap)
                    deltaPhiPos=np.abs(boxPhi-phi)
                    deltaRadPos=np.sqrt(deltaPseudoPos**2+deltaPhiPos**2)
                    if minRad > deltaRadPos:
                        minRad=deltaRadPos
                        minPhi=deltaPhiPos
                        minPseudo=deltaPseudoPos
                #Add one to hash in delta{} because they need to account for the barrel doublet in the first spot
                print("output: " +str(1+pixel))
                deltaPseudo[1+pixel].append(minPseudo)
                deltaPhi[1+pixel].append(minPhi)
                deltaR[1+pixel].append(minRad)
            
            else:
                print("box: " +str(1+pixel))
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