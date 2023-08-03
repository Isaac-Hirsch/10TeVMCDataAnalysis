#Imports are meant to be run on snowmass21 server
from array import array
from pyLCIO import IOIMPL, EVENT, UTIL
import ROOT
from ROOT import TH1D, TH2D, TFile, TLorentzVector, TTree, TMath
import glob
from optparse import OptionParser

#Code needed to download data to a root file
Bfield = 3.56  # T
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_hits_SiTracksNOBIB.root',
                  type=str, default='ntup_hits_SiTracksNOBIB.root')
(options, args) = parser.parse_args()

tree = TTree("tracks_tree", "tracks_tree")

fnames = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/reco/*100.slcio")

collections=[
    "AllTracks",
    "ECALBarrel",
    "ECALEndcap",
    "ECalBarrelCollection",
    "ECalEndcapCollection",
    "HCALBarrel",
    "HCALEndcap",
    "HCALOther",
    "HCalBarrelCollection",
    "HCalEndcapCollection",
    "HCalRingCollection",
    "IBTrackerHits",
    "IBTrackerHitsRelations",
    "IETrackerHits",
    "IETrackerHitsRelations",
    "InnerTrackerBarrelCollection",
    "InnerTrackerEndcapCollection",
    "JetOut",
    "LE_LooseSelectedPandoraPFOs",
    "LE_SelectedPandoraPFOs",
    "LE_TightSelectedPandoraPFOs",
    "LooseSelectedPandoraPFOs",
    "MCParticle",
    "MCPhysicsParticles",
    "MUON",
    "OBTrackerHits",
    "OBTrackerHitsRelations",
    "OETrackerHits",
    "OETrackerHitsRelations",
    "OuterTrackerBarrelCollection",
    "OuterTrackerEndcapCollection",
    "PandoraClusters",
    "PandoraPFOs",
    "PandoraStartVertices",
    "RelationCaloHit",
    "RelationMuonHit",
    "SeedTracks",
    "SelectedPandoraPFOs",
    "SiTracks",
    "SiTracks_Refitted",
    "TightSelectedPandoraPFOs",
    "VBTrackerHits",
    "VBTrackerHitsRelations",
    "VETrackerHits",
    "VETrackerHitsRelations",
    "VertexBarrelCollection",
    "VertexEndcapCollection",
    "YokeBarrelCollection",
    "YokeEndcapCollection"
    ]
methods=[]

for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)
    events=[event for event in reader]
    for collection in collections:
        print(collection + ":")
        print(dir(events[0].getCollection(collection)))
        methods.append(dir(events[0].getCollection(collection)))
    print("Parameters:")
    for i in range(len(collections)):
        if ("getParameters" in methods[i]):
            print(collections[i])
