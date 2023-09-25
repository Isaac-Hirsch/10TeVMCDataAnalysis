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

fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_0_50/muonGun_pT_0_50_reco_2?00.slcio")

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
    "YokeBarrelCollection",
    "YokeEndcapCollection"
    ]

maxPseudo=0
minPseudo=0
i=0
for f in fnames:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    for event in reader:

        vertexHitsCollection = event.getCollection('VBTrackerHits')
        encoding = vertexHitsCollection.getParameters(
        ).getStringVal(EVENT.LCIO.CellIDEncoding)
        decoder = UTIL.BitField64(encoding)

        relationCollection=event.getCollection('VBTrackerHitsRelations')
        relation = UTIL.LCRelationNavigator(relationCollection)

        for hit in vertexHitsCollection:
            particle = relation.getRelatedToObjects(hit)[0]
            print(f"From objects: {relation.getRelatedFromObjects(hit)}")
            print(f"To objects: {relation.getRelatedToObjects(hit)}")
            print(f"From types: {relation.getFromType(hit)}")
            print(f"To types: {relation.getToType(hit)}")
            mcp=particle.getMCParticle()
            pdg=mcp.getPDG()
            sim=mcp.getSimulatorStatus()
            gen=mcp.getGeneratorStatus()
            decay=mcp.isCreatedInSimulation() 
            overlay=mcp.isOverlay()
            print(f"Pdg: {pdg}")
            print(f"Sim: {sim}")
            print(f"Gen: {gen}")
            print(f"Decay: {decay}")
            print(f"Overlay: {overlay}")
            print("")