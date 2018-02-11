import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIKK')

input_file = "file:FABC2662-9AC8-E711-BF94-02163E019BB9.root "

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v11')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(input_file)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootuple-PsiTrakTrakRootupler.root'),
)

hltList = [
#Phi
'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi',
'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi',
'HLT_Mu20_TkMu0_Phi',
'HLT_Dimuon14_Phi_Barrel_Seagulls',
'HLT_Mu25_TkMu0_Phi',
'HLT_Dimuon24_Phi_noCorrL1',
#JPsi
'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
'HLT_DoubleMu4_JpsiTrk_Displaced',
'HLT_DoubleMu4_Jpsi_Displaced',
'HLT_DoubleMu4_3_Jpsi_Displaced',
'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
'HLT_Dimuon25_Jpsi',
'HLT_Dimuon0_Jpsi'
]

hltpaths = cms.vstring(hltList)

hltpathsV = cms.vstring([h + '_v*' for h in hltList])

filters = cms.vstring(
                                #PHI TRIGGERS FILTER
                                #HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi
                                'hltDiMuonGlbOrTrkFiltered0v2', #Phi
                                #'hltDiMuonGlbOrTrk0zFiltered0p2v2',
                                'hltDoubleMu2JpsiL3Filtered', ##JPsi
                                #'hltMumuVtxProducerDoubleMu2Jpsi',
                                'hltMumuFilterDoubleMu2Jpsi',
                                #HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi
                                #'hltDoubleMu2JpsiDoubleTrkL3Filtered',
                                #'hltDoubleTrkmumuVtxProducerDoubleMu2Jpsi',
                                'hltDoubleTrkmumuFilterDoubleMu2Jpsi',
                                #'hltJpsiTkAllConeTracksIterDoubleTrk',
                                #'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v2',
                                'hltJpsiTkTkVertexFilterPhiDoubleTrk1v2',
                                #HLT_Mu20_TkMu0_Phi
                                #'hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered20',
                                #'hltDiMuonGlbFiltered20TrkFiltered0',
                                'hltDiMuonGlb20Trk0DzFiltered0p2',
                                #HLT_Dimuon14_Phi_Barrel_Seagulls
                                #'hltDimuon14PhiBarrelnoCowL3Filtered',
                                #'hltDisplacedmumuVtxProducerDimuon14PhiBarrelnoCow',
                                'hltDisplacedmumuFilterDimuon14PhiBarrelnoCow',
                                #HLT_Mu25_TkMu0_Phi
                                #'hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered20',
                                #'hltDiMuonGlbFiltered25TrkFiltered0',
                                'hltDiMuonGlb25Trk0DzFiltered0p2',
                                #HLT_Dimuon24_Phi_noCorrL1
                                #'hltDisplacedmumuFilterDimuon24PhiBarrelNoCorrL1',
                                #'hltDisplacedmumuVtxProducerDimuon24PhiNoCorrL1',
                                'hltDimuon24PhiNoCorrL1L3fL3Filtered',
                                #JPSI Trigger Filters
                                #HLT_DoubleMu4_JpsiTrkTrk_Displaced_v4
                                #'hltDoubleMu4JpsiDisplacedL3Filtered'
                                #'hltDisplacedmumuVtxProducerDoubleMu4Jpsi',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                #'hltJpsiTkAllConeTracksIter',
                                #'hltJpsiTrkTrkVertexProducerPhiKstar',
                                #'hltJpsiTkTkVertexFilterPhiKstar',
                                #HLT_DoubleMu4_JpsiTrk_Displaced_v12
                                #'hltDoubleMu4JpsiDisplacedL3Filtered',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                #'hltJpsiTkVertexProducer',
                                #'hltJpsiTkVertexFilter',
                                #HLT_DoubleMu4_Jpsi_Displaced
                                #'hltDoubleMu4JpsiDisplacedL3Filtered',
                                #'hltDisplacedmumuVtxProducerDoubleMu4Jpsi',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                #HLT_DoubleMu4_3_Jpsi_Displaced
                                #'hltDoubleMu43JpsiDisplacedL3Filtered',
                                #'hltDisplacedmumuVtxProducerDoubleMu43Jpsi',
                                'hltDisplacedmumuFilterDoubleMu43Jpsi',
                                #HLT_Dimuon20_Jpsi_Barrel_Seagulls
                                #'hltDimuon20JpsiBarrelnoCowL3Filtered',
                                #'hltDisplacedmumuVtxProducerDimuon20JpsiBarrelnoCow',
                                'hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow',
                                #HLT_Dimuon25_Jpsi
                                'hltDisplacedmumuFilterDimuon25Jpsis'
                                #HLT_Dimuon0_Jpsi
                                #'hltDimuon0JpsiL3Filtered',
                                #'hltDisplacedmumuVtxProducerDimuon0Jpsi',
                                'hltDisplacedmumuFilterDimuon0Jpsi'
                                )

# process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
# process.CandidateSelectedTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
#                 src=cms.InputTag("oniaSelectedTracks"),
#                 particleType=cms.string('K+')
#                 )
#
# from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
# process.patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.JPsi2MuMuPAT = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('slimmedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.9 < mass && mass < 3.3 && charge==0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.PsiPhiProducer = cms.EDProducer('OniaPFPFProducer',
    Onia = cms.InputTag('JPsi2MuMuPAT'),
    PFCandidates = cms.InputTag('packedPFCandidates'),
    OniaMassCuts = cms.vdouble(2.946916,3.246916),      # J/psi mass window 3.096916 +/- 0.150
    TrakTrakMassCuts = cms.vdouble(0.919461,1.119461),  # phi mass window 1.019461 +/- .015
    OniaTrakTrakMassCuts = cms.vdouble(4.0,6.0),            # b-hadron mass window
    MassTraks = cms.vdouble(0.493677,0.493677),         # traks masses
    OnlyBest  = cms.bool(False)
)

process.PsiPhiFitter = cms.EDProducer('PsiPFPFKinematicFit',
    PsiPFPF     = cms.InputTag('PsiPhiProducer','OniaPFPFCandidates'),
    mass_constraint = cms.double(3.096916),              # J/psi mass in GeV
    OniaTrakTrakMassCuts = cms.vdouble(4.0,6.0),            # b-hadron mass window
    MassTraks = cms.vdouble(0.493677,0.493677),         # traks masses
    product_name    = cms.string('PsiPhiCandidates')
)

# process.rootuple = cms.EDAnalyzer('PsiTrakTrakRootupler',
#     jpsitrktrk_cand = cms.InputTag('PsiPhiProducer','OniaTrakTrakCandidates'),
#     jpsitrktrk_rf_cand = cms.InputTag("PsiPhiFitter","PsiPhiCandidates"),
#     beamSpotTag = cms.InputTag("offlineBeamSpot"),
#     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
#     TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#     isMC = cms.bool(False),
#     OnlyBest = cms.bool(False),
#     HLTs = hltpaths,
#     filters = filters
# )

# process.Phi2KKPAT = cms.EDProducer('Phi2KKPAT',
#   kaons = cms.InputTag("patSelectedTracks"),
#   beamSpotTag = cms.InputTag("offlineBeamSpot"),
#   primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
#   OniaTag = cms.InputTag("onia2MuMuPAT"),                      ## Use Onia2MuMu as seed for PV, only tracks in this PV are used, PV=0 is used otherwise
#   higherPuritySelection = cms.string(""),                      ## At least one kaon must pass this selection
#   lowerPuritySelection  = cms.string(""),                      ## BOTH kaons must pass this selection
#   dikaonSelection  = cms.string("0.85 < mass && mass < 1.2 && charge==0 && userFloat('deltar') < 0.7"),  ## The dikaon must pass this selection before vertexing
#   addCommonVertex = cms.bool(True),                            ## Embed the full reco::Vertex out of the common vertex fit
#   resolvePileUpAmbiguity = cms.bool(True)                      ## Order PVs by their vicinity to the Phi vertex, not by sumPt
# )

# process.rootupleKK = cms.EDAnalyzer('Phi2KKRootupler',
#                           dikaons = cms.InputTag("Phi2KKPAT"),
#                           kaons = cms.InputTag("patSelectedTracks"),
#                           primaryVertices = cms.InputTag("offlinePrimaryVertices"),
#                           TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
# 			              TestFilterNames =  filters,
#                           kk_mass_cuts = cms.vdouble(0.85,1.2),
#                           isMC = cms.bool(False),
#                           OnlyBest = cms.bool(False),
#                           OnlyGen = cms.bool(False)
#                           )

# process.rootupleMuMu = cms.EDAnalyzer('Onia2MuMuRootupler',
#                           dimuons = cms.InputTag("onia2MuMuPAT"),
#                           muons = cms.InputTag("replaceme"),
#                           primaryVertices = cms.InputTag("offlinePrimaryVertices"),
#                           TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#                           onia_pdgid = cms.uint32(443),
#                           onia_mass_cuts = cms.vdouble(2.5,3.5),
#                           isMC = cms.bool(False),
#                           OnlyBest = cms.bool(False),
#                           OnlyGen = cms.bool(False),
#                           HLTs = hltpaths
#                           )



process.p = cms.Path(process.triggerSelection * process.JPsi2MuMuPAT * process.PsiPhiProducer * process.PsiPhiFitter) # * process.rootuple * process.rootupleMuMu)# * process.Phi2KKPAT * process.patSelectedTracks *process.rootupleKK)
