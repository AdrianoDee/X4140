#input_filename = '/store/data/Run2017B/MuOnia/MINIAOD/PromptReco-v1/000/297/723/00000/9040368C-DE5E-E711-ACFF-02163E0134FF.root'
ouput_filename = 'rootuple.root'
input_filename = 'file:FABC2662-9AC8-E711-BF94-02163E019BB9.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v11', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.load("mmkk.mmkk.slimmedMuonsTriggerMatcher2017_cfi")
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")

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


process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.JPsi2MuMuPAT = cms.EDProducer('Onia2MuMuPAT',
  muons = cms.InputTag("patMuons"),
  beamSpotTag = cms.InputTag("offlineBeamSpot"),
  primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
  higherPuritySelection = cms.string(""), ## At least one muon must pass this selection
  lowerPuritySelection  = cms.string(""), ## BOTH muons must pass this selection
  dimuonSelection  = cms.string("2.9 < mass < 3.2"), ## The dimuon must pass this selection before vertexing
  addCommonVertex = cms.bool(True), ## Embed the full reco::Vertex out of the common vertex fit
  addMuonlessPrimaryVertex = cms.bool(True), ## Embed the primary vertex re-made from all the tracks except the two muons
  addMCTruth = cms.bool(True),      ## Add the common MC mother of the two muons, if any
  resolvePileUpAmbiguity = cms.bool(True)   ## Order PVs by their vicinity to the J/psi vertex, not by sumPt
)

process.Phi2MuMuPAT = cms.EDProducer('Onia2MuMuPAT',
  muons = cms.InputTag("patMuons"),
  beamSpotTag = cms.InputTag("offlineBeamSpot"),
  primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
  higherPuritySelection = cms.string(""), ## At least one muon must pass this selection
  lowerPuritySelection  = cms.string(""), ## BOTH muons must pass this selection
  dimuonSelection  = cms.string("0.55 < mass < 1.25"), ## The dimuon must pass this selection before vertexing
  addCommonVertex = cms.bool(True), ## Embed the full reco::Vertex out of the common vertex fit
  addMuonlessPrimaryVertex = cms.bool(True), ## Embed the primary vertex re-made from all the tracks except the two muons
  addMCTruth = cms.bool(True),      ## Add the common MC mother of the two muons, if any
  resolvePileUpAmbiguity = cms.bool(True)   ## Order PVs by their vicinity to the J/psi vertex, not by sumPt
)

process.Onia2MuMuFilteredJpsi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("JPsi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.9 < mass && mass < 3.3"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.Onia2MuMuFilteredPhi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("Phi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("0.6 < mass && mass < 1.2"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters

)

process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFilteredJpsi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.DiMuonCounterPhi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFilteredPhi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.PsiPhiProducer = cms.EDProducer('PsiPhiFourMuonsProducer',
    PsiCollection_ = cms.InputTag('onia2MuMuPAT'),
    PhiCollection_ = cms.InputTag('onia2MuMuPAT'),
    JPsiMassCuts_ = cms.vdouble(2.946916,3.246916),      # J/psi mass window 3.096916 +/- 0.150
    PhiMassCuts_ = cms.vdouble(1.004461,1.034461),  # phi mass window 1.019461 +/- .015
    OniaTrakTrakMassCuts = cms.vdouble(4.0,6.0),            # b-hadron mass window
    OnlyBest  = cms.bool(False)
)

process.xCandSequence = cms.Sequence(
                   process.triggerSelection *
                   process.slimmedMuonsWithTriggerSequence *
				   process.oniaSelectedMuons *
                   process.Onia2MuMuFilteredPhi *
                   process.Onia2MuMuFilteredPhi *
                  # process.DiMuonCounterPhi *
                   process.Onia2MuMuFilteredJpsi *
                  # process.DiMuonCounterJPsi *
                   process.xProducer*
                   #process.BkgProducer
				   )

process.rootuple = cms.EDAnalyzer('x4MuRootupler',
                          phidimuons = cms.InputTag("FourOnia2MuMuPhi"),
                          jpsidimuons = cms.InputTag("FourOnia2MuMuJPsi"),
                          HLTs = hltpaths,
			              x_cand = cms.InputTag("xProducer"),
                          bkg_cand = cms.InputTag("BkgProducer"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC = cms.bool(False)
                         )
process.p = cms.Path(process.xCandSequence * process.rootuple)