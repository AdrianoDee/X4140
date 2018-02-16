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
process.load("mmkk.mmkk.slimmedTracksTriggerMatcher2017_cfi")

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
'HLT_Dimuon0_Jpsi',
'HLT_Dimuon16_Jpsi', #v4.2 2016
'HLT_Dimuon20_Jpsi' #v4.2 2016
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

process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuonsWithTrigger'),
   cut = cms.string(''),
   filter = cms.bool(True)
)

process.FourOnia2MuMuPhi = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('oniaSelectedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("0.6 < mass && mass < 1.2 && charge==0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.FourOnia2MuMuJPsi = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('oniaSelectedMuons'),
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

process.Onia2MuMuFilteredJpsi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("FourOnia2MuMuJPsi"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.9 < mass && mass < 3.3"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.Onia2MuMuFilteredPhi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("FourOnia2MuMuPhi"),
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

process.xProducer = cms.EDProducer('FourOniaProducer',
    phidimuons                  = cms.InputTag("FourOnia2MuMuPhi"),
    jpsidimuons                 = cms.InputTag("FourOnia2MuMuJPsi"),
    primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
    dzmax                       = cms.double(20.0),
    triggerMatch                = cms.bool(False),
    addCommonVertex             = cms.bool(True),
    addMuonlessPrimaryVertex    = cms.bool(True),
    resolvePileUpAmbiguity      = cms.bool(True),
    quadmuonSelection           = cms.string("4.0 < mass && mass < 6.0 && charge==0")
)

process.BkgOnia2MuMuPhi = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('oniaSelectedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("0.6 < mass && mass < 1.2 && charge != 0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.BkgOnia2MuMuJPsi = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('oniaSelectedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.9 < mass && mass < 3.3 && charge != 0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.BkgProducer = cms.EDProducer('FourOniaProducer',
    phidimuons                  = cms.InputTag("BkgOnia2MuMuPhi"),
    jpsidimuons                 = cms.InputTag("BkgOnia2MuMuJPsi"),
    primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
    dzmax                       = cms.double(20.0),
    triggerMatch                = cms.bool(False),
    addCommonVertex             = cms.bool(True),
    addMuonlessPrimaryVertex    = cms.bool(True),
    resolvePileUpAmbiguity      = cms.bool(True),
    quadmuonSelection           = cms.string("4.0 < mass && mass < 6.0")
)


process.xCandSequence = cms.Sequence(
                   process.triggerSelection *
                   proces.slimmedPFCandsWithTriggerSequence *
                   process.slimmedMuonsWithTriggerSequence *
				   process.oniaSelectedMuons *
                   process.FourOnia2MuMuPhi *
                  # process.Onia2MuMuFilteredPhi *
                  # process.DiMuonCounterPhi *
				   process.FourOnia2MuMuJPsi *
                   process.BkgOnia2MuMuPhi *
                   process.BkgOnia2MuMuJPsi *
                  # process.Onia2MuMuFilteredJpsi *
                  # process.DiMuonCounterJPsi *
                   process.xProducer*
                   process.BkgProducer
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
