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
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v4', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.load("mmkk.mmkk.slimmedMuonsTriggerMatcher2017_cfi")

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(
                                                                        'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v*',
                                                                        # 'HLT_Mu20_TkMu0_Phi_v*',
                                                                        # 'HLT_Dimuon14_Phi_Barrel_Seagulls_v*',
                                                                        # 'HLT_Mu25_TkMu0_Phi_v*',
                                                                        # 'HLT_Mu20_TkMu0_Phi_v*',
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

#make patTracks
# from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
# makeTrackCandidates(process,
#                        label        = 'KaonPCands',                  # output collection
#                        tracks       = cms.InputTag('generalTracks'), # input track collection
#                        particleType = 'k+',                           # particle type (for assigning a mass)
#                        preselection = 'pt > 0.1',                    # preselection cut on candidates
#                        selection    = 'pt > 0.1',                    # selection on PAT Layer 1 objects
#                        isolation    = {},                            # isolations to use (set to {} for None)
#                        isoDeposits  = [],
#                        mcAs         = None                           # replicate MC match as the one used for Muons
#    )
#
from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
makeTrackCandidates(process,
                       label        = 'KaonTracks',                  # output collection
                       tracks       = cms.InputTag('generalTracks'), # input track collection
                       particleType = 'K+',                           # particle type (for assigning a mass)
                       preselection = 'pt > 0.2',                    # preselection cut on candidates
                       selection    = 'pt > 0.4',                    # selection on PAT Layer 1 objects
                       isolation    = {},                            # isolations to use (set to {} for None)
                       isoDeposits  = [],
                       mcAs         = None                           # replicate MC match as the one used for Muons
   )

l1cands = getattr(process, 'cleanPatKaonTracks')

process.cleanPatKaonTracks.embedTrack = cms.bool(True)

process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuonsWithTrigger'),
   cut = cms.string(' ( abs(eta) <= 3 )'
   ),
   filter = cms.bool(True)
)

filters = cms.vstring(          #HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi
                                'hltDiMuonGlbOrTrkFiltered0v2',
                                'hltDiMuonGlbOrTrk0zFiltered0p2v2',
                                'hltDoubleMu2JpsiL3Filtered',
                                'hltMumuVtxProducerDoubleMu2Jpsi',
                                'hltMumuFilterDoubleMu2Jpsi',
                                #HLT_Dimuon14_Phi_Barrel_Seagulls
                                'hltDimuon14PhiBarrelnoCowL3Filtered',
                                'hltDisplacedmumuVtxProducerDimuon14PhiBarrelnoCow',
                                'hltDisplacedmumuFilterDimuon14PhiBarrelnoCow',
                                #HLT_Mu25_TkMu0_Phi
                                'hltDimuon14PhiBarrelnoCowL3Filtered',
                                'hltDisplacedmumuVtxProducerDimuon14PhiBarrelnoCow',
                                'hltDisplacedmumuFilterDimuon14PhiBarrelnoCow',
                                #HLT_Mu20_TkMu0_Phi
                                'hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered20',
                                'hltDiMuonGlbFiltered20TrkFiltered0',
                                'hltDiMuonGlb20Trk0DzFiltered0p2',

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
      dimuonSelection     = cms.string("2.9 < mass && mass < 3.3 && userFloat('vProb') > 0.0 "),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.FourOnia2KKPhi = cms.EDProducer('FourOnia2KKPAT',
    tracks                      = cms.InputTag('cleanPatKaonTracks'),
    primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
    higherPuritySelection       = cms.string(""),
    lowerPuritySelection        = cms.string(""),
    dimuonSelection             = cms.string("0.5 < mass && mass < 1.5 && charge==0 "),
    addCommonVertex             = cms.bool(True),
    addMuonlessPrimaryVertex    = cms.bool(False),
    addMCTruth                  = cms.bool(False),
    resolvePileUpAmbiguity      = cms.bool(True),
    HLTFilters                  = filters
)

# process.xFitter = cms.EDProducer('FourOniaKinFit',
#                           x_cand = cms.InputTag("xProducer"),
#                           x_mass = cms.double(5.36679), # GeV   1S = 9.46030   2S = 10.02326    3S = 10.35520  J/psi=3.0969
#                           product_name = cms.string("xCand"),
#                           pdgID = cms.int32(531)
#                          )


process.xCandSequence = cms.Sequence(
                   process.triggerSelection *
                   process.slimmedMuonsWithTriggerSequence *
				   process.oniaSelectedMuons *
				   process.FourOnia2MuMuJPsi *
                   process.Onia2MuMuFilteredJpsi *
                   #process.DiMuonCounterJPsi *
                   process.content *
                   process.FourOnia2KKPhi
                   #process.xFitter
				   )

# process.rootuple = cms.EDAnalyzer('x4MuRootupler',
#                           phidimuons = cms.InputTag("Onia2MuMuFilteredPhi"),
#                           jpsidimuons = cms.InputTag("Onia2MuMuFilteredJpsi"),
#                           #chi_cand = cms.InputTag("chiProducer"),
# 			              x_cand = cms.InputTag("xProducer"),
#                           #xrefit = cms.InputTag("xFitter","xCand"),
# 			              # refit2S  = cms.InputTag("chiFitter2S","y2S"),
# 			              # refit3S  = cms.InputTag("chiFitter3S","y3S"),
#                           primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                           TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
#                           isMC = cms.bool(False)
#                          )
process.p = cms.Path(process.triggerSelection * process.xCandSequence)# * process.rootuple)
