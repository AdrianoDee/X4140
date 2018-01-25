import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIKK')

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputfile')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("CompactSkim.Examples.PsiTrakTrakRootupler_cfi")

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


process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

PsiTrakProducer = cms.EDProducer('OniaTrakProducer',
    Onia = cms.InputTag('onia2MuMuPAT'),
    Trak = cms.InputTag('cleanPatTrackCands'),
    OniaMassCuts = cms.vdouble(2.947,3.247),
    OniaTrakMassCuts = cms.vdouble(5.0,5.7),
    OnlyBest = cms.bool(True)
)

PsiTrakFitter = cms.EDProducer('PsiTrakKinematicFit',
    PsiTrak         = cms.InputTag('PsiTrakProducer','OniaTrakCandidates'),
    mass_constraint = cms.double(3.0969),              # J/psi mass in GeV
    product_name    = cms.string('PsiTrakCandidates')
)

rootuple = cms.EDAnalyzer('PsiTrakRootupler',
    oniat_cand = cms.InputTag("PsiTrakProducer","OniaTrakCandidates"),
    oniat_rf_cand = cms.InputTag("PsiTrakFitter","PsiTrakCandidates"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    oniat_pdgid = cms.uint32(521),
    onia_pdgid = cms.uint32(443),
    trak_pdgid = cms.uint32(321),
    isMC = cms.bool(False),
    OnlyBest = cms.bool(True)
)

PsiTrakSequence = cms.Sequence(PsiTrakProducer*PsiTrakFitter*rootuple)

process.p = cms.Path(process.PsiTrakTrakSequence)
