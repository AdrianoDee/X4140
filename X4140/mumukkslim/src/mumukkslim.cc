// -*- C++ -*-
//
// Package:    mumukk
// Class:      mumukk
//
/**\class mumukk mumukk.cc myAnalyzers/mumukk/src/mumukk.cc

Description: <one line class summary>
Make rootTuple for JPsiKK reconstruction

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//


/// system include files
#include <memory>

/// user include files
#include "../interface/mumukkslim.h"
#include "../interface/VertexReProducer.h"
//#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

/// for 53x
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "TMath.h"
#include "Math/VectorUtil.h"

/// useless so far
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "HepMC/GenVertex.h"
//#include <HepMC/GenVertex.h>
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

///
/// constants, enums and typedefs
///

typedef math::Error<3>::type CovarianceMatrix;

const ParticleMass muon_mass = 0.10565837; //pdg mass
const ParticleMass kaon_mass = 0.493667; //pdg mass
ParticleMass JPsi_mass = 3.096916;
const ParticleMass Phi_mass = 1.0194;


/// Setting insignificant mass sigma to avoid singularities in the covariance matrix.
float small_sigma = muon_mass*1.e-6;
//float small_sigma = kaon_mass*1.e-6; /// SEMRA

///
/// static data member definitions
///

///
/// constructors and destructor
///
mumukk::mumukk(const edm::ParameterSet& iConfig) :
hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN",edm::InputTag("genParticles"))),
vtxSample(iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices"))),

doData( iConfig.getUntrackedParameter<bool>("DoDataAnalysis", true) ),
doMC( iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", true) ),
MCParticle( iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443) ), /// 20443 X, 100443 Psi(2S), 9120443 X from B / decide later for X(4140)
MCExclusiveDecay( iConfig.getUntrackedParameter<bool>("MonteCarloExclusiveDecay", true) ),
MCMother( iConfig.getUntrackedParameter<int>("MonteCarloMotherId", 511) ), /// 511 B0 (=anti-B0), 531 B0 / decide later MCMotherId for X(4140)
MCDaughtersN( iConfig.getUntrackedParameter<int>(" MonteCarloDaughtersN", 3) ), /// will be same
doMuMuMassConst( iConfig.getUntrackedParameter<bool>("DoMuMuMassConstraint", true) ),
skipJPsi(iConfig.getUntrackedParameter<bool>("SkipJPsi", false)),

MuMinPixHits(iConfig.getUntrackedParameter<int>("MinNumMuPixHits", 0)),
MuMinSiHits(iConfig.getUntrackedParameter<int>("MinNumMuSiHits", 0)),
MuMaxNormChi(iConfig.getUntrackedParameter<double>("MaxMuNormChi2", 1000)),
MuMaxD0(iConfig.getUntrackedParameter<double>("MaxMuD0", 1000)),
sharedFraction(iConfig.getUntrackedParameter<double>("sharedFraction", 0.5)),

TrMinSiHits(iConfig.getUntrackedParameter<int>("MinNumTrSiHits", 0)),
TrMinPt(iConfig.getUntrackedParameter<double>("MinTrPt", 0)),
TrMaxNormChi2(iConfig.getUntrackedParameter<double>("MaxTrChi2NDF", 10)),
TriggersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
FiltersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching")),
resolveAmbiguity_(iConfig.getUntrackedParameter<bool>("resolvePileUpAmbiguity",true)),
addMuMulessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addMuMulessPrimaryVertex", true)),

JPsiMinMass(iConfig.getUntrackedParameter<double>("MinJPsiMass", 2.8)),
JPsiMaxMass(iConfig.getUntrackedParameter<double>("MaxJPsiMass", 3.4)),
PhiMinMass(iConfig.getUntrackedParameter<double>("MinPhiMass", 0.97)),
PhiMaxMass(iConfig.getUntrackedParameter<double>("MaxPhiMass", 1.07)),
JPsiPhiMaxXMass(iConfig.getUntrackedParameter<double>("MaxJPsiPhiXMass", 4.35)),
JPsiPhiMinB0Mass(iConfig.getUntrackedParameter<double>("MinJPsiPhiB0Mass", 5.1)),
JPsiPhiMaxB0Mass(iConfig.getUntrackedParameter<double>("MaxJPsiPhiB0Mass", 5.6)),
MuMuTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxMuMuTrackDR", 1)),

XTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxXCandTrackDR", 1.1)),
UseXDR(iConfig.getUntrackedParameter<bool>("UseXDR", false)),
MuMuKKMinB0Mass(iConfig.getUntrackedParameter<double>("MinMuMuKKB0Mass", 0)),
MuMuKKMaxB0Mass(iConfig.getUntrackedParameter<double>("MaxMuMuKKB0Mass", 10)),
MuMuKKMaxXMass(iConfig.getUntrackedParameter<double>("MaxMuMuKKXMass", 10)),
addXlessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addXlessPrimaryVertex", true)),

Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",true)),
DeDxEstimator_(iConfig.getUntrackedParameter<std::string>("DeDxEstimator", std::string("dedxHarmonic2"))),
m_dEdxDiscrimTag(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),
m_dEdxDiscrimTag_kaon(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),

xTree(0),
runNum(0), evtNum(0), lumiNum(0),
trigRes(0), trigNames(0), L1TT(0), MatchTriggerNames(0),
/// counters for X
nMu(0), nMuMu(0), nX(0), nKK(0),
nX_pre0(0), nX_pre1(0), nX_pre2(0), nX_pre3(0), nX_pre4(0), nX_pre5(0), nX_pre6(0), nX_pre7(0), nX_pre8(0), nX_pre9(0), nX_pre10(0), nX_pre11(0), nX_pre12(0), nX_pre13(0), nX_pre14(0), nX_pre15(0),

priVtx_n(0), priVtx_X(0), priVtx_Y(0), priVtx_Z(0), priVtx_XE(0), priVtx_YE(0), priVtx_ZE(0), priVtx_NormChi2(0), priVtx_Chi2(0), priVtx_CL(0), priVtx_tracks(0), priVtx_tracksPtSq(0),
/// indices
mu1Idx(0), mu2Idx(0), MuMuType(0), ka1Idx(0), ka2Idx(0),
X_MuMuIdx(0), X_ka1Idx(0), X_ka2Idx(0),
/// MC Analysis /// n_B0Ancestors & no for X
n_genEvtVtx(0), genEvtVtx_X(0), genEvtVtx_Y(0), genEvtVtx_Z(0), genEvtVtx_particles(0), n_XAncestors(0),
nMCAll(0), nMCX(0), /*nMCXVtx(0),*/ MCPdgIdAll(0), MCDanNumAll(0),
// Gen Primary Vertex
PriVtxGen_X(0), PriVtxGen_Y(0), PriVtxGen_Z(0), PriVtxGen_EX(0), PriVtxGen_EY(0), PriVtxGen_EZ(0),
PriVtxGen_Chi2(0), PriVtxGen_CL(0), PriVtxGen_Ndof(0), PriVtxGen_tracks(0),
MCJPsiPx(0), MCJPsiPy(0), MCJPsiPz(0),
MCmupPx(0), MCmupPy(0), MCmupPz(0),
MCmumPx(0), MCmumPy(0), MCmumPz(0),
MCPhiPx(0), MCPhiPy(0), MCPhiPz(0),
MCkpPx(0), MCkpPy(0), MCkpPz(0),
MCkmPx(0), MCkmPy(0), MCkmPz(0),
//MCpionPx(0), MCpionPy(0), MCpionPz(0),
//MCkaonPx(0), MCkaonPy(0), MCkaonPz(0),
//MCpionCh(0), MCkaonCh(0),
MCPx(0), MCPy(0), MCPz(0),
/// generic muons
muPx(0), muPy(0), muPz(0), muCharge(0),
muPhits(0), muShits(0), muLayersTr(0), muLayersPix(0),
muD0(0),  muD0E(0), muDz(0), muChi2(0), muNDF(0),
mufHits(0), muFirstBarrel(0), muFirstEndCap(0),
muDzVtx(0), muDxyVtx(0), muDzVtxErr(0), muKey(0),
muIsGlobal(0), muIsPF(0),
muGlMuHits(0), muGlChi2(0), muGlNDF(0), muGlMatchedStation(0),
muGlDzVtx(0), muGlDxyVtx(0),
nMatchedStations(0),
muType(0), muQual(0), muTrack(0), muNOverlap(0), muNSharingSegWith(0),
/// generic tracks
trNotRef(0), trRef(0),
trPx(0), trPy(0), trPz(0), trE(0),
trNDF(0), trPhits(0), trShits(0), trChi2(0),
trD0(0), trD0E(0), trCharge(0),
trfHits(0), trFirstBarrel(0), trFirstEndCap(0),
trDzVtx(0), trDxyVtx(0),
trQualityHighPurity(0), trQualityTight(0),
tr_nsigdedx(0), tr_dedx(0), tr_dedxMass(0), tr_theo(0), tr_sigma(0),
tr_dedx_byHits(0), tr_dedxErr_byHits(0), tr_saturMeas_byHits(0), tr_Meas_byHits(0),
/// MuMu cand & KK cand
MuMuMass(0), MuMuPx(0), MuMuPy(0), MuMuPz(0),
MuMuVtx_CL(0), MuMuVtx_Chi2(0),
MuMuDecayVtx_X(0), MuMuDecayVtx_Y(0), MuMuDecayVtx_Z(0),
MuMuDecayVtx_XE(0), MuMuDecayVtx_YE(0), MuMuDecayVtx_ZE(0),
MuMuMuonTrigMatch(0),
KKMass(0), KKPx(0), KKPy(0), KKPz(0),
KKVtx_CL(0), KKVtx_Chi2(0),
KKDecayVtx_X(0), KKDecayVtx_Y(0), KKDecayVtx_Z(0),
KKDecayVtx_XE(0), KKDecayVtx_YE(0), KKDecayVtx_ZE(0),
/// muons after JPsi (MuMu) fit &kaons after Phi (KK) fit
mu1_MuMu_Px(0), mu1_MuMu_Py(0), mu1_MuMu_Pz(0), mu1_MuMu_Chi2(0), mu1_MuMu_NDF(0),
mu2_MuMu_Px(0), mu2_MuMu_Py(0), mu2_MuMu_Pz(0), mu2_MuMu_Chi2(0), mu2_MuMu_NDF(0),
ka1_KK_Px(0), ka1_KK_Py(0), ka1_KK_Pz(0), ka1_KK_Chi2(0), ka1_KK_NDF(0),
ka2_KK_Px(0), ka2_KK_Py(0), ka2_KK_Pz(0), ka2_KK_Chi2(0), ka2_KK_NDF(0),
DR_MuMu_K1(0), DR_MuMu_K2(0), DR_MuMuKK_K1(0), DR_MuMuKK_K2(0),
/// Primary Vertex with "MuMu correction"
PriVtxMuMuCorr_n(0),
PriVtxMuMuCorr_X(0), PriVtxMuMuCorr_Y(0), PriVtxMuMuCorr_Z(0), PriVtxMuMuCorr_EX(0), PriVtxMuMuCorr_EY(0), PriVtxMuMuCorr_EZ(0),
PriVtxMuMuCorr_Chi2(0), PriVtxMuMuCorr_CL(0), PriVtxMuMuCorr_tracks(0),
nTrk(0),
/// X candidates
xMass(0), xVtx_CL(0), xVtx_Chi2(0),
xPx(0), xPy(0), xPz(0), xPxE(0), xPyE(0), xPzE(0),
xDecayVtx_X(0), xDecayVtx_Y(0), xDecayVtx_Z(0), xDecayVtx_XE(0), xDecayVtx_YE(0), xDecayVtx_ZE(0),
/// Muons and tracks after X candidates fit
mu1Px_MuMuKK(0), mu1Py_MuMuKK(0), mu1Pz_MuMuKK(0), mu1E_MuMuKK(0),
mu2Px_MuMuKK(0), mu2Py_MuMuKK(0), mu2Pz_MuMuKK(0), mu2E_MuMuKK(0),
k1Px_MuMuKK(0), k1Py_MuMuKK(0), k1Pz_MuMuKK(0), k1E_MuMuKK(0),
kaon1_nsigdedx(0), kaon1_dedx(0), kaon1_dedxMass(0), kaon1_theo(0), kaon1_sigma(0),
kaon1_dedx_byHits(0), kaon1_dedxErr_byHits(0), kaon1_saturMeas_byHits(0), kaon1_Meas_byHits(0),
k2Px_MuMuKK(0), k2Py_MuMuKK(0), k2Pz_MuMuKK(0), k2E_MuMuKK(0),
kaon2_nsigdedx(0), kaon2_dedx(0), kaon2_dedxMass(0), kaon2_theo(0), kaon2_sigma(0),
kaon2_dedx_byHits(0), kaon2_dedxErr_byHits(0), kaon2_saturMeas_byHits(0), kaon2_Meas_byHits(0),
/// Primary Vertex with largest B0_cos(alpha) no less values for X
PriVtx_XCosAlpha_n(0),
PriVtx_XCosAlpha_X(0), PriVtx_XCosAlpha_Y(0), PriVtx_XCosAlpha_Z(0), PriVtx_XCosAlpha_EX(0), PriVtx_XCosAlpha_EY(0), PriVtx_XCosAlpha_EZ(0),
PriVtx_XCosAlpha_Chi2(0), PriVtx_XCosAlpha_CL(0), PriVtx_XCosAlpha_tracks(0),
PriVtx_XCosAlpha3D_n(0),
PriVtx_XCosAlpha3D_X(0), PriVtx_XCosAlpha3D_Y(0), PriVtx_XCosAlpha3D_Z(0), PriVtx_XCosAlpha3D_EX(0), PriVtx_XCosAlpha3D_EY(0), PriVtx_XCosAlpha3D_EZ(0),
PriVtx_XCosAlpha3D_Chi2(0), PriVtx_XCosAlpha3D_CL(0), PriVtx_XCosAlpha3D_tracks(0),
XLessPV_tracksPtSq(0), XLessPV_4tracksPtSq(0),
PriVtxXLess_n(0),
PriVtxXLess_X(0), PriVtxXLess_Y(0), PriVtxXLess_Z(0), PriVtxXLess_EX(0), PriVtxXLess_EY(0), PriVtxXLess_EZ(0),
PriVtxXLess_Chi2(0), PriVtxXLess_CL(0), PriVtxXLess_tracks(0),
PriVtxXLess_XCosAlpha_n(0),
PriVtxXLess_XCosAlpha_X(0), PriVtxXLess_XCosAlpha_Y(0), PriVtxXLess_XCosAlpha_Z(0), PriVtxXLess_XCosAlpha_EX(0), PriVtxXLess_XCosAlpha_EY(0), PriVtxXLess_XCosAlpha_EZ(0),
PriVtxXLess_XCosAlpha_Chi2(0), PriVtxXLess_XCosAlpha_CL(0), PriVtxXLess_XCosAlpha_tracks(0),
PriVtxXLess_XCosAlpha3D_n(0),
PriVtxXLess_XCosAlpha3D_X(0), PriVtxXLess_XCosAlpha3D_Y(0), PriVtxXLess_XCosAlpha3D_Z(0), PriVtxXLess_XCosAlpha3D_EX(0), PriVtxXLess_XCosAlpha3D_EY(0), PriVtxXLess_XCosAlpha3D_EZ(0),
PriVtxXLess_XCosAlpha3D_Chi2(0), PriVtxXLess_XCosAlpha3D_CL(0), PriVtxXLess_XCosAlpha3D_tracks(0),
/// Primary Vertex with "B0 correction"
PriVtxXCorr_n(0),
PriVtxXCorr_X(0), PriVtxXCorr_Y(0), PriVtxXCorr_Z(0), PriVtxXCorr_EX(0), PriVtxXCorr_EY(0), PriVtxXCorr_EZ(0),
PriVtxXCorr_Chi2(0), PriVtxXCorr_CL(0), PriVtxXCorr_tracks(0),
/// Lifetime variables for B0
xCosAlphaBS(0), xCosAlpha3DBS(0), xCTauBS(0), xCTauBSE(0), xLxyBS(0), xLxyBSE(0), xLxyzBS(0), xLxyzBSE(0),
xCosAlphaPV(0), xCosAlpha3DPV(0), xCTauPV(0), xCTauPVE(0), xLxyPV(0), xLxyPVE(0), xLxyzPV(0), xLxyzPVE(0),
xCosAlphaPVCosAlpha(0), xCosAlpha3DPVCosAlpha(0), xCTauPVCosAlpha(0), xCTauPVCosAlphaE(0), xLxyPVCosAlpha(0), xLxyPVCosAlphaE(0), xLxyzPVCosAlpha(0), xLxyzPVCosAlphaE(0),
xCosAlphaPVCosAlpha3D(0), xCosAlpha3DPVCosAlpha3D(0), xCTauPVCosAlpha3D(0), xCTauPVCosAlpha3DE(0), xLxyPVCosAlpha3D(0), xLxyPVCosAlpha3DE(0), xLxyzPVCosAlpha3D(0), xLxyzPVCosAlpha3DE(0),
xCosAlphaXLessPV(0), xCosAlpha3DXLessPV(0),xCTauXLessPV(0),
xCTauXLessPVE(0), xLxyXLessPV(0), xLxyXLessPVE(0), xLxyzXLessPV(0), xLxyzXLessPVE(0),
xCosAlphaXLessPVCosAlpha(0), xCosAlpha3DXLessPVCosAlpha(0), xCTauXLessPVCosAlpha(0), xCTauXLessPVCosAlphaE(0), xLxyXLessPVCosAlpha(0), xLxyXLessPVCosAlphaE(0), xLxyzXLessPVCosAlpha(0), xLxyzXLessPVCosAlphaE(0),
xCosAlphaXLessPVCosAlpha3D(0), xCosAlpha3DXLessPVCosAlpha3D(0), xCTauXLessPVCosAlpha3D(0), xCTauXLessPVCosAlpha3DE(0), xLxyXLessPVCosAlpha3D(0), xLxyXLessPVCosAlpha3DE(0), xLxyzXLessPVCosAlpha3D(0), xLxyzXLessPVCosAlpha3DE(0),
xCosAlphaPVX(0), xCTauPVX(0), xCTauPVXE(0), xLxyPVX(0), xLxyzPVX(0), xLxyzPVXE(0),
xCTauPVX_3D(0), xCTauPVX_3D_err(0),
/// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
kaon1_dxy_PV(0), kaon1_dz_PV(0), kaon2_dxy_PV(0), kaon2_dz_PV(0),
kaon1_dxy_BS(0), kaon1_dz_BS(0), kaon2_dxy_BS(0), kaon2_dz_BS(0),
kaon1_dxy_XLessPV(0), kaon1_dz_XLessPV(0), kaon2_dxy_XLessPV(0), kaon2_dz_XLessPV(0),
kaon1_dxyE(0), kaon1_dzE(0), kaon2_dxyE(0), kaon2_dzE(0),

Kaon1FromPV(0), Kaon2FromPV(0)

{
  /// now do what ever initialization is needed
  MuMuMinMass = JPsiMinMass;
  MuMuMaxMass = JPsiMaxMass;
  KKMinMass = PhiMinMass;
  KKMaxMass = PhiMaxMass;
  MuMuKKMinB0Mass = JPsiPhiMinB0Mass;
  MuMuKKMaxB0Mass = JPsiPhiMaxB0Mass;
  MuMuKKMaxXMass = JPsiPhiMaxXMass;
}

mumukk::~mumukk()
{
  /// do anything here that needs to be done at desctruction time
  /// (e.g. close files, deallocate resources etc.)

}


///
/// member functions
///

/// ------------ method called to for each event  ------------
void mumukk::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool decayChainOK = false;
  /// Event content information
  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiNum = iEvent.id().luminosityBlock();


  bool hasRequestedTrigger = false;
  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  /// HLT results
  map<string,int> HLTPreScaleMap;
  edm::Handle<edm::TriggerResults> hltresults;
  try {
    iEvent.getByLabel(hlTriggerResults_, hltresults);
  }
  catch ( ... ) {
    cout << "Couldn't get handle on HLT Trigger!" << endl;
  }
  if (!hltresults.isValid()) {
    cout << "No Trigger Results!" << endl;
  }
  else {
    int ntrigs = hltresults->size();
    if (ntrigs==0){
      cout << "No trigger name given in TriggerResults of the input " << endl;
    }

    /// get hold of trigger names - based on TriggerResults object!
    edm::TriggerNames triggerNames_;
    triggerNames_ = iEvent.triggerNames(*hltresults);
    int ntriggers = TriggersForMatching_.size();
    // search for matched triggers
    //initialize
    for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) { // initialize MatchingTriggerResult array
      MatchingTriggerResult[MatchTrig] = 0;
    }

    //cycle on all triggers for the event
    for (int itrig = 0; itrig < ntrigs; itrig++) {

      //get the event hlt flag
      string trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();

      if (Debug_) if (hltflag) cout << trigName << " " <<hltflag <<endl;

      trigRes->push_back(hltflag);
      trigNames->push_back(trigName);

      int ntriggers = TriggersForMatching_.size();
      //cycle on trigger we want
      for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) {
        if (TriggersForMatching_[MatchTrig] == triggerNames_.triggerName(itrig)){

          MatchingTriggerResult[MatchTrig] = hltflag;
          //true if at least one trigger is on
          if (hltflag==1) hasRequestedTrigger = true;
          break;
        }
      }
    }
    for (int MatchTrig = 0; MatchTrig<ntriggers; MatchTrig++){
      if (Debug_) cout << TriggersForMatching_[MatchTrig]<<endl;
      MatchTriggerNames->push_back(TriggersForMatching_[MatchTrig]);
    }

    ///
    /// Get HLT map : triggername associated with its prescale, saved only for accepted trigger
    ///
    for (unsigned int i=0; i<triggerNames_.size(); i++){
      if ( hltresults->accept(i) ) { //  save trigger info only for accepted paths
        /// get the prescale from the HLTConfiguration, initialized at beginRun
        int prescale = hltConfig_.prescaleValue(iEvent,iSetup,triggerNames_.triggerNames().at(i));
        if (Debug_) std::cout<<" HLT===> "<<triggerNames_.triggerNames().at(i)<<" prescale ="<<prescale<<std::endl;
        HLTPreScaleMap[triggerNames_.triggerNames().at(i)] = prescale;
      }
    }
    HLTTrig = &HLTPreScaleMap; // store in the branch

  } /// end valid trigger
  ////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  /// L1 trigger infos
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
  const DecisionWord dWord = gtRecord->decisionWord();
  const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
  for(unsigned int l1i = 0; l1i != ttWord.size(); ++l1i){
    L1TT->push_back(ttWord.at(l1i));
  }

  ////////////////////////////////////////////////////////////////////////
  //Primary Vertex
  Vertex thePrimaryVtx, theBeamSpotVtx;
  math::XYZPoint RefVtx;
  Int_t thePrimaryVtx_multiplicity = -1 ;

  //Beam Spot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) {
    beamSpot = *beamSpotHandle;
    theBeamSpotVtx = Vertex(beamSpot.position(), beamSpot.covariance3D());
  }
  else cout << "No beam spot available from EventSetup" << endl;

  //Vertex collection
  Handle<VertexCollection> recVtxs;
  iEvent.getByLabel(vtxSample, recVtxs);
  unsigned int nVtxTrks = 0;
  //Multiple verteces -> primary = first
  if ( recVtxs->begin() != recVtxs->end() ) {
    thePrimaryVtx_multiplicity = recVtxs->size() ;

    if (addMuMulessPrimaryVertex_ || addXlessPrimaryVertex_ || resolveAmbiguity_) {
      //thePrimaryVtx = Vertex(*(recVtxs->begin()));
      //cout <<"here" <<endl;
      thePrimaryVtx = *(recVtxs->begin());
    }
    else {
      for ( reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx) {
        //take only the vertex with more tracks
        if (nVtxTrks < vtx->tracksSize() ) {
          nVtxTrks = vtx->tracksSize();
          thePrimaryVtx = Vertex(*vtx);
        }
      }
    }
  } else {
    //NO vtx multiplicity
    thePrimaryVtx = Vertex(beamSpot.position(), beamSpot.covariance3D());
    thePrimaryVtx_multiplicity = 1 ;
  }

  ////////////////////////////////////////////////////////////////////////
  //Track Builder
  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  ///Take primary vertex properties
  RefVtx = thePrimaryVtx.position(); /// reference primary vertex choosen
  priVtx_n = thePrimaryVtx_multiplicity ;
  priVtx_X = (thePrimaryVtx.position().x()) ;
  priVtx_Y = (thePrimaryVtx.position().y()) ;
  priVtx_Z = (thePrimaryVtx.position().z()) ;
  priVtx_XE = (thePrimaryVtx.xError()) ;
  priVtx_YE = (thePrimaryVtx.yError()) ;
  priVtx_ZE = (thePrimaryVtx.zError()) ;
  priVtx_NormChi2 = (thePrimaryVtx.normalizedChi2()) ;
  priVtx_Chi2 = thePrimaryVtx.chi2() ;
  priVtx_CL = ChiSquaredProbability( (double)(thePrimaryVtx.chi2()), (double)(thePrimaryVtx.ndof())) ;
  priVtx_tracks = thePrimaryVtx.tracksSize() ;
  VertexHigherPtSquared vertexHigherPtSquared ;
  priVtx_tracksPtSq = vertexHigherPtSquared.sumPtSquared(thePrimaryVtx) ;
  ///

  Handle< vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);
  /// ^ container of tracks not refitted
  Handle< vector<pat::GenericParticle> > theKaonRefittedPATTrackHandle;
  iEvent.getByLabel("cleanPatTrackKaonCands", theKaonRefittedPATTrackHandle);
  /// ^ container of tracks with kaon mass hypothesis

  //vectors of refitted and not refitted tracks
  for ( vector<pat::GenericParticle>::const_iterator TrackNotRefitted = thePATTrackHandle->begin(); TrackNotRefitted != thePATTrackHandle->end(); ++TrackNotRefitted ) {
    for ( vector<pat::GenericParticle>::const_iterator TrackRefitted = theKaonRefittedPATTrackHandle->begin(); TrackRefitted != theKaonRefittedPATTrackHandle->end(); ++TrackRefitted ) {
      if ( TrackNotRefitted->track().key() == TrackRefitted->track().key() ) {
        trNotRef->push_back( TrackNotRefitted->p() ) ;
        trRef->push_back( TrackRefitted->p() ) ;
        break ;
      }
    }
    break ;
  }

  ////////////////////////////////////////////////////////////////////////
  //Muons
  Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel("patMuonsWithTrigger", thePATMuonHandle);

  //de/dx values collection
  Handle<reco::DeDxDataValueMap> elossCollection;
  energyLoss = 0;
  iexception_dedx = 0;
  try {
    iEvent.getByLabel(DeDxEstimator_, elossCollection);
    energyLoss = elossCollection.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: elossCollection";
    iexception_dedx = 1;
  }

  /// dE/dx hits
  Handle<edm::ValueMap<reco::DeDxData> > dEdxTrackHandle;
  try {
    iEvent.getByLabel(m_dEdxDiscrimTag, dEdxTrackHandle);
    dEdxTrack = *dEdxTrackHandle.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: dEdxTrackHandle";
    iexception_dedx = 1;
  }
  /// dE/dx hits for kaons
  Handle<edm::ValueMap<reco::DeDxData> > dEdxTrackHandle_Kaon;
  try {
    iEvent.getByLabel(m_dEdxDiscrimTag_kaon, dEdxTrackHandle_Kaon);
    dEdxTrack_Kaon = *dEdxTrackHandle_Kaon.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: dEdxTrackHandle_Kaon";
    iexception_dedx = 1;
  }


  ////////////////// check MC truth //////////////////
  if (doMC) {
    /*
    // Get generated event
    //Handle<edm::HepMCProduct> hepEv;
    //iEvent.getByLabel("generator", hepEv);
    Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator", genEvtInfo);

    //const HepMC::GenEvent *myGenEvent = hepEv->GetEvent();
    const HepMC::GenEvent *myGenEvent = genEvtInfo->GetEvent();
    n_genEvtVtx = myGenEvent->vertices_size() ;

    HepMC::GenVertex* primaryGenVtx = *(myGenEvent->vertices_begin()) ;

    genEvtVtx_X->push_back( primaryGenVtx->point3d().x() );
    genEvtVtx_Y->push_back( primaryGenVtx->point3d().y() );
    genEvtVtx_Z->push_back( primaryGenVtx->point3d().z() );
    //genEvtVtx_XE = (primaryGenVtx->xError()) ;
    //genEvtVtx_YE = (primaryGenVtx->yError()) ;
    //genEvtVtx_ZE = (primaryGenVtx->zError()) ;
    //genEvtVtx_NormChi2 = (primaryGenVtx->normalizedChi2()) ;
    //genEvtVtx_Chi2 = primaryGenVtx->chi2() ;
    //genEvtVtx_CL = ChiSquaredProbability( (double)(primaryGenVtx.chi2()), (double)(primaryGenVtx.ndof())) ;
    genEvtVtx_particles->push_back( primaryGenVtx->particles_out_size() );
    */

    Handle< vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
    vector<PileupSummaryInfo>::const_iterator PVI;
    if (Debug_) cout <<"\nBunchXing multiplicity = " <<PupInfo->size() <<endl ;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
    if (Debug_) cout <<"Pileup Information: bunchXing, nvtx: " <<PVI->getBunchCrossing() <<" " <<PVI->getPU_NumInteractions() <<endl;

    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    if (Debug_) cout << "############### GenParticles Analysis ###############" << endl;
    float jpsiPx=0., jpsiPy=0., jpsiPz=0.;
    float  mupPx=0., mupPy=0., mupPz=0., mumPx=0., mumPy=0., mumPz=0.;
    float phiPx=0., phiPy=0., phiPz=0.;
    float  kpPx=0., kpPy=0., kpPz=0., kmPx=0., kmPy=0., kmPz=0.;
    //float pionPx=0., pionPy=0., pionPz=0., kaonPx=0., kaonPy=0., kaonPz=0.;
    //int pionCh=0, kaonCh=0 ;

    for (size_t i = 0; i < genParticles->size(); ++ i) {
      nMCAll++;
      const reco::GenParticle &p = (*genParticles)[i];
      int pdgid = p.pdgId() ;
      int dauNum = p.numberOfDaughters();
      MCPdgIdAll->push_back( pdgid );
      MCDanNumAll->push_back( dauNum );

      if ( MCExclusiveDecay ) {
        /// check if there is a MCMother which has MCDaughtersN daughters
        if ( abs(pdgid) == MCMother  &&  dauNum == MCDaughtersN ) {
          bool mumuOK = false;
          bool kkOK = false;
          //bool pionOK = false, kaonOK = false;

          for (int j=0; j<dauNum; ++j) {
            const Candidate *dau = p.daughter(j);
            if (Debug_) cout << "dauPdgId = " << dau->pdgId() << endl;

            /// check if one of B0 daughters is a psi(nS) whitch has 2 muons as daughters /// SEMRA ask again !!!
            int mumuId = 0 ;
            if (skipJPsi) /// SEMRA cleaned skipPsi2S
            if (Debug_) cout <<"Skipping J/psi!" <<endl ; /// SEMRA cleaned skipPsi2S
            //else if (skipPsi2S) /// SEMRA
            //  mumuId = 443 ; /// SEMRA (JPsi ID)

            if ( ((skipJPsi) && (dau->pdgId() == mumuId)) ||
            ((!skipJPsi) && (dau->pdgId()%1000 == 443)) ) {
              jpsiPx = dau->px(); jpsiPy = dau->py(); jpsiPz = dau->pz();
              int jpsiDauNum = dau->numberOfDaughters();
              if (Debug_) cout << "jpsiDauNum = " << jpsiDauNum << endl;
              int muNum = 0;
              for (int k=0; k<jpsiDauNum; ++k) {
                const Candidate *grandDau = dau->daughter(k);
                if (Debug_)  cout << "grandDauPdgId = " << grandDau->pdgId() << endl;
                if ( abs(grandDau->pdgId()) == 13 ) {
                  muNum++;
                  if (grandDau->pdgId() < 0) {
                    mupPx = grandDau->px(); mupPy = grandDau->py(); mupPz = grandDau->pz();
                  } else {
                    mumPx = grandDau->px(); mumPy = grandDau->py(); mumPz = grandDau->pz();
                  }
                }
              }
              if ( muNum == 2 ) mumuOK = true ;

            } /// end check if one of the MCMother daughters is a J/Psi or psi'

            /// for Phi
            phiPx = dau->px(); phiPy = dau->py(); phiPz = dau->pz();
            int phiDauNum = dau->numberOfDaughters();
            if (Debug_) cout << "phiDauNum = " << phiDauNum << endl;
            int kNum = 0;
            for (int n=0; n<phiDauNum; ++n) {
              const Candidate *grandDau = dau->daughter(n);
              if (Debug_)  cout << "grandDauPdgId = " << grandDau->pdgId() << endl;
              if ( abs(grandDau->pdgId()) == 321 ) {
                kNum++;
                if (grandDau->pdgId() < 0) {
                  kpPx = grandDau->px(); kpPy = grandDau->py(); kpPz = grandDau->pz();
                } else {
                  kmPx = grandDau->px(); kmPy = grandDau->py(); kmPz = grandDau->pz();
                }
              }
            }
            if ( kNum == 2 ) kkOK = true ;


            /*else if ( abs(dau->pdgId()) == 211 ) { // check if one of B0 daughters is a pion /// SEMRA ask again !!!
            pionPx = dau->px(); pionPy = dau->py(); pionPz = dau->pz();
            pionCh = (dau->pdgId() == 211)? 1 : -1;
            pionOK = true; /// SEMRA pions change with kaons for B0 ?
          } else if ( abs(dau->pdgId()) == 321 ) { // check if one of B0 daughters is a kaon /// SEMRA ask again !!!
          kaonPx = dau->px(); kaonPy=dau->py(); kaonPz=dau->pz();
          kaonCh = (dau->pdgId() == 321)? 1 : -1;
          kaonOK = true;
        }*/

      } /// end loop on MCMother daughters

      if (Debug_) cout << "mumuOK = " << mumuOK << ", kkOK = " << kkOK << endl;
      if ( mumuOK && kkOK ) {
        if (Debug_) {
          cout <<"\nnumber of X mothers = " <<p.numberOfMothers() <<endl ;
          cout <<"X mother pdgID = " <<p.mother(0)->pdgId() <<endl ;
        }
        ++nMCX ;
        PriVtxGen_X->push_back( p.vx() ) ;
        PriVtxGen_Y->push_back( p.vy() ) ;
        PriVtxGen_Z->push_back( p.vz() ) ;
        PriVtxGen_CL->push_back( p.vertexNormalizedChi2() ) ;
        PriVtxGen_Chi2->push_back( p.vertexChi2() ) ;
        PriVtxGen_Ndof->push_back( p.vertexNdof() ) ;

        Bool_t status = kTRUE ;
        const Candidate *x_ancestor = p.mother(0) ; /// a particle can have several mothers
        Int_t n_ancestors = 1 ;
        while ( status ) {
          if ( abs(x_ancestor->pdgId()) <= 8 || x_ancestor->pdgId() == 21 || x_ancestor->status() == 3 ) {
            status = kFALSE ;
            if (Debug_) cout <<"X ancestor ID = " <<x_ancestor->pdgId() <<endl ;
            genEvtVtx_X->push_back( x_ancestor->daughter(0)->vx() ) ;
            genEvtVtx_Y->push_back( x_ancestor->daughter(0)->vy() ) ;
            genEvtVtx_Z->push_back( x_ancestor->daughter(0)->vz() ) ;
            genEvtVtx_particles->push_back( x_ancestor->numberOfDaughters() ) ;
            n_XAncestors->push_back( n_ancestors ) ;
          }
          else {
            x_ancestor = x_ancestor->mother(0) ;
            n_ancestors++ ;
          }
        }

        MCJPsiPx->push_back(jpsiPx); MCJPsiPy->push_back(jpsiPy); MCJPsiPz->push_back(jpsiPz);
        MCmupPx->push_back(mupPx); MCmupPy->push_back(mupPy); MCmupPz->push_back(mupPz);
        MCmumPx->push_back(mumPx); MCmumPy->push_back(mumPy); MCmumPz->push_back(mumPz);
        MCPhiPx->push_back(phiPx); MCPhiPy->push_back(phiPy); MCPhiPz->push_back(phiPz);
        MCkpPx->push_back(kpPx); MCkpPy->push_back(kpPy); MCkpPz->push_back(kpPz);
        MCkmPx->push_back(kmPx); MCkmPy->push_back(kmPy); MCkmPz->push_back(kmPz);
        //MCpionPx->push_back(pionPx); MCpionPy->push_back(pionPy); MCpionPz->push_back(pionPz);
        //MCkaonPx->push_back(kaonPx); MCkaonPy->push_back(kaonPy); MCkaonPz->push_back(kaonPz);
        //MCpionCh->push_back(pionCh) ; MCkaonCh->push_back(kaonCh) ;
        decayChainOK = true;
        MCPx->push_back( p.px() );
        MCPy->push_back( p.py() );
        MCPz->push_back( p.pz() );
      }
      if (Debug_) cout << "decayChainOK = " << decayChainOK << endl;
    } // if ( abs(pdgid) == MCMother  &&  dauNum == 3 )
  } // if ( !MCExclusiveDecay )

} // for (size_t i = 0; i < genParticles->size(); ++ i)
} // if (doMC)


/// reconstruction only for events with B decaying in psi(nS)+Pi+K /// SEMRA JPsiPhi !!!
if ( (doMC && !MCExclusiveDecay) || (doMC && (MCExclusiveDecay && decayChainOK)) || doData ) {

  bool isEventWithInvalidMu = false;

  if (Debug_) cout << "starting event with " << thePATTrackHandle->size() << " tracks, and " << thePATMuonHandle->size() << " muons" << endl;

  if ((thePATMuonHandle->size()) * (thePATTrackHandle->size()) > 1E6) {
    cout << "Too many Muons: " << thePATMuonHandle->size() << ", and Tracks: " << thePATTrackHandle->size() << endl;
  } else //if (thePATMuonHandle->size() >= 2) { // check
    if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
      //There are at least 2 muons and 1 trigger firing
      if (Debug_) cout <<"============================  evt: " <<evtNum <<" Accept event with 2 mu and TRIGGER ==============================================" <<endl;

      ////////////////////////////////////////////////////////////////////////
      //Muon candidates
      for ( std::vector<pat::Muon>::const_iterator Muon1 = thePATMuonHandle->begin(); Muon1 != thePATMuonHandle->end(); ++Muon1 ) {

        /// push back all muon information
        ++nMu;
        const reco::Muon* rmu1 = dynamic_cast<const reco::Muon * >(Muon1->originalObject());

        ////////////////// check for muon1 //////////////////
        TrackRef muTrack1 = Muon1->track();
        if ( muTrack1.isNull() )
        continue;
        /// cuts on muon1
        if (rmu1->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
        || rmu1->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
        || rmu1->track()->chi2()/rmu1->track()->ndof() > MuMaxNormChi
        || fabs(rmu1->track()->dxy(RefVtx)) > MuMaxD0) {
          continue ;
        }

        //////////////////////////////////////////////////////
        //MUON 2 Loop
        for ( std::vector<pat::Muon>::const_iterator Muon2 = Muon1+1; Muon2 != thePATMuonHandle->end(); ++Muon2) {
          if(Muon2->charge() * Muon1->charge() > 0)
          continue ;
          const reco::Muon* rmu2 = dynamic_cast<const reco::Muon *>(Muon2->originalObject()) ;
          if (muon::overlap(*rmu1, *rmu2) )
          continue ;
          TrackRef muTrack2 = Muon2->track() ;
          if ( muTrack2.isNull() )
          continue ;
          /// cuts on muon2
          if (rmu2->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
          || rmu2->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
          || rmu2->track()->chi2()/rmu1->track()->ndof() > MuMaxNormChi
          || fabs(rmu2->track()->dxy(RefVtx)) > MuMaxD0) {
            continue ;
          }


          //take muon tracks
          TransientTrack muon1TT( muTrack1, &(*bFieldHandle) );
          TransientTrack muon2TT( muTrack2, &(*bFieldHandle) );
          //particle factory
          KinematicParticleFactoryFromTransientTrack pFactory;

          float chi = 0., ndf = 0.;
          vector<RefCountedKinematicParticle> muons; /// the final state muons produced by the KinematicParticleFactory
          //refitted muons particle
          muons.push_back( pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
          muons.push_back( pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
          //creating the vertex fitter for JPsi
          KinematicParticleVertexFitter MuMuFitter;
          RefCountedKinematicTree MuMuVertexFitTree;
          //fit the common vertex
          MuMuVertexFitTree = MuMuFitter.fit(muons);
          if (!MuMuVertexFitTree->isValid())
          continue ;
          //go to the top mother (J/Psi)
          MuMuVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle MuMuCand_fromFit = MuMuVertexFitTree->currentParticle();
          RefCountedKinematicVertex MuMuCand_vertex_fromFit = MuMuVertexFitTree->currentDecayVertex();

          //Muon one
          MuMuVertexFitTree->movePointerToTheFirstChild();
          RefCountedKinematicParticle Mu1Cand_fromFit = MuMuVertexFitTree->currentParticle();
          KinematicParameters Mu1Cand_KP = Mu1Cand_fromFit->currentState().kinematicParameters();
          //Muon two
          MuMuVertexFitTree->movePointerToTheNextChild();
          RefCountedKinematicParticle Mu2Cand_fromFit = MuMuVertexFitTree->currentParticle();
          KinematicParameters Mu2Cand_KP = Mu2Cand_fromFit->currentState().kinematicParameters();

          //////////////////////////////////////////////////////
          //Now I have a MuMu candidate

          if (MuMuCand_fromFit->currentState().mass() < MuMuMinMass  ||  MuMuCand_fromFit->currentState().mass() > MuMuMaxMass)
          continue ;


          //MUMU type
          Int_t dimuonType = 0;   //0 nothing,  1 J/psi  , 2 psi(2S)
          if ( MuMuCand_fromFit->currentState().mass() > JPsiMinMass  &&  MuMuCand_fromFit->currentState().mass() < JPsiMaxMass ) {
            dimuonType = 1 ;
          }
          if (Debug_) cout <<dimuonType <<endl;

          if (Debug_) cout <<"evt:" <<evtNum <<" MuMu with diMuonType = " <<dimuonType <<endl;

          muons.clear();


          //////////////////////////////////////////////////////////////////////
          /// for B0
          if (Debug_) cout <<"evt:"<<evtNum<< " is Invalid Muon ?  " <<isEventWithInvalidMu << endl;

          nTrk->push_back( thePATTrackHandle->size() ) ;

          if (thePATTrackHandle->size() < 2) {
            nX_pre0++;
          }
          if (Debug_) cout<<"nmumu : "<<nMuMu<<endl;

          ////////////////////////////
          //JPsi MASS WINDOW
          ++nMuMu;
          if (MuMuMass->at(nMuMu-1) < MuMuMinMass  ||  MuMuMass->at(nMuMu-1) > MuMuMaxMass)
            continue ;


          //cout<< "POINT 14" <<endl;


          ////////////////////////////////////////////////////////////////////////
          //Tracks cycle refitted tracks!
          for ( vector<pat::GenericParticle>::const_iterator Track1 = theKaonRefittedPATTrackHandle->begin(); Track1 != theKaonRefittedPATTrackHandle->end(); ++Track1 ) {

            //SKIP muon tracks
            if (Track1->track().key() == rmu1->track().key()  ||  Track1->track().key() == rmu2->track().key())
            continue ;
            //Track pt and track quality
            if (( Track1->track()->chi2()/Track1->track()->ndof() > TrMaxNormChi2 )  ||  Track1->pt() < TrMinPt)
            continue ;

            //TRACK 2 LOOP (always refitted)
            for ( vector<pat::GenericParticle>::const_iterator Track2 = Track1+1; Track2 != theKaonRefittedPATTrackHandle->end(); ++Track2 ){

              //Check that this second track doesn't overlap with the the first track candidate
              if (Track2->track().key() == Track1->track().key())
              continue ;

              //Track pt and track quality
              if ((Track2->track()->chi2() / Track2->track()->ndof() > TrMaxNormChi2)  ||  Track2->pt() < TrMinPt)
              continue;

              //SKIP muon tracks
              if (Track2->track().key() == rmu1->track().key()  ||  Track2->track().key() == rmu2->track().key())
              continue ;
              //Same charge tracks! (phi->kk)
              if (Track1->charge() * Track2->charge() > 0)
              continue ;


              //KAON TRACKs
              TransientTrack kaon1TT( Track1->track(), &(*bFieldHandle) );
              TransientTrack kaon2TT( Track2->track(), &(*bFieldHandle) );

              KinematicParticleFactoryFromTransientTrack pFactory;

              //Let's create the kaons from these tracks
              float chi = 0., ndf = 0.;
              vector<RefCountedKinematicParticle> kaons;
              kaons.push_back( pFactory.particle( kaon1TT, kaon_mass, chi, ndf, small_sigma));
              kaons.push_back( pFactory.particle( kaon2TT, kaon_mass, chi, ndf, small_sigma));

              //Kaon vertex fit
              KinematicParticleVertexFitter KKFitter;
              RefCountedKinematicTree KKVertexFitTree;
              KKVertexFitTree = KKFitter.fit(kaons);
              if (!KKVertexFitTree->isValid())
              continue ;

              //KK candidate
              KKVertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle KKCand_fromFit = KKVertexFitTree->currentParticle();
              RefCountedKinematicVertex KKCand_vertex_fromFit = KKVertexFitTree->currentDecayVertex();

              //Kaon 1
              KKVertexFitTree->movePointerToTheFirstChild();
              RefCountedKinematicParticle Ka1Cand_fromFit = KKVertexFitTree->currentParticle();
              KinematicParameters Ka1Cand_KP = Ka1Cand_fromFit->currentState().kinematicParameters();

              //Kaon 2
              KKVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle Ka2Cand_fromFit = KKVertexFitTree->currentParticle();
              KinematicParameters Ka2Cand_KP = Ka2Cand_fromFit->currentState().kinematicParameters();

              //KK mass window
              if (KKCand_fromFit->currentState().mass() < KKMinMass  ||  KKCand_fromFit->currentState().mass() > KKMaxMass)
              continue ;
              //KK paramters
              //mass
              KKMass->push_back( KKCand_fromFit->currentState().mass() );
              //KK vertex
              KKDecayVtx_X->push_back( KKCand_vertex_fromFit->position().x() );
              KKDecayVtx_Y->push_back( KKCand_vertex_fromFit->position().y() );
              KKDecayVtx_Z->push_back( KKCand_vertex_fromFit->position().z() );
              KKDecayVtx_XE->push_back( sqrt( KKCand_vertex_fromFit->error().cxx()) );
              KKDecayVtx_YE->push_back( sqrt( KKCand_vertex_fromFit->error().cyy()) );
              KKDecayVtx_ZE->push_back( sqrt( KKCand_vertex_fromFit->error().czz()) );
              KKVtx_CL->push_back( ChiSquaredProbability((double)( KKCand_vertex_fromFit->chiSquared()),(double)( KKCand_vertex_fromFit->degreesOfFreedom())) );
              KKVtx_Chi2->push_back( MuMuCand_vertex_fromFit->chiSquared() ) ;
              //KK momentum
              KKPx->push_back( Ka1Cand_KP.momentum().x() + Ka2Cand_KP.momentum().x() );
              KKPy->push_back( Ka1Cand_KP.momentum().y() + Ka2Cand_KP.momentum().y() );
              KKPz->push_back( Ka1Cand_KP.momentum().z() + Ka2Cand_KP.momentum().z() );
              //kaon ids
              ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track1));
              ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track2));

              //Single kaons from refitted phi
              //One
              ka1_KK_Px->push_back( Ka1Cand_KP.momentum().x());
              ka1_KK_Py->push_back( Ka1Cand_KP.momentum().y());
              ka1_KK_Pz->push_back( Ka1Cand_KP.momentum().z());
              ka1_KK_Chi2->push_back( Ka1Cand_fromFit->chiSquared());
              ka1_KK_NDF->push_back( Ka1Cand_fromFit->degreesOfFreedom());
              //two
              ka2_KK_Px->push_back( Ka2Cand_KP.momentum().x());
              ka2_KK_Py->push_back( Ka2Cand_KP.momentum().y());
              ka2_KK_Pz->push_back( Ka2Cand_KP.momentum().z());
              ka2_KK_Chi2->push_back( Ka2Cand_fromFit->chiSquared());
              ka2_KK_NDF->push_back( Ka2Cand_fromFit->degreesOfFreedom());

              ++nKK;
              kaons.clear();

              ////////////////// cuts on tracks' delta R for B0 //////////////////
              math::XYZTLorentzVector MuMu = (rmu1->p4() + rmu2->p4());
              math::XYZTLorentzVector MuMuKK = (MuMu + Track1->p4() + Track2->p4());
              float MuMu_K1_DR = sqrt( pow(MuMu.eta() - Track1->p4().eta(),2) + pow(MuMu.phi() - Track1->p4().phi(), 2) );
              float MuMu_K2_DR = sqrt( pow(MuMu.eta() - Track2->p4().eta(),2) + pow(MuMu.phi() - Track2->p4().phi(), 2) );
              float MuMuKK_K1_DR = sqrt( pow(MuMuKK.eta() - Track1->p4().eta(),2) + pow(MuMuKK.phi() - Track1->p4().phi(), 2) );
              float MuMuKK_K2_DR = sqrt( pow(MuMuKK.eta() - Track2->p4().eta(),2) + pow(MuMuKK.phi() - Track2->p4().phi(), 2) );

              DR_MuMu_K1->push_back(MuMu_K1_DR);
              DR_MuMu_K2->push_back(MuMu_K2_DR);
              DR_MuMuKK_K1->push_back(MuMuKK_K1_DR);
              DR_MuMuKK_K2->push_back(MuMuKK_K2_DR);


              if (UseXDR) {
                if (MuMuKK_K1_DR > XTrackMaxDR || MuMuKK_K2_DR > XTrackMaxDR)
                XTrackMaxDR = 2;
              } else {
                if (MuMu_K1_DR > MuMuTrackMaxDR || MuMu_K2_DR > MuMuTrackMaxDR)
                MuMuTrackMaxDR = 3.5;
              }



              //Control region mass
              if (((Track1->p4() + Track2->p4() + MuMu).M() > MuMuKKMaxB0Mass  ||  (Track1->p4() + Track2->p4() + MuMu).M() < MuMuKKMinB0Mass) && ((Track1->p4() + Track2->p4() + MuMu).M() >  MuMuKKMaxXMass))
              continue ;

              /// having two oppositely charged muons, and two oppositely charged tracks: try to vertex them
              //TransientTrack kaon1TT( Track1->track(), &(*bFieldHandle) );
              //TransientTrack kaon2TT( Track2->track(), &(*bFieldHandle) );

              //seek for the non-refitted track partner
              TransientTrack kaon2TT_notRefit ;
              Bool_t notRefittedPartner = false ;

              for ( vector<pat::GenericParticle>::const_iterator Track2_notRefit = thePATTrackHandle->begin(); Track2_notRefit != thePATTrackHandle->end(); ++Track2_notRefit )
              if ( Track2_notRefit->track().key() == Track2->track().key() ) {
                notRefittedPartner = true ;
                kaon2TT_notRefit = TransientTrack( Track2_notRefit->track(), &(*bFieldHandle) ) ;
                break ;
              }

              ////////////////////////////////////////////////
              //Search for the X cand

              vector<RefCountedKinematicParticle> xDaughters;
              xDaughters.push_back(pFactory.particle( muon1TT, muon_mass, chi, ndf, small_sigma));
              xDaughters.push_back(pFactory.particle( muon2TT, muon_mass, chi, ndf, small_sigma));
              xDaughters.push_back(pFactory.particle( kaon1TT, kaon_mass, chi, ndf, small_sigma));
              xDaughters.push_back(pFactory.particle( kaon2TT, kaon_mass, chi, ndf, small_sigma));

              RefCountedKinematicTree XVertexFitTree, XVertexFitTree_noKrefit ;
              KinematicConstrainedVertexFitter XFitter ;

              if (doMuMuMassConst) { // MassConst = 'MC' in the following
              MultiTrackKinematicConstraint *MuMu = 0;
              if (dimuonType == 1) { // constrain to JPsi mass
                MuMu = new TwoTrackMassKinematicConstraint(JPsi_mass);
              } //else if (dimuonType == 2) { // constrain to Psi(2S) mass /// SEMRA will we use this or not ?
                //MuMu = new TwoTrackMassKinematicConstraint(psi2S_mass);
                //} // already asked for: if (dimuonType == 0) continue ;

                XVertexFitTree = XFitter.fit( xDaughters, MuMu );
                if (notRefittedPartner) { // use not refitted kaons
                  xDaughters.pop_back() ;
                  xDaughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma));
                  XVertexFitTree_noKrefit = XFitter.fit( xDaughters, MuMu );
                }
              }
              else {
                XVertexFitTree = XFitter.fit( xDaughters );
                if (notRefittedPartner) { // use not refitted kaons
                  xDaughters.pop_back() ;
                  xDaughters.push_back(pFactory.particle( kaon2TT_notRefit, kaon_mass, chi, ndf, small_sigma));
                  XVertexFitTree_noKrefit = XFitter.fit( xDaughters );
                }
              }


              if ( !XVertexFitTree->isValid() )
              continue ;


              XVertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle XCand_fromMCFit = XVertexFitTree->currentParticle();
              RefCountedKinematicVertex XCand_vertex_fromMCFit = XVertexFitTree->currentDecayVertex();

              if ( !XCand_vertex_fromMCFit->vertexIsValid() )
              continue ;
              if ( XCand_vertex_fromMCFit->chiSquared() < 0  ||  XCand_vertex_fromMCFit->chiSquared() > 10000 )
              continue ;

              if (XCand_vertex_fromMCFit->chiSquared() / XCand_vertex_fromMCFit->degreesOfFreedom() > 7 )
              continue ;

              if ( XCand_fromMCFit->currentState().mass() > 100 )
              continue ;

              double xVtxProb = ChiSquaredProbability((double)(XCand_vertex_fromMCFit->chiSquared()), (double)(XCand_vertex_fromMCFit->degreesOfFreedom()));
              if ( xVtxProb < 0.005 ) //0.0001 )
              continue ;


              //////////////////// Lifetimes calculations ////////////////////
              TVector3 X_vtx((*XCand_vertex_fromMCFit).position().x(), (*XCand_vertex_fromMCFit).position().y(), 0) ;
              TVector3 X_pperp(XCand_fromMCFit->currentState().globalMomentum().x(), XCand_fromMCFit->currentState().globalMomentum().y(), 0);
              TVector3 X_vtx3D((*XCand_vertex_fromMCFit).position().x(), (*XCand_vertex_fromMCFit).position().y(), (*XCand_vertex_fromMCFit).position().z()) ;
              TVector3 X_pperp3D(XCand_fromMCFit->currentState().globalMomentum().x(),XCand_fromMCFit->currentState().globalMomentum().y(), XCand_fromMCFit->currentState().globalMomentum().z());

              AlgebraicVector3 X_v3pperp ;
              X_v3pperp[0] = X_pperp.x(); X_v3pperp[1] = X_pperp.y(); X_v3pperp[2] = 0.;
              TVector3 X_pvtx, X_pvtx3D, X_vdiff, X_vdiff3D ;
              double X_cosAlpha, X_cosAlpha3D, X_ctau ;
              VertexDistanceXY X_vdistXY ;
              Measurement1D X_distXY ;
              GlobalError X_v1e = (Vertex(*XCand_vertex_fromMCFit)).error();
              GlobalError X_v2e ;
              AlgebraicSymMatrix33 X_vXYe ;
              double X_ctauErr ;
              float X_lxy, X_lxyErr, X_lxyz, X_lxyzErr ;
              ROOT::Math::SVector<double, 3> X_vDiff, X_vDiff3D ; // needed by Similarity method


              ////////////////// Lifetime wrt PV for B0 //////////////////
              X_v2e = thePrimaryVtx.error();
              X_vXYe = X_v1e.matrix() + X_v2e.matrix() ;
              /// 2D
              X_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ;
              X_vdiff = X_vtx - X_pvtx ;
              X_cosAlpha = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp()) ;
              X_lxy = X_vdiff.Perp();
              X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
              X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
              X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(thePrimaryVtx));
              X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
              X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2()) ;
              /// 3D
              X_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z());
              X_vdiff3D = X_vtx3D - X_pvtx3D;
              X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/(X_vdiff3D.Mag()*X_pperp3D.Mag());
              X_lxyz = X_vdiff3D.Mag();
              X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
              X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();


              ////////////////// Last cuts for B0 //////////////////
              //if ( !(X_ctau/X_ctauErr > 2.8) || !(X_cosAlpha > 0.8) ) /// Alexis, BC_CTau_CTauErr plot has a cut which we don't want so we closed it.
              //continue ;


              ////////////////// fill X candidate variables //////////////////
              xMass->push_back( XCand_fromMCFit->currentState().mass()) ;
              xPx->push_back( XCand_fromMCFit->currentState().globalMomentum().x()) ;
              xPy->push_back( XCand_fromMCFit->currentState().globalMomentum().y()) ;
              xPz->push_back( XCand_fromMCFit->currentState().globalMomentum().z()) ;
              xPxE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(3,3) ) ) ;
              xPyE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(4,4) ) ) ;
              xPzE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(5,5) ) ) ;
              xVtx_CL->push_back( xVtxProb );
              xVtx_Chi2->push_back( XCand_vertex_fromMCFit->chiSquared() ) ;
              xDecayVtx_X->push_back((*XCand_vertex_fromMCFit).position().x());
              xDecayVtx_Y->push_back((*XCand_vertex_fromMCFit).position().y());
              xDecayVtx_Z->push_back((*XCand_vertex_fromMCFit).position().z());
              xDecayVtx_XE->push_back(sqrt((*XCand_vertex_fromMCFit).error().cxx()));
              xDecayVtx_YE->push_back(sqrt((*XCand_vertex_fromMCFit).error().cyy()));
              xDecayVtx_ZE->push_back(sqrt((*XCand_vertex_fromMCFit).error().czz()));
              XVertexFitTree->movePointerToTheFirstChild();
              RefCountedKinematicParticle mu1_MuMuKK = XVertexFitTree->currentParticle();
              XVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle mu2_MuMuKK = XVertexFitTree->currentParticle();
              XVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle k1_MuMuKK = XVertexFitTree->currentParticle();
              XVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle k2_MuMuKK = XVertexFitTree->currentParticle();

              /// muon1 & muon2
              mu1Px_MuMuKK->push_back( mu1_MuMuKK->currentState().globalMomentum().x() );
              mu1Py_MuMuKK->push_back( mu1_MuMuKK->currentState().globalMomentum().y() );
              mu1Pz_MuMuKK->push_back( mu1_MuMuKK->currentState().globalMomentum().z() );
              mu1E_MuMuKK->push_back( mu1_MuMuKK->currentState().kinematicParameters().energy() );
              mu2Px_MuMuKK->push_back( mu2_MuMuKK->currentState().globalMomentum().x() );
              mu2Py_MuMuKK->push_back( mu2_MuMuKK->currentState().globalMomentum().y() );
              mu2Pz_MuMuKK->push_back( mu2_MuMuKK->currentState().globalMomentum().z() );
              mu2E_MuMuKK->push_back( mu2_MuMuKK->currentState().kinematicParameters().energy() );
              /// kaon1 & kaon2
              k1Px_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().x() );
              k1Py_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().y() );
              k1Pz_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().z() );
              k1E_MuMuKK->push_back( k1_MuMuKK->currentState().kinematicParameters().energy() );
              Double_t theo = 0., sigma = 0. ;
              kaon1_nsigdedx->push_back( nsigmaofdedx(Track1->track(),theo,sigma) );
              kaon1_dedx->push_back( getEnergyLoss(Track1->track()) );
              kaon1_dedxMass->push_back( GetMass(Track1->track()) );
              kaon1_theo->push_back( theo );
              kaon1_sigma->push_back( sigma );
              kaon1_dedx_byHits->push_back( (dEdxTrack)[Track1->track()].dEdx() );
              kaon1_dedxErr_byHits->push_back( (dEdxTrack)[Track1->track()].dEdxError() );
              kaon1_saturMeas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfSaturatedMeasurements() );
              kaon1_Meas_byHits->push_back( (dEdxTrack)[Track1->track()].numberOfMeasurements() );
              k2Px_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().x() );
              k2Py_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().y() );
              k2Pz_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().z() );
              k2E_MuMuKK->push_back( k2_MuMuKK->currentState().kinematicParameters().energy() );
              theo = 0.; sigma = 0. ;
              kaon2_nsigdedx->push_back(nsigmaofdedx(Track2->track(),theo,sigma));
              kaon2_dedx->push_back(getEnergyLoss(Track2->track()));
              kaon2_dedxMass->push_back(GetMass(Track2->track()));
              kaon2_theo->push_back(theo);
              kaon2_sigma->push_back(sigma);
              kaon2_dedx_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdx() );
              kaon2_dedxErr_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].dEdxError() );
              kaon2_saturMeas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfSaturatedMeasurements() );
              kaon2_Meas_byHits->push_back( (dEdxTrack_Kaon)[Track2->track()].numberOfMeasurements() );
              /// PV
              xCosAlphaPV->push_back( X_cosAlpha ); xCosAlpha3DPV->push_back( X_cosAlpha3D );
              xCTauPV->push_back( X_ctau ); xCTauPVE->push_back( X_ctauErr );
              xLxyPV->push_back( X_lxy ); xLxyPVE->push_back( X_lxyErr );
              xLxyzPV->push_back( X_lxyz ); xLxyzPVE->push_back( X_lxyzErr );
              /// dxy, dz, dxyE, dzE for kaons from PV
              kaon1_dxy_PV->push_back( Track1->track()->dxy(RefVtx) );
              kaon1_dz_PV->push_back( Track1->track()->dz(RefVtx) );
              kaon2_dxy_PV->push_back( Track2->track()->dxy(RefVtx) );
              kaon2_dz_PV->push_back( Track2->track()->dz(RefVtx) );

              //MuMu

              //mass
              MuMuMass->push_back( MuMuCand_fromFit->currentState().mass() );
              //Vertex
              MuMuDecayVtx_X->push_back( MuMuCand_vertex_fromFit->position().x() );
              MuMuDecayVtx_Y->push_back( MuMuCand_vertex_fromFit->position().y() );
              MuMuDecayVtx_Z->push_back( MuMuCand_vertex_fromFit->position().z() );
              MuMuDecayVtx_XE->push_back( sqrt( MuMuCand_vertex_fromFit->error().cxx()) );
              MuMuDecayVtx_YE->push_back( sqrt( MuMuCand_vertex_fromFit->error().cyy()) );
              MuMuDecayVtx_ZE->push_back( sqrt( MuMuCand_vertex_fromFit->error().czz()) );
              MuMuVtx_CL->push_back( ChiSquaredProbability((double)( MuMuCand_vertex_fromFit->chiSquared()),(double)( MuMuCand_vertex_fromFit->degreesOfFreedom())) );
              MuMuVtx_Chi2->push_back( MuMuCand_vertex_fromFit->chiSquared() ) ;
              //P
              MuMuPx->push_back( Mu1Cand_KP.momentum().x() + Mu2Cand_KP.momentum().x() );
              MuMuPy->push_back( Mu1Cand_KP.momentum().y() + Mu2Cand_KP.momentum().y() );
              MuMuPz->push_back( Mu1Cand_KP.momentum().z() + Mu2Cand_KP.momentum().z() );

              mu1Idx->push_back(std::distance(thePATMuonHandle->begin(), Muon1));
              mu2Idx->push_back(std::distance(thePATMuonHandle->begin(), Muon2));

              //Refitted muons parameters
              //one
              mu1_MuMu_Px->push_back( Mu1Cand_KP.momentum().x());
              mu1_MuMu_Py->push_back( Mu1Cand_KP.momentum().y());
              mu1_MuMu_Pz->push_back( Mu1Cand_KP.momentum().z());
              mu1_MuMu_Chi2->push_back( Mu1Cand_fromFit->chiSquared());
              mu1_MuMu_NDF->push_back( Mu1Cand_fromFit->degreesOfFreedom());
              //two
              mu2_MuMu_Px->push_back( Mu2Cand_KP.momentum().x());
              mu2_MuMu_Py->push_back( Mu2Cand_KP.momentum().y());
              mu2_MuMu_Pz->push_back( Mu2Cand_KP.momentum().z());
              mu2_MuMu_Chi2->push_back( Mu2Cand_fromFit->chiSquared());
              mu2_MuMu_NDF->push_back( Mu2Cand_fromFit->degreesOfFreedom());

              //cout << "POINT 0" << endl;
              MuMuType->push_back(dimuonType);
              //cout << "POINT 2" << endl;

              int ntriggers = TriggersForMatching_.size();
              if (Debug_) cout << "ntriggers: " << ntriggers << endl;
              for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++)
              {
                if (Debug_) cout << "MatchingTriggerResult[" << MatchTrig << "]: " << MatchingTriggerResult[MatchTrig] << endl;
                if ( MatchingTriggerResult[MatchTrig]!=0 )
                {
                  if (Debug_) cout << "CHECKING FiltersForMatching_[" << MatchTrig << "]: " << FiltersForMatching_[MatchTrig] << endl;
                  pat::TriggerObjectStandAloneCollection mu1HLTMatches = Muon1->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
                  pat::TriggerObjectStandAloneCollection mu2HLTMatches = Muon2->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
                  bool pass1 = mu1HLTMatches.size() > 0;
                  bool pass2 = mu2HLTMatches.size() > 0;

                  if ((pass1) && (pass2))
                  {
                    //cout << "POINT 8" << endl;
                    MuMuMuonTrigMatch->push_back(true);
                    if (Debug_) cout <<"Matched MuMu" <<endl ;
                  } else
                  //cout << "POINT 9" << endl;
                  MuMuMuonTrigMatch->push_back(false);
                }
                else
                //cout << "POINT 10" << endl;
                MuMuMuonTrigMatch->push_back(false);
              }


              /// vertex without matched muons
              vector<TransientVertex> pvs ;
              Vertex MuMuLessPV = thePrimaryVtx ;

              if (addMuMulessPrimaryVertex_)
              {
                VertexReProducer revertex(recVtxs, iEvent);
                Handle<TrackCollection> pvtracks;
                iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
                Handle<BeamSpot>        pvbeamspot;
                iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

                if ( pvbeamspot.isValid() < 0 )
                continue ;
                if (pvbeamspot.id() != beamSpotHandle.id()) {
                  edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
                }
                const reco::Muon *rmu_1 = dynamic_cast<const reco::Muon*>( Muon1->originalObject() ) ;
                const reco::Muon *rmu_2 = dynamic_cast<const reco::Muon*>( Muon2->originalObject() ) ;

                if (rmu_1 != 0  &&  rmu_2 != 0  &&  rmu_1->track().id() == pvtracks.id()  &&  rmu_2->track().id() == pvtracks.id() ) {
                  TrackCollection MuMuLess;
                  MuMuLess.reserve(pvtracks->size());
                  for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
                    if (i == rmu_1->track().key()) continue;
                    if (i == rmu_2->track().key()) continue;
                    MuMuLess.push_back((*pvtracks)[i]);
                  }
                  if (Debug_) cout <<"pvbeamspot.isValid() = " <<pvbeamspot.isValid() <<endl ;
                  pvs = revertex.makeVertices(MuMuLess, *pvbeamspot, iSetup) ;
                  if (!pvs.empty()) {
                    MuMuLessPV = Vertex(pvs.front());
                  }
                }
              }
              PriVtxMuMuCorr_n->push_back( pvs.size() ) ;
              PriVtxMuMuCorr_X->push_back( MuMuLessPV.position().x() ) ;
              PriVtxMuMuCorr_Y->push_back( MuMuLessPV.position().y() ) ;
              PriVtxMuMuCorr_Z->push_back( MuMuLessPV.position().z() ) ;
              PriVtxMuMuCorr_EX->push_back( MuMuLessPV.xError() ) ;
              PriVtxMuMuCorr_EY->push_back( MuMuLessPV.yError() ) ;
              PriVtxMuMuCorr_EZ->push_back( MuMuLessPV.zError() ) ;
              PriVtxMuMuCorr_CL->push_back( ChiSquaredProbability( (double)(MuMuLessPV.chi2()), (double)(MuMuLessPV.ndof())) ) ;
              PriVtxMuMuCorr_Chi2->push_back( MuMuLessPV.chi2() ) ;
              PriVtxMuMuCorr_tracks->push_back( MuMuLessPV.tracksSize() ) ;




              ////////////////// Lifetime wrt BS for B0 //////////////////
              X_v2e = theBeamSpotVtx.error();
              X_vXYe = X_v1e.matrix() + X_v2e.matrix();
              /// 2D
              X_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
              X_vdiff = X_vtx - X_pvtx;
              X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp());
              X_lxy = X_vdiff.Perp();
              X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
              X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
              X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(theBeamSpotVtx));
              X_ctau = X_distXY.value() * X_cosAlpha * (XCand_fromMCFit->currentState().mass() / X_pperp.Perp()) ;
              X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass()/X_pperp.Perp2();
              /// 3D
              X_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
              X_vdiff3D = X_vtx3D - X_pvtx3D;
              X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/(X_vdiff3D.Mag()*X_pperp3D.Mag());
              X_lxyz = X_vdiff3D.Mag();
              X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
              X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();


              ////////////////// BS (beam spot) for B0 //////////////////
              xCosAlphaBS->push_back( X_cosAlpha ); xCosAlpha3DBS->push_back( X_cosAlpha3D );
              xCTauBS->push_back( X_ctau ); xCTauBSE->push_back( X_ctauErr );
              xLxyBS->push_back( X_lxy ); xLxyBSE->push_back( X_lxyErr );
              xLxyzBS->push_back( X_lxyz ); xLxyzBSE->push_back( X_lxyzErr );

              vector<TransientVertex> X_pvs ;
              Vertex XLessPV = thePrimaryVtx ;


              if (addXlessPrimaryVertex_)
              {
                VertexReProducer revertex(recVtxs, iEvent);
                Handle<TrackCollection> pvtracks;
                iEvent.getByLabel(revertex.inputTracks(), pvtracks);
                Handle<BeamSpot>        pvbeamspot;
                iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

                if (pvbeamspot.id() != beamSpotHandle.id() )
                edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

                const reco::Muon *Xrmu_1 = dynamic_cast<const reco::Muon *>(Muon1->originalObject());
                const reco::Muon *Xrmu_2 = dynamic_cast<const reco::Muon *>(Muon2->originalObject());

                if (Xrmu_1 != 0  &&  Xrmu_2 != 0  &&  Xrmu_1->track().id() == pvtracks.id()  &&  Xrmu_2->track().id() == pvtracks.id()
                &&  Track1->track().id() == pvtracks.id()  &&  Track2->track().id() ==  pvtracks.id()) {
                  vector<TransientTrack> XLess; // need TransientTrack to keep the TrackRef
                  XLess.reserve( pvtracks->size() );
                  Double_t removedTrksPtSq = 0. ;
                  for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
                    if (i == Xrmu_1->track().key()) { removedTrksPtSq += (Xrmu_1->track()->pt())*(Xrmu_1->track()->pt()) ;
                      continue; }
                      if (i == Xrmu_2->track().key()) { removedTrksPtSq += (Xrmu_2->track()->pt())*(Xrmu_2->track()->pt()) ;
                        continue; }
                        if (i == Track1->track().key()) { removedTrksPtSq += (Track1->track()->pt())*(Track1->track()->pt()) ;
                          continue; }
                          if (i == Track2->track().key()) { removedTrksPtSq += (Track2->track()->pt())*(Track2->track()->pt()) ;
                            continue; }

                            reco::TrackRef trk_now(pvtracks, i) ;
                            TransientTrack transientTrack = theTTBuilder->build( trk_now );
                            transientTrack.setBeamSpot( beamSpot );
                            XLess.push_back( transientTrack );
                          }
                          if ( removedTrksPtSq > 0. ) {
                            X_pvs = revertex.makeVertices(XLess, *pvbeamspot, iSetup) ; // list of PV
                          } else
                          if (Debug_) cout <<"\n\\\\\\\\\\\\\\\\\\\\ excluded tracks pT^2 = 0 \\\\\\\\\\\\\\\\\\\\\n" <<endl ;
                          if ( !X_pvs.empty() ) {
                            XLessPV = Vertex(X_pvs.front());
                            XLessPV_tracksPtSq->push_back( vertexHigherPtSquared.sumPtSquared(XLessPV) ) ;
                            XLessPV_4tracksPtSq->push_back( removedTrksPtSq ) ;
                            if (Debug_) {
                              cout <<"\nXLessPV_z = " <<XLessPV.position().z() <<endl ;
                              cout <<"XLessPV_tracks = " <<XLessPV.tracksSize() <<endl ;
                              cout <<"XLessPV_tracksPtSq = " <<vertexHigherPtSquared.sumPtSquared(XLessPV) <<endl ;
                              cout <<"XLessPV_removedTracksPtSq = " <<removedTrksPtSq <<endl ;
                              cout <<"X_pvs->size() = " <<X_pvs.size() <<endl ;
                              cout <<"priVtx_z = " << priVtx_Z <<endl ;
                              cout <<"priVtx_tracks = " <<priVtx_tracks <<endl ;
                              cout <<"priVtx_tracksPtSq = " <<priVtx_tracksPtSq <<endl ;
                              cout <<"recVtxs->size() = " <<recVtxs->size() <<endl ;
                            }
                          }
                        }
                      }


                      PriVtxXLess_n->push_back( X_pvs.size() ) ;
                      PriVtxXLess_X->push_back( XLessPV.position().x() ) ;
                      PriVtxXLess_Y->push_back( XLessPV.position().y() ) ;
                      PriVtxXLess_Z->push_back( XLessPV.position().z() ) ;
                      PriVtxXLess_EX->push_back( XLessPV.xError() ) ;
                      PriVtxXLess_EY->push_back( XLessPV.yError() ) ;
                      PriVtxXLess_EZ->push_back( XLessPV.zError() ) ;
                      PriVtxXLess_CL->push_back( ChiSquaredProbability( (double)(XLessPV.chi2()), (double)(XLessPV.ndof())) );
                      PriVtxXLess_Chi2->push_back( XLessPV.chi2() ) ;
                      PriVtxXLess_tracks->push_back( XLessPV.tracksSize() ) ;

                      /// dxy, dz, dxyE, dzE for kaons from BS
                      math::XYZPoint BSVtx;
                      BSVtx = theBeamSpotVtx.position();
                      kaon1_dxy_BS->push_back( Track1->track()->dxy(BSVtx) );
                      kaon1_dz_BS->push_back( Track1->track()->dz(BSVtx) );
                      kaon2_dxy_BS->push_back( Track2->track()->dxy(BSVtx) );
                      kaon2_dz_BS->push_back( Track2->track()->dz(BSVtx) );


                      ////////////////// Lifetime wrt B0LessPV for B0 //////////////////
                      X_v2e = XLessPV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(XLessPV.position().x(), XLessPV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp());
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(XLessPV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      /// 3D
                      X_pvtx3D.SetXYZ(XLessPV.position().x(), XLessPV.position().y(), XLessPV.position().z());
                      X_vdiff3D = X_vtx3D - X_pvtx3D;
                      X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/( X_vdiff3D.Mag()*X_pperp3D.Mag() );
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaXLessPV->push_back( X_cosAlpha ) ; xCosAlpha3DXLessPV->push_back( X_cosAlpha3D ) ;
                      xCTauXLessPV->push_back( X_ctau ) ; xCTauXLessPVE->push_back( X_ctauErr ) ;
                      xLxyXLessPV->push_back( X_lxy ) ; xLxyXLessPVE->push_back( X_lxyErr ) ;
                      xLxyzXLessPV->push_back( X_lxyz ) ; xLxyzXLessPVE->push_back( X_lxyzErr ) ;

                      /// dxy, dz, dxyE, dzE for kaons from B0LessPV
                      math::XYZPoint XLessPVvtx;
                      XLessPVvtx = XLessPV.position();
                      kaon1_dxy_XLessPV->push_back( Track1->track()->dxy(XLessPVvtx) );
                      kaon1_dz_XLessPV->push_back( Track1->track()->dz(XLessPVvtx) );
                      kaon2_dxy_XLessPV->push_back( Track2->track()->dxy(XLessPVvtx) );
                      kaon2_dz_XLessPV->push_back( Track2->track()->dz(XLessPVvtx) );

                      kaon1_dxyE->push_back( Track1->track()->dxyError() );
                      kaon1_dzE->push_back( Track1->track()->dzError() );
                      kaon2_dxyE->push_back( Track2->track()->dxyError() );
                      kaon2_dzE->push_back( Track2->track()->dzError() );


                      /// Find the PV among the original offlinePV with the largest B0_cos(alpha)
                      Vertex theCosAlphaV = thePrimaryVtx ;
                      float maxCosAlpha = -1. ;

                      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
                        X_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;
                        X_vdiff = X_vtx - X_pvtx ;
                        float cosAlpha_temp = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp()) ; // Perp() == Mag() when z = 0

                        if ( cosAlpha_temp > maxCosAlpha ) {
                          maxCosAlpha = cosAlpha_temp ;
                          theCosAlphaV = Vertex(*itv) ;
                        }
                      }

                      PriVtx_XCosAlpha_n->push_back( recVtxs->size() ) ;
                      PriVtx_XCosAlpha_X->push_back( theCosAlphaV.position().x() ) ;
                      PriVtx_XCosAlpha_Y->push_back( theCosAlphaV.position().y() ) ;
                      PriVtx_XCosAlpha_Z->push_back( theCosAlphaV.position().z() ) ;
                      PriVtx_XCosAlpha_EX->push_back( theCosAlphaV.xError() ) ;
                      PriVtx_XCosAlpha_EY->push_back( theCosAlphaV.yError() ) ;
                      PriVtx_XCosAlpha_EZ->push_back( theCosAlphaV.zError() ) ;
                      PriVtx_XCosAlpha_CL->push_back( ChiSquaredProbability((double)(theCosAlphaV.chi2()), (double)(theCosAlphaV.ndof())) ) ;
                      PriVtx_XCosAlpha_Chi2->push_back( theCosAlphaV.chi2() ) ;
                      PriVtx_XCosAlpha_tracks->push_back( theCosAlphaV.tracksSize() ) ;


                      /// Find the PV among the original offlinePV with the largest B0_cos(alpha) 3D
                      Vertex theCosAlpha3DV = thePrimaryVtx ;
                      float maxCosAlpha3D = -1. ;

                      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
                        X_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;
                        X_vdiff3D = X_vtx3D - X_pvtx3D ;
                        float cosAlpha_temp3D = X_vdiff3D.Dot(X_pperp3D) / (X_vdiff3D.Mag()*X_pperp3D.Mag()) ;

                        if ( cosAlpha_temp3D > maxCosAlpha3D ) {
                          maxCosAlpha3D = cosAlpha_temp3D ;
                          theCosAlpha3DV = Vertex(*itv) ;
                        }
                      }


                      ////////////////// Lifetime wrt PV with largest B0_cos(alpha) 3D candidate //////////////////
                      X_v2e = theCosAlpha3DV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp()); ;
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), theCosAlpha3DV.position().z()) ;
                      X_vdiff3D = X_vtx3D - X_pvtx3D ;
                      X_cosAlpha3D =  maxCosAlpha3D ;
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaPVCosAlpha3D->push_back( X_cosAlpha ) ; xCosAlpha3DPVCosAlpha3D->push_back( X_cosAlpha3D ) ;
                      xCTauPVCosAlpha3D->push_back( X_ctau ) ; xCTauPVCosAlpha3DE->push_back( X_ctauErr ) ;
                      xLxyPVCosAlpha3D->push_back( X_lxy ) ; xLxyPVCosAlpha3DE->push_back( X_lxyErr ) ;
                      xLxyzPVCosAlpha3D->push_back( X_lxyz ) ; xLxyzPVCosAlpha3DE->push_back( X_lxyzErr ) ;


                      /// Find the PV among the B0lessPV with the largest B0_cos(alpha)
                      Vertex theXLessCosAlphaV = thePrimaryVtx ;
                      maxCosAlpha = -1. ;

                      for (vector<TransientVertex>::iterator itv = X_pvs.begin(), itvend = X_pvs.end(); itv != itvend; ++itv) {
                        X_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;
                        X_vdiff = X_vtx - X_pvtx ;
                        float cosAlpha_temp = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp()) ; // Perp() == Mag() when z = 0

                        if ( cosAlpha_temp > maxCosAlpha ) {
                          maxCosAlpha = cosAlpha_temp ;
                          theXLessCosAlphaV = Vertex(*itv) ;
                        }
                      }

                      PriVtxXLess_XCosAlpha_n->push_back( X_pvs.size() ) ;
                      PriVtxXLess_XCosAlpha_X->push_back( theXLessCosAlphaV.position().x() ) ;
                      PriVtxXLess_XCosAlpha_Y->push_back( theXLessCosAlphaV.position().y() ) ;
                      PriVtxXLess_XCosAlpha_Z->push_back( theXLessCosAlphaV.position().z() ) ;
                      PriVtxXLess_XCosAlpha_EX->push_back( theXLessCosAlphaV.xError() ) ;
                      PriVtxXLess_XCosAlpha_EY->push_back( theXLessCosAlphaV.yError() ) ;
                      PriVtxXLess_XCosAlpha_EZ->push_back( theXLessCosAlphaV.zError() ) ;
                      PriVtxXLess_XCosAlpha_CL->push_back( ChiSquaredProbability((double)(theXLessCosAlphaV.chi2()), (double)(theXLessCosAlphaV.ndof())) ) ;
                      PriVtxXLess_XCosAlpha_Chi2->push_back( theXLessCosAlphaV.chi2() ) ;
                      PriVtxXLess_XCosAlpha_tracks->push_back( theXLessCosAlphaV.tracksSize() ) ;


                      ////////////////// Lifetime wrt PV with largest B0_cos(alpha) candidate //////////////////
                      X_v2e = theCosAlphaV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha =  maxCosAlpha ;
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), theCosAlphaV.position().z());
                      X_vdiff3D = X_vtx3D - X_pvtx3D;
                      X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/( X_vdiff3D.Mag()*X_pperp3D.Mag() );
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaPVCosAlpha->push_back( X_cosAlpha ) ; xCosAlpha3DPVCosAlpha->push_back( X_cosAlpha3D ) ;
                      xCTauPVCosAlpha->push_back( X_ctau ) ; xCTauPVCosAlphaE->push_back( X_ctauErr ) ;
                      xLxyPVCosAlpha->push_back( X_lxy ) ; xLxyPVCosAlphaE->push_back( X_lxyErr ) ;
                      xLxyzPVCosAlpha->push_back( X_lxyz ) ; xLxyzPVCosAlphaE->push_back( X_lxyzErr ) ;

                      PriVtx_XCosAlpha3D_n->push_back( recVtxs->size() ) ;
                      PriVtx_XCosAlpha3D_X->push_back( theCosAlpha3DV.position().x() ) ;
                      PriVtx_XCosAlpha3D_Y->push_back( theCosAlpha3DV.position().y() ) ;
                      PriVtx_XCosAlpha3D_Z->push_back( theCosAlpha3DV.position().z() ) ;
                      PriVtx_XCosAlpha3D_EX->push_back( theCosAlpha3DV.xError() ) ;
                      PriVtx_XCosAlpha3D_EY->push_back( theCosAlpha3DV.yError() ) ;
                      PriVtx_XCosAlpha3D_EZ->push_back( theCosAlpha3DV.zError() ) ;
                      PriVtx_XCosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theCosAlpha3DV.chi2()), (double)(theCosAlpha3DV.ndof())) ) ;
                      PriVtx_XCosAlpha3D_Chi2->push_back( theCosAlpha3DV.chi2() ) ;
                      PriVtx_XCosAlpha3D_tracks->push_back( theCosAlpha3DV.tracksSize() ) ;


                      ////////////////// Lifetime wrt B0LessPV with largest B0_cos(alpha) candidate
                      X_v2e = theXLessCosAlphaV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theXLessCosAlphaV.position().x(), theXLessCosAlphaV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha =  maxCosAlpha ;
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp() ;
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theXLessCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theXLessCosAlphaV.position().x(), theXLessCosAlphaV.position().y(), theXLessCosAlphaV.position().z());
                      X_vdiff3D = X_vtx3D - X_pvtx3D;
                      X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/( X_vdiff3D.Mag()*X_pperp3D.Mag() );
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaXLessPVCosAlpha->push_back( X_cosAlpha ) ; xCosAlpha3DXLessPVCosAlpha->push_back( X_cosAlpha3D ) ;
                      xCTauXLessPVCosAlpha->push_back( X_ctau ) ; xCTauXLessPVCosAlphaE->push_back( X_ctauErr ) ;
                      xLxyXLessPVCosAlpha->push_back( X_lxy ) ; xLxyXLessPVCosAlphaE->push_back( X_lxyErr ) ;
                      xLxyzXLessPVCosAlpha->push_back( X_lxyz ) ; xLxyzXLessPVCosAlphaE->push_back( X_lxyzErr ) ;


                      /// Find the PV among the B0lessPV with the largest B0_cos(alpha) 3D
                      Vertex theXLessCosAlpha3DV = thePrimaryVtx ;
                      maxCosAlpha3D = -1. ;

                      for (vector<TransientVertex>::iterator itv = X_pvs.begin(), itvend = X_pvs.end(); itv != itvend; ++itv) {
                        X_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;
                        X_vdiff3D = X_vtx3D - X_pvtx3D ;
                        float cosAlpha_temp3D = X_vdiff3D.Dot(X_pperp3D) / (X_vdiff3D.Mag()*X_pperp3D.Mag()) ;

                        if ( cosAlpha_temp3D > maxCosAlpha3D ) {
                          maxCosAlpha3D = cosAlpha_temp3D ;
                          theXLessCosAlpha3DV = Vertex(*itv) ;
                        }
                      }

                      PriVtxXLess_XCosAlpha3D_n->push_back( X_pvs.size() ) ;
                      PriVtxXLess_XCosAlpha3D_X->push_back( theXLessCosAlpha3DV.position().x() ) ;
                      PriVtxXLess_XCosAlpha3D_Y->push_back( theXLessCosAlpha3DV.position().y() ) ;
                      PriVtxXLess_XCosAlpha3D_Z->push_back( theXLessCosAlpha3DV.position().z() ) ;
                      PriVtxXLess_XCosAlpha3D_EX->push_back( theXLessCosAlpha3DV.xError() ) ;
                      PriVtxXLess_XCosAlpha3D_EY->push_back( theXLessCosAlpha3DV.yError() ) ;
                      PriVtxXLess_XCosAlpha3D_EZ->push_back( theXLessCosAlpha3DV.zError() ) ;
                      PriVtxXLess_XCosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theXLessCosAlpha3DV.chi2()), (double)(theXLessCosAlpha3DV.ndof())) ) ;
                      PriVtxXLess_XCosAlpha3D_Chi2->push_back( theXLessCosAlpha3DV.chi2() ) ;
                      PriVtxXLess_XCosAlpha3D_tracks->push_back( theXLessCosAlpha3DV.tracksSize() ) ;


                      ////////////////// Lifetime wrt B0LessPV with largest B0_cos(alpha) 3D candidate
                      X_v2e = theXLessCosAlpha3DV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theXLessCosAlpha3DV.position().x(), theXLessCosAlpha3DV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp());
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theXLessCosAlpha3DV.position().x(), theXLessCosAlpha3DV.position().y(), theXLessCosAlpha3DV.position().z()) ;
                      X_vdiff3D = X_vtx3D - X_pvtx3D ;
                      X_cosAlpha3D =  maxCosAlpha3D ;
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();
                      X_lxy = X_vdiff3D.Dot(X_pperp) / X_pperp.Mag() ;

                      xCosAlphaXLessPVCosAlpha3D->push_back( X_cosAlpha ) ; xCosAlpha3DXLessPVCosAlpha3D->push_back( X_cosAlpha3D ) ;
                      xCTauXLessPVCosAlpha3D->push_back( X_ctau ) ; xCTauXLessPVCosAlpha3DE->push_back( X_ctauErr ) ;
                      xLxyXLessPVCosAlpha3D->push_back( X_lxy ) ; xLxyXLessPVCosAlpha3DE->push_back( X_lxyErr ) ;
                      xLxyzXLessPVCosAlpha3D->push_back( X_lxyz ) ; xLxyzXLessPVCosAlpha3DE->push_back( X_lxyzErr ) ;

                      Vertex theOtherV = thePrimaryVtx;
                      if (resolveAmbiguity_) {
                        float minDz = 999999. ;
                        if (!addXlessPrimaryVertex_) {
                          for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv)
                          {
                            float deltaZ = fabs((*XCand_vertex_fromMCFit).position().z() - itv->position().z()) ;
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              thePrimaryVtx = Vertex(*itv);
                              theOtherV = thePrimaryVtx;
                            }
                          }
                        } else {
                          for (vector<TransientVertex>::iterator itv2 = X_pvs.begin(), itvend2 = X_pvs.end(); itv2 != itvend2; ++itv2)
                          {
                            float deltaZ = fabs((*XCand_vertex_fromMCFit).position().z() - itv2->position().z()) ;
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              Vertex XLessPV = Vertex(*itv2);
                              thePrimaryVtx = XLessPV;
                              theOtherV = XLessPV;
                            }
                          }
                        }
                      }

                      Vertex TheOtherVertex3D = thePrimaryVtx;
                      if (Debug_) cout<<" choose PV ="<< endl;
                      Int_t theXCorrPV_multiplicity = -1 ;
                      if (resolveAmbiguity_) {
                        float minDz = 999999.;
                        if (!addXlessPrimaryVertex_) {
                          theXCorrPV_multiplicity = recVtxs->size() ;
                          for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
                            float deltaZ = fabs((*XCand_vertex_fromMCFit).position().z() - itv->position().z()) ;
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              TheOtherVertex3D = Vertex(*itv);
                            }
                          }
                        } else {
                          theXCorrPV_multiplicity = X_pvs.size() ;
                          for (vector<TransientVertex>::iterator itv2 = X_pvs.begin(), itvend2 = X_pvs.end(); itv2 != itvend2; ++itv2) {
                            VertexDistance3D a3d;
                            float deltaZ   = a3d.distance(Vertex(*itv2), Vertex(*XCand_vertex_fromMCFit)).value();
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              Vertex XLessPV = Vertex(*itv2);
                              TheOtherVertex3D = XLessPV;
                              //cout<<" z(X) - z(vtx) min="<<minDz<<endl;
                            }

                          }
                        }
                      }

                      PriVtxXCorr_n->push_back( theXCorrPV_multiplicity ) ;
                      PriVtxXCorr_X->push_back( thePrimaryVtx.position().x() ) ;
                      PriVtxXCorr_Y->push_back( thePrimaryVtx.position().y() ) ;
                      PriVtxXCorr_Z->push_back( thePrimaryVtx.position().z() ) ;
                      PriVtxXCorr_EX->push_back( thePrimaryVtx.xError() ) ;
                      PriVtxXCorr_EY->push_back( thePrimaryVtx.yError() ) ;
                      PriVtxXCorr_EZ->push_back( thePrimaryVtx.zError() ) ;
                      PriVtxXCorr_CL->push_back( ChiSquaredProbability( (double)(thePrimaryVtx.chi2()), (double)(thePrimaryVtx.ndof())) );
                      PriVtxXCorr_Chi2->push_back( thePrimaryVtx.chi2() ) ;
                      PriVtxXCorr_tracks->push_back( thePrimaryVtx.tracksSize() ) ;


                      ////////////////// Lifetime wrt PV with smaller longitudinal X impact parameter for B0  //////////////////
                      X_pvtx.SetXYZ(theOtherV.position().x(), theOtherV.position().y(), 0);
                      X_vdiff = X_vtx - X_pvtx;
                      X_cosAlpha = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp());
                      X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(theOtherV));
                      double X_ctauPVX = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      GlobalError X_v1eX = (Vertex(*XCand_vertex_fromMCFit)).error();
                      GlobalError X_v2eX = theOtherV.error();
                      AlgebraicSymMatrix33 X_vXYeX = X_v1eX.matrix() + X_v2eX.matrix();
                      double ctauErrPVX = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYeX)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      float lxyPVX = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      float lxyzPVX = X_vdiff3D.Dot(X_pperp3D) / X_pperp3D.Mag() ;
                      xCosAlphaPVX->push_back(X_cosAlpha);
                      xCTauPVX->push_back(X_ctauPVX); xCTauPVXE->push_back(ctauErrPVX);
                      xLxyPVX->push_back(lxyPVX);
                      xLxyzPVX->push_back(lxyzPVX);
                      VertexDistance3D a3d;
                      float Dist3DPV     = a3d.distance(TheOtherVertex3D, Vertex(*XCand_vertex_fromMCFit)).value() ;
                      float Dist3DPV_err = a3d.distance(TheOtherVertex3D, Vertex(*XCand_vertex_fromMCFit)).error() ;
                      xCTauPVX_3D->push_back(Dist3DPV);
                      xCTauPVX_3D_err->push_back(Dist3DPV_err);
                      //cout << Dist3DPV << " " << Dist3DPV_err << endl;
                      X_MuMuIdx->push_back(nMuMu-1);
                      X_ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track1));
                      X_ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), Track2));
                      nX++;
                      xDaughters.clear();


                      ////////////////// flag for checking the Kaons from PV or not PV //////////////////
                      /// flag for kaon1
                      vector<TransientTrack> vertexTracksKaon1;
                      //cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << endl;
                      //cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << endl;
                      for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {

                        TrackRef trackRefkaon1 = iTrack->castTo<TrackRef>();
                        //cout << "\ntrackRefkaon1 = " << trackRefkaon1 << endl;
                        //cout <<"before match" ;
                        if ( (Track1->track().key() == trackRefkaon1.key()) ) {
                          cout << "\ninside match" << endl;
                          TransientTrack kaon1TT(trackRefkaon1, &(*bFieldHandle) );
                          vertexTracksKaon1.push_back(kaon1TT);
                        }
                      }
                      //cout << "\nvertexTracksKaon1.size() = " << vertexTracksKaon1.size() << endl;
                      if (vertexTracksKaon1.size()==0)
                      Kaon1FromPV->push_back(false);
                      else
                      Kaon1FromPV->push_back(true);

                      /// flag for kaon2
                      vector<TransientTrack> vertexTracksKaon2;
                      //cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << endl;
                      //cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << endl;
                      for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {

                        TrackRef trackRefkaon2 = iTrack->castTo<TrackRef>();
                        //cout << "\ntrackRefkaon2 = " << trackRefkaon2 << endl;
                        //cout <<"before match" ;
                        if (  (Track2->track().key() == trackRefkaon2.key()) ) {
                          TransientTrack kaon2TT(trackRefkaon2, &(*bFieldHandle) );
                          vertexTracksKaon2.push_back(kaon2TT);
                        }
                      }
                      //cout << "\nvertexTracksKaon1.size() = " << vertexTracksKaon1.size() << endl;
                      if (vertexTracksKaon2.size()==0)
                      Kaon2FromPV->push_back(false);
                      else
                      Kaon2FromPV->push_back(true);

                    } // 2nd loop over track (look for k2)
                  } // 1st loop over track (look for k1)
                } // 2nd loop over muons (look for mu-)
              } //first loop over muons (look for mu+)
            } // if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
            } // if (doMC || doData)
            // AT THE END OF THE EVENT fill the tree and clear the vectors
            // ===========================================================

            if (nX > 0)
            xTree->Fill() ;

            /// trigger stuff
            trigRes->clear(); trigNames->clear(); L1TT->clear(); MatchTriggerNames->clear();
            /// event numbers
            runNum = 0; evtNum = 0; lumiNum = 0;
            /// counters for x(4140)
            nMu = 0; nMuMu = 0; nX = 0; nKK = 0;
            nX_pre0 = 0; nX_pre1 = 0; nX_pre2 = 0; nX_pre3 = 0; nX_pre4 = 0; nX_pre5 = 0; nX_pre6 = 0; nX_pre7 = 0; nX_pre8 = 0; nX_pre9 = 0; nX_pre10 = 0; nX_pre11 = 0; nX_pre12 = 0; nX_pre13 = 0; nX_pre14 = 0; nX_pre15 = 0;
            //nX = 0;
            /// indices
            mu1Idx->clear(); mu2Idx->clear();
            ka1Idx->clear(); ka2Idx->clear();
            X_MuMuIdx->clear(); X_ka1Idx->clear(); X_ka2Idx->clear();

            /// MC Analysis
            if (doMC) {
              // Gen Primary Vertex
              n_genEvtVtx = 0;
              genEvtVtx_X->clear(); genEvtVtx_Y->clear(); genEvtVtx_Z->clear();
              genEvtVtx_particles->clear();
              n_XAncestors->clear();
              nMCAll = 0, nMCX = 0; //nMCXVtx = 0;
              // Gen Primary Vertex
              PriVtxGen_X->clear(); PriVtxGen_Y->clear(); PriVtxGen_Z->clear();
              PriVtxGen_EX->clear(); PriVtxGen_EY->clear(); PriVtxGen_EZ->clear();
              PriVtxGen_Chi2->clear(); PriVtxGen_CL->clear(); PriVtxGen_Ndof->clear();
              PriVtxGen_tracks->clear();

              MCPdgIdAll->clear(); MCDanNumAll->clear();
              MCJPsiPx->clear(); MCJPsiPy->clear(); MCJPsiPz->clear();
              MCmupPx->clear(); MCmupPy->clear(); MCmupPz->clear();
              MCmumPx->clear(); MCmumPy->clear(); MCmumPz->clear();
              MCPhiPx->clear(); MCPhiPy->clear(); MCPhiPz->clear();
              MCkpPx->clear(); MCkpPy->clear(); MCkpPz->clear();
              MCkmPx->clear(); MCkmPy->clear(); MCkmPz->clear();
              //MCpionPx->clear(); MCpionPy->clear(); MCpionPz->clear();
              //MCkaonPx->clear(); MCkaonPy->clear(); MCkaonPz->clear();
              //MCpionCh->clear(); MCkaonCh->clear();
              MCPx->clear(); MCPy->clear(); MCPz->clear();
            }
            if (Debug_) cout <<"after MC stuff clear" <<endl ;
            /// Primary Vertex
            priVtx_n = 0;
            priVtx_X = 0; priVtx_Y = 0; priVtx_Z = 0 ;
            priVtx_XE = 0; priVtx_YE = 0; priVtx_ZE = 0 ;
            priVtx_NormChi2 = 0; priVtx_Chi2 = 0; priVtx_CL = 0; priVtx_tracks = 0; priVtx_tracksPtSq = 0 ;
            /// MuMu cand & KK cand
            MuMuMass->clear(); MuMuVtx_CL->clear(); MuMuVtx_Chi2->clear();
            MuMuPx->clear(); MuMuPy->clear(); MuMuPz->clear();
            MuMuDecayVtx_X->clear(); MuMuDecayVtx_Y->clear(); MuMuDecayVtx_Z->clear();
            MuMuDecayVtx_XE->clear(); MuMuDecayVtx_YE->clear(); MuMuDecayVtx_ZE->clear();
            MuMuMuonTrigMatch->clear();
            KKMass->clear(); KKPx->clear(); KKPy->clear(); KKPz->clear();
            KKVtx_CL->clear(); KKVtx_Chi2->clear();
            KKDecayVtx_X->clear(); KKDecayVtx_Y->clear(); KKDecayVtx_Z->clear();
            KKDecayVtx_XE->clear(); KKDecayVtx_YE->clear(); KKDecayVtx_ZE->clear();
            /// muons from JPsi (MuMu) fit & kaons from Phi (KK) fit
            mu1_MuMu_Px->clear(); mu1_MuMu_Py->clear(); mu1_MuMu_Pz->clear(); mu1_MuMu_Chi2->clear(); mu1_MuMu_NDF->clear();
            mu2_MuMu_Px->clear(); mu2_MuMu_Py->clear(); mu2_MuMu_Pz->clear(); mu2_MuMu_Chi2->clear(); mu2_MuMu_NDF->clear();
            MuMuType->clear();
            ka1_KK_Px->clear(); ka1_KK_Py->clear(); ka1_KK_Pz->clear(); ka1_KK_Chi2->clear(); ka1_KK_NDF->clear();
            ka2_KK_Px->clear(); ka2_KK_Py->clear();  ka2_KK_Pz->clear(); ka2_KK_Chi2->clear(); ka2_KK_NDF->clear();
            DR_MuMu_K1->clear(); DR_MuMu_K2->clear(); DR_MuMuKK_K1->clear(); DR_MuMuKK_K2->clear();
            /// Primary Vertex with "MuMu correction"
            PriVtxMuMuCorr_n->clear();
            PriVtxMuMuCorr_X->clear(); PriVtxMuMuCorr_Y->clear(); PriVtxMuMuCorr_Z->clear();
            PriVtxMuMuCorr_EX->clear(); PriVtxMuMuCorr_EY->clear(); PriVtxMuMuCorr_EZ->clear();
            PriVtxMuMuCorr_Chi2->clear(); PriVtxMuMuCorr_CL->clear(); PriVtxMuMuCorr_tracks->clear();
            nTrk->clear();
            /// X candidates
            xMass->clear(); xVtx_CL->clear(); xVtx_Chi2->clear();
            xPx->clear(); xPy->clear(); xPz->clear();
            xPxE->clear(); xPyE->clear(); xPzE->clear();
            xDecayVtx_X->clear(); xDecayVtx_Y->clear(); xDecayVtx_Z->clear();
            xDecayVtx_XE->clear(); xDecayVtx_YE->clear(); xDecayVtx_ZE->clear();
            /// Muons and tracks after X candidates fit
            mu1Px_MuMuKK->clear(); mu1Py_MuMuKK->clear(); mu1Pz_MuMuKK->clear(); mu1E_MuMuKK->clear();
            mu2Px_MuMuKK->clear(); mu2Py_MuMuKK->clear(); mu2Pz_MuMuKK->clear(); mu2E_MuMuKK->clear();
            k1Px_MuMuKK->clear(); k1Py_MuMuKK->clear(); k1Pz_MuMuKK->clear(); k1E_MuMuKK->clear();
            kaon1_nsigdedx->clear(); kaon1_dedx->clear(); kaon1_dedxMass->clear(); kaon1_theo->clear(); kaon1_sigma->clear();
            kaon1_dedx_byHits->clear(); kaon1_dedxErr_byHits->clear(); kaon1_saturMeas_byHits->clear(); kaon1_Meas_byHits->clear();
            k2Px_MuMuKK->clear(); k2Py_MuMuKK->clear(); k2Pz_MuMuKK->clear(); k2E_MuMuKK->clear();
            kaon2_nsigdedx->clear(); kaon2_dedx->clear(); kaon2_dedxMass->clear(); kaon2_theo->clear(); kaon2_sigma->clear();
            kaon2_dedx_byHits->clear(); kaon2_dedxErr_byHits->clear(); kaon2_saturMeas_byHits->clear(); kaon2_Meas_byHits->clear();
            /// Primary Vertex with largest B0_cos(alpha)
            PriVtxXLess_n->clear();
            PriVtxXLess_X->clear(); PriVtxXLess_Y->clear(); PriVtxXLess_Z->clear();
            PriVtxXLess_EX->clear(); PriVtxXLess_EY->clear(); PriVtxXLess_EZ->clear();
            PriVtxXLess_Chi2->clear(); PriVtxXLess_CL->clear(); PriVtxXLess_tracks->clear();

            XLessPV_tracksPtSq->clear(); XLessPV_4tracksPtSq->clear();
            PriVtx_XCosAlpha_n->clear();
            PriVtx_XCosAlpha_X->clear(); PriVtx_XCosAlpha_Y->clear(); PriVtx_XCosAlpha_Z->clear();
            PriVtx_XCosAlpha_EX->clear(); PriVtx_XCosAlpha_EY->clear(); PriVtx_XCosAlpha_EZ->clear();
            PriVtx_XCosAlpha_Chi2->clear(); PriVtx_XCosAlpha_CL->clear(); PriVtx_XCosAlpha_tracks->clear();
            PriVtxXLess_XCosAlpha_n->clear();
            PriVtxXLess_XCosAlpha_X->clear(); PriVtxXLess_XCosAlpha_Y->clear(); PriVtxXLess_XCosAlpha_Z->clear();
            PriVtxXLess_XCosAlpha_EX->clear(); PriVtxXLess_XCosAlpha_EY->clear(); PriVtxXLess_XCosAlpha_EZ->clear();
            PriVtxXLess_XCosAlpha_Chi2->clear(); PriVtxXLess_XCosAlpha_CL->clear(); PriVtxXLess_XCosAlpha_tracks->clear();
            /// Primary Vertex with "B0 correction"
            PriVtxXCorr_n->clear();
            PriVtxXCorr_X->clear(); PriVtxXCorr_Y->clear(); PriVtxXCorr_Z->clear();
            PriVtxXCorr_EX->clear(); PriVtxXCorr_EY->clear(); PriVtxXCorr_EZ->clear();
            PriVtxXCorr_Chi2->clear(); PriVtxXCorr_CL->clear(); PriVtxXCorr_tracks->clear();
            /// Lifetime variables for B0
            xCosAlphaBS->clear(); xCosAlpha3DBS->clear(); xCTauBS->clear(); xCTauBSE->clear(); xLxyBS->clear(); xLxyBSE->clear(); xLxyzBS->clear(); xLxyzBSE->clear();
            xCosAlphaPV->clear(); xCosAlpha3DPV->clear(); xCTauPV->clear(); xCTauPVE->clear(); xLxyPV->clear(); xLxyPVE->clear(); xLxyzPV->clear(); xLxyzPVE->clear();
            xCosAlphaPVCosAlpha->clear(); xCosAlpha3DPVCosAlpha->clear(); xCTauPVCosAlpha->clear(); xCTauPVCosAlphaE->clear(); xLxyPVCosAlpha->clear(); xLxyPVCosAlphaE->clear(); xLxyzPVCosAlpha->clear(); xLxyzPVCosAlphaE->clear();
            xCosAlphaPVCosAlpha3D->clear(); xCosAlpha3DPVCosAlpha3D->clear(); xCTauPVCosAlpha3D->clear(); xCTauPVCosAlpha3DE->clear(); xLxyPVCosAlpha3D->clear(); xLxyPVCosAlpha3DE->clear(); xLxyzPVCosAlpha3D->clear(); xLxyzPVCosAlpha3DE->clear();
            xCosAlphaXLessPV->clear(); xCosAlpha3DXLessPV->clear(); xCTauXLessPV->clear() ; xCTauXLessPVE->clear() ; xLxyXLessPV->clear() ; xLxyXLessPVE->clear() ; xLxyzXLessPV->clear() ; xLxyzXLessPVE->clear() ;
            xCosAlphaXLessPVCosAlpha->clear(); xCosAlpha3DXLessPVCosAlpha->clear(); xCTauXLessPVCosAlpha->clear() ; xCTauXLessPVCosAlphaE->clear() ; xLxyXLessPVCosAlpha->clear() ; xLxyXLessPVCosAlphaE->clear() ; xLxyzXLessPVCosAlpha->clear() ; xLxyzXLessPVCosAlphaE->clear() ;
            xCosAlphaXLessPVCosAlpha3D->clear(); xCosAlpha3DXLessPVCosAlpha3D->clear(); xCTauXLessPVCosAlpha3D->clear() ; xCTauXLessPVCosAlpha3DE->clear() ; xLxyXLessPVCosAlpha3D->clear() ; xLxyXLessPVCosAlpha3DE->clear() ; xLxyzXLessPVCosAlpha3D->clear() ; xLxyzXLessPVCosAlpha3DE->clear() ;
            xCosAlphaPVX->clear(); xCTauPVX->clear(); xCTauPVXE->clear(); xLxyPVX->clear(); xLxyzPVX->clear();
            xCTauPVX_3D->clear(); xCTauPVX_3D_err->clear();
            /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
            kaon1_dxy_PV->clear(); kaon1_dz_PV->clear(); kaon2_dxy_PV->clear(); kaon2_dz_PV->clear();
            kaon1_dxy_BS->clear(); kaon1_dz_BS->clear(); kaon2_dxy_BS->clear(); kaon2_dz_BS->clear();
            kaon1_dxy_XLessPV->clear(); kaon1_dz_XLessPV->clear(); kaon2_dxy_XLessPV->clear(); kaon2_dz_XLessPV->clear();
            kaon1_dxyE->clear(); kaon1_dzE->clear(); kaon2_dxyE->clear(); kaon2_dzE->clear();

            Kaon1FromPV->clear(); Kaon2FromPV->clear();

            if (Debug_) cout <<"before muon stuff clear" <<endl ;
            /// muons
            muPx->clear(); muPy->clear(); muPz->clear(); muCharge->clear();
            muD0->clear(); muDz->clear(); muChi2->clear(); muGlChi2->clear();
            mufHits->clear(); muFirstBarrel->clear(); muFirstEndCap->clear(); muD0E->clear() ;  muDzVtxErr->clear() ; muKey->clear() ;
            muIsGlobal->clear(); muIsPF->clear();
            muDzVtx->clear(); muDxyVtx->clear(); muGlMatchedStation->clear(); muGlDzVtx->clear(); muGlDxyVtx->clear();
            nMatchedStations->clear();
            muNDF->clear(); muGlNDF->clear(); muPhits->clear(); muShits->clear(); muGlMuHits->clear(); muType->clear();
            muQual->clear(); muTrack->clear(); muNOverlap->clear(); muNSharingSegWith->clear();

            if (Debug_) cout <<"after muon stuff clear" <<endl ;
            /// tracks
            trNotRef->clear(); trRef->clear();
            trPx->clear(); trPy->clear(); trPz->clear(); trE->clear();
            trNDF->clear(); trPhits->clear(); trShits->clear(); trChi2->clear();
            trD0->clear(); trD0E->clear(); trCharge->clear();
            trQualityHighPurity->clear(); trQualityTight->clear();
            trfHits->clear(); trFirstBarrel->clear(); trFirstEndCap->clear();
            trDzVtx->clear(); trDxyVtx->clear();
            tr_nsigdedx->clear(); tr_dedx->clear(); tr_dedxMass->clear(); tr_theo->clear(); tr_sigma->clear();
            tr_dedx_byHits->clear(); tr_dedxErr_byHits->clear(); tr_saturMeas_byHits->clear(); tr_Meas_byHits->clear();

            if (Debug_) cout <<"end of branches clear" <<endl ;
          }
          //}/// analyze
          /// ------------ method called once each job just before starting event loop  ------------
          void mumukk::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
          {
          }
          void mumukk::beginJob()
          {
            edm::Service<TFileService> fs;

            xTree = fs->make<TTree>("X_data", "X(4140) Data");

            xTree->Branch("TrigRes", &trigRes);
            xTree->Branch("TrigNames", &trigNames);
            xTree->Branch("MatchTriggerNames", &MatchTriggerNames);
            xTree->Branch("L1TrigRes", &L1TT);
            xTree->Branch("evtNum", &evtNum,"evtNum/i");
            xTree->Branch("runNum", &runNum,"runNum/i");
            xTree->Branch("lumiNum", &lumiNum, "lumiNum/i");
            xTree->Branch("priVtx_n", &priVtx_n, "priVtx_n/i");
            xTree->Branch("priVtx_X", &priVtx_X, "priVtx_X/f");
            xTree->Branch("priVtx_Y", &priVtx_Y, "priVtx_Y/f");
            xTree->Branch("priVtx_Z", &priVtx_Z, "priVtx_Z/f");
            xTree->Branch("priVtx_XE", &priVtx_XE, "priVtx_XE/f");
            xTree->Branch("priVtx_YE", &priVtx_YE, "priVtx_YE/f");
            xTree->Branch("priVtx_ZE", &priVtx_ZE, "priVtx_ZE/f");
            xTree->Branch("priVtx_NormChi2",&priVtx_NormChi2, "priVtx_NormChi2/f");
            xTree->Branch("priVtx_Chi2",&priVtx_Chi2, "priVtx_Chi2/f");
            xTree->Branch("priVtx_CL",&priVtx_CL, "priVtx_CL/f");
            xTree->Branch("priVtx_tracks", &priVtx_tracks, "priVtx_tracks/i");
            xTree->Branch("priVtx_tracksPtSq", &priVtx_tracksPtSq, "priVtx_tracksPtSq/f");
            /// MC Analysis
            if (doMC) {
              // Gen Primary Vertex
              xTree->Branch("genEvtVtx_X", &genEvtVtx_X);
              xTree->Branch("genEvtVtx_Y", &genEvtVtx_Y);
              xTree->Branch("genEvtVtx_Z", &genEvtVtx_Z);
              xTree->Branch("genEvtVtx_particles", &genEvtVtx_particles);
              xTree->Branch("n_XAncestors", &n_XAncestors);
              xTree->Branch("nMCAll", &nMCAll, "nMCAll/i");
              xTree->Branch("MCPdgIdAll", &MCPdgIdAll);
              xTree->Branch("MCDanNumAll", &MCDanNumAll);
              xTree->Branch("nMCX",&nMCX,"nMCX/i");
              // Gen Primary Vertex
              xTree->Branch("PriVtxGen_X",&PriVtxGen_X);
              xTree->Branch("PriVtxGen_Y",&PriVtxGen_Y);
              xTree->Branch("PriVtxGen_Z",&PriVtxGen_Z);
              xTree->Branch("PriVtxGen_EX",&PriVtxGen_EX);
              xTree->Branch("PriVtxGen_EY",&PriVtxGen_EY);
              xTree->Branch("PriVtxGen_EZ",&PriVtxGen_EZ);
              xTree->Branch("PriVtxGen_Chi2",&PriVtxGen_Chi2);
              xTree->Branch("PriVtxGen_CL",&PriVtxGen_CL);
              xTree->Branch("PriVtxGen_Ndof",&PriVtxGen_Ndof);
              xTree->Branch("PriVtxGen_tracks",&PriVtxGen_tracks);
              xTree->Branch("MCJPsiPx",&MCJPsiPx);
              xTree->Branch("MCJPsiPy",&MCJPsiPy);
              xTree->Branch("MCJPsiPz",&MCJPsiPz);
              xTree->Branch("MCmupPx",&MCmupPx);
              xTree->Branch("MCmupPy",&MCmupPy);
              xTree->Branch("MCmupPz",&MCmupPz);
              xTree->Branch("MCmumPx",&MCmumPx);
              xTree->Branch("MCmumPy",&MCmumPy);
              xTree->Branch("MCmumPz",&MCmumPz);
              xTree->Branch("MCPhiPx",&MCPhiPx);
              xTree->Branch("MCPhiPy",&MCPhiPy);
              xTree->Branch("MCPhiPz",&MCPhiPz);
              xTree->Branch("MCkpPx",&MCkpPx);
              xTree->Branch("MCkpPy",&MCkpPy);
              xTree->Branch("MCkpPz",&MCkpPz);
              xTree->Branch("MCkmPx",&MCkmPx);
              xTree->Branch("MCkmPy",&MCkmPy);
              xTree->Branch("MCkmPz",&MCkmPz);
              //xTree->Branch("MCpionPx",&MCpionPx);
              //xTree->Branch("MCpionPy",&MCpionPy);
              //xTree->Branch("MCpionPz",&MCpionPz);
              //xTree->Branch("MCpionCh",&MCpionCh);
              //xTree->Branch("MCkaonPx",&MCkaonPx);
              //xTree->Branch("MCkaonPy",&MCkaonPy);
              //xTree->Branch("MCkaonPz",&MCkaonPz);
              //xTree->Branch("MCkaonCh",&MCkaonCh);
              xTree->Branch("MCPx", &MCPx);
              xTree->Branch("MCPy", &MCPy);
              xTree->Branch("MCPz", &MCPz);
            }
            /// generic tracks
            xTree->Branch("trNotRef", &trNotRef);
            xTree->Branch("trRef", &trRef);
            xTree->Branch("trackPx", &trPx);
            xTree->Branch("trackPy", &trPy);
            xTree->Branch("trackPz", &trPz);
            xTree->Branch("trackEnergy", &trE);
            xTree->Branch("trackNDF", &trNDF);
            xTree->Branch("trackPhits", &trPhits);
            xTree->Branch("trackShits", &trShits);
            xTree->Branch("trackChi2", &trChi2);
            xTree->Branch("trackD0", &trD0);
            xTree->Branch("trackD0Err", &trD0E);
            xTree->Branch("trackCharge", &trCharge);
            xTree->Branch("TrackHighPurity", &trQualityHighPurity);
            xTree->Branch("TrackTight", &trQualityTight);
            xTree->Branch("trackfHits", &trfHits);
            xTree->Branch("trackFirstBarrel", &trFirstBarrel);
            xTree->Branch("trackFirstEndCap", &trFirstEndCap);
            xTree->Branch("trackDzVtx", &trDzVtx);
            xTree->Branch("trackDxyVtx", &trDxyVtx);
            xTree->Branch("tr_nsigdedx", &tr_nsigdedx);
            xTree->Branch("tr_dedx", &tr_dedx);
            xTree->Branch("tr_dedxMass", &tr_dedxMass);
            xTree->Branch("tr_theo", &tr_theo);
            xTree->Branch("tr_sigma", &tr_sigma);
            xTree->Branch("tr_dedx_byHits", &tr_dedx_byHits );
            xTree->Branch("tr_dedxErr_byHits", &tr_dedxErr_byHits );
            xTree->Branch("tr_saturMeas_byHits", &tr_saturMeas_byHits );
            xTree->Branch("tr_Meas_byHits", &tr_Meas_byHits );
            /// Generic muons
            xTree->Branch("nMu", &nMu, "nMu/i");
            xTree->Branch("muPx",&muPx);
            xTree->Branch("muPy",&muPy);
            xTree->Branch("muPz",&muPz);
            xTree->Branch("muCharge", &muCharge);
            xTree->Branch("muD0",&muD0);
            xTree->Branch("muDz",&muDz);
            xTree->Branch("muChi2",&muChi2);
            xTree->Branch("muNDF",&muNDF);
            xTree->Branch("muPhits",&muPhits);
            xTree->Branch("muShits",&muShits);
            xTree->Branch("muLayersTr",&muLayersTr);
            xTree->Branch("muLayersPix",&muLayersPix);
            xTree->Branch("muD0E",&muD0E);
            xTree->Branch("muDzVtxErr",&muDzVtxErr);
            xTree->Branch("muKey",&muKey);
            xTree->Branch("muIsGlobal",&muIsGlobal);
            xTree->Branch("muIsPF",&muIsPF);
            xTree->Branch("muGlMuHits",&muGlMuHits);
            xTree->Branch("muGlChi2",&muGlChi2);
            xTree->Branch("muGlNDF",&muGlNDF);
            xTree->Branch("muGlMatchedStation",&muGlMatchedStation);
            xTree->Branch("muGlDzVtx", &muGlDzVtx);
            xTree->Branch("muGlDxyVtx", &muGlDxyVtx);
            xTree->Branch("nMatchedStations", &nMatchedStations);
            xTree->Branch("muType",&muType);
            xTree->Branch("muQual",&muQual);
            xTree->Branch("muTrack",&muTrack);
            xTree->Branch("muNOverlap",&muNOverlap);
            xTree->Branch("muNSharingSegWith",&muNSharingSegWith);
            xTree->Branch("mufHits", &mufHits);
            xTree->Branch("muFirstBarrel", &muFirstBarrel);
            xTree->Branch("muFirstEndCap", &muFirstEndCap);
            xTree->Branch("muDzVtx", &muDzVtx);
            xTree->Branch("muDxyVtx", &muDxyVtx);
            /// MuMu cand
            xTree->Branch("nMuMu",&nMuMu,"nMuMu/i");
            xTree->Branch("MuMuMass",&MuMuMass);
            xTree->Branch("MuMuPx",&MuMuPx);
            xTree->Branch("MuMuPy",&MuMuPy);
            xTree->Branch("MuMuPz",&MuMuPz);
            xTree->Branch("MuMuVtx_CL",&MuMuVtx_CL);
            xTree->Branch("MuMuVtx_Chi2",&MuMuVtx_Chi2);
            xTree->Branch("MuMuDecayVtx_X",&MuMuDecayVtx_X);
            xTree->Branch("MuMuDecayVtx_Y",&MuMuDecayVtx_Y);
            xTree->Branch("MuMuDecayVtx_Z",&MuMuDecayVtx_Z);
            xTree->Branch("MuMuDecayVtx_XE",&MuMuDecayVtx_XE);
            xTree->Branch("MuMuDecayVtx_YE",&MuMuDecayVtx_YE);
            xTree->Branch("MuMuDecayVtx_ZE",&MuMuDecayVtx_ZE);
            /// muons from JPsi (MuMu) fit
            xTree->Branch("mu1Idx",&mu1Idx);
            xTree->Branch("mu2Idx",&mu2Idx);
            xTree->Branch("mu1Px_MuMu",&mu1_MuMu_Px);
            xTree->Branch("mu1Py_MuMu",&mu1_MuMu_Py);
            xTree->Branch("mu1Pz_MuMu",&mu1_MuMu_Pz);
            xTree->Branch("mu1Chi2_MuMu",&mu1_MuMu_Chi2);
            xTree->Branch("mu1NDF_MuMu",&mu1_MuMu_NDF);
            xTree->Branch("mu2Px_MuMu",&mu2_MuMu_Px);
            xTree->Branch("mu2Py_MuMu",&mu2_MuMu_Py);
            xTree->Branch("mu2Pz_MuMu",&mu2_MuMu_Pz);
            xTree->Branch("mu2Chi2_MuMu",&mu2_MuMu_Chi2);
            xTree->Branch("mu2NDF_MuMu",&mu2_MuMu_NDF);
            xTree->Branch("MuMuType",&MuMuType);
            xTree->Branch("MuMuMuonTrigMatch",&MuMuMuonTrigMatch);
            /// Primary Vertex with "MuMu correction"
            xTree->Branch("PriVtxMuMuCorr_n", &PriVtxMuMuCorr_n);
            xTree->Branch("PriVtxMuMuCorr_X", &PriVtxMuMuCorr_X);
            xTree->Branch("PriVtxMuMuCorr_Y", &PriVtxMuMuCorr_Y);
            xTree->Branch("PriVtxMuMuCorr_Z", &PriVtxMuMuCorr_Z);
            xTree->Branch("PriVtxMuMuCorr_EX", &PriVtxMuMuCorr_EX);
            xTree->Branch("PriVtxMuMuCorr_EY", &PriVtxMuMuCorr_EY);
            xTree->Branch("PriVtxMuMuCorr_EZ", &PriVtxMuMuCorr_EZ);
            xTree->Branch("PriVtxMuMuCorr_Chi2", &PriVtxMuMuCorr_Chi2);
            xTree->Branch("PriVtxMuMuCorr_CL", &PriVtxMuMuCorr_CL);
            xTree->Branch("PriVtxMuMuCorr_tracks", &PriVtxMuMuCorr_tracks);
            xTree->Branch("nTrk_afterMuMu", &nTrk);
            /// KK cand
            xTree->Branch("nKK",&nKK,"nKK/i");
            xTree->Branch("KKMass",&KKMass);
            xTree->Branch("KKPx",&KKPx);
            xTree->Branch("KKPy",&KKPy);
            xTree->Branch("KKPz",&KKPz);
            xTree->Branch("KKVtx_CL",&KKVtx_CL);
            xTree->Branch("KKVtx_Chi2",&KKVtx_Chi2);
            xTree->Branch("KKDecayVtx_X",&KKDecayVtx_X);
            xTree->Branch("KKDecayVtx_Y",&KKDecayVtx_Y);
            xTree->Branch("KKDecayVtx_Z",&KKDecayVtx_Z);
            xTree->Branch("KKDecayVtx_XE",&KKDecayVtx_XE);
            xTree->Branch("KKDecayVtx_YE",&KKDecayVtx_YE);
            xTree->Branch("KKDecayVtx_ZE",&KKDecayVtx_ZE);
            /// kaons from Phi (KK) fit
            xTree->Branch("ka1Idx",&ka1Idx);
            xTree->Branch("ka2Idx",&ka2Idx);
            xTree->Branch("ka1Px_KK",&ka1_KK_Px);
            xTree->Branch("ka1Py_KK",&ka1_KK_Py);
            xTree->Branch("ka1Pz_KK",&ka1_KK_Pz);
            xTree->Branch("ka1Chi2_KK",&ka1_KK_Chi2);
            xTree->Branch("ka1NDF_KK",&ka1_KK_NDF);
            xTree->Branch("ka2Px_KK",&ka2_KK_Px);
            xTree->Branch("ka2Py_KK",&ka2_KK_Py);
            xTree->Branch("ka2Pz_KK",&ka2_KK_Pz);
            xTree->Branch("ka2Chi2_KK",&ka2_KK_Chi2);
            xTree->Branch("ka2NDF_KK",&ka2_KK_NDF);
            xTree->Branch("DR_MuMu_K1",&DR_MuMu_K1);
            xTree->Branch("DR_MuMu_K2",&DR_MuMu_K2);
            xTree->Branch("DR_MuMuKK_K1",&DR_MuMuKK_K1);
            xTree->Branch("DR_MuMuKK_K2",&DR_MuMuKK_K2);
            /// counters for X
            xTree->Branch("nX",&nX,"nX/i");
            xTree->Branch("nX_pre0",&nX_pre0,"nX_pre0/i");
            xTree->Branch("nX_pre1",&nX_pre1,"nX_pre1/i");
            xTree->Branch("nX_pre2",&nX_pre2,"nX_pre2/i");
            xTree->Branch("nX_pre3",&nX_pre3,"nX_pre3/i");
            xTree->Branch("nX_pre4",&nX_pre4,"nX_pre4/i");
            xTree->Branch("nX_pre5",&nX_pre5,"nX_pre5/i");
            xTree->Branch("nX_pre6",&nX_pre6,"nX_pre6/i");
            xTree->Branch("nX_pre7",&nX_pre7,"nX_pre7/i");
            xTree->Branch("nX_pre8",&nX_pre8,"nX_pre8/i");
            xTree->Branch("nX_pre9",&nX_pre9,"nX_pre9/i");
            xTree->Branch("nX_pre10",&nX_pre10,"nX_pre10/i");
            xTree->Branch("nX_pre11",&nX_pre11,"nX_pre11/i");
            xTree->Branch("nX_pre12",&nX_pre12,"nX_pre12/i");
            xTree->Branch("nX_pre13",&nX_pre13,"nX_pre13/i");
            xTree->Branch("nX_pre14",&nX_pre14,"nX_pre14/i");
            xTree->Branch("nX_pre15",&nX_pre15,"nX_pre15/i");
            /// B0 cand
            xTree->Branch("XMass",&xMass);
            xTree->Branch("XPx",&xPx);
            xTree->Branch("XPy",&xPy);
            xTree->Branch("XPz",&xPz);
            xTree->Branch("XPxE",&xPxE);
            xTree->Branch("XPyE",&xPyE);
            xTree->Branch("XPzE",&xPzE);
            xTree->Branch("XVtx_CL",&xVtx_CL);
            xTree->Branch("XVtx_Chi2",&xVtx_Chi2);
            xTree->Branch("XDecayVtx_X",&xDecayVtx_X);
            xTree->Branch("XDecayVtx_Y",&xDecayVtx_Y);
            xTree->Branch("XDecayVtx_Z",&xDecayVtx_Z);
            xTree->Branch("XDecayVtx_XE",&xDecayVtx_XE);
            xTree->Branch("XDecayVtx_YE",&xDecayVtx_YE);
            xTree->Branch("XDecayVtx_ZE",&xDecayVtx_ZE);
            xTree->Branch("XCosAlphaBS", &xCosAlphaBS);
            xTree->Branch("XCosAlpha3DBS", &xCosAlpha3DBS);
            xTree->Branch("XCTauBS", &xCTauBS);
            xTree->Branch("XCTauBSE", &xCTauBSE);
            xTree->Branch("XLxyBS", &xLxyBS);
            xTree->Branch("XLxyBSE", &xLxyBSE);
            xTree->Branch("XLxyzBS", &xLxyzBS);
            xTree->Branch("XLxyzBSE", &xLxyzBSE);
            xTree->Branch("XCosAlphaPV", &xCosAlphaPV);
            xTree->Branch("XCosAlpha3DPV", &xCosAlpha3DPV);
            xTree->Branch("XCTauPV", &xCTauPV);
            xTree->Branch("XCTauPVE", &xCTauPVE);
            xTree->Branch("XLxyPV", &xLxyPV);
            xTree->Branch("XLxyPVE", &xLxyPVE);
            xTree->Branch("XLxyzPV", &xLxyzPV);
            xTree->Branch("XLxyzPVE", &xLxyzPVE);
            /// Primary Vertex with largest B0_cos(alpha)
            xTree->Branch("PriVtx_XCosAlpha_n",&PriVtx_XCosAlpha_n);
            xTree->Branch("PriVtx_XCosAlpha_X",&PriVtx_XCosAlpha_X);
            xTree->Branch("PriVtx_XCosAlpha_Y",&PriVtx_XCosAlpha_Y);
            xTree->Branch("PriVtx_XCosAlpha_Z",&PriVtx_XCosAlpha_Z);
            xTree->Branch("PriVtx_XCosAlpha_EX",&PriVtx_XCosAlpha_EX);
            xTree->Branch("PriVtx_XCosAlpha_EY",&PriVtx_XCosAlpha_EY);
            xTree->Branch("PriVtx_XCosAlpha_EZ",&PriVtx_XCosAlpha_EZ);
            xTree->Branch("PriVtx_XCosAlpha_Chi2",&PriVtx_XCosAlpha_Chi2);
            xTree->Branch("PriVtx_XCosAlpha_CL",&PriVtx_XCosAlpha_CL);
            xTree->Branch("PriVtx_XCosAlpha_tracks",&PriVtx_XCosAlpha_tracks);
            xTree->Branch("XCosAlphaPVCosAlpha", &xCosAlphaPVCosAlpha);
            xTree->Branch("XCosAlpha3DPVCosAlpha", &xCosAlpha3DPVCosAlpha);
            xTree->Branch("XCTauPVCosAlpha", &xCTauPVCosAlpha);
            xTree->Branch("XCTauPVCosAlphaE", &xCTauPVCosAlphaE);
            xTree->Branch("XLxyPVCosAlpha", &xLxyPVCosAlpha);
            xTree->Branch("XLxyPVCosAlphaE", &xLxyPVCosAlphaE);
            xTree->Branch("XLxyzPVCosAlpha", &xLxyzPVCosAlpha);
            xTree->Branch("XLxyzPVCosAlphaE", &xLxyzPVCosAlphaE);
            xTree->Branch("PriVtx_XCosAlpha3D_n",&PriVtx_XCosAlpha3D_n);
            xTree->Branch("PriVtx_XCosAlpha3D_X",&PriVtx_XCosAlpha3D_X);
            xTree->Branch("PriVtx_XCosAlpha3D_Y",&PriVtx_XCosAlpha3D_Y);
            xTree->Branch("PriVtx_XCosAlpha3D_Z",&PriVtx_XCosAlpha3D_Z);
            xTree->Branch("PriVtx_XCosAlpha3D_EX",&PriVtx_XCosAlpha3D_EX);
            xTree->Branch("PriVtx_XCosAlpha3D_EY",&PriVtx_XCosAlpha3D_EY);
            xTree->Branch("PriVtx_XCosAlpha3D_EZ",&PriVtx_XCosAlpha3D_EZ);
            xTree->Branch("PriVtx_XCosAlpha3D_Chi2",&PriVtx_XCosAlpha3D_Chi2);
            xTree->Branch("PriVtx_XCosAlpha3D_CL",&PriVtx_XCosAlpha3D_CL);
            xTree->Branch("PriVtx_XCosAlpha3D_tracks",&PriVtx_XCosAlpha3D_tracks);
            xTree->Branch("XCosAlphaPVCosAlpha3D", &xCosAlphaPVCosAlpha3D);
            xTree->Branch("XCosAlpha3DPVCosAlpha3D", &xCosAlpha3DPVCosAlpha3D);
            xTree->Branch("XCTauPVCosAlpha3D", &xCTauPVCosAlpha3D);
            xTree->Branch("XCTauPVCosAlpha3DE", &xCTauPVCosAlpha3DE);
            xTree->Branch("XLxyPVCosAlpha3D", &xLxyPVCosAlpha3D);
            xTree->Branch("XLxyPVCosAlpha3DE", &xLxyPVCosAlpha3DE);
            xTree->Branch("XLxyzPVCosAlpha3D", &xLxyzPVCosAlpha3D);
            xTree->Branch("XLxyzPVCosAlpha3DE", &xLxyzPVCosAlpha3DE);

            xTree->Branch("XLessPV_tracksPtSq",&XLessPV_tracksPtSq);
            xTree->Branch("XLessPV_4tracksPtSq",&XLessPV_4tracksPtSq);
            xTree->Branch("PriVtxXLess_n",&PriVtxXLess_n);
            xTree->Branch("PriVtxXLess_X",&PriVtxXLess_X);
            xTree->Branch("PriVtxXLess_Y",&PriVtxXLess_Y);
            xTree->Branch("PriVtxXLess_Z",&PriVtxXLess_Z);
            xTree->Branch("PriVtxXLess_EX",&PriVtxXLess_EX);
            xTree->Branch("PriVtxXLess_EY",&PriVtxXLess_EY);
            xTree->Branch("PriVtxXLess_EZ",&PriVtxXLess_EZ);
            xTree->Branch("PriVtxXLess_Chi2",&PriVtxXLess_Chi2);
            xTree->Branch("PriVtxXLess_CL",&PriVtxXLess_CL);
            xTree->Branch("PriVtxXLess_tracks",&PriVtxXLess_tracks);
            xTree->Branch("XCosAlphaXLessPV", &xCosAlphaXLessPV);
            xTree->Branch("XCosAlpha3DXLessPV", &xCosAlpha3DXLessPV);
            xTree->Branch("XCTauXLessPV", &xCTauXLessPV);
            xTree->Branch("XCTauXLessPVE", &xCTauXLessPVE);
            xTree->Branch("XLxyXLessPV", &xLxyXLessPV);
            xTree->Branch("XLxyXLessPVE", &xLxyXLessPVE);
            xTree->Branch("XLxyzXLessPV", &xLxyzXLessPV);
            xTree->Branch("XLxyzXLessPVE", &xLxyzXLessPVE);
            xTree->Branch("PriVtxXLess_XCosAlpha_n",&PriVtxXLess_XCosAlpha_n);
            xTree->Branch("PriVtxXLess_XCosAlpha_X",&PriVtxXLess_XCosAlpha_X);
            xTree->Branch("PriVtxXLess_XCosAlpha_Y",&PriVtxXLess_XCosAlpha_Y);
            xTree->Branch("PriVtxXLess_XCosAlpha_Z",&PriVtxXLess_XCosAlpha_Z);
            xTree->Branch("PriVtxXLess_XCosAlpha_EX",&PriVtxXLess_XCosAlpha_EX);
            xTree->Branch("PriVtxXLess_XCosAlpha_EY",&PriVtxXLess_XCosAlpha_EY);
            xTree->Branch("PriVtxXLess_XCosAlpha_EZ",&PriVtxXLess_XCosAlpha_EZ);
            xTree->Branch("PriVtxXLess_XCosAlpha_Chi2",&PriVtxXLess_XCosAlpha_Chi2);
            xTree->Branch("PriVtxXLess_XCosAlpha_CL",&PriVtxXLess_XCosAlpha_CL);
            xTree->Branch("PriVtxXLess_XCosAlpha_tracks",&PriVtxXLess_XCosAlpha_tracks);
            xTree->Branch("XCosAlphaXLessPVCosAlpha", &xCosAlphaXLessPVCosAlpha);
            xTree->Branch("XCosAlpha3DXLessPVCosAlpha", &xCosAlpha3DXLessPVCosAlpha);
            xTree->Branch("XCTauXLessPVCosAlpha", &xCTauXLessPVCosAlpha);
            xTree->Branch("XCTauXLessPVCosAlphaE", &xCTauXLessPVCosAlphaE);
            xTree->Branch("XLxyXLessPVCosAlpha", &xLxyXLessPVCosAlpha);
            xTree->Branch("XLxyXLessPVCosAlphaE", &xLxyXLessPVCosAlphaE);
            xTree->Branch("XLxyzXLessPVCosAlpha", &xLxyzXLessPVCosAlpha);
            xTree->Branch("XLxyzXLessPVCosAlphaE", &xLxyzXLessPVCosAlphaE);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_n",&PriVtxXLess_XCosAlpha3D_n);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_X",&PriVtxXLess_XCosAlpha3D_X);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_Y",&PriVtxXLess_XCosAlpha3D_Y);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_Z",&PriVtxXLess_XCosAlpha3D_Z);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_EX",&PriVtxXLess_XCosAlpha3D_EX);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_EY",&PriVtxXLess_XCosAlpha3D_EY);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_EZ",&PriVtxXLess_XCosAlpha3D_EZ);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_Chi2",&PriVtxXLess_XCosAlpha3D_Chi2);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_CL",&PriVtxXLess_XCosAlpha3D_CL);
            xTree->Branch("PriVtxXLess_XCosAlpha3D_tracks",&PriVtxXLess_XCosAlpha3D_tracks);
            xTree->Branch("XCosAlphaXLessPVCosAlpha3D", &xCosAlphaXLessPVCosAlpha3D);
            xTree->Branch("XCosAlpha3DXLessPVCosAlpha3D", &xCosAlpha3DXLessPVCosAlpha3D);
            xTree->Branch("XCTauXLessPVCosAlpha3D", &xCTauXLessPVCosAlpha3D);
            xTree->Branch("XCTauXLessPVCosAlpha3DE", &xCTauXLessPVCosAlpha3DE);
            xTree->Branch("XLxyXLessPVCosAlpha3D", &xLxyXLessPVCosAlpha3D);
            xTree->Branch("XLxyXLessPVCosAlpha3DE", &xLxyXLessPVCosAlpha3DE);
            xTree->Branch("XLxyzXLessPVCosAlpha3D", &xLxyzXLessPVCosAlpha3D);
            xTree->Branch("XLxyzXLessPVCosAlpha3DE", &xLxyzXLessPVCosAlpha3DE);
            /// Primary Vertex with "B0 correction"
            xTree->Branch("PriVtxXCorr_n",&PriVtxXCorr_n);
            xTree->Branch("PriVtxXCorr_X",&PriVtxXCorr_X);
            xTree->Branch("PriVtxXCorr_Y",&PriVtxXCorr_Y);
            xTree->Branch("PriVtxXCorr_Z",&PriVtxXCorr_Z);
            xTree->Branch("PriVtxXCorr_EX",&PriVtxXCorr_EX);
            xTree->Branch("PriVtxXCorr_EY",&PriVtxXCorr_EY);
            xTree->Branch("PriVtxXCorr_EZ",&PriVtxXCorr_EZ);
            xTree->Branch("PriVtxXCorr_Chi2",&PriVtxXCorr_Chi2);
            xTree->Branch("PriVtxXCorr_CL",&PriVtxXCorr_CL);
            xTree->Branch("PriVtxXCorr_tracks",&PriVtxXCorr_tracks);
            /// Lifetime variables for B0
            xTree->Branch("XCosAlphaPVX", &xCosAlphaPVX);
            xTree->Branch("XCTauPVX", &xCTauPVX);
            xTree->Branch("XCTauPVXE", &xCTauPVXE);
            xTree->Branch("XLxyPVX", &xLxyPVX);
            xTree->Branch("XLxyzPVX", &xLxyzPVX);
            xTree->Branch("XCTauPVX_3D", &xCTauPVX_3D);
            xTree->Branch("XCTauPVX_3D_err", &xCTauPVX_3D_err);
            /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
            xTree->Branch("kaon1_dxy_PV", &kaon1_dxy_PV);
            xTree->Branch("kaon1_dz_PV", &kaon1_dz_PV);
            xTree->Branch("kaon2_dxy_PV", &kaon2_dxy_PV);
            xTree->Branch("kaon2_dz_PV", &kaon2_dz_PV);
            xTree->Branch("kaon1_dxy_BS", &kaon1_dxy_BS);
            xTree->Branch("kaon1_dz_BS", &kaon1_dz_BS);
            xTree->Branch("kaon2_dxy_BS", &kaon2_dxy_BS);
            xTree->Branch("kaon2_dz_BS", &kaon2_dz_BS);
            xTree->Branch("kaon1_dxy_XLessPV", &kaon1_dxy_XLessPV);
            xTree->Branch("kaon1_dz_XLessPV", &kaon1_dz_XLessPV);
            xTree->Branch("kaon2_dxy_XLessPV", &kaon2_dxy_XLessPV);
            xTree->Branch("kaon2_dz_XLessPV", &kaon2_dz_XLessPV);
            xTree->Branch("kaon1_dxyE", &kaon1_dxyE);
            xTree->Branch("kaon1_dzE", &kaon1_dzE);
            xTree->Branch("kaon2_dxyE", &kaon2_dxyE);
            xTree->Branch("kaon2_dzE", &kaon2_dzE);

            xTree->Branch("XMuMuIdx", &X_MuMuIdx);
            xTree->Branch("XKaon1Idx", &X_ka1Idx);
            xTree->Branch("XKaon2Idx", &X_ka2Idx);

            xTree->Branch("Kaon1FromPV",&Kaon1FromPV);
            xTree->Branch("Kaon2FromPV",&Kaon2FromPV );

            /// Muons and tracks after X candidates fit
            xTree->Branch("Muon1Px_MuMuKK", &mu1Px_MuMuKK);
            xTree->Branch("Muon1Py_MuMuKK", &mu1Py_MuMuKK);
            xTree->Branch("Muon1Pz_MuMuKK", &mu1Pz_MuMuKK);
            xTree->Branch("Muon1E_MuMuKK", &mu1E_MuMuKK);
            xTree->Branch("Muon2Px_MuMuKK", &mu2Px_MuMuKK);
            xTree->Branch("Muon2Py_MuMuKK", &mu2Py_MuMuKK);
            xTree->Branch("Muon2Pz_MuMuKK", &mu2Pz_MuMuKK);
            xTree->Branch("Muon2E_MuMuKK", &mu2E_MuMuKK);
            xTree->Branch("Kaon1Px_MuMuKK", &k1Px_MuMuKK);
            xTree->Branch("Kaon1Py_MuMuKK", &k1Py_MuMuKK);
            xTree->Branch("Kaon1Pz_MuMuKK", &k1Pz_MuMuKK);
            xTree->Branch("Kaon1E_MuMuKK", &k1E_MuMuKK);
            xTree->Branch("kaon1_nsigdedx", &kaon1_nsigdedx);
            xTree->Branch("kaon1_dedx", &kaon1_dedx);
            xTree->Branch("kaon1_dedxMass", &kaon1_dedxMass);
            xTree->Branch("kaon1_theo", &kaon1_theo);
            xTree->Branch("kaon1_sigma", &kaon1_sigma);
            xTree->Branch("kaon1_dedx_byHits", &kaon1_dedx_byHits);
            xTree->Branch("kaon1_dedxErr_byHits", &kaon1_dedxErr_byHits);
            xTree->Branch("kaon1_saturMeas_byHits", &kaon1_saturMeas_byHits);
            xTree->Branch("kaon1_Meas_byHits", &kaon1_Meas_byHits);
            xTree->Branch("Kaon2Px_MuMuKK", &k2Px_MuMuKK);
            xTree->Branch("Kaon2Py_MuMuKK", &k2Py_MuMuKK);
            xTree->Branch("Kaon2Pz_MuMuKK", &k2Pz_MuMuKK);
            xTree->Branch("Kaon2E_MuMuKK", &k2E_MuMuKK);
            xTree->Branch("kaon2_nsigdedx", &kaon2_nsigdedx);
            xTree->Branch("kaon2_dedx", &kaon2_dedx);
            xTree->Branch("kaon2_dedxMass", &kaon2_dedxMass);
            xTree->Branch("kaon2_theo", &kaon2_theo);
            xTree->Branch("kaon2_sigma", &kaon2_sigma);
            xTree->Branch("kaon2_dedx_byHits", &kaon2_dedx_byHits);
            xTree->Branch("kaon2_dedxErr_byHits", &kaon2_dedxErr_byHits);
            xTree->Branch("kaon2_saturMeas_byHits", &kaon2_saturMeas_byHits);
            xTree->Branch("kaon2_Meas_byHits", &kaon2_Meas_byHits);

          }/// begin Job

          /// ------------ method called once each job just after ending the event loop  ------------
          void mumukk::endJob() {
            xTree->GetDirectory()->cd();
            xTree->Write();
          }/// endjob


          bool mumukk::isAbHadron(int pdgID) {

            if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
            return false;

          }

          bool mumukk::isAMixedbHadron(int pdgID, int momPdgID) {

            if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
            (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
            return true;
            return false;

          }

          std::pair<int, float> mumukk::findCandMCInfo(reco::GenParticleRef genCand) {

            int momJpsiID = 0;
            float trueLife = -99.;
            //cout <<"externalmodule"<<endl;

            if (genCand->numberOfMothers()>0) {

              TVector3 trueVtx(0.0,0.0,0.0);
              TVector3 trueP(0.0,0.0,0.0);
              TVector3 trueVtxMom(0.0,0.0,0.0);

              trueVtx.SetXYZ(genCand->vertex().x(),genCand->vertex().y(),genCand->vertex().z());
              trueP.SetXYZ(genCand->momentum().x(),genCand->momentum().y(),genCand->momentum().z());

              bool aBhadron = false;
              reco::GenParticleRef Candmom = genCand->motherRef();       // find mothers
              if (Candmom.isNull()) {
                std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
                return result;
              } else {
                reco::GenParticleRef CandGrandMom = Candmom->motherRef();
                if (isAbHadron(Candmom->pdgId())) {
                  if (CandGrandMom.isNonnull() && isAMixedbHadron(Candmom->pdgId(),CandGrandMom->pdgId())) {
                    momJpsiID = CandGrandMom->pdgId();
                    trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
                  } else {
                    momJpsiID = Candmom->pdgId();
                    trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z());
                  }
                  aBhadron = true;
                } else {
                  if (CandGrandMom.isNonnull() && isAbHadron(CandGrandMom->pdgId())) {
                    reco::GenParticleRef JpsiGrandgrandmom = CandGrandMom->motherRef();
                    if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(CandGrandMom->pdgId(),JpsiGrandgrandmom->pdgId())) {
                      momJpsiID = JpsiGrandgrandmom->pdgId();
                      trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
                    } else {
                      momJpsiID = CandGrandMom->pdgId();
                      trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
                    }
                    aBhadron = true;
                  }
                }
                if (!aBhadron) {
                  momJpsiID = Candmom->pdgId();
                  trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z());
                }
              }

              TVector3 vdiff = trueVtx - trueVtxMom;
              trueLife = vdiff.Perp()*genCand->mass()/trueP.Perp();
            }
            std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
            return result;

          }

          double mumukk::getSigmaOfLogdEdx(double logde)
          {
            return 0.3;
          }

          float mumukk::getEnergyLoss(const reco::TrackRef & track)
          {
            if (iexception_dedx==1) return 9999.;
            const reco::DeDxDataValueMap & eloss = *energyLoss;
            return eloss[track].dEdx();
          }

          double mumukk::nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma)
          {

            // no usable dE/dx if p > 2
            double nsigma = 99 ;
            if (iexception_dedx==1) return nsigma ;

            double m  = 0.13957;
            double bg = track->p() / m;

            theo = getLogdEdx(bg);


            int nhitr = track->numberOfValidHits();
            double meas = log(getEnergyLoss(track));
            sigma = getSigmaOfLogdEdx(theo) * pow(nhitr,-0.65);
            if (sigma>0)
            nsigma = (meas-theo) / sigma ;
            return nsigma;
          }


          double mumukk::getLogdEdx(double bg)
          {
            const double a =  3.25 ;
            const double b =  0.288;
            const double c = -0.852;

            double beta = bg/sqrt(bg*bg + 1);
            double dedx = log( a/(beta*beta) + b * log(bg) + c );

            return dedx;

          }


          double mumukk::GetMass(const reco::TrackRef & track){
            double P = track->p();
            double C = 2.625;
            double K = 2.495;
            double I = getEnergyLoss(track);
            return sqrt((I-C)/K)*P;
          }


          template<typename T>
          bool mumukk::isBetterMuon(const T &mu1, const T &mu2) const {
            if (mu2.track().isNull()) return true;
            if (mu1.track().isNull()) return false;
            if (mu1.isPFMuon() != mu2.isPFMuon()) return mu1.isPFMuon();
            if (mu1.isGlobalMuon() != mu2.isGlobalMuon()) return mu1.isGlobalMuon();
            if (mu1.charge() == mu2.charge() && deltaR2(mu1,mu2) < 0.0009) {
              return mu1.track()->ptError()/mu1.track()->pt() < mu2.track()->ptError()/mu2.track()->pt();
            } else {
              int nm1 = mu1.numberOfMatches(reco::Muon::SegmentArbitration);
              int nm2 = mu2.numberOfMatches(reco::Muon::SegmentArbitration);
              return (nm1 != nm2 ? nm1 > nm2 : mu1.pt() > mu2.pt());
            }
          }

          bool mumukk::isSameMuon(const reco::Muon &mu1, const reco::Muon &mu2) const {
            return (& mu1 == & mu2) ||
            //(mu1.originalObjectRef() == mu2.originalObjectRef()) ||
            (mu1.reco::Muon::innerTrack().isNonnull() ?
            mu1.reco::Muon::innerTrack() == mu2.reco::Muon::innerTrack() :
            mu1.reco::Muon::outerTrack() == mu2.reco::Muon::outerTrack());
          }


          /// define this as a plug-in
          DEFINE_FWK_MODULE(mumukk);

          // rsync -vut --existing src/mumukk.cc semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/mumukk/src/mumukk.cc