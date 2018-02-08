/*
   Package:    PsiTrakTrakRootupler
   Class:      PsiTrakTrakRootupler

   Description: make rootuple of J/psi Track combination

   Original Author:  Alberto Sanchez Hernandez
   Created:  based on Alessandro Degano RootupleChic

*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//
// class declaration
//

class PsiTrakTrakRootupler : public edm::EDAnalyzer {
   public:
      explicit PsiTrakTrakRootupler(const edm::ParameterSet&);
      ~PsiTrakTrakRootupler() override;

      bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginJob() override ;
      void analyze(const edm::Event&, const edm::EventSetup&) override;
      void endJob() override ;

      void beginRun(edm::Run const&, edm::EventSetup const&) override;
      void endRun(edm::Run const&, edm::EventSetup const&) override;
      void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      UInt_t PsiTrakTrakRootupler::isTriggerMatched(pat::CompositeCandidate *diMuon_cand);

  // ----------member data ---------------------------
  std::string file_name;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsitrktrk_cand_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsitrktrk_rf_cand_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  jpsitrktrk_pdgid_, onia_pdgid_, trak_pdgid_;
  bool isMC_,OnlyBest_;
  std::vector<std::string>  HLTs_;
  std::vector<std::string>  HLTFilters_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector jpsitrktrk_p4;
  TLorentzVector jpsi_p4;
  TLorentzVector phi_p4;
  TLorentzVector muonp_p4;
  TLorentzVector muonn_p4;
  TLorentzVector kaonp_p4;
  TLorentzVector kaonn_p4;

  TLorentzVector jpsitrktrk_rf_p4;
  TLorentzVector jpsitrktrk_not_rf_p4;
  TLorentzVector jpsi_rf_p4, jpsi_not_rf_p4;
  TLorentzVector phi_rf_p4, phi_not_rf_p4;
  TLorentzVector muonp_rf_p4;
  TLorentzVector muonn_rf_p4;
  TLorentzVector kaonp_rf_p4;
  TLorentzVector kaonn_rf_p4;

  Int_t    jpsitrktrk_charge;

  UInt_t jpsi_triggerMatch, jpsi_triggerMatch_rf;

  Double_t jpsitrktrk_vProb,  jpsitrktrk_vChi2, jpsitrktrk_cosAlpha, jpsitrktrk_ctauPV, jpsitrktrk_ctauErrPV;
  Double_t jpsitrktrk_rf_vProb,  jpsitrktrk_rf_vChi2, jpsitrktrk_rf_cosAlpha, jpsitrktrk_rf_ctauPV, jpsitrktrk_rf_ctauErrPV;
  Double_t track_d0, track_d0Err, track_dz, track_dxy;

  Double_t jpsi_vProb, jpsi_vChi2, jpsi_DCA, jpsi_ctauPV, jpsi_ctauErrPV, jpsi_cosAlpha;
  Double_t jpsi_vProb_rf, jpsi_vChi2_rf, jpsi_DCA_rf, jpsi_ctauPV_rf, jpsi_ctauErrPV_rf, jpsi_cosAlpha_rf;

  Bool_t muonP_isLoose, muonP_isSoft, muonP_isMedium, muonP_isHighPt;
  Bool_t muonN_isLoose, muonN_isSoft, muonN_isMedium, muonN_isHighPt;

  Bool_t muonP_isTracker, muonP_isGlobal, muonN_isTracker, muonN_isGlobal;
  UInt_t muonP_type, muonN_type;

  Bool_t muonP_rf_isLoose, muonP_rf_isSoft, muonP_rf_isMedium, muonP_rf_isHighPt;
  Bool_t muonN_rf_isLoose, muonN_rf_isSoft, muonN_rf_isMedium, muonN_rf_isHighPt;

  Bool_t muonP_rf_isTracker, muonP_rf_isGlobal, muonN_rf_isTracker, muonN_rf_isGlobal;
  UInt_t muonP_rf_type, muonN_rf_type;

  Double_t track_KP_d0, track_KP_d0Err, track_KP_dz, track_KP_dxy;
  Int_t track_KP_nvsh, track_KP_nvph;

  Double_t track_KN_d0, track_KN_d0Err, track_KN_dz, track_KN_dxy;
  Int_t track_KN_nvsh, track_KN_nvph;

  Int_t jpsitrktrk_rf_bindx;

  Int_t          gen_jpsitrktrk_pdgId;
  TLorentzVector gen_jpsitrktrk_p4;
  TLorentzVector gen_jpsi_p4;
  TLorentzVector gen_phi_p4;
  TLorentzVector gen_muonp_p4;
  TLorentzVector gen_muonn_p4;
  TLorentzVector gen_kaonp_p4;
  TLorentzVector gen_kaonn_p4;

  TTree* jpsitrktrk_tree, *jpsitrktrk_tree_rf;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

UInt_t PsiTrakTrakRootupler::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
    // std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
    // if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) std::cout << std::endl << HLTFilters_[iTr] << std::endl;
  }

  return matched;
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const Double_t psi1SMass =  3.09691;


static const double Y_sig_par_C = -20.77;

//
// constructors and destructor
//
PsiTrakTrakRootupler::PsiTrakTrakRootupler(const edm::ParameterSet& iConfig):
        jpsitrktrk_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("jpsitrktrk_cand"))),
        jpsitrktrk_rf_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("jpsitrktrk_rf_cand"))),
        thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
        HLTFilters_(iConfig.getParameter<std::vector<std::string>>("filters"))
{
	      edm::Service<TFileService> fs;
        jpsitrktrk_tree = fs->make<TTree>("OniaTrkTrkTree","Tree of Onia and Phi");

        jpsitrktrk_tree->Branch("run",                &run,                "run/I");
        jpsitrktrk_tree->Branch("event",              &event,              "event/I");
        jpsitrktrk_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        jpsitrktrk_tree->Branch("trigger",            &trigger,            "trigger/I");

        //p4s
        jpsitrktrk_tree->Branch("jpsitrktrk_p4",   "TLorentzVector", &jpsitrktrk_p4);
        jpsitrktrk_tree->Branch("phi_p4",     "TLorentzVector", &phi_p4);
        jpsitrktrk_tree->Branch("jpsi_p4",     "TLorentzVector", &jpsi_p4);
        jpsitrktrk_tree->Branch("muonp_p4",   "TLorentzVector", &muonp_p4);
        jpsitrktrk_tree->Branch("muonn_p4",   "TLorentzVector", &muonn_p4);
        jpsitrktrk_tree->Branch("kaonp_p4",   "TLorentzVector", &kaonp_p4);
        jpsitrktrk_tree->Branch("kaonn_p4",   "TLorentzVector", &kaonn_p4);

        //refitted p4s
        jpsitrktrk_tree->Branch("jpsitrktrk_rf_p4",   "TLorentzVector", &jpsitrktrk_rf_p4);
        jpsitrktrk_tree->Branch("phi_rf_p4",     "TLorentzVector", &phi_rf_p4);
        jpsitrktrk_tree->Branch("jpsi_rf_p4",     "TLorentzVector", &jpsi_rf_p4);
        jpsitrktrk_tree->Branch("muonp_rf_p4",   "TLorentzVector", &muonp_rf_p4);
        jpsitrktrk_tree->Branch("muonn_rf_p4",   "TLorentzVector", &muonn_rf_p4);
        jpsitrktrk_tree->Branch("kaonp_rf_p4",   "TLorentzVector", &kaonp_rf_p4);
        jpsitrktrk_tree->Branch("kaonn_rf_p4",   "TLorentzVector", &kaonn_rf_p4);

        //2mu vertexing
        jpsitrktrk_tree->Branch("jpsi_vProb",        &jpsi_vProb,        "jpsi_vProb/D");
        jpsitrktrk_tree->Branch("jpsi_vNChi2",       &jpsi_vChi2,        "jpsi_vNChi2/D");
        jpsitrktrk_tree->Branch("jpsi_DCA",          &jpsi_DCA,          "jpsi_DCA/D");
        jpsitrktrk_tree->Branch("jpsi_ctauPV",       &jpsi_ctauPV,       "jpsi_ctauPV/D");
        jpsitrktrk_tree->Branch("jpsi_ctauErrPV",    &jpsi_ctauErrPV,    "jpsi_ctauErrPV/D");
        jpsitrktrk_tree->Branch("jpsi_cosAlpha",     &jpsi_cosAlpha,     "jpsi_cosAlpha/D");
        jpsitrktrk_tree->Branch("jpsi_triggerMatch", &jpsi_triggerMatch, "jpsi_triggerMatch/I");

        //2mu+2Trk vertexing
        jpsitrktrk_tree->Branch("jpsitrktrk_vProb",      &jpsitrktrk_vProb,        "jpsitrktrk_vProb/D");
        jpsitrktrk_tree->Branch("jpsitrktrk_vChi2",      &jpsitrktrk_vChi2,        "jpsitrktrk_vChi2/D");
        jpsitrktrk_tree->Branch("jpsitrktrk_cosAlpha",   &jpsitrktrk_cosAlpha,     "jpsitrktrk_cosAlpha/D");
        jpsitrktrk_tree->Branch("jpsitrktrk_ctauPV",     &jpsitrktrk_ctauPV,       "jpsitrktrk_ctauPV/D");
        jpsitrktrk_tree->Branch("jpsitrktrk_ctauErrPV",  &jpsitrktrk_ctauErrPV,    "jpsitrktrk_ctauErrPV/D");
        jpsitrktrk_tree->Branch("jpsitrktrk_charge",     &jpsitrktrk_charge,       "jpsitrktrk_charge/I");

        //Muon flags
        jpsitrktrk_tree->Branch("muonP_isLoose",        &muonP_isLoose,        "muonP_isLoose/O");
        jpsitrktrk_tree->Branch("muonP_isSoft",        &muonP_isSoft,        "muonP_isSoft/O");
        jpsitrktrk_tree->Branch("muonP_isMedium",        &muonP_isMedium,        "muonP_isMedium/O");
        jpsitrktrk_tree->Branch("muonP_isHighPt",        &muonP_isHighPt,        "muonP_isHighPt/O");

        jpsitrktrk_tree->Branch("muonP_isTracker",        &muonP_isTracker,        "muonP_isTracker/O");
        jpsitrktrk_tree->Branch("muonP_isGlobal",        &muonP_isGlobal,        "muonP_isGlobal/O");

        jpsitrktrk_tree->Branch("muonN_isLoose",        &muonN_isLoose,        "muonN_isLoose/O");
        jpsitrktrk_tree->Branch("muonN_isSoft",        &muonN_isSoft,        "muonN_isSoft/O");
        jpsitrktrk_tree->Branch("muonN_isMedium",        &muonN_isMedium,        "muonN_isMedium/O");
        jpsitrktrk_tree->Branch("muonN_isHighPt",        &muonN_isHighPt,        "muonN_isHighPt/O");

        jpsitrktrk_tree->Branch("muonN_isTracker",        &muonN_isTracker,        "muonN_isTracker/O");
        jpsitrktrk_tree->Branch("muonN_isGlobal",        &muonN_isGlobal,        "muonN_isGlobal/O");

        jpsitrktrk_tree->Branch("muonP_type",     &muonP_type,       "muonP_type/i");
        jpsitrktrk_tree->Branch("muonN_type",     &muonN_type,       "muonN_type/i");

        //Track flags

        jpsitrktrk_tree->Branch("track_KN_d0",    &track_KN_d0,    "track_KN_d0/D");
        jpsitrktrk_tree->Branch("track_KN_d0Err", &track_KN_d0Err, "track_KN_d0Err/D");
        jpsitrktrk_tree->Branch("track_KN_dz",    &track_KN_dz,    "track_KN_dz/D");
        jpsitrktrk_tree->Branch("track_KN_dxy",   &track_KN_dxy,   "track_KN_dxy/D");
        jpsitrktrk_tree->Branch("track_KN_nvsh",  &track_KN_nvsh,  "track_KN_nvsh/I");
        jpsitrktrk_tree->Branch("track_KN_nvph",  &track_KN_nvph,  "track_KN_nvph/I");

        jpsitrktrk_tree->Branch("track_KP_d0",    &track_KP_d0,    "track_KP_d0/D");
        jpsitrktrk_tree->Branch("track_KP_d0Err", &track_KP_d0Err, "track_KP_d0Err/D");
        jpsitrktrk_tree->Branch("track_KP_dz",    &track_KP_dz,    "track_KP_dz/D");
        jpsitrktrk_tree->Branch("track_KP_dxy",   &track_KP_dxy,   "track_KP_dxy/D");
        jpsitrktrk_tree->Branch("track_KP_nvsh",  &track_KP_nvsh,  "track_KP_nvsh/I");
        jpsitrktrk_tree->Branch("track_KP_nvph",  &track_KP_nvph,  "track_KP_nvph/I");

}

PsiTrakTrakRootupler::~PsiTrakTrakRootupler() {}

//
// member functions
//

bool PsiTrakTrakRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void PsiTrakTrakRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> jpsitrktrk_cand_handle;
  iEvent.getByToken(jpsitrktrk_cand_Label, jpsitrktrk_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> jpsitrktrk_rf_cand_handle;
  iEvent.getByToken(jpsitrktrk_rf_cand_Label, jpsitrktrk_rf_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  reco::Vertex thePrimaryV;
  reco::Vertex theBeamSpotV;

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;

  if ( primaryVertices_handle->begin() != primaryVertices_handle->end() ) {
    thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
  }
  else {
    thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
  }

  trigger = 0;

  if (triggerResults_handle.isValid()) {
     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
     unsigned int NTRIGGERS = HLTs_.size();

     for (unsigned int i = 0; i < NTRIGGERS; i++) {
        for (int version = 1; version < 20; version++) {
           std::stringstream ss;
           ss << HLTs_[i] << "_v" << version;
           unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
           if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
              trigger += (1<<i);
              break;
           }
        }
     }
   } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

// grabbing oniat information
  if (!jpsitrktrk_cand_handle.isValid()) std::cout<< "No oniat information " << run << "," << event <<std::endl;
  if (!jpsitrktrk_rf_cand_handle.isValid()) std::cout<< "No jpsitrktrk_rf information " << run << "," << event <<std::endl;
// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if (jpsitrktrk_rf_cand_handle.isValid() && jpsitrktrk_cand_handle.isValid()) {

    pat::CompositeCandidate jpsitrktrk_rf_cand, jpsitrktrk_cand, *onia_cand, *phi_cand, *onia_cand_rf, *phi_cand_rf;

    //Refitted Handle
    for (unsigned int i=0; i< jpsitrktrk_rf_cand_handle->size(); i++){

      jpsitrktrk_rf_cand   = jpsitrktrk_rf_cand_handle->at(i);
      jpsitrktrk_rf_bindx = jpsitrktrk_rf_cand.userInt("bIndex");

      if (jpsitrktrk_rf_bindx<0 || jpsitrktrk_rf_bindx>(int) jpsitrktrk_cand_handle->size()) {
        std::cout << "Incorrect index for oniat combination " << run << "," << event <<"," << jpsitrktrk_rf_bindx << std::endl;
        continue;
      }

      jpsitrktrk_vProb     = jpsitrktrk_rf_cand.userFloat("vProb");
      jpsitrktrk_vChi2     = jpsitrktrk_rf_cand.userFloat("vChi2");
      jpsitrktrk_cosAlpha  = jpsitrktrk_rf_cand.userFloat("cosAlpha");
      jpsitrktrk_ctauPV    = jpsitrktrk_rf_cand.userFloat("ctauPV");
      jpsitrktrk_ctauErrPV = jpsitrktrk_rf_cand.userFloat("ctauErrPV");
      jpsitrktrk_charge    = jpsitrktrk_cand.charge();

      jpsitrktrk_rf_p4.SetPtEtaPhiM(jpsitrktrk_rf_cand.pt(),jpsitrktrk_rf_cand.eta(),jpsitrktrk_rf_cand.phi(),jpsitrktrk_rf_cand.mass());
      jpsi_rf_p4.SetPtEtaPhiM(jpsitrktrk_rf_cand.daughter("onia")->pt(),jpsitrktrk_rf_cand.daughter("onia")->eta(),
                              jpsitrktrk_rf_cand.daughter("onia")->phi(),jpsitrktrk_rf_cand.daughter("onia")->mass());
      phi_rf_p4.SetPtEtaPhiM(jpsitrktrk_rf_cand.daughter("ditrak")->pt(),jpsitrktrk_rf_cand.daughter("ditrak")->eta(),
                              jpsitrktrk_rf_cand.daughter("ditrak")->phi(),jpsitrktrk_rf_cand.daughter("ditrak")->mass());

      onia_cand_rf = dynamic_cast <pat::CompositeCandidate *>(jpsitrktrk_rf_cand.daughter("onia"));
      phi_cand_rf = dynamic_cast <pat::CompositeCandidate *>(jpsitrktrk_rf_cand.daughter("ditrak"));

      reco::Candidate::LorentzVector vP = onia_cand_rf->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vM = onia_cand_rf->daughter("muon2")->p4();

      const pat::Muon *muon_rf_P, *muon_rf_N;

      if (onia_cand_rf->daughter("muon1")->charge() < 0) {
         vP = onia_cand_rf->daughter("muon2")->p4();
         vM = onia_cand_rf->daughter("muon1")->p4();
         muon_rf_N = dynamic_cast<const pat::Muon*>(onia_cand_rf->daughter("muon1"));
         muon_rf_P = dynamic_cast<const pat::Muon*>(onia_cand_rf->daughter("muon2"));
      } else
      {
        muon_rf_P = dynamic_cast<const pat::Muon*>(onia_cand_rf->daughter("muon1"));
        muon_rf_N = dynamic_cast<const pat::Muon*>(onia_cand_rf->daughter("muon2"));
      }


      muonp_rf_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
      muonn_rf_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

      kP = phi_cand_rf->daughter("trak1")->p4();
      kM = phi_cand_rf->daughter("trak2")->p4();

      if (phi_cand_rf->daughter("trak1")->charge() < 0) {
         kP = phi_cand_rf->daughter("trak2")->p4();
         kM = phi_cand_rf->daughter("trak1")->p4();
      }

      kaonp_rf_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
      kaonn_rf_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

      //unref corresponding

      jpsitrktrk_cand = jpsitrktrk_cand_handle->at(jpsitrktrk_rf_bindx);
      onia_cand = dynamic_cast <pat::CompositeCandidate *>(jpsitrktrk_cand.daughter("onia"));
      phi_cand = dynamic_cast <pat::CompositeCandidate *>(jpsitrktrk_cand.daughter("ditrak"));

      const pat::Muon *muonP, *muonN;

      if (onia_cand_rf->daughter("muon1")->charge() < 0) {
         vP = onia_cand->daughter("muon2")->p4();
         vM = onia_cand->daughter("muon1")->p4();
         muonN = dynamic_cast<const pat::Muon*>(onia_cand->daughter("muon1"));
         muonP = dynamic_cast<const pat::Muon*>(onia_cand->daughter("muon2"));
      } else
      {
        muonP = dynamic_cast<const pat::Muon*>(onia_cand->daughter("muon1"));
        muonN = dynamic_cast<const pat::Muon*>(onia_cand->daughter("muon2"));
      }

      muonP_isLoose    =  muonP->isLooseMuon();
      muonP_isSoft     =  muonP->isSoftMuon(thePrimaryV);
      muonP_isMedium   = muonP->isMediumMuon();
      muonP_isHighPt   = muonP->isHighPtMuon(thePrimaryV);
      muonP_isTracker  = muonP->isTrackerMuon();
      muonP_isGlobal   = muonP->isGlobalMuon();
      muonN_isLoose    = muonN->isLooseMuon();
      muonN_isSoft     = muonN->isSoftMuon(thePrimaryV);
      muonN_isMedium   = muonN->isMediumMuon();
      muonN_isHighPt   = muonN->isHighPtMuon(thePrimaryV);
      muonN_isTracker  = muonN->isTrackerMuon();
      muonN_isGlobal   = muonN->isGlobalMuon();
      muonP_type       = muonP->type();
      muonN_type       = muonN->type();

      muonp_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
      muonn_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

      reco::Candidate::LorentzVector kP = phi_cand->daughter("trak1")->p4();
      reco::Candidate::LorentzVector kM = phi_cand->daughter("trak2")->p4();

      track_KN_d0      = phi_cand->userFloat("trak_d0_2");
      track_KN_d0Err   = phi_cand->userFloat("trak_d0err_2");
      track_KN_dz      = phi_cand->userFloat("trak_dz_2");
      track_KN_dxy     = phi_cand->userFloat("trak_dxy_2");
      track_KN_nvsh    = phi_cand->userInt("trak_nvsh_2");
      track_KN_nvph    = phi_cand->userInt("trak_nvph_2");

      track_KP_d0      = phi_cand->userFloat("trak_d0_1");
      track_KP_d0Err   = phi_cand->userFloat("trak_d0err_1");
      track_KP_dz      = phi_cand->userFloat("trak_dz_1");
      track_KP_dxy     = phi_cand->userFloat("trak_dxy_1");
      track_KP_nvsh    = phi_cand->userInt("trak_nvsh_1");
      track_KP_nvph    = phi_cand->userInt("trak_nvph_1");

      if (phi_cand->daughter("trak1")->charge() < 0) {
         kP = phi_cand->daughter("trak2")->p4();
         kM = phi_cand->daughter("trak1")->p4();
         track_KN_d0      = phi_cand->userFloat("trak_d0_1");
         track_KN_d0Err   = phi_cand->userFloat("trak_d0err_1");
         track_KN_dz      = phi_cand->userFloat("trak_dz_1");
         track_KN_dxy     = phi_cand->userFloat("trak_dxy_1");
         track_KN_nvsh    = phi_cand->userInt("trak_nvsh_1");
         track_KN_nvph    = phi_cand->userInt("trak_nvph_1");

         track_KP_d0      = phi_cand->userFloat("trak_d0_2");
         track_KP_d0Err   = phi_cand->userFloat("trak_d0err_2");
         track_KP_dz      = phi_cand->userFloat("trak_dz_2");
         track_KP_dxy     = phi_cand->userFloat("trak_dxy_2");
         track_KP_nvsh    = phi_cand->userInt("trak_nvsh_2");
         track_KP_nvph    = phi_cand->userInt("trak_nvph_2");
      }

      kaonp_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
      kaonn_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

      //double kmass = 0.4936770;
      jpsitrktrk_p4.SetPtEtaPhiM(jpsitrktrk_cand.pt(),jpsitrktrk_cand.eta(),jpsitrktrk_cand.phi(),jpsitrktrk_cand.mass());
      jpsi_p4.SetPtEtaPhiM(onia_cand->pt(),onia_cand->eta(),onia_cand->phi(),onia_cand->mass());
      phi_p4.SetPtEtaPhiM(phi_cand->pt(), phi_cand->eta(), phi_cand->phi(), phi_cand->mass());

      jpsi_vProb_rf        = onia_cand->userFloat("vProb");
      jpsi_vChi2_rf        = onia_cand->userFloat("vNChi2");
      jpsi_DCA_rf          = onia_cand->userFloat("DCA");
      jpsi_ctauPV_rf       = onia_cand->userFloat("ppdlPV");
      jpsi_ctauErrPV_rf    = onia_cand->userFloat("ppdlErrPV");
      jpsi_cosAlpha_rf     = onia_cand->userFloat("cosAlpha");
      jpsi_triggerMatch = PsiTrakTrakRootupler::isTriggerMatched(onia_cand);

      jpsitrktrk_tree->Fill();

      if (OnlyBest_) break;  // oniat candidates are sorted by vProb
    }

  }

}

// ------------ method called once each job just before starting event loop  ------------
void PsiTrakTrakRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PsiTrakTrakRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void PsiTrakTrakRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void PsiTrakTrakRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void PsiTrakTrakRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void PsiTrakTrakRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PsiTrakTrakRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PsiTrakTrakRootupler);
