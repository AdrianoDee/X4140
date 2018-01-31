/*
   Package:    PsiPhiFourMuonsRootupler
   Class:      PsiPhiFourMuonsRootupler

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

class PsiPhiFourMuonsRootupler : public edm::EDAnalyzer {
   public:
      explicit PsiPhiFourMuonsRootupler(const edm::ParameterSet&);
      ~PsiPhiFourMuonsRootupler() override;

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

  // ----------member data ---------------------------
  std::string file_name;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsiphi_cand_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsiphi_rf_cand_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  jpsiphi_pdgid_, jpsi_pdgid_, phi_pdgid_;
  bool isMC_,OnlyBest_;
  std::vector<std::string>                            HLTs_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector jpsiphi_p4;
  TLorentzVector psi_p4;
  TLorentzVector phi_p4;
  TLorentzVector muonPhiP_p4;
  TLorentzVector muonJpsiN_p4;
  TLorentzVector muonJpsiP_p4;
  TLorentzVector muonPhiN_p4;

  TLorentzVector jpsiphi_rf_p4;
  TLorentzVector psi_rf_p4;
  TLorentzVector phi_rf_p4;

  Int_t    jpsiphi_charge, psi_triggerMatch, track_nvsh, track_nvph;
  Double_t jpsiphi_vProb,  jpsiphi_vChi2, jpsiphi_cosAlpha, jpsiphi_ctauPV, jpsiphi_ctauErrPV;
  Double_t track_d0, track_d0Err, track_dz, track_dxy;
  Double_t psi_vProb, psi_vChi2, psi_DCA, psi_ctauPV, psi_ctauErrPV, psi_cosAlpha, psi_nSigma;

  Bool_t muonJpsiP_isLoose, muonJpsiP_isSoft, muonJpsiP_isMedium, muonJpsiP_isHighPt;
  Bool_t muonJpsiN_isLoose, muonJpsiN_isSoft, muonJpsiN_isMedium, muonJpsiN_isHighPt;
  Bool_t muonPhiP_isLoose, muonPhiP_isSoft, muonPhiP_isMedium, muonPhiP_isHighPt;
  Bool_t muonPhiN_isLoose, muonPhiN_isSoft, muonPhiN_isMedium, muonPhiN_isHighPt;

  Bool_t muonJpsiP_isTracker, muonJpsiP_isGlobal, muonJpsiN_isTracker, muonJpsiN_isGlobal;
  Bool_t muonPhiP_isTracker, muonPhiP_isGlobal, muonPhiN_isTracker, muonPhiN_isGlobal;

  UInt_t muonJpsiP_type, muonJpsiM_type, muonPhiP_type, muonPhiM_type;
  Int_t          gen_jpsiphi_pdgId;
  TLorentzVector gen_jpsiphi_p4;
  TLorentzVector gen_psi_p4;
  TLorentzVector gen_phi_p4;
  TLorentzVector gen_muonPhiP_p4;
  TLorentzVector gen_muonJpsiN_p4;
  TLorentzVector gen_muonJpsiP_p4;
  TLorentzVector gen_muonPhiN_p4;

  TTree* jpsiphi_tree;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constants, enums and typedefs
//

//
// constructors and destructor
//
PsiPhiFourMuonsRootupler::PsiPhiFourMuonsRootupler(const edm::ParameterSet& iConfig):
        jpsiphi_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("jpsiphi_cand"))),
        jpsiphi_rf_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("jpsiphi_rf_cand"))),
        thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs"))
{
	edm::Service<TFileService> fs;
        jpsiphi_tree = fs->make<TTree>("OniaPhiTree","Tree of Onia and Phi");

        jpsiphi_tree->Branch("run",                &run,                "run/I");
        jpsiphi_tree->Branch("event",              &event,              "event/I");
        jpsiphi_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        jpsiphi_tree->Branch("trigger",            &trigger,            "trigger/I");

        jpsiphi_tree->Branch("jpsiphi_p4",   "TLorentzVector", &jpsiphi_p4);
        jpsiphi_tree->Branch("phi_p4",     "TLorentzVector", &phi_p4);
        jpsiphi_tree->Branch("psi_p4",     "TLorentzVector", &psi_p4);
        jpsiphi_tree->Branch("muonPhiN_p4",   "TLorentzVector", &muonPhiN_p4);
        jpsiphi_tree->Branch("muonJpsiN_p4",   "TLorentzVector", &muonJpsiN_p4);
        jpsiphi_tree->Branch("muonJpsiP_p4",   "TLorentzVector", &muonJpsiP_p4);
        jpsiphi_tree->Branch("muonPhiN_p4",   "TLorentzVector", &muonPhiN_p4);

        jpsiphi_tree->Branch("jpsiphi_rf_p4", "TLorentzVector", &jpsiphi_rf_p4);
        jpsiphi_tree->Branch("psi_rf_p4",   "TLorentzVector", &psi_rf_p4);
        jpsiphi_tree->Branch("phi_rf_p4",   "TLorentzVector", &phi_rf_p4);

        jpsiphi_tree->Branch("psi_vProb",        &psi_vProb,        "psi_vProb/D");
        jpsiphi_tree->Branch("psi_vNChi2",       &psi_vChi2,        "psi_vNChi2/D");
        jpsiphi_tree->Branch("psi_DCA",          &psi_DCA,          "psi_DCA/D");
        jpsiphi_tree->Branch("psi_ctauPV",       &psi_ctauPV,       "psi_ctauPV/D");
        jpsiphi_tree->Branch("psi_ctauErrPV",    &psi_ctauErrPV,    "psi_ctauErrPV/D");
        jpsiphi_tree->Branch("psi_cosAlpha",     &psi_cosAlpha,     "psi_cosAlpha/D");
        jpsiphi_tree->Branch("psi_triggerMatch", &psi_triggerMatch, "psi_triggerMatch/I");

        jpsiphi_tree->Branch("jpsiphi_vProb",      &jpsiphi_vProb,        "jpsiphi_vProb/D");
        jpsiphi_tree->Branch("jpsiphi_vChi2",      &jpsiphi_vChi2,        "jpsiphi_vChi2/D");
        jpsiphi_tree->Branch("jpsiphi_cosAlpha",   &jpsiphi_cosAlpha,     "jpsiphi_cosAlpha/D");
        jpsiphi_tree->Branch("jpsiphi_ctauPV",     &jpsiphi_ctauPV,       "jpsiphi_ctauPV/D");
        jpsiphi_tree->Branch("jpsiphi_ctauErrPV",  &jpsiphi_ctauErrPV,    "jpsiphi_ctauErrPV/D");
        jpsiphi_tree->Branch("jpsiphi_charge",     &jpsiphi_charge,       "jpsiphi_charge/I");

        jpsiphi_tree->Branch("muonJpsi",        &muonJpsi,        "muonJpsi/O");

        jpsiphi_tree->Branch("muonJpsiP_isLoose",        &muonJpsiP_isLoose,        "muonJpsiP_isLoose/O");
        jpsiphi_tree->Branch("muonJpsiP_isSoft",        &muonJpsiP_isSoft,        "muonJpsiP_isSoft/O");
        jpsiphi_tree->Branch("muonJpsiP_isMedium",        &muonJpsiP_isMedium,        "muonJpsiP_isMedium/O");
        jpsiphi_tree->Branch("muonJpsiP_isHighPt",        &muonJpsiP_isHighPt,        "muonJpsiP_isHighPt/O");

        jpsiphi_tree->Branch("muonJpsiP_isTracker",        &muonJpsiP_isTracker,        "muonJpsiP_isTracker/O");
        jpsiphi_tree->Branch("muonJpsiP_isGlobal",        &muonJpsiP_isGlobal,        "muonJpsiP_isGlobal/O");

        jpsiphi_tree->Branch("muonJpsiN_isLoose",        &muonJpsiN_isLoose,        "muonJpsiN_isLoose/O");
        jpsiphi_tree->Branch("muonJpsiN_isSoft",        &muonJpsiN_isSoft,        "muonJpsiN_isSoft/O");
        jpsiphi_tree->Branch("muonJpsiN_isMedium",        &muonJpsiN_isMedium,        "muonJpsiN_isMedium/O");
        jpsiphi_tree->Branch("muonJpsiN_isHighPt",        &muonJpsiN_isHighPt,        "muonJpsiN_isHighPt/O");

        jpsiphi_tree->Branch("muonJpsiN_isTracker",        &muonJpsiN_isTracker,        "muonJpsiN_isTracker/O");
        jpsiphi_tree->Branch("muonJpsiN_isGlobal",        &muonJpsiN_isGlobal,        "muonJpsiN_isGlobal/O");

        jpsiphi_tree->Branch("muonJpsiP_type",     &muonJpsiP_type,       "muonJpsiP_type/i");
        jpsiphi_tree->Branch("muonJpsiN_type",     &muonJpsiN_type,       "muonJpsiN_type/i");

        jpsiphi_tree->Branch("muonPhiP_isLoose",        &muonPhiP_isLoose,        "muonPhiP_isLoose/O");
        jpsiphi_tree->Branch("muonPhiP_isSoft",        &muonPhiP_isSoft,        "muonPhiP_isSoft/O");
        jpsiphi_tree->Branch("muonPhiP_isMedium",        &muonPhiP_isMedium,        "muonPhiP_isMedium/O");
        jpsiphi_tree->Branch("muonPhiP_isHighPt",        &muonPhiP_isHighPt,        "muonPhiP_isHighPt/O");

        jpsiphi_tree->Branch("muonPhiP_isTracker",        &muonPhiP_isTracker,        "muonPhiP_isTracker/O");
        jpsiphi_tree->Branch("muonPhiP_isGlobal",        &muonPhiP_isGlobal,        "muonPhiP_isGlobal/O");

        jpsiphi_tree->Branch("muonPhiN_isLoose",        &muonPhiN_isLoose,        "muonPhiN_isLoose/O");
        jpsiphi_tree->Branch("muonPhiN_isSoft",        &muonPhiN_isSoft,        "muonPhiN_isSoft/O");
        jpsiphi_tree->Branch("muonPhiN_isMedium",        &muonPhiN_isMedium,        "muonPhiN_isMedium/O");
        jpsiphi_tree->Branch("muonPhiN_isHighPt",        &muonPhiN_isHighPt,        "muonPhiN_isHighPt/O");

        jpsiphi_tree->Branch("muonPhiN_isTracker",        &muonPhiN_isTracker,        "muonPhiN_isTracker/O");
        jpsiphi_tree->Branch("muonPhiN_isGlobal",        &muonPhiN_isGlobal,        "muonPhiN_isGlobal/O");

        jpsiphi_tree->Branch("muonPhiP_type",     &muonPhiP_type,       "muonPhiP_type/i");
        jpsiphi_tree->Branch("muonPhiN_type",     &muonPhiN_type,       "muonPhiN_type/i");


	if(isMC_)
	  {
            jpsiphi_tree->Branch("gen_jpsiphi_pdgId", &gen_jpsiphi_pdgId, "gen_jpsiphi_pdgId/I");
      	    jpsiphi_tree->Branch("gen_jpsiphi_p4",    "TLorentzVector", &gen_jpsiphi_p4);
      	    jpsiphi_tree->Branch("gen_psi_p4",      "TLorentzVector", &gen_psi_p4);
      	    jpsiphi_tree->Branch("gen_phi_p4",      "TLorentzVector", &gen_phi_p4);
            jpsiphi_tree->Branch("gen_muonPhiP_p4",    "TLorentzVector", &gen_muonPhiP_p4);
            jpsiphi_tree->Branch("gen_muonJpsiN_p4",    "TLorentzVector", &gen_muonJpsiN_p4);
            jpsiphi_tree->Branch("gen_muonJpsiP_p4",    "TLorentzVector", &gen_muonJpsiP_p4);
            jpsiphi_tree->Branch("gen_muonPhiN_p4",    "TLorentzVector", &gen_muonPhiN_p4);
	  }
        genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
        packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

PsiPhiFourMuonsRootupler::~PsiPhiFourMuonsRootupler() {}

//
// member functions
//

bool PsiPhiFourMuonsRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void PsiPhiFourMuonsRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> jpsiphi_cand_handle;
  iEvent.getByToken(jpsiphi_cand_Label, jpsiphi_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> jpsiphi_rf_cand_handle;
  iEvent.getByToken(jpsiphi_rf_cand_Label, jpsiphi_rf_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  Vertex thePrimaryV;
  Vertex theBeamSpotV;

  Handle<BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  if ( primaryVertices_handle->begin() != primaryVertices_handle->end() ) {
    thePrimaryV = Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
  }

	// grab Trigger information
	// save it in variable trigger, trigger is an int between 0 and 7, in binary it is:
	// (pass 10)(pass 8)(pass 0)
	// ex. 7 = pass 0, 8 and 10
	// ex. 6 = pass 8, 10
        // ex. 1 = pass 0

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

// grabbing jpsiphi information
  if (!jpsiphi_cand_handle.isValid()) std::cout<< "No jpsiphi information " << run << "," << event <<std::endl;
  if (!jpsiphi_rf_cand_handle.isValid()) std::cout<< "No jpsiphi_rf information " << run << "," << event <<std::endl;
// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if (jpsiphi_rf_cand_handle.isValid() && jpsiphi_cand_handle.isValid()) {
    pat::CompositeCandidate jpsiphi_rf_cand, jpsiphi_cand, *jpsi_cand, *phi_cand;
    for (unsigned int i=0; i< jpsiphi_rf_cand_handle->size(); i++){
      jpsiphi_rf_cand   = jpsiphi_rf_cand_handle->at(i);
      int  bindx        = jpsiphi_rf_cand.userInt("bIndex");
      jpsiphi_vProb     = jpsiphi_rf_cand.userFloat("vProb");
      jpsiphi_vChi2     = jpsiphi_rf_cand.userFloat("vChi2");
      jpsiphi_cosAlpha  = jpsiphi_rf_cand.userFloat("cosAlpha");
      jpsiphi_ctauPV    = jpsiphi_rf_cand.userFloat("ctauPV");
      jpsiphi_ctauErrPV = jpsiphi_rf_cand.userFloat("ctauErrPV");
      jpsiphi_charge    = jpsiphi_rf_cand.charge();
      jpsiphi_rf_p4.SetPtEtaPhiM(jpsiphi_rf_cand.pt(),jpsiphi_rf_cand.eta(),jpsiphi_rf_cand.phi(),jpsiphi_rf_cand.mass());
      psi_rf_p4.SetPtEtaPhiM(jpsiphi_rf_cand.daughter("jpsi")->pt(),jpsiphi_rf_cand.daughter("jpsi")->eta(),
                              jpsiphi_rf_cand.daughter("jpsi")->phi(),jpsiphi_rf_cand.daughter("jpsi")->mass());
      phi_rf_p4.SetPtEtaPhiM(jpsiphi_rf_cand.daughter("phi")->pt(),jpsiphi_rf_cand.daughter("phi")->eta(),
                              jpsiphi_rf_cand.daughter("phi")->phi(),jpsiphi_rf_cand.daughter("phi")->mass());
      if (bindx<0 || bindx>(int) jpsiphi_cand_handle->size()) {
        std::cout << "Incorrect index for jpsiphi combination " << run << "," << event <<"," << bindx << std::endl;
        continue;
      }
      jpsiphi_cand = jpsiphi_cand_handle->at(bindx);
      jpsi_cand = dynamic_cast <pat::CompositeCandidate *>(jpsiphi_cand.daughter("jpsi"));
      phi_cand = dynamic_cast <pat::CompositeCandidate *>(jpsiphi_cand.daughter("phi"));

      psi_vProb        = jpsi_cand->userFloat("vProb");
      psi_vChi2        = jpsi_cand->userFloat("vNChi2");
      psi_DCA          = jpsi_cand->userFloat("DCA");
      psi_ctauPV       = jpsi_cand->userFloat("ppdlPV");
      psi_ctauErrPV    = jpsi_cand->userFloat("ppdlErrPV");
      psi_cosAlpha     = jpsi_cand->userFloat("cosAlpha");
      psi_triggerMatch = 0; //jpsi_cand->userInt("isTriggerMatched");

      reco::Candidate::LorentzVector vJpsiP = jpsi_cand->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vJpsiM = jpsi_cand->daughter("muon2")->p4();

      pat::Muon *jpsiPatMuonP,  *jpsiPatMuonN, *phiPatMuonP, *phiPatMuonN;

      if (jpsi_cand->daughter("muon1")->charge() < 0) {
         vJpsiP = jpsi_cand->daughter("muon2")->p4();
         vJpsiM = jpsi_cand->daughter("muon1")->p4();
         jpsiPatMuonN = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"));
         jpsiPatMuonP = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"));
      } else
      {
        jpsiPatMuonP = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"));
        jpsiPatMuonN = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"));
      }
      muonPhiN_p4.SetPtEtaPhiM(vJpsiP.pt(), vJpsiP.eta(), vJpsiP.phi(), vJpsiP.mass());
      muonJpsiN_p4.SetPtEtaPhiM(vJpsiM.pt(), vJpsiM.eta(), vJpsiM.phi(), vJpsiM.mass());

      muonJpsiP_isLoose   =  jpsiPatMuonP->isLooseMuon();
      muonJpsiP_isSoft   =  jpsiPatMuonP->isSoftMuon(thePrimaryV);
      muonJpsiP_isMedium   = jpsiPatMuonP->isMediumMuon();
      muonJpsiP_isHighPt   = jpsiPatMuonP->isHighPtMuon(thePrimaryV);
      muonJpsiP_isTracker   = jpsiPatMuonP->isTrackerMuon();
      muonJpsiP_isGlobal   = jpsiPatMuonP->isGlobalMuon();
      muonJpsiN_isLoose   = jpsiPatMuonN->isLooseMuon();
      muonJpsiN_isSoft   = jpsiPatMuonN->isSoftMuon(thePrimaryV);
      muonJpsiN_isMedium   = jpsiPatMuonN->isMediumMuon();
      muonJpsiN_isHighPt   = jpsiPatMuonN->isHighPtMuon(thePrimaryV);
      muonJpsiN_isTracker   = jpsiPatMuonN->isTrackerMuon();
      muonJpsiN_isGlobal   = jpsiPatMuonN->isGlobalMuon();
      muonJpsiP_type   =  jpsiPatMuonP->type();
      muonJpsiN_type   = jpsiPatMuonN->type();

      //double kmass = 0.4936770;
      jpsiphi_p4.SetPtEtaPhiM(jpsiphi_cand.pt(),jpsiphi_cand.eta(),jpsiphi_cand.phi(),jpsiphi_cand.mass());
      psi_p4.SetPtEtaPhiM(jpsi_cand->pt(),jpsi_cand->eta(),jpsi_cand->phi(),jpsi_cand->mass());
      phi_p4.SetPtEtaPhiM(phi_cand->pt(), phi_cand->eta(), phi_cand->phi(), phi_cand->mass());

      reco::Candidate::LorentzVector vPhiP = phi_cand->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vPhiM = phi_cand->daughter("muon2")->p4();
      if (phi_cand->daughter("muon1")->charge() < 0) {
         vPhiP = phi_cand->daughter("muon2")->p4();
         vPhiM = phi_cand->daughter("muon1")->p4();
         phiPatMuonN = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"));
         phiPatMuonP = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"));
      } else
      {
        phiPatMuonP = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"));
        phiPatMuonN = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"));
      }

      muonPhiP_isLoose   =  phiPatMuonP->isLooseMuon();
      muonPhiP_isSoft   =  phiPatMuonP->isSoftMuon(thePrimaryV);
      muonPhiP_isMedium   = phiPatMuonP->isMediumMuon();
      muonPhiP_isHighPt   = phiPatMuonP->isHighPtMuon(thePrimaryV);
      muonPhiP_isTracker   = phiPatMuonP->isTrackerMuon();
      muonPhiP_isGlobal   = phiPatMuonP->isGlobalMuon();
      muonPhiN_isLoose   = phiPatMuonN->isLooseMuon();
      muonPhiN_isSoft   = phiPatMuonN->isSoftMuon(thePrimaryV);
      muonPhiN_isMedium   = phiPatMuonN->isMediumMuon();
      muonPhiN_isHighPt   = phiPatMuonN->isHighPtMuon(thePrimaryV);
      muonPhiN_isTracker   = phiPatMuonN->isTrackerMuon();
      muonPhiN_isGlobal   = phiPatMuonN->isGlobalMuon();
      muonPhiP_type   =  phiPatMuonP->type();
      muonPhiN_type   = phiPatMuonN->type();

      muonJpsiP_p4.SetPtEtaPhiM(vPhiP.pt(), vPhiP.eta(), vPhiP.phi(), vPhiP.mass());
      muonPhiN_p4.SetPtEtaPhiM(vPhiM.pt(), vPhiM.eta(), vPhiM.phi(), vPhiM.mass());

      jpsiphi_tree->Fill();
      if (OnlyBest_) break;  // jpsiphi candidates are sorted by vProb
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void PsiPhiFourMuonsRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PsiPhiFourMuonsRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void PsiPhiFourMuonsRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void PsiPhiFourMuonsRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void PsiPhiFourMuonsRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void PsiPhiFourMuonsRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PsiPhiFourMuonsRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PsiPhiFourMuonsRootupler);