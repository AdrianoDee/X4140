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

  // ----------member data ---------------------------
  std::string file_name;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> oniat_cand_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> oniat_rf_cand_Label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  oniat_pdgid_, onia_pdgid_, trak_pdgid_;
  bool isMC_,OnlyBest_;
  std::vector<std::string>                            HLTs_;

  UInt_t run,event,numPrimaryVertices, trigger;

  TLorentzVector oniat_p4;
  TLorentzVector psi_p4;
  TLorentzVector phi_p4;
  TLorentzVector muonp_p4;
  TLorentzVector muonn_p4;
  TLorentzVector kaonp_p4;
  TLorentzVector kaonn_p4;

  TLorentzVector oniat_rf_p4;
  TLorentzVector psi_rf_p4;
  TLorentzVector phi_rf_p4;

  Int_t    oniat_charge, psi_triggerMatch, track_nvsh, track_nvph;
  Double_t oniat_vProb,  oniat_vChi2, oniat_cosAlpha, oniat_ctauPV, oniat_ctauErrPV;
  Double_t track_d0, track_d0Err, track_dz, track_dxy;
  Double_t psi_vProb, psi_vChi2, psi_DCA, psi_ctauPV, psi_ctauErrPV, psi_cosAlpha, psi_nSigma;

  Int_t          gen_oniat_pdgId;
  TLorentzVector gen_oniat_p4;
  TLorentzVector gen_psi_p4;
  TLorentzVector gen_phi_p4;
  TLorentzVector gen_muonp_p4;
  TLorentzVector gen_muonn_p4;
  TLorentzVector gen_kaonp_p4;
  TLorentzVector gen_kaonn_p4;

  TTree* oniat_tree;
  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const double pi0_mass    =  0.1349766;
static const Double_t psi1SMass =  3.09691;
static const Double_t psi2SMass =  3.68610;
static const Double_t ups1SMass =  9.46030;

// 2011 par
//static const double Y_sig_par_A = 0.058;
//static const double Y_sig_par_B = 0.047;
//static const double Y_sig_par_C = 0.22;

// 2012 par
static const double Y_sig_par_A = 62.62;
static const double Y_sig_par_B = 56.3;
static const double Y_sig_par_C = -20.77;

//
// constructors and destructor
//
PsiTrakTrakRootupler::PsiTrakTrakRootupler(const edm::ParameterSet& iConfig):
        oniat_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("oniat_cand"))),
        oniat_rf_cand_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("oniat_rf_cand"))),
        primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
        triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	      isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
        HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs"))
{
	edm::Service<TFileService> fs;
        oniat_tree = fs->make<TTree>("OniaPhiTree","Tree of Onia and Phi");

        oniat_tree->Branch("run",                &run,                "run/I");
        oniat_tree->Branch("event",              &event,              "event/I");
        oniat_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        oniat_tree->Branch("trigger",            &trigger,            "trigger/I");

        oniat_tree->Branch("oniat_p4",   "TLorentzVector", &oniat_p4);
        oniat_tree->Branch("phi_p4",     "TLorentzVector", &phi_p4);
        oniat_tree->Branch("psi_p4",     "TLorentzVector", &psi_p4);
        oniat_tree->Branch("muonp_p4",   "TLorentzVector", &muonp_p4);
        oniat_tree->Branch("muonn_p4",   "TLorentzVector", &muonn_p4);
        oniat_tree->Branch("kaonp_p4",   "TLorentzVector", &kaonp_p4);
        oniat_tree->Branch("kaonn_p4",   "TLorentzVector", &kaonn_p4);

        oniat_tree->Branch("oniat_rf_p4", "TLorentzVector", &oniat_rf_p4);
        oniat_tree->Branch("psi_rf_p4",   "TLorentzVector", &psi_rf_p4);
        oniat_tree->Branch("phi_rf_p4",   "TLorentzVector", &phi_rf_p4);

        oniat_tree->Branch("psi_vProb",        &psi_vProb,        "psi_vProb/D");
        oniat_tree->Branch("psi_vNChi2",       &psi_vChi2,        "psi_vNChi2/D");
        oniat_tree->Branch("psi_DCA",          &psi_DCA,          "psi_DCA/D");
        oniat_tree->Branch("psi_ctauPV",       &psi_ctauPV,       "psi_ctauPV/D");
        oniat_tree->Branch("psi_ctauErrPV",    &psi_ctauErrPV,    "psi_ctauErrPV/D");
        oniat_tree->Branch("psi_cosAlpha",     &psi_cosAlpha,     "psi_cosAlpha/D");
        oniat_tree->Branch("psi_nSigma",       &psi_nSigma,       "psi_nSigma/D");
        oniat_tree->Branch("psi_triggerMatch", &psi_triggerMatch, "psi_triggerMatch/I");

        oniat_tree->Branch("oniat_vProb",      &oniat_vProb,        "oniat_vProb/D");
        oniat_tree->Branch("oniat_vChi2",      &oniat_vChi2,        "oniat_vChi2/D");
        oniat_tree->Branch("oniat_cosAlpha",   &oniat_cosAlpha,     "oniat_cosAlpha/D");
        oniat_tree->Branch("oniat_ctauPV",     &oniat_ctauPV,       "oniat_ctauPV/D");
        oniat_tree->Branch("oniat_ctauErrPV",  &oniat_ctauErrPV,    "oniat_ctauErrPV/D");
        oniat_tree->Branch("oniat_charge",     &oniat_charge,       "oniat_charge/I");

	if(isMC_)
	  {
            oniat_tree->Branch("gen_oniat_pdgId", &gen_oniat_pdgId, "gen_oniat_pdgId/I");
	    oniat_tree->Branch("gen_oniat_p4",    "TLorentzVector", &gen_oniat_p4);
	    oniat_tree->Branch("gen_psi_p4",      "TLorentzVector", &gen_psi_p4);
	    oniat_tree->Branch("gen_phi_p4",      "TLorentzVector", &gen_phi_p4);
            oniat_tree->Branch("gen_muonp_p4",    "TLorentzVector", &gen_muonp_p4);
            oniat_tree->Branch("gen_muonn_p4",    "TLorentzVector", &gen_muonn_p4);
            oniat_tree->Branch("gen_kaonp_p4",    "TLorentzVector", &gen_kaonp_p4);
            oniat_tree->Branch("gen_kaonn_p4",    "TLorentzVector", &gen_kaonn_p4);
	  }
        genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
        packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
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

  edm::Handle<std::vector<pat::CompositeCandidate>> oniat_cand_handle;
  iEvent.getByToken(oniat_cand_Label, oniat_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> oniat_rf_cand_handle;
  iEvent.getByToken(oniat_rf_cand_Label, oniat_rf_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

// Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packCands_,  packed);

  if (isMC_ && packed.isValid() && pruned.isValid()) {
     gen_oniat_pdgId  = 0;
     gen_oniat_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     int foundit   = 0;
     for (size_t i=0; i<pruned->size(); i++) {
         int p_id = (*pruned)[i].pdgId();
         const reco::Candidate *aonia = &(*pruned)[i];
         if ( (abs(p_id) ==  531) && aonia->status() == 2) {
            gen_oniat_p4.SetPtEtaPhiM(aonia->pt(),aonia->eta(),aonia->phi(),aonia->mass());
            gen_oniat_pdgId = p_id;
            foundit++;
            for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
               const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
               const reco::Candidate * d = &(*packed)[j];
               /*if ( motherInPrunedCollection != nullptr && (abs(d->pdgId()) == onia_pdgid_ && d->status() == 2) && isAncestor(aonia , motherInPrunedCollection) ){
                  gen_psi_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }*/
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == 321)  && isAncestor(aonia , motherInPrunedCollection) ){
                  gen_kaonp_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == -321)  && isAncestor(aonia , motherInPrunedCollection) ){
                  gen_kaonn_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
                  gen_muonn_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
                  gen_muonp_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( foundit == 5 ) break;
            }
            if ( foundit == 5 ) {
               gen_psi_p4 = gen_muonn_p4 + gen_muonp_p4;   // this should take into account FSR
               gen_phi_p4 = gen_kaonp_p4 + gen_kaonn_p4;
               break;
            } else {
               foundit = 0;
               gen_oniat_pdgId = 0;
               gen_oniat_p4.SetPtEtaPhiM(0.,0.,0.,0.);
            }
         }  // if ( p_id
     } // for (size

//... sanity check
     if (!gen_oniat_pdgId) std::cout << "PsiTrakTrakRootupler: didn't find the given decay " << run << "," << event << std::endl;
  } // end if isMC

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

// grabbing oniat information
  if (!oniat_cand_handle.isValid()) std::cout<< "No oniat information " << run << "," << event <<std::endl;
  if (!oniat_rf_cand_handle.isValid()) std::cout<< "No oniat_rf information " << run << "," << event <<std::endl;
// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if (oniat_rf_cand_handle.isValid() && oniat_cand_handle.isValid()) {
    pat::CompositeCandidate oniat_rf_cand, oniat_cand, *onia_cand, *phi_cand;
    for (unsigned int i=0; i< oniat_rf_cand_handle->size(); i++){
      oniat_rf_cand   = oniat_rf_cand_handle->at(i);
      int   bindx     = oniat_rf_cand.userInt("bIndex");
      oniat_vProb     = oniat_rf_cand.userFloat("vProb");
      oniat_vChi2     = oniat_rf_cand.userFloat("vChi2");
      oniat_cosAlpha  = oniat_rf_cand.userFloat("cosAlpha");
      oniat_ctauPV    = oniat_rf_cand.userFloat("ctauPV");
      oniat_ctauErrPV = oniat_rf_cand.userFloat("ctauErrPV");
      oniat_charge    = oniat_rf_cand.charge();
      oniat_rf_p4.SetPtEtaPhiM(oniat_rf_cand.pt(),oniat_rf_cand.eta(),oniat_rf_cand.phi(),oniat_rf_cand.mass());
      psi_rf_p4.SetPtEtaPhiM(oniat_rf_cand.daughter("onia")->pt(),oniat_rf_cand.daughter("onia")->eta(),
                              oniat_rf_cand.daughter("onia")->phi(),oniat_rf_cand.daughter("onia")->mass());
      phi_rf_p4.SetPtEtaPhiM(oniat_rf_cand.daughter("ditrak")->pt(),oniat_rf_cand.daughter("ditrak")->eta(),
                              oniat_rf_cand.daughter("ditrak")->phi(),oniat_rf_cand.daughter("ditrak")->mass());
      if (bindx<0 || bindx>(int) oniat_cand_handle->size()) {
        std::cout << "Incorrect index for oniat combination " << run << "," << event <<"," << bindx << std::endl;
        continue;
      }
      oniat_cand = oniat_cand_handle->at(bindx);
      onia_cand = dynamic_cast <pat::CompositeCandidate *>(oniat_cand.daughter("onia"));
      phi_cand = dynamic_cast <pat::CompositeCandidate *>(oniat_cand.daughter("ditrak"));

      psi_vProb        = onia_cand->userFloat("vProb");
      psi_vChi2        = onia_cand->userFloat("vNChi2");
      psi_DCA          = onia_cand->userFloat("DCA");
      psi_ctauPV       = onia_cand->userFloat("ppdlPV");
      psi_ctauErrPV    = onia_cand->userFloat("ppdlErrPV");
      psi_cosAlpha     = onia_cand->userFloat("cosAlpha");
      psi_triggerMatch = 0; //onia_cand->userInt("isTriggerMatched");

      reco::Candidate::LorentzVector vP = onia_cand->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vM = onia_cand->daughter("muon2")->p4();
      if (onia_cand->daughter("muon1")->charge() < 0) {
         vP = onia_cand->daughter("muon2")->p4();
         vM = onia_cand->daughter("muon1")->p4();
      }
      muonp_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
      muonn_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

      //double kmass = 0.4936770;
      oniat_p4.SetPtEtaPhiM(oniat_cand.pt(),oniat_cand.eta(),oniat_cand.phi(),oniat_cand.mass());
      psi_p4.SetPtEtaPhiM(onia_cand->pt(),onia_cand->eta(),onia_cand->phi(),onia_cand->mass());
      phi_p4.SetPtEtaPhiM(phi_cand->pt(), phi_cand->eta(), phi_cand->phi(), phi_cand->mass());

      reco::Candidate::LorentzVector kP = phi_cand->daughter("trak1")->p4();
      reco::Candidate::LorentzVector kM = phi_cand->daughter("trak2")->p4();
      if (phi_cand->daughter("trak1")->charge() < 0) {
         kP = phi_cand->daughter("trak2")->p4();
         kM = phi_cand->daughter("trak1")->p4();
      }

      kaonp_p4.SetPtEtaPhiM(kP.pt(), kP.eta(), kP.phi(), kP.mass());
      kaonn_p4.SetPtEtaPhiM(kM.pt(), kM.eta(), kM.phi(), kM.mass());

      double sigma = Y_sig_par_A + Y_sig_par_B*pow(fabs(psi_p4.Rapidity()),2) + Y_sig_par_C*pow(fabs(psi_p4.Rapidity()),3);
      sigma /= 1000.;
      sigma *= psi1SMass/ups1SMass;
      psi_nSigma = fabs(psi_p4.M() - psi1SMass) / sigma;

      oniat_tree->Fill();
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
