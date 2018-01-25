// -*- C++ -*-
//
// Package:    Phi2KKRootupler
// Class:      Phi2KKRootupler
//
// Description: Dump  Phi(k+ k-)  decays
//
// Author:  Alberto Sanchez Hernandez
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class Phi2KKRootupler:public edm::EDAnalyzer {
      public:
	explicit Phi2KKRootupler(const edm::ParameterSet &);
	~Phi2KKRootupler() override;

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        UInt_t getTriggerBits(const edm::Event &);
        bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
        const  reco::Candidate* GetAncestor(const reco::Candidate *);
	void analyze(const edm::Event &, const edm::EventSetup &) override;

	// ----------member data ---------------------------
	std::string file_name;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dikaon_Label;
        edm::EDGetTokenT<std::vector<pat::GenericParticle>> kaon_Label;
        edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
        edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
	std::vector<std::string> TestFilterNames_;
        std::vector<double> KKMassCuts_;
	bool isMC_;
        bool OnlyBest_;
        bool OnlyGen_;

	UInt_t    run;
	ULong64_t event;
        UInt_t    lumiblock;
        UInt_t    nkk;
        UInt_t    nkaons;
        UInt_t    trigger;
        Int_t     charge;
	Int_t     ipv;
	Int_t     ntracks, npvtracks;
	Int_t     ntracks_pv, opv;
	Float_t   dzpv;
	Float_t   dphi,deta,dr;

	TLorentzVector dikaon_p4;
	TLorentzVector kaonP_p4;
	TLorentzVector kaonN_p4;

        Float_t MassErr;
        Float_t vProb;
        Float_t DCA;
        Float_t ppdlPV;
        Float_t ppdlErrPV;
        Float_t ppdlBS;
        Float_t ppdlErrBS;
        Float_t cosAlpha;
        Float_t lxyPV;
        Float_t lxyBS;

	UInt_t numPrimaryVertices;

	TTree *kk_tree;

        Int_t mother_pdgId;
        Int_t dikaon_pdgId;
	TLorentzVector gen_dikaon_p4;
	TLorentzVector gen_kaonP_p4;
	TLorentzVector gen_kaonM_p4;

        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

};

//
// constructors and destructor
//

Phi2KKRootupler::Phi2KKRootupler(const edm::ParameterSet & iConfig):
dikaon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dikaons"))),
kaon_Label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter< edm::InputTag>("kaons"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	TestFilterNames_(iConfig.getParameter<std::vector<std::string>>("TestFilterNames")),
KKMassCuts_(iConfig.getParameter<std::vector<double>>("kk_mass_cuts")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen"))
{
  edm::Service < TFileService > fs;
  kk_tree = fs->make < TTree > ("kkTree", "Tree of Phi2KK");

  kk_tree->Branch("run",      &run,      "run/i");
  kk_tree->Branch("event",    &event,    "event/l");
  kk_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  if (!OnlyGen_) {
    kk_tree->Branch("nkk",    &nkk,    "nkk/i");
    kk_tree->Branch("nkaons",   &nkaons,   "nkaons/i");
    kk_tree->Branch("trigger",  &trigger,  "trigger/i");
    kk_tree->Branch("charge",   &charge,   "charge/I");
    kk_tree->Branch("ipv",      &ipv,      "ipv/I");
    // kk_tree->Branch("ntracks",  &ntracks,  "ntracks/I");
    kk_tree->Branch("npvtracks",&npvtracks,"npvtracks/I");
    kk_tree->Branch("ntracks_pv",&ntracks_pv,"ntracks_pv/I");
    kk_tree->Branch("opv",       &opv,       "opv/I");
    kk_tree->Branch("dzpv",     &dzpv,     "dzpv/F");
    kk_tree->Branch("dphi",     &dphi,     "dphi/F");
    kk_tree->Branch("deta",     &deta,     "deta/F");
    kk_tree->Branch("dr",       &dr,       "dr/F");

    kk_tree->Branch("dikaon_p4", "TLorentzVector", &dikaon_p4);
    kk_tree->Branch("kaonP_p4",  "TLorentzVector", &kaonP_p4);
    kk_tree->Branch("kaonN_p4",  "TLorentzVector", &kaonN_p4);

    kk_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
    kk_tree->Branch("vProb",     &vProb,      "vProb/F");
    kk_tree->Branch("DCA",       &DCA,        "DCA/F");
    kk_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
    kk_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
    kk_tree->Branch("ppdlBS",    &ppdlBS,     "ppdlBS/F");
    kk_tree->Branch("ppdlErrBS", &ppdlErrBS,  "ppdlErrBS/F");
    kk_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
    kk_tree->Branch("lxyPV",     &lxyPV,      "lxyPV/F");
    kk_tree->Branch("lxyBS",     &lxyBS,      "lxyBS/F");

    kk_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
  }

  if (isMC_ || OnlyGen_) {
     kk_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
     kk_tree->Branch("dikaon_pdgId",  &dikaon_pdgId,     "dikaon_pdgId/I");
     kk_tree->Branch("gen_dikaon_p4", "TLorentzVector",  &gen_dikaon_p4);
     kk_tree->Branch("gen_kaonP_p4",  "TLorentzVector",  &gen_kaonP_p4);
     kk_tree->Branch("gen_kaonN_p4",  "TLorentzVector",  &gen_kaonM_p4);
  }
  genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

Phi2KKRootupler::~Phi2KKRootupler() {}

//
// member functions
//

const reco::Candidate* Phi2KKRootupler::GetAncestor(const reco::Candidate* p) {
   if (p->numberOfMothers()) {
      if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
      else return p->mother(0);
   }
   return p;
}

//Check recursively if any ancestor of particle is the given one
bool Phi2KKRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 1, 2
   ex. 1 = pass 0
*/

UInt_t Phi2KKRootupler::getTriggerBits(const edm::Event& iEvent ) {
   UInt_t trigger = 0;
   edm::Handle<edm::TriggerResults> triggerResults_handle;
   iEvent.getByToken(triggerResults_Label, triggerResults_handle);
   if (triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      for (unsigned int i = 0; i < TestFilterNames_.size(); i++) {
         for (int version = 1; version < 9; version++) {
            std::stringstream ss;
            ss << TestFilterNames_[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
               trigger += (1<<i);
               break;
            }
         }
      }
   } else std::cout << "Phi2KKRootupler::getTriggerBits: *** NO triggerResults found *** " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
   return trigger;
}

// ------------ method called for each event  ------------
void Phi2KKRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<std::vector<pat::GenericParticle>> kaons;
  iEvent.getByToken(kaon_Label,kaons);

  edm::Handle<pat::CompositeCandidateCollection> dikaons;
  iEvent.getByToken(dikaon_Label,dikaons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  dikaon_pdgId = 0;
  mother_pdgId = 0;
  nkk  = 0;
  nkaons = 0;

  dikaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  kaonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  kaonN_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_dikaon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_kaonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_kaonM_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packCands_,  packed);

  if ( (isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid() ) {
    for (size_t i=0; i<pruned->size(); i++) {
      const reco::Candidate *aonia = &(*pruned)[i];
      if ( (abs(aonia->pdgId()) == 333) && (aonia->status() == 2) ) {
        int foundit = 1;
        dikaon_pdgId = aonia->pdgId();
        for ( size_t j=0; j<packed->size(); j++ ) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
          const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
          const reco::Candidate * d = &(*packed)[j];
          if ( motherInPrunedCollection != nullptr && (d->pdgId() ==  321 ) && isAncestor(aonia , motherInPrunedCollection) ) {
            gen_kaonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
            foundit++;
          }
          if ( motherInPrunedCollection != nullptr && (d->pdgId() == -321 ) && isAncestor(aonia , motherInPrunedCollection) ) {
            gen_kaonP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
            foundit++;
          }
          if ( foundit == 3 ) break;
        }
        if ( foundit == 3 ) {
          gen_dikaon_p4 = gen_kaonM_p4 + gen_kaonP_p4;   // this should take into account FSR
          mother_pdgId  = GetAncestor(aonia)->pdgId();
          break;
        } else dikaon_pdgId = 0;
      }  // if ( p_id
    } // for (size
    if ( ! dikaon_pdgId ) std::cout << "Phi2KKRootupler: does not found the given decay " << run << "," << event << std::endl; // sanity check
  }  // end if isMC

  float KKMassMax_ = KKMassCuts_[1];
  float KKMassMin_ = KKMassCuts_[0];

  bool already_stored = false;
  if ( ! OnlyGen_ ) { // we will look for dikaons, then for kaons
    if ( dikaons.isValid() && !dikaons->empty()) {
      for ( pat::CompositeCandidateCollection::const_iterator dikaonCand = dikaons->begin(); dikaonCand != dikaons->end(); ++dikaonCand ) {
        if (dikaonCand->mass() > KKMassMin_ && dikaonCand->mass() < KKMassMax_ && dikaonCand->charge() == 0) {
          dikaon_p4.SetPtEtaPhiM(dikaonCand->pt(),dikaonCand->eta(),dikaonCand->phi(),dikaonCand->mass());
          reco::Candidate::LorentzVector vP = dikaonCand->daughter("kaon1")->p4();
          reco::Candidate::LorentzVector vM = dikaonCand->daughter("kaon2")->p4();
          if ( dikaonCand->daughter("kaon1")->charge() < 0 ) {
              vP = dikaonCand->daughter("kaon2")->p4();
              vM = dikaonCand->daughter("kaon1")->p4();
          }
          kaonP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
          kaonN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());
          MassErr = dikaonCand->userFloat("MassErr");
          vProb = dikaonCand->userFloat("vProb");
          DCA = -1.;
          if (dikaonCand->hasUserFloat("DCA"))  DCA = dikaonCand->userFloat("DCA");
          ppdlPV = dikaonCand->userFloat("ppdlPV");
          ppdlErrPV = dikaonCand->userFloat("ppdlErrPV");
          ppdlBS = dikaonCand->userFloat("ppdlBS");
          ppdlErrBS = dikaonCand->userFloat("ppdlErrBS");
          cosAlpha = dikaonCand->userFloat("cosAlpha");
          charge = dikaonCand->charge();
          TVector3 pperp(dikaonCand->px(),dikaonCand->py(),0);
          lxyPV = ppdlPV * pperp.Perp() / dikaonCand->mass();
          lxyBS = ppdlBS * pperp.Perp() / dikaonCand->mass();
	  ipv   = dikaonCand->userInt("iPV");
	  dzpv  = dikaonCand->userFloat("dzPV");
	  dphi  = vP.Phi() - vM.Phi();
	  while (dphi < -M_PI) dphi += 2*M_PI;
	  while (dphi >  M_PI) dphi -= 2*M_PI;
	  deta  = vP.Eta() - vM.Eta();
	  //dr    = sqrt(pow(dphi,2) + pow(deta,2));
	  dr      = dikaonCand->userFloat("deltar");
	  // ntracks = dikaonCand->userInt("ntracks");
	  npvtracks = dikaonCand->userInt("npvtracks");
	  ntracks_pv = dikaonCand->userInt("ntracks_pv");
	  opv   = dikaonCand->userInt("oniaPV");
          nkk++;
          if (OnlyBest_) break;
          else {
            kk_tree->Fill();   // be aware, we are storing all combinations
            already_stored = true;
          }
        }
      }
    } //..else {
      if ( nkk == 0 && kaons.isValid() && !kaons->empty() ) {
        int mcharge1 = 0, mcharge2 = 0;
        reco::Candidate::LorentzVector v1, v2;
        for ( std::vector<pat::GenericParticle>::const_iterator kaonCand = kaons->begin(); kaonCand!= kaons->end(); ++kaonCand ) {
          nkaons++;
          if (nkaons == 1) {
            mcharge1 = kaonCand->charge();
            v1 = kaonCand->p4();
          } else {
            if ( mcharge1*kaonCand->charge() < 0  && mcharge2 == 0 ) {
              mcharge2 = kaonCand->charge();
              v2 = kaonCand->p4();
              nkaons = 2;
              break;    // we store only 2 kaons
            }
          }
        }
        if ( mcharge1 > 0 ) {
          kaonP_p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
          if (mcharge2 < 0 ) kaonN_p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
        } else {
          kaonN_p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
          if (mcharge2 > 0 ) kaonP_p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
        }
      }
  }  // !OnlyGen_

  if ( !already_stored ) {  // we have to make sure, we are not double storing an combination
    if ( !isMC_ ) {
      if ( nkk > 0 ) kk_tree->Fill();   // if not MC filter out
    } else kk_tree->Fill();
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phi2KKRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phi2KKRootupler);
