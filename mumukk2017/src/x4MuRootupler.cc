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
#include <DataFormats/PatCandidates/interface/UserData.h>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TLorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TVector3.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

typedef math::XYZPoint Point;

class x4MuRootupler:public edm::EDAnalyzer {
      public:
	explicit x4MuRootupler(const edm::ParameterSet &);
	~x4MuRootupler() override {};
	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

      private:
	void analyze(const edm::Event &, const edm::EventSetup &) override;

	std::string file_name;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> xcand_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> refit1_;
        edm::EDGetTokenT<reco::VertexCollection>            primaryVertices_;
        edm::EDGetTokenT<edm::TriggerResults>               triggerResults_;

	bool isMC_;

	UInt_t    run;
        ULong64_t event;
        UInt_t    lumiblock;

	TLorentzVector x_p4;
	TLorentzVector jpsi_p4;
	TLorentzVector muonP_jpsi_p4;
	TLorentzVector muonM_jpsi_p4;
  TLorentzVector phi_p4;
  TLorentzVector muonP_phi_p4;
  TLorentzVector muonM_phi_p4;


	Double_t invm1S;
  Double_t probFit1S;
	Double_t y1S_nsigma;

	Double_t ele_lowerPt_pt;
	Double_t ele_higherPt_pt;
	Double_t ctpv;
	Double_t ctpv_error;
	Double_t conv_vertex;
	Double_t dz;
  Double_t dz_jpsi;
  Double_t dz_phi;

	UInt_t photon_flags;
	UInt_t numPrimaryVertices;
	UInt_t trigger;
	UInt_t rf1S_rank;

	TTree *x_tree;

	Int_t chi_pdgId;
	Int_t yns_pdgId;
	TLorentzVector gen_chi_p4;
	TLorentzVector gen_yns_p4;
  TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_photon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;
  Point xVertex;
  Point jpsVertex;
  Point phiVertex;

  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;

  TTree *upsilon_tree;
  TLorentzVector mumu_p4, muP_p4, muM_p4;
  UInt_t x_rank;

};

static const double pi0_mass =  0.134977;
static const double y1SMass  =  3.0969;

/*
// 2011 par
static const double Y_sig_par_A = 0.058;
static const double Y_sig_par_B = 0.047;
static const double Y_sig_par_C = 0.22;
*/

// 2012 par
static const double Y_sig_par_A = 62.62;
static const double Y_sig_par_B = 56.3;
static const double Y_sig_par_C = -20.77;

x4MuRootupler::x4MuRootupler(const edm::ParameterSet & iConfig):
// chi_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("chi_cand"))),
xcand_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("x_cand"))),
// refit1_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refit1S"))),
primaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("primaryVertices"))),
triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))),
isMC_(iConfig.getParameter < bool > ("isMC"))
{
    edm::Service < TFileService > fs;
    x_tree = fs->make < TTree > ("chiTree", "Tree of chic");

    x_tree->Branch("run",      &run,      "run/i");
    x_tree->Branch("event",    &event,    "event/l");
    x_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

    x_tree->Branch("x_p4",    "TLorentzVector", &x_p4);
    x_tree->Branch("trigger", &trigger,            "trigger/i");

    x_tree->Branch("jpsi_p4", "TLorentzVector", &jpsi_p4);
    x_tree->Branch("muonP_jpsi_p4",  "TLorentzVector", &muonP_jpsi_p4);
    x_tree->Branch("muonM_jpsi_p4",  "TLorentzVector", &muonM_jpsi_p4);

    x_tree->Branch("phi_p4", "TLorentzVector", &phi_p4);
    x_tree->Branch("muonP_phi_p4",  "TLorentzVector", &muonP_phi_p4);
    x_tree->Branch("muonM_phi_p4",  "TLorentzVector", &muonM_phi_p4);

    x_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");

    x_tree->Branch("dz",           &dz,           "dz/D");
    x_tree->Branch("dzjpsi",           &dz_jpsi,           "dz_jpsi/D");
    x_tree->Branch("dzphi",           &dz_phi,           "dz_phi/D");


    x_tree->Branch("xVertex",  "Point", &xVertex);
    x_tree->Branch("jpsVertex",  "Point", &jpsVertex);
    x_tree->Branch("phiVertex",  "Point", &phiVertex);

    //
    // x_tree->Branch("rf1S_chi_p4", "TLorentzVector", &rf1S_chi_p4);
    // x_tree->Branch("invm1S",      &invm1S,          "invm1S/D");
    // x_tree->Branch("probFit1S",   &probFit1S,       "probFit1S/D");
    // x_tree->Branch("y1S_nsigma",  &y1S_nsigma,      "y1S_nsigma/D");
    //
    // x_tree->Branch("ele_lowerPt_pt",  &ele_lowerPt_pt,  "ele_lowerPt_pt/D");
    // x_tree->Branch("ele_higherPt_pt", &ele_higherPt_pt, "ele_higherPt_pt/D");
    //
    // x_tree->Branch("ctpv",         &ctpv,         "ctpv/D");
    // x_tree->Branch("ctpv_error",   &ctpv_error,   "ctpv_error/D");
    // x_tree->Branch("conv_vertex",  &conv_vertex,  "conv_vertex/D");
    // x_tree->Branch("dz",           &dz,           "dz/D");
    //
    // x_tree->Branch("photon_flags", &photon_flags, "photon_flags/i");
    //
    //
    // x_tree->Branch("trigger",            &trigger,            "trigger/i");
    // x_tree->Branch("rf1S_rank",          &rf1S_rank,          "rf1S_rank/i");
    //
    // if (isMC_) {
    //    x_tree->Branch("chi_pdgId",     &chi_pdgId,        "chi_pdgId/I");
    //    x_tree->Branch("yns_pdgId",     &yns_pdgId,        "yns_pdgId/I");
    //    x_tree->Branch("gen_chi_p4",    "TLorentzVector",  &gen_chi_p4);
    //    x_tree->Branch("gen_yns_p4",    "TLorentzVector",  &gen_yns_p4);
    //    x_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
    //    x_tree->Branch("gen_photon_p4", "TLorentzVector",  &gen_photon_p4);
    //    x_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
    //    x_tree->Branch("gen_muonM_p4",  "TLorentzVector",  &gen_muonM_p4);
    // }
    // genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
    //
    // upsilon_tree = fs->make<TTree>("psiTree","Tree of Jpsi");
    // upsilon_tree->Branch("mumu_p4",  "TLorentzVector", &mumu_p4);
    // upsilon_tree->Branch("muP_p4",   "TLorentzVector", &muP_p4);
    // upsilon_tree->Branch("muM_p4",   "TLorentzVector", &muM_p4);
    // upsilon_tree->Branch("trigger",  &trigger,         "trigger/i");
    // upsilon_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
    // upsilon_tree->Branch("x_rank",&x_rank,       "x_rank/i");
}

//Check recursively if any ancestor of particle is the given one
bool x4MuRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

// ------------ method called for each event  ------------
void x4MuRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {


  edm::Handle < pat::CompositeCandidateCollection >xcand_hand;
  iEvent.getByToken(xcand_, xcand_hand);

  // edm::Handle < pat::CompositeCandidateCollection >refit1S_handle;
  // iEvent.getByToken(refit1_, refit1S_handle);

  edm::Handle < reco::VertexCollection  >primaryVertices_handle;
  iEvent.getByToken(primaryVertices_, primaryVertices_handle);

  edm::Handle < edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken(triggerResults_, triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  pat::CompositeCandidate chi_cand;
  pat::CompositeCandidate refit1S;

  edm::Handle<reco::GenParticleCollection> pruned;
  // iEvent.getByToken(genCands_,pruned);

  if (false && isMC_) {
   gen_chi_p4.SetPtEtaPhiM(0, 0, 0, 0);
   gen_yns_p4.SetPtEtaPhiM(0, 0, 0, 0);
   gen_dimuon_p4.SetPtEtaPhiM(0, 0, 0, 0);
   chi_pdgId = 0;
   for (size_t i=0; i<pruned->size(); i++) {
      int p_id = abs((*pruned)[i].pdgId());
      int p_status = (*pruned)[i].status();
      yns_pdgId = 0;
      int foundit = 0;
      if ( ( p_id == 20443 || p_id == 445 || p_id == 10441) && p_status == 2)  yns_pdgId = 443;
      if (yns_pdgId > 0) {
         chi_pdgId = p_id;
         foundit++;
         const reco::Candidate * pwave = &(*pruned)[i];
         gen_chi_p4.SetPtEtaPhiM(pwave->pt(),pwave->eta(),pwave->phi(),pwave->mass());
         for (size_t j=0; j<pwave->numberOfDaughters(); j++) {
            const reco::Candidate *dau = pwave->daughter(j);
            if (dau->pdgId() == yns_pdgId && dau->status() == 2) {
               gen_yns_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
               uint nmuons = 0;
               for (size_t k=0; k<dau->numberOfDaughters(); k++) {
                  const reco::Candidate *gdau = dau->daughter(k);
                  if (gdau->pdgId() == 13 && gdau->status()==1) {
                     nmuons++;
                     gen_muonM_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                  } else {
                     if (gdau->pdgId() == -13 && gdau->status()==1) {
                        nmuons++;
                        gen_muonP_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                     }
                  }
               }
               if (nmuons == 2 ) {
                  foundit += 3;                                  // found complete dimuon decay
                  gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;   // will account fsr
               }
            } else {
               if (dau->pdgId() == 22 && dau->status() ==1) {
                  foundit++;
                  gen_photon_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
               }  else std::cout << "Rootupler: unexpected pdg_id " << dau->pdgId() << " (" << run << "," << event << ")" << std::endl;
            }
            if (foundit == 5 ) break;                             // decay found !
         }
      }
      if (chi_pdgId && yns_pdgId && foundit==5) break;        // just one decay of this kind is expected
      else chi_pdgId = 0;
   }
   if (!chi_pdgId)  std::cout << "Rootupler does not found the given decay " << run << "," << event << std::endl;
  }

   //grab Trigger informations
   // save it in variable trigger, trigger is an int between 0 and 15, in binary it is:
   // (pass 11)(pass 8)(pass 7)(pass 5)
   // es. 11 = pass 5, 7 and 11
   // es. 4 = pass only 8

   trigger = 0;
   if (triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      unsigned int NTRIGGERS = 7;
      std::string TriggersToTest[NTRIGGERS] = {
	      "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi","HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi"};

      for (unsigned int i = 0; i < NTRIGGERS; i++) {
         for (int version = 1; version < 19; version++) {
            std::stringstream ss;
            ss << TriggersToTest[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
               trigger += (1<<i);
               break;
            }
         }
      }
    } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

    bool bestCandidateOnly_ = false;
    rf1S_rank = 0;
    photon_flags = 0; //else std::cout << "no valid chi handle" << std::endl;

    x_rank = 0;
    if (!xcand_hand.isValid()) std::cout << "Invalid" << std::endl;
    if (xcand_hand->empty()) std::cout << "Empty" << std::endl;
    if (xcand_hand.isValid() && !xcand_hand->empty()) {
      for (unsigned int i=0; i< xcand_hand->size(); i++) {
        pat::CompositeCandidate x_ = xcand_hand->at(i);

        xVertex  = x_.vertex();
        phiVertex = x_.daughter("phi")->vertex();
        jpsVertex = x_.daughter("jpsi")->vertex();

        dz = Getdz(x_,xVertex);
        dz_jpsi = x_.userFloat("dzJpsi");
        dz_phi = x_.userFloat("dzPhi");

        std::cout<<"Cycling on ups"<<std::endl;
        x_p4.SetPtEtaPhiM(x_.pt(), x_.eta(), x_.phi(), x_.mass());

        x_tree->Fill();
        x_rank++;
      }
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void x4MuRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

float x4MuRootupler::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {

  const reco::Candidate::LorentzVector& mom = c.p4();
  const reco::Candidate::Point& vtx = c.vertex();

  double dz = (vtx.Z()-p.Z()) - ((vtx.X()-p.X())*mom.X()+(vtx.Y()-p.Y())*mom.Y())/mom.Rho() * mom.Z()/mom.Rho();
  return (float) dz;

}

//define this as a plug-in
DEFINE_FWK_MODULE(x4MuRootupler);
