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

class diMuonRootupler:public edm::EDAnalyzer {
      public:
	explicit diMuonRootupler(const edm::ParameterSet &);
	~diMuonRootupler() override {};
	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

      private:
	void analyze(const edm::Event &, const edm::EventSetup &) override;

	    std::string file_name;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> jcand_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> refit1_;
        edm::EDGetTokenT<reco::VertexCollection>            primaryVertices_;
        edm::EDGetTokenT<edm::TriggerResults>               triggerResults_;
        std::vector<std::string>                            HLTs_;
	      bool isMC_;

	UInt_t run;
  ULong64_t event;
  UInt_t lumiblock;
  UInt_t trigger;
  UInt_t numPrimaryVertices;
  UInt_t countTksOfPV;

	TLorentzVector p4;

  Double_t cosAlpha, cosAlphaMuLess, ctauErrPV, ctauPV, ctauPVMuLess, ctauErrPVMuLess;
  Double_t ctauErrBS, ctauBS, vNChi2, vProb, sumPTPV;
  Double_t vertexWeight, dz, dz_diMuon, dz_phi;
  Double_t MassErr;

	TTree *mumu_tree;

  Point vertex;
  Point jpsVertex;
  Point phiVertex;
  reco::Vertex commonVertex;
  Point PVwithmuons;
  reco::Vertex muLessVertex;

  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;

  TTree *upsilon_tree;
  TLorentzVector muP_p4, muM_p4;
  UInt_t j_rank;

};

diMuonRootupler::diMuonRootupler(const edm::ParameterSet & iConfig):
// chi_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("chi_cand"))),
jcand_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("dimuon_cand"))),
primaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("primaryVertices"))),
triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
isMC_(iConfig.getParameter < bool > ("isMC"))
{

    edm::Service < TFileService > fs;
    mumu_tree = fs->make < TTree > ("diMuonTree", "Tree of diMuons");

    mumu_tree->Branch("run", &run, "run/i");
    mumu_tree->Branch("event", &event, "event/l");
    mumu_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

    mumu_tree->Branch("p4", "TLorentzVector", &p4);

    mumu_tree->Branch("muP_p4",  "TLorentzVector", &muP_p4);
    mumu_tree->Branch("muM_p4",  "TLorentzVector", &muM_p4);

    mumu_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");

    mumu_tree->Branch("dz", &dz, "dz/D");

    mumu_tree->Branch("vertex",  "Point", &vertex);

}

//Check recursively if any ancestor of particle is the given one
bool diMuonRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

// ------------ method called for each event  ------------
void diMuonRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  int debug = 0;

  edm::Handle < pat::CompositeCandidateCollection >jcand_hand;
  iEvent.getByToken(jcand_, jcand_hand);

  edm::Handle < reco::VertexCollection  >primaryVertices_handle;
  iEvent.getByToken(primaryVertices_, primaryVertices_handle);

  edm::Handle < edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken(triggerResults_, triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  edm::Handle<reco::GenParticleCollection> pruned;
  // iEvent.getByToken(genCands_,pruned);

  // if (false && isMC_) {
  //  gen_chi_p4.SetPtEtaPhiM(0, 0, 0, 0);
  //  gen_yns_p4.SetPtEtaPhiM(0, 0, 0, 0);
  //  gen_dimuon_p4.SetPtEtaPhiM(0, 0, 0, 0);
  //  chi_pdgId = 0;
  //  for (size_t i=0; i<pruned->size(); i++) {
  //     int p_id = abs((*pruned)[i].pdgId());
  //     int p_status = (*pruned)[i].status();
  //     yns_pdgId = 0;
  //     int foundit = 0;
  //     if ( ( p_id == 20443 || p_id == 445 || p_id == 10441) && p_status == 2)  yns_pdgId = 443;
  //     if (yns_pdgId > 0) {
  //        chi_pdgId = p_id;
  //        foundit++;
  //        const reco::Candidate * pwave = &(*pruned)[i];
  //        gen_chi_p4.SetPtEtaPhiM(pwave->pt(),pwave->eta(),pwave->phi(),pwave->mass());
  //        for (size_t j=0; j<pwave->numberOfDaughters(); j++) {
  //           const reco::Candidate *dau = pwave->daughter(j);
  //           if (dau->pdgId() == yns_pdgId && dau->status() == 2) {
  //              gen_yns_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
  //              uint nmuons = 0;
  //              for (size_t k=0; k<dau->numberOfDaughters(); k++) {
  //                 const reco::Candidate *gdau = dau->daughter(k);
  //                 if (gdau->pdgId() == 13 && gdau->status()==1) {
  //                    nmuons++;
  //                    gen_muonP_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
  //                 } else {
  //                    if (gdau->pdgId() == -13 && gdau->status()==1) {
  //                       nmuons++;
  //                       gen_muonM_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
  //                    }
  //                 }
  //              }
  //              if (nmuons == 2 ) {
  //                 foundit += 3;                                  // found complete dimuon decay
  //                 gen_dimuon_p4 = gen_muonP_p4 + gen_muonM_p4;   // will account fsr
  //              }
  //           } else {
  //              if (dau->pdgId() == 22 && dau->status() ==1) {
  //                 foundit++;
  //                 gen_photon_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
  //              }  else std::cout << "Rootupler: unexpected pdg_id " << dau->pdgId() << " (" << run << "," << event << ")" << std::endl;
  //           }
  //           if (foundit == 5 ) break;                             // decay found !
  //        }
  //     }
  //     if (chi_pdgId && yns_pdgId && foundit==5) break;        // just one decay of this kind is expected
  //     else chi_pdgId = 0;
  //  }
  //  if (!chi_pdgId)  std::cout << "Rootupler does not found the given decay " << run << "," << event << std::endl;
  // }

   //grab Trigger informations
   // save it in variable trigger, trigger is an int between 0 and 15, in binary it is:
   // (pass 11)(pass 8)(pass 7)(pass 5)
   // es. 11 = pass 5, 7 and 11
   // es. 4 = pass only 8

   trigger = 0;
   if (triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      unsigned int NTRIGGERS = HLTs_.size();

      for (unsigned int i = 0; i < NTRIGGERS; i++) {
         for (int version = 1; version < 19; version++) {
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

    bool bestCandidateOnly_ = false;

    j_rank = 0;
    // std::string getdata = "";
    if (jcand_hand.isValid() && !jcand_hand->empty()) {
      for (unsigned int i=0; i< jcand_hand->size(); i++) {
        pat::CompositeCandidate diMuon = jcand_hand->at(i);

        vertex  = diMuon.vertex();

        // PVwithmuons = (diMuon.userData<reco::Vertex>("PVwithmuons"))->Point();
        // muLessVertex = (diMuon.userData<reco::Vertex>("muonlessPV"));
        // commonVertex = (diMuon.userData<reco::Vertex>("commonVertex"));

        // filter = diMuon.userInt("isTriggerMatched");
        countTksOfPV = diMuon.userInt("countTksOfPV");
        vertexWeight = diMuon.userFloat("vertexWeight");
        sumPTPV       = diMuon.userFloat("sumPTPV");

        vProb           = diMuon.userFloat("vProb");
        vNChi2          = diMuon.userFloat("vNChi2");

        cosAlpha = diMuon.userFloat("cosAlpha");

        MassErr = diMuon.userFloat("MassErr");

        MassErr = diMuon.userFloat("MassErr");

        p4.SetPtEtaPhiM(diMuon.pt(), diMuon.eta(), diMuon.phi(), diMuon.mass());

        if ((diMuon.daughter("muon1")->charge()) > 0 )
        {
          muM_p4.SetPtEtaPhiM(diMuon.daughter("muon2")->pt(), diMuon.eta(), diMuon.daughter("muon2")->phi(), diMuon.daughter("muon2")->mass());
          muP_p4.SetPtEtaPhiM(diMuon.daughter("muon1")->pt(), diMuon.daughter("muon1")->eta(), diMuon.daughter("muon1")->phi(), diMuon.daughter("muon1")->mass());
        } else
        {
          muP_p4.SetPtEtaPhiM(diMuon.daughter("muon2")->pt(), diMuon.eta(), diMuon.daughter("muon2")->phi(), diMuon.daughter("muon2")->mass());
          muM_p4.SetPtEtaPhiM(diMuon.daughter("muon1")->pt(), diMuon.daughter("muon1")->eta(), diMuon.daughter("muon1")->phi(), diMuon.daughter("muon1")->mass());
        }

        mumu_tree->Fill();

        j_rank++;
      }
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void diMuonRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(diMuonRootupler);
