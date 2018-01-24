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

#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef math::XYZPoint Point;

class xMuMuKKRootupler:public edm::EDAnalyzer {
      public:
	explicit xMuMuKKRootupler(const edm::ParameterSet &);
	~xMuMuKKRootupler() override {};
	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

      private:
	void analyze(const edm::Event &, const edm::EventSetup &) override;

	      std::string file_name;
        // edm::EDGetTokenT<pat::CompositeCandidateCollection> xcand_;
        // edm::EDGetTokenT<pat::CompositeCandidateCollection> bkgcand_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> phi_dimuon_Label;
        // edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsi_dimuon_Label;
        edm::EDGetTokenT<reco::VertexCollection>            primaryVertices_;
        edm::EDGetTokenT<edm::TriggerResults>               triggerResults_;
        std::vector<std::string>                            HLTs_;

	bool isMC_;

  TTree *x_tree;

  Double_t kkM;


};

xMuMuKKRootupler::xMuMuKKRootupler(const edm::ParameterSet & iConfig):
// chi_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("chi_cand"))),
//xcand_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("x_cand"))),
//bkgcand_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("bkg_cand"))),
phi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("phidimuons"))),
// jpsi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("jpsidimuons"))),
// refit1_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refit1S"))),
primaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("primaryVertices"))),
triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
isMC_(iConfig.getParameter < bool > ("isMC"))
{

    edm::Service < TFileService > fs;
    x_tree = fs->make < TTree > ("xTree", "Tree of xs");

    //bkg tree
    x_tree->Branch("kkM", &kkM, "kkM/i");


}

//Check recursively if any ancestor of particle is the given one
bool xMuMuKKRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

// ------------ method called for each event  ------------
void xMuMuKKRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  int debug = 0;

  edm::Handle<pat::CompositeCandidateCollection> dimuonsPhi;
  iEvent.getByToken(phi_dimuon_Label,dimuonsPhi);

  edm::Handle < reco::VertexCollection  >primaryVertices_handle;
  iEvent.getByToken(primaryVertices_, primaryVertices_handle);

  edm::Handle < edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken(triggerResults_, triggerResults_handle);


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

    // bool bestCandidateOnly_ = false;

    // std::cout << "JPsi " << dimuonsJPsi->size() << std::endl;

    if (dimuonsPhi.isValid() && !dimuonsPhi->empty())
    {


      for (unsigned int i=0; i< dimuonsPhi->size(); i++)
      {

        // std::cout << i << std::endl;

        pat::CompositeCandidate p_ = dimuonsPhi->at(i);

        TLorentzVector kkP4;
        kkP4.SetPtEtaPhiM(p_.pt(), p_.eta(), p_.phi(), p_.mass());
        kkM = kkP4.M();


        x_tree->Fill();

      }
    }

    // std::cout << "X " <<  xcand_hand->size() << std::endl;

    x_rank = 0;
    // std::string getdata = "";
    if (xcand_hand.isValid() && !xcand_hand->empty()) {
      for (unsigned int i=0; i< xcand_hand->size(); i++) {

        // std::cout << i << std::endl;

        pat::CompositeCandidate x_ = xcand_hand->at(i);

        xVertex  = x_.vertex();
        phiVertex = x_.daughter("phi")->vertex();
        jpsiVertex = x_.daughter("jpsi")->vertex();

        jpsi_trigger = x_.userInt("jpsi_isTriggerMatched");
        phi_trigger = x_.userInt("phi_isTriggerMatched");

        jpsi_deltaR = x_.userFloat("jpsi_deltaR");
        phi_deltaR = x_.userFloat("phi_deltaR");
        // PVwithmuons = (x_.userData<reco::Vertex>("PVwithmuons"))->Point();
        // muLessVertex = (x_.userData<reco::Vertex>("muonlessPV"));
        // commonVertex = (x_.userData<reco::Vertex>("commonVertex"));

        countTksOfPV = x_.userInt("countTksOfPV");
        vertexWeight = x_.userFloat("vertexWeight");
        sumPTPV       = x_.userFloat("sumPTPV");

        vProb           = x_.userFloat("vProb");
        vNChi2          = x_.userFloat("vNChi2");

        ctauBS          = x_.userFloat("ctauBS");
        ctauErrBS       = x_.userFloat("ctauErrBS");

        ctauPV          = x_.userFloat("ctauPV");
        ctauErrPV       = x_.userFloat("ctauErrPV");

        ctauPVMuLess    = x_.userFloat("ctauPVMuLess");
        ctauErrPVMuLess = x_.userFloat("ctauErrPVMuLess");

        cosAlpha = x_.userFloat("cosAlpha");
        cosAlphaMuLess = x_.userFloat("cosAlphaMuLess");
        cosAlphaBS = x_.userFloat("cosAlphaBS");
        cosAlpha3D = x_.userFloat("cosAlpha3D");
        cosAlphaBS3D = x_.userFloat("cosAlphaBS3D");

        l_xy = x_.userFloat("l_xy");
        l_xyBS = x_.userFloat("l_xyBS");
        l_xyz = x_.userFloat("l_xyz");
        l_xyzBS = x_.userFloat("l_xyzBS");

        lErr_xy = x_.userFloat("lErr_xy");
        lErr_xyBS = x_.userFloat("lErr_xyBS");
        lErr_xyz = x_.userFloat("lErr_xyz");
        lErr_xyzBS = x_.userFloat("lErr_xyzBS");

        MassErr = x_.userFloat("MassErr");

        dz = x_.userFloat("dzFourMuons");
        dz_jpsi = x_.userFloat("dzJpsi");
        dz_phi = x_.userFloat("dzPhi");


        x_p4.SetPtEtaPhiM(x_.pt(), x_.eta(), x_.phi(), x_.mass());

        xM = x_p4.M();

        jpsi_p4.SetPtEtaPhiM(x_.daughter("jpsi")->pt(), x_.daughter("jpsi")->eta(), x_.daughter("jpsi")->phi(), x_.daughter("jpsi")->mass());
        if ((x_.daughter("jpsi")->daughter("muon1")->charge()) > 0 )
        {
          muonP_jpsi_p4.SetPtEtaPhiM(x_.daughter("jpsi")->daughter("muon1")->pt(), x_.daughter("jpsi")->daughter("muon1")->eta(), x_.daughter("jpsi")->daughter("muon1")->phi(), x_.daughter("jpsi")->daughter("muon1")->mass());
          muonM_jpsi_p4.SetPtEtaPhiM(x_.daughter("jpsi")->daughter("muon2")->pt(), x_.daughter("jpsi")->daughter("muon2")->eta(), x_.daughter("jpsi")->daughter("muon2")->phi(), x_.daughter("jpsi")->daughter("muon2")->mass());

          jpsi_muonP_type = dynamic_cast<const pat::Muon*>(x_.daughter("jpsi")->daughter("muon1"))->type();
          jpsi_muonM_type = dynamic_cast<const pat::Muon*>(x_.daughter("jpsi")->daughter("muon2"))->type();

        } else
        {
          muonP_jpsi_p4.SetPtEtaPhiM(x_.daughter("jpsi")->daughter("muon2")->pt(), x_.daughter("jpsi")->daughter("muon2")->eta(), x_.daughter("jpsi")->daughter("muon2")->phi(), x_.daughter("jpsi")->daughter("muon2")->mass());
          muonM_jpsi_p4.SetPtEtaPhiM(x_.daughter("jpsi")->daughter("muon1")->pt(), x_.daughter("jpsi")->daughter("muon1")->eta(), x_.daughter("jpsi")->daughter("muon1")->phi(), x_.daughter("jpsi")->daughter("muon1")->mass());

          jpsi_muonP_type = dynamic_cast<const pat::Muon*>(x_.daughter("jpsi")->daughter("muon2"))->type();
          jpsi_muonM_type = dynamic_cast<const pat::Muon*>(x_.daughter("jpsi")->daughter("muon1"))->type();

        }

        phi_p4.SetPtEtaPhiM(x_.daughter("phi")->pt(), x_.daughter("phi")->eta(), x_.daughter("phi")->phi(), x_.daughter("phi")->mass());
        if((x_.daughter("phi")->daughter("muon1")->charge()) > 0 )
        {
          muonM_phi_p4.SetPtEtaPhiM(x_.daughter("phi")->daughter("muon2")->pt(), x_.daughter("phi")->daughter("muon2")->eta(), x_.daughter("phi")->daughter("muon2")->phi(), x_.daughter("phi")->daughter("muon2")->mass());
          muonP_phi_p4.SetPtEtaPhiM(x_.daughter("phi")->daughter("muon1")->pt(), x_.daughter("phi")->daughter("muon1")->eta(), x_.daughter("phi")->daughter("muon1")->phi(), x_.daughter("phi")->daughter("muon1")->mass());

          phi_muonP_type = dynamic_cast<const pat::Muon*>(x_.daughter("phi")->daughter("muon1"))->type();
          phi_muonM_type = dynamic_cast<const pat::Muon*>(x_.daughter("phi")->daughter("muon2"))->type();

        }
        else
        {
          muonP_phi_p4.SetPtEtaPhiM(x_.daughter("phi")->daughter("muon2")->pt(), x_.daughter("phi")->daughter("muon2")->eta(), x_.daughter("phi")->daughter("muon2")->phi(), x_.daughter("phi")->daughter("muon2")->mass());
          muonM_phi_p4.SetPtEtaPhiM(x_.daughter("phi")->daughter("muon1")->pt(), x_.daughter("phi")->daughter("muon1")->eta(), x_.daughter("phi")->daughter("muon1")->phi(), x_.daughter("phi")->daughter("muon1")->mass());

          phi_muonP_type = dynamic_cast<const pat::Muon*>(x_.daughter("phi")->daughter("muon2"))->type();
          phi_muonM_type = dynamic_cast<const pat::Muon*>(x_.daughter("phi")->daughter("muon1"))->type();

        }

        phi_M = phi_p4.M();
        jpsi_M = jpsi_p4.M();

        x_tree->Fill();

        x_rank++;
      }
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void xMuMuKKRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(xMuMuKKRootupler);
