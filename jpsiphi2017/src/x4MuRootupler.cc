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
        edm::EDGetTokenT<pat::CompositeCandidateCollection> bkgcand_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> phi_dimuon_Label;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsi_dimuon_Label;
        edm::EDGetTokenT<reco::VertexCollection>            primaryVertices_;
        edm::EDGetTokenT<edm::TriggerResults>               triggerResults_;
        std::vector<std::string>                            HLTs_;

	bool isMC_;

  TTree *x_tree;
  TTree *j_tree;
  TTree *p_tree;
  TTree *b_tree;

	UInt_t run;
  ULong64_t event;
  UInt_t lumiblock;
  UInt_t trigger;
  UInt_t numPrimaryVertices;
  UInt_t countTksOfPV;

  //x tree variables
	TLorentzVector x_p4;
	TLorentzVector jpsi_p4;
	TLorentzVector muonM_jpsi_p4;
	TLorentzVector muonP_jpsi_p4;
  TLorentzVector phi_p4;
  TLorentzVector muonM_phi_p4;
  TLorentzVector muonP_phi_p4;

  Double_t xM,jpsi_M,phi_M;
  Double_t cosAlpha, cosAlpha3D, cosAlphaMuLess, ctauErrPV, ctauPV, ctauPVMuLess, ctauErrPVMuLess;
  Double_t cosAlphaBS,cosAlphaBS3D, ctauErrBS, ctauBS, vNChi2, vProb, sumPTPV;
  Double_t phi_deltaR, jpsi_deltaR;
  Double_t l_xy, l_xyBS, l_xyz, l_xyzBS;
  Double_t lErr_xy, lErr_xyBS, lErr_xyz, lErr_xyzBS;
  Double_t vertexWeight, dz, dz_jpsi, dz_phi;
  Double_t MassErr;

  UInt_t jpsi_i,phi_i;
  UInt_t x_rank, jpsi_muonM_type, jpsi_muonP_type, phi_muonP_type, phi_muonM_type;
  UInt_t jpsi_trigger, phi_trigger;

  //bkg tree variables
	TLorentzVector x_p4Bkg;
  Double_t xMBkg, jpsi_deltaRBkg, phi_deltaRBkg;

  //jpsi tree variables

  TLorentzVector j_muonM_p4, j_muonP_p4, j_p4;

  Double_t jM;
  Double_t j_cosAlpha,j_cosAlphaBS, j_cosAlphaBS3D, j_cosAlpha3D, j_vNChi2, j_vProb, j_dz;
  Double_t j_l_xy, j_l_xyBS, j_l_xyz, j_l_xyzBS, j_deltaR;
  Double_t j_lErr_xy, j_lErr_xyBS, j_lErr_xyz, j_lErr_xyzBS;
  Double_t j_ctauErrPV, j_ctauPV, j_ctauErrBS, j_ctauBS;
  UInt_t j_rank, j_muonM_type, j_muonP_type, j_triggerMatch;
  UInt_t j_muonM_isGlobal,j_muonM_isTracker, j_muonP_isGlobal,j_muonP_isTracker;

  //phi tree variables

  TLorentzVector p_muonM_p4, p_muonP_p4, p_p4;

  Double_t pM;
  Double_t p_cosAlpha, p_cosAlphaBS, p_cosAlphaBS3D, p_cosAlpha3D;
  Double_t p_ctauErrPV, p_ctauPV, p_ctauErrBS, p_ctauBS, p_vNChi2, p_vProb, p_dz, p_refit;
  Double_t p_l_xy, p_l_xyBS, p_l_xyz, p_l_xyzBS, p_deltaR;
  Double_t p_lErr_xy, p_lErr_xyBS, p_lErr_xyz, p_lErr_xyzBS;
  UInt_t p_rank, p_muonM_type, p_muonP_type, p_triggerMatch;
  UInt_t p_muonM_isGlobal,p_muonM_isTracker, p_muonP_isGlobal,p_muonP_isTracker ;

  //Verteces
  Point jVertex, pVertex;
  Point xVertex, jpsiVertex, phiVertex;
  Point xVertexBkg, jpsiVertexBkg, phiVertexBkg;

  reco::Vertex commonVertex;
  Point PVwithmuons;
  reco::Vertex muLessVertex;

  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;

  TTree *upsilon_tree;
  TLorentzVector mumu_p4, muP_p4, muM_p4;


};

x4MuRootupler::x4MuRootupler(const edm::ParameterSet & iConfig):
// chi_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("chi_cand"))),
xcand_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("x_cand"))),
bkgcand_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("bkg_cand"))),
phi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("phidimuons"))),
jpsi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("jpsidimuons"))),
// refit1_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refit1S"))),
primaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("primaryVertices"))),
triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
isMC_(iConfig.getParameter < bool > ("isMC"))
{

    edm::Service < TFileService > fs;
    x_tree = fs->make < TTree > ("xTree", "Tree of xs");
    j_tree = fs->make < TTree > ("jTree", "Tree of jpsis");
    p_tree = fs->make < TTree > ("pTree", "Tree of phis");
    b_tree = fs->make < TTree > ("bTree", "Tree of bkg");

    //bkg tree
    b_tree->Branch("run", &run, "run/i");
    b_tree->Branch("event", &event, "event/l");
    x_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

    b_tree->Branch("x_p4Bkg", "TLorentzVector", &x_p4Bkg);
    b_tree->Branch("xMBkg",  &xMBkg, "xMBkg/D");

    b_tree->Branch("jpsi_deltaRBkg", &jpsi_deltaRBkg, "jpsi_deltaRBkg/D");
    b_tree->Branch("phi_deltaRBkg", &phi_deltaRBkg, "phi_deltaRBkg/D");

    //signal tree

    x_tree->Branch("run", &run, "run/i");
    x_tree->Branch("event", &event, "event/l");
    x_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

    x_tree->Branch("x_p4", "TLorentzVector", &x_p4);
    x_tree->Branch("xM",  &xM, "xM/D");
    x_tree->Branch("trigger", &trigger, "trigger/i");
    // x_tree->Branch("filter", &filter, "filter/i");

    x_tree->Branch("jpsi_i",&jpsi_i, "jpsi_i/i");
    x_tree->Branch("phi_i", &phi_i, "phi_i/i");

    x_tree->Branch("jpsi_p4", "TLorentzVector", &jpsi_p4);
    x_tree->Branch("jpsi_M", &jpsi_M, "jpsi_M/D");
    x_tree->Branch("muonM_jpsi_p4",  "TLorentzVector", &muonM_jpsi_p4);
    x_tree->Branch("muonP_jpsi_p4",  "TLorentzVector", &muonP_jpsi_p4);
    x_tree->Branch("jpsi_muonM_type", &jpsi_muonM_type, "jpsi_muonM_type/I");
    x_tree->Branch("jpsi_muonP_type", &jpsi_muonP_type, "jpsi_muonP_type/I");

    x_tree->Branch("phi_p4", "TLorentzVector", &phi_p4);
    x_tree->Branch("phi_M", &phi_M, "phi_M/D");
    x_tree->Branch("muonM_phi_p4",  "TLorentzVector", &muonM_phi_p4);
    x_tree->Branch("muonP_phi_p4",  "TLorentzVector", &muonP_phi_p4);
    x_tree->Branch("phi_muonM_type", &phi_muonM_type, "phi_muonM_type/I");
    x_tree->Branch("phi_muonP_type", &phi_muonP_type, "phi_muonP_type/I");

    x_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");

    x_tree->Branch("dz", &dz, "dz/D");
    x_tree->Branch("dzjpsi", &dz_jpsi, "dz_jpsi/D");
    x_tree->Branch("dzphi", &dz_phi, "dz_phi/D");

    x_tree->Branch("xVertex",  "Point", &xVertex);
    x_tree->Branch("muLessVertex",  "reco::Vertex", &muLessVertex);
    x_tree->Branch("PVwithmuons",  "Point", &PVwithmuons);
    x_tree->Branch("jpsiVertex",  "Point", &jpsiVertex);
    x_tree->Branch("phiVertex",  "Point", &phiVertex);
    x_tree->Branch("commonVertex",  "reco::Vertex", &commonVertex);

    x_tree->Branch("jpsi_trigger", &jpsi_trigger, "jpsi_trigger/i");
    x_tree->Branch("phi_trigger", &phi_trigger, "phi_trigger/i");

    x_tree->Branch("jpsi_deltaR", &jpsi_deltaR, "jpsi_deltaR/i");
    x_tree->Branch("phi_deltaR", &phi_deltaR, "phi_deltaR/i");

    x_tree->Branch("countTksOfPV", &countTksOfPV, "countTksOfPV/i");
    x_tree->Branch("vertexWeight", &vertexWeight, "vertexWeight/D");
    x_tree->Branch("sumPTPV", &sumPTPV, "sumPTPV/D");

    x_tree->Branch("vProb", &vProb, "vProb/D");
    x_tree->Branch("vNChi2", &vNChi2, "vNChi2/D");

    x_tree->Branch("ctauBS", &ctauBS, "ctauBS/D");
    x_tree->Branch("ctauErrBS", &ctauErrBS, "ctauErrBS/D");

    x_tree->Branch("ctauPV", &ctauPV, "ctauPV/D");
    x_tree->Branch("ctauErrPV", &ctauErrPV, "ctauErrPV/D");

    x_tree->Branch("ctauPVMuLess", &ctauPVMuLess, "ctauPVMuLess/D");
    x_tree->Branch("ctauErrPVMuLess", &ctauErrPVMuLess, "ctauErrPVMuLess/D");

    x_tree->Branch("cosAlpha", &cosAlpha, "cosAlpha/D");
    x_tree->Branch("cosAlphaMuLess", &cosAlphaMuLess, "cosAlphaMuLess/D");

    x_tree->Branch("cosAlphaBS", &cosAlphaBS, "cosAlphaBS/D");
    x_tree->Branch("cosAlpha3D", &cosAlpha3D, "cosAlpha3D/D");
    x_tree->Branch("cosAlphaBS3D", &cosAlphaBS3D, "cosAlphaBS3D/D");

    x_tree->Branch("l_xy", &l_xy, "l_xy/D");
    x_tree->Branch("l_xyBS", &l_xyBS, "l_xyBS/D");
    x_tree->Branch("l_xyz", &l_xyz, "l_xyz/D");
    x_tree->Branch("l_xyzBS", &l_xyzBS, "l_xyzBS/D");

    x_tree->Branch("lErr_xy", &lErr_xy, "lErr_xy/D");
    x_tree->Branch("lErr_xyBS", &lErr_xyBS, "lErr_xyBS/D");
    x_tree->Branch("lErr_xyz", &lErr_xyz, "lErr_xyz/D");
    x_tree->Branch("lErr_xyzBS", &lErr_xyzBS, "lErr_xyzBS/D");

    x_tree->Branch("x_rank", &x_rank, "x_rank/I");

    //jpsi tree
    j_tree->Branch("run", &run, "run/i");
    j_tree->Branch("event", &event, "event/l");
    j_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

    j_tree->Branch("j_p4", "TLorentzVector", &x_p4);
    j_tree->Branch("jM",  &jM, "jM/D");
    j_tree->Branch("trigger", &trigger, "trigger/i");
    j_tree->Branch("jVertex",  "Point", &xVertex);

    j_tree->Branch("j_muonM_p4",  "TLorentzVector", &j_muonM_p4);
    j_tree->Branch("j_muonP_p4",  "TLorentzVector", &j_muonP_p4);
    j_tree->Branch("j_muonM_type", &j_muonM_type, "j_muonM_type/I");
    j_tree->Branch("j_muonP_type", &j_muonP_type, "j_muonP_type/I");

    j_tree->Branch("j_muonM_isGlobal", &j_muonM_isGlobal, "j_muonM_isGlobal/I");
    j_tree->Branch("j_muonP_isGlobal", &j_muonP_isGlobal, "j_muonP_isGlobal/I");
    j_tree->Branch("j_muonM_isTracker", &j_muonM_isTracker, "j_muonM_isTracker/I");
    j_tree->Branch("j_muonP_isTracker", &j_muonP_isTracker, "j_muonP_isTracker/I");

    j_tree->Branch("j_vProb", &j_vProb, "j_vProb/D");
    j_tree->Branch("j_triggerMatch", &j_triggerMatch, "j_triggerMatch/I");
    // j_tree->Branch("j_dz", &j_dz, "j_dz/D");
    j_tree->Branch("j_vNChi2", &j_vNChi2, "j_vNChi2/D");
    j_tree->Branch("j_refit",&j_refit,"j_refit/D");

    j_tree->Branch("j_cosAlpha", &j_cosAlpha, "j_cosAlpha/D");
    j_tree->Branch("j_cosAlphaBS", &j_cosAlphaBS, "j_cosAlphaBS/D");
    j_tree->Branch("j_cosAlpha3D", &j_cosAlpha3D, "j_cosAlpha3D/D");
    j_tree->Branch("j_cosAlphaBS3D", &j_cosAlphaBS3D, "j_cosAlphaBS3D/D");

    j_tree->Branch("j_deltaR", &j_deltaR, "j_deltaR/D");

    j_tree->Branch("j_l_xy", &j_l_xy, "j_l_xy/D");
    j_tree->Branch("j_l_xyBS", &j_l_xyBS, "j_l_xyBS/D");
    j_tree->Branch("j_l_xyz", &j_l_xyz, "j_l_xyz/D");
    j_tree->Branch("j_l_xyzBS", &j_l_xyzBS, "j_l_xyzBS/D");

    j_tree->Branch("j_lErr_xy", &j_lErr_xy, "j_lErr_xy/D");
    j_tree->Branch("j_lErr_xyBS", &j_lErr_xyBS, "j_lErr_xyBS/D");
    j_tree->Branch("j_lErr_xyz", &j_lErr_xyz, "j_lErr_xyz/D");
    j_tree->Branch("j_lErr_xyzBS", &j_lErr_xyzBS, "j_lErr_xyzBS/D");

    j_tree->Branch("j_ctauPV", &j_ctauPV, "j_ctauPV/D");
    j_tree->Branch("j_ctauErrPV", &j_ctauErrPV, "j_ctauErrPV/D");
    j_tree->Branch("j_ctauBS", &j_ctauBS, "j_ctauBS/D");
    j_tree->Branch("j_ctauErrBS", &j_ctauErrBS, "j_ctauErrBS/D");

    j_tree->Branch("j_rank", &j_rank, "j_rank/I");

    //phi tree
    p_tree->Branch("run", &run, "run/i");
    p_tree->Branch("event", &event, "event/l");
    p_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

    p_tree->Branch("p_p4", "TLorentzVector", &x_p4);
    p_tree->Branch("pM",  &pM, "pM/D");
    p_tree->Branch("trigger", &trigger, "trigger/i");
    p_tree->Branch("jVertex",  "Point", &xVertex);

    p_tree->Branch("p_muonM_p4",  "TLorentzVector", &p_muonM_p4);
    p_tree->Branch("p_muonP_p4",  "TLorentzVector", &p_muonP_p4);
    p_tree->Branch("p_muonM_type", &p_muonM_type, "p_muonM_type/I");
    p_tree->Branch("p_muonP_type", &p_muonP_type, "p_muonP_type/I");
    p_tree->Branch("p_muonM_isGlobal", &p_muonM_isGlobal, "p_muonM_isGlobal/I");
    p_tree->Branch("p_muonP_isGlobal", &p_muonP_isGlobal, "p_muonP_isGlobal/I");
    p_tree->Branch("p_muonM_isTracker", &p_muonM_isTracker, "p_muonM_isTracker/I");
    p_tree->Branch("p_muonP_isTracker", &p_muonP_isTracker, "p_muonP_isTracker/I");

    p_tree->Branch("p_refit",&p_refit,"p_refit/D");
    p_tree->Branch("p_vProb", &p_vProb, "p_vProb/D");
    p_tree->Branch("p_triggerMatch", &p_triggerMatch, "p_triggerMatch/I");
    // p_tree->Branch("p_dz", &p_dz, "p_dz/D");
    p_tree->Branch("p_vNChi2", &p_vNChi2, "p_vNChi2/D");

    p_tree->Branch("p_cosAlpha", &p_cosAlpha, "p_cosAlpha/D");
    p_tree->Branch("p_cosAlphaBS", &p_cosAlphaBS, "p_cosAlphaBS/D");
    p_tree->Branch("p_cosAlpha3D", &p_cosAlpha3D, "p_cosAlpha3D/D");
    p_tree->Branch("p_cosAlphaBS3D", &p_cosAlphaBS3D, "p_cosAlphaBS3D/D");

    p_tree->Branch("p_deltaR", &p_deltaR, "p_deltaR/D");

    p_tree->Branch("p_l_xy", &p_l_xy, "p_l_xy/D");
    p_tree->Branch("p_l_xyBS", &p_l_xyBS, "p_l_xyBS/D");
    p_tree->Branch("p_l_xyz", &p_l_xyz, "p_l_xyz/D");
    p_tree->Branch("p_l_xyzBS", &p_l_xyzBS, "p_l_xyzBS/D");

    p_tree->Branch("p_lErr_xy", &p_lErr_xy, "p_lErr_xy/D");
    p_tree->Branch("p_lErr_xyBS", &p_lErr_xyBS, "p_lErr_xyBS/D");
    p_tree->Branch("p_lErr_xyz", &p_lErr_xyz, "p_lErr_xyz/D");
    p_tree->Branch("p_lErr_xyzBS", &p_lErr_xyzBS, "p_lErr_xyzBS/D");

    p_tree->Branch("p_ctauPV", &p_ctauPV, "p_ctauPV/D");
    p_tree->Branch("p_ctauErrPV", &p_ctauErrPV, "p_ctauErrPV/D");
    p_tree->Branch("p_ctauBS", &p_ctauBS, "p_ctauBS/D");
    p_tree->Branch("p_ctauErrBS", &p_ctauErrBS, "p_ctauErrBS/D");

    p_tree->Branch("p_rank", &p_rank, "p_rank/I");


}

//Check recursively if any ancestor of particle is the given one
bool x4MuRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

// ------------ method called for each event  ------------
void x4MuRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  int debug = 0;

  edm::Handle < pat::CompositeCandidateCollection >xcand_hand;
  iEvent.getByToken(xcand_, xcand_hand);

  edm::Handle < pat::CompositeCandidateCollection >bkg_hand;
  iEvent.getByToken(bkgcand_, bkg_hand);

  edm::Handle<pat::CompositeCandidateCollection> dimuonsPhi;
  iEvent.getByToken(phi_dimuon_Label,dimuonsPhi);

  edm::Handle<pat::CompositeCandidateCollection> dimuonsJPsi;
  iEvent.getByToken(jpsi_dimuon_Label,dimuonsJPsi);

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

    if (bkg_hand.isValid() && !bkg_hand->empty())
    {

      for (unsigned int i=0; i< bkg_hand->size(); i++)
      {

        pat::CompositeCandidate b_ = bkg_hand->at(i);

        x_p4Bkg.SetPtEtaPhiM(b_.pt(), b_.eta(), b_.phi(), b_.mass());
        xMBkg = x_p4Bkg.M();

        jpsi_deltaRBkg = b_.userFloat("jpsi_deltaR");
        phi_deltaRBkg  = b_.userFloat("phi_deltaR");

        b_tree->Fill();

      }

    }

    // std::cout << "JPsi " << dimuonsJPsi->size() << std::endl;

    if (dimuonsJPsi.isValid() && !dimuonsJPsi->empty())
    {
      j_rank = dimuonsJPsi->size();

      for (unsigned int i=0; i< dimuonsJPsi->size(); i++)
      {
        // std::cout << i << std::endl;

        pat::CompositeCandidate j_ = dimuonsJPsi->at(i);

        j_p4.SetPtEtaPhiM(j_.pt(), j_.eta(), j_.phi(), j_.mass());
        jM = j_p4.M();
        jVertex = j_.vertex();

        if ((j_.daughter("muon1")->charge()) > 0 )
        {
          j_muonM_p4.SetPtEtaPhiM(j_.daughter("muon2")->pt(), j_.eta(), j_.daughter("muon2")->phi(), j_.daughter("muon2")->mass());
          j_muonP_p4.SetPtEtaPhiM(j_.daughter("muon1")->pt(), j_.daughter("muon1")->eta(), j_.daughter("muon1")->phi(), j_.daughter("muon1")->mass());

          j_muonP_type = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"))->type();
          j_muonM_type = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"))->type();

          j_muonP_isTracker = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"))->isTrackerMuon();
          j_muonM_isTracker = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"))->isTrackerMuon();

          j_muonP_isGlobal = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"))->isGlobalMuon();
          j_muonM_isGlobal = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"))->isGlobalMuon();

        } else
        {
          j_muonP_p4.SetPtEtaPhiM(j_.daughter("muon2")->pt(), j_.eta(), j_.daughter("muon2")->phi(), j_.daughter("muon2")->mass());
          j_muonM_p4.SetPtEtaPhiM(j_.daughter("muon1")->pt(), j_.daughter("muon1")->eta(), j_.daughter("muon1")->phi(), j_.daughter("muon1")->mass());

          j_muonP_type = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"))->type();
          j_muonM_type = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"))->type();

          j_muonP_isTracker = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"))->isTrackerMuon();
          j_muonM_isTracker = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"))->isTrackerMuon();

          j_muonP_isGlobal = dynamic_cast<const pat::Muon*>(j_.daughter("muon2"))->isGlobalMuon();
          j_muonM_isGlobal = dynamic_cast<const pat::Muon*>(j_.daughter("muon1"))->isGlobalMuon();

        }

        j_vProb           = j_.userFloat("vProb");
        j_vNChi2          = j_.userFloat("vNChi2");

        j_refit           = j_.userFloat("refittedMass");

        j_ctauBS          = j_.userFloat("ctauBS");
        j_ctauErrBS       = j_.userFloat("ctauErrBS");

        j_ctauPV          = j_.userFloat("ctauPV");
        j_ctauErrPV       = j_.userFloat("ctauErrPV");

        j_cosAlpha = j_.userFloat("cosAlpha");
        j_cosAlphaBS = j_.userFloat("cosAlphaBS");
        j_cosAlpha3D = j_.userFloat("cosAlpha3D");
        j_cosAlphaBS3D = j_.userFloat("cosAlphaBS3D");

        j_l_xy = j_.userFloat("l_xy");
        j_l_xyBS = j_.userFloat("l_xyBS");
        j_l_xyz = j_.userFloat("l_xyz");
        j_l_xyzBS = j_.userFloat("l_xyzBS");

        j_lErr_xy = j_.userFloat("lErr_xy");
        j_lErr_xyBS = j_.userFloat("lErr_xyBS");
        j_lErr_xyz = j_.userFloat("lErr_xyz");
        j_lErr_xyzBS = j_.userFloat("lErr_xyzBS");

        j_triggerMatch = j_.userInt("isTriggerMatched");

        j_deltaR = j_.userFloat("deltaR");

        j_tree->Fill();

      }
    }

    // std::cout << "Phi " <<  dimuonsPhi->size() << std::endl;

    if (dimuonsPhi.isValid() && !dimuonsPhi->empty())
    {
      p_rank = dimuonsPhi->size();

      for (unsigned int i=0; i< dimuonsPhi->size(); i++)
      {

        // std::cout << i << std::endl;

        pat::CompositeCandidate p_ = dimuonsPhi->at(i);

        p_p4.SetPtEtaPhiM(p_.pt(), p_.eta(), p_.phi(), p_.mass());
        pM = p_p4.M();
        pVertex = p_.vertex();

        if ((p_.daughter("muon1")->charge()) > 0 )
        {
          p_muonP_p4.SetPtEtaPhiM(p_.daughter("muon1")->pt(), p_.daughter("muon1")->eta(), p_.daughter("muon1")->phi(), p_.daughter("muon1")->mass());
          p_muonM_p4.SetPtEtaPhiM(p_.daughter("muon2")->pt(), p_.daughter("muon2")->eta(), p_.daughter("muon2")->phi(), p_.daughter("muon2")->mass());

          p_muonP_type = dynamic_cast<const pat::Muon*>(p_.daughter("muon1"))->type();
          p_muonM_type = dynamic_cast<const pat::Muon*>(p_.daughter("muon2"))->type();

          p_muonP_isTracker = dynamic_cast<const pat::Muon*>(p_.daughter("muon1"))->isTrackerMuon();
          p_muonM_isTracker = dynamic_cast<const pat::Muon*>(p_.daughter("muon2"))->isTrackerMuon();

          p_muonP_isGlobal = dynamic_cast<const pat::Muon*>(p_.daughter("muon1"))->isGlobalMuon();
          p_muonM_isGlobal = dynamic_cast<const pat::Muon*>(p_.daughter("muon2"))->isGlobalMuon();


        } else
        {
          p_muonP_p4.SetPtEtaPhiM(p_.daughter("muon2")->pt(), p_.daughter("muon2")->eta(), p_.daughter("muon2")->phi(), p_.daughter("muon2")->mass());
          p_muonM_p4.SetPtEtaPhiM(p_.daughter("muon1")->pt(), p_.daughter("muon1")->eta(), p_.daughter("muon1")->phi(), p_.daughter("muon1")->mass());

          p_muonP_type = dynamic_cast<const pat::Muon*>(p_.daughter("muon2"))->type();
          p_muonM_type = dynamic_cast<const pat::Muon*>(p_.daughter("muon1"))->type();

          p_muonP_isTracker = dynamic_cast<const pat::Muon*>(p_.daughter("muon2"))->isTrackerMuon();
          p_muonM_isTracker = dynamic_cast<const pat::Muon*>(p_.daughter("muon1"))->isTrackerMuon();

          p_muonP_isGlobal = dynamic_cast<const pat::Muon*>(p_.daughter("muon2"))->isGlobalMuon();
          p_muonM_isGlobal = dynamic_cast<const pat::Muon*>(p_.daughter("muon1"))->isGlobalMuon();

        }

        p_vNChi2          = p_.userFloat("vNChi2");
        p_vProb           = p_.userFloat("vProb");

        p_refit           = p_.userFloat("refittedMass");
        p_ctauBS          = p_.userFloat("ctauBS");
        p_ctauErrBS       = p_.userFloat("ctauErrBS");

        p_ctauPV          = p_.userFloat("ctauPV");
        p_ctauErrPV       = p_.userFloat("ctauErrPV");

        p_cosAlpha = p_.userFloat("cosAlpha");
        p_cosAlphaBS = p_.userFloat("cosAlphaBS");
        p_cosAlpha3D = p_.userFloat("cosAlpha3D");
        p_cosAlphaBS3D = p_.userFloat("cosAlphaBS3D");

        p_l_xy = p_.userFloat("l_xy");
        p_l_xyBS = p_.userFloat("l_xyBS");
        p_l_xyz = p_.userFloat("l_xyz");
        p_l_xyzBS = p_.userFloat("l_xyzBS");

        p_lErr_xy = p_.userFloat("lErr_xy");
        p_lErr_xyBS = p_.userFloat("lErr_xyBS");
        p_lErr_xyz = p_.userFloat("lErr_xyz");
        p_lErr_xyzBS = p_.userFloat("lErr_xyzBS");

        p_triggerMatch = p_.userInt("isTriggerMatched");

        p_deltaR = p_.userFloat("deltaR");

        p_tree->Fill();

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
void x4MuRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(x4MuRootupler);
