// -*- C++ -*-
//
// Package:    DiMuonDiTrackKF
// Class:      DiMuonDiTrackKF
//
/**\class DiMuonDiTrackKF Ponia/OniaTrak/src/DiMuonDiTrackKF.cc

 Description: performs vertex kinematical fit for (MuMu) (TrackTrack)

**/
//
// Original Author:  Alberto Sanchez-Hernandez
//         Created:  Octubre 2016

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//

class DiMuonDiTrackKF : public edm::EDProducer {
   public:
      explicit DiMuonDiTrackKF(const edm::ParameterSet&);
      ~DiMuonDiTrackKF() override;
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void produce(edm::Event&, const edm::EventSetup&) override;

// ----------member data ---------------------------
  edm::EDGetTokenT<pat::CompositeCandidateCollection> combination_cand_tag_;
  double dimuon_mass_constraint_;
  double ditrack_mass_constraint_;
  std::vector<double> combination_mass_cut_;
  std::vector<double> track_mass_;
  bool make_ditrack_vertex_;
  bool make_ditrack_mass_constraint_;
  std::string product_name_;

  template<typename T>
  struct GreaterByVProb {
     typedef T first_argument_type;
     typedef T second_argument_type;
     bool operator()( const T & t1, const T & t2 ) const {
        return t1.userFloat("vProb") > t2.userFloat("vProb");
     }
  };
};

//
// constructors and destructor
//
DiMuonDiTrackKF::DiMuonDiTrackKF(const edm::ParameterSet& iConfig) {
  combination_cand_tag_         = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("combination_cand_tag"));
  dimuon_mass_constraint_       = iConfig.getParameter<double>("dimuon_mass_constraint");
  ditrack_mass_constraint_      = iConfig.getParameter<double>("ditrack_mass_constraint");
  combination_mass_cut_         = iConfig.getParameter<std::vector<double>>("combination_mass_cut");
  track_mass_                   = iConfig.getParameter<std::vector<double>>("track_mass");
  make_ditrack_vertex_          = iConfig.getParameter<bool>("make_ditrack_vertex");
  make_ditrack_mass_constraint_ = iConfig.getParameter<bool>("make_ditrack_mass_constraint");
  product_name_                 = iConfig.getParameter<std::string>("product_name");

// kinematic refit collections
  produces<pat::CompositeCandidateCollection>(product_name_);
  if (make_ditrack_mass_constraint_) make_ditrack_vertex_ = make_ditrack_mass_constraint_;
}

DiMuonDiTrackKF::~DiMuonDiTrackKF() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void DiMuonDiTrackKF::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> combination_cand_handle;
  iEvent.getByToken(combination_cand_tag_, combination_cand_handle);

// Kinematic refit collection
  std::unique_ptr< pat::CompositeCandidateCollection > CombinationCollection(new pat::CompositeCandidateCollection);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> tt_builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",tt_builder);

  int index_combination=-1;
  if ( combination_cand_handle.isValid() && !combination_cand_handle->empty() ) {
    for (pat::CompositeCandidateCollection::const_iterator i_combination=combination_cand_handle->begin(); i_combination!=combination_cand_handle->end(); ++i_combination) {
    index_combination++;

    reco::TrackRef tk_combination[4]={
      ( dynamic_cast<const pat::Muon*>(i_combination->daughter("onia")->daughter("muon1")) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(i_combination->daughter("onia")->daughter("muon2")) )->innerTrack(),
      ( dynamic_cast<const pat::GenericParticle*>(i_combination->daughter("ditrak")->daughter("trak1")) )->track(),
      ( dynamic_cast<const pat::GenericParticle*>(i_combination->daughter("ditrak")->daughter("trak2")) )->track()
    };

    //const pat::CompositeCandidate *dimuonC = dynamic_cast<const pat::CompositeCandidate *>(i_combination->daughter("onia"));
    //const reco::Vertex the_pv = *dimuonC->userData<reco::Vertex>("PVwithmuons");

    const reco::Vertex the_pv = *( dynamic_cast<const pat::CompositeCandidate *>(i_combination->daughter("onia")) )->userData<reco::Vertex>("PVwithmuons");

    std::vector<reco::TransientTrack> tt_combination;
    tt_combination.push_back((*tt_builder).build(&tk_combination[0])); // muon
    tt_combination.push_back((*tt_builder).build(&tk_combination[1])); // muon
    tt_combination.push_back((*tt_builder).build(&tk_combination[2])); // track
    tt_combination.push_back((*tt_builder).build(&tk_combination[3])); // track

    KinematicParticleFactoryFromTransientTrack pFactory;

    const ParticleMass muonMass(0.1056583);
    float muonSigma = muonMass*1E-6;
    const ParticleMass trackMass1(track_mass_[0]);
    float trackSigma1 = trackMass1*1E-6;
    const ParticleMass trackMass2(track_mass_[1]);
    float trackSigma2 = trackMass2*1E-6;

    std::vector<RefCountedKinematicParticle> daughters_ditrack;
    daughters_ditrack.push_back(pFactory.particle(tt_combination[2], trackMass1, float(0), float(0), trackSigma1));
    daughters_ditrack.push_back(pFactory.particle(tt_combination[3], trackMass2, float(0), float(0), trackSigma2));

    KinematicParticleVertexFitter ditrack_fitter;
    RefCountedKinematicTree ditrack_tree = nullptr;
    ditrack_tree = ditrack_fitter.fit(daughters_ditrack);

    if (make_ditrack_vertex_  && !ditrack_tree->isValid()) continue;

    // apply mass constarint if any
    KinematicParticleFitter ditrack_cs_fitter;
    KinematicConstraint * ditrack_mtc = new MassKinematicConstraint(ditrack_mass_constraint_,float(0));
    if (make_ditrack_mass_constraint_ && ditrack_tree->isValid()) {
       ditrack_tree->movePointerToTheTop();
       ditrack_tree = ditrack_cs_fitter.fit(ditrack_mtc,ditrack_tree);
       if (!ditrack_tree->isValid()) continue;
    }

    std::vector<RefCountedKinematicParticle> daughters_combination;
    daughters_combination.push_back(pFactory.particle (tt_combination[0], muonMass, float(0), float(0), muonSigma));
    daughters_combination.push_back(pFactory.particle (tt_combination[1], muonMass, float(0), float(0), muonSigma));

    if (make_ditrack_vertex_) {
       ditrack_tree->movePointerToTheTop();
       daughters_combination.push_back(ditrack_tree->currentParticle());
    } else {
       daughters_combination.push_back(pFactory.particle (tt_combination[2], trackMass1, float(0), float(0), trackSigma1));
       daughters_combination.push_back(pFactory.particle (tt_combination[3], trackMass2, float(0), float(0), trackSigma2));
    }

    KinematicConstrainedVertexFitter constVertexFitter;
    MultiTrackKinematicConstraint *dimuon_mtc = new  TwoTrackMassKinematicConstraint(dimuon_mass_constraint_);
    RefCountedKinematicTree tree_combination = constVertexFitter.fit(daughters_combination,dimuon_mtc);

    if (tree_combination->isValid() && !tree_combination->isEmpty()) {
       tree_combination->movePointerToTheTop();
       RefCountedKinematicParticle combination_particle = tree_combination->currentParticle();
       RefCountedKinematicVertex combination_vertex = tree_combination->currentDecayVertex();

// Store combination reffited
       double combination_ma_fit = 14000.;
       double combination_vp_fit = -9999.;
       double combination_x2_fit = 10000.;
       if (combination_particle && combination_vertex && combination_particle->currentState().isValid()) {
         combination_ma_fit = combination_particle->currentState().mass();
         combination_x2_fit = combination_vertex->chiSquared();
	 if (combination_x2_fit < 100.) combination_vp_fit = ChiSquaredProbability(combination_x2_fit,
                                              (double)(combination_vertex->degreesOfFreedom()));
       }
       if ( combination_ma_fit > combination_mass_cut_[0] && combination_ma_fit < combination_mass_cut_[1] && combination_vp_fit > .0001) {
            TVector3 vtx;
            TVector3 pvtx;
            VertexDistanceXY vdistXY;
            int    combination_ch_fit = i_combination->charge();
            double combination_px_fit = combination_particle->currentState().kinematicParameters().momentum().x();
            double combination_py_fit = combination_particle->currentState().kinematicParameters().momentum().y();
            double combination_pz_fit = combination_particle->currentState().kinematicParameters().momentum().z();
            double combination_en_fit = sqrt(combination_ma_fit*combination_ma_fit+combination_px_fit*combination_px_fit+
                                      combination_py_fit*combination_py_fit+combination_pz_fit*combination_pz_fit);
            double combination_vx_fit = combination_vertex->position().x();
	    double combination_vy_fit = combination_vertex->position().y();
            double combination_vz_fit = combination_vertex->position().z();

            vtx.SetXYZ(combination_vx_fit,combination_vy_fit,0);
            TVector3 pperp(combination_px_fit, combination_py_fit, 0);
            AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
            pvtx.SetXYZ(the_pv.position().x(),the_pv.position().y(),0);
            TVector3 vdiff = vtx - pvtx;
            double cosAlpha = -100.;
	    if ( TMath::Abs(vdiff.Perp()*pperp.Perp())>0.00001 ) cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
            Measurement1D distXY = vdistXY.distance(reco::Vertex(*combination_vertex), the_pv);
            double ctauPV = -100.;
            if ( TMath::Abs(pperp.Perp()) > 0.00001) ctauPV = distXY.value()*cosAlpha * combination_ma_fit/pperp.Perp();
            GlobalError v1e = (reco::Vertex(*combination_vertex)).error();
            GlobalError v2e = the_pv.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
            double ctauErrPV = 100000.;
            if ( TMath::Abs(pperp.Perp2())>0.00001) ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*combination_ma_fit/(pperp.Perp2());

	    reco::CompositeCandidate recoCombo(combination_ch_fit,math::XYZTLorentzVector(combination_px_fit,combination_py_fit,combination_pz_fit,combination_en_fit),
                                              math::XYZPoint(combination_vx_fit,combination_vy_fit,combination_vz_fit),441);
	    pat::CompositeCandidate patCombo(recoCombo);
            patCombo.addUserFloat("vProb",combination_vp_fit);
            patCombo.addUserFloat("vChi2",combination_x2_fit);
            patCombo.addUserFloat("cosAlpha",cosAlpha);
            patCombo.addUserFloat("ctauPV",ctauPV);
            patCombo.addUserFloat("ctauErrPV",ctauErrPV);
            patCombo.addUserInt("comboIndex",index_combination);

// get first muon
            bool child = tree_combination->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitMu1 = tree_combination->currentParticle();
            if (!child) continue;
            float m1_ma_fit = fitMu1->currentState().mass();
            int   m1_ch_fit = fitMu1->currentState().particleCharge();
            float m1_px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
            float m1_py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
            float m1_pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
            float m1_en_fit = sqrt(m1_ma_fit*m1_ma_fit+m1_px_fit*m1_px_fit+m1_py_fit*m1_py_fit+m1_pz_fit*m1_pz_fit);
            reco::CompositeCandidate recoMu1(m1_ch_fit,math::XYZTLorentzVector(m1_px_fit,m1_py_fit,m1_pz_fit,m1_en_fit),
                                             math::XYZPoint(combination_vx_fit,combination_vy_fit,combination_vz_fit),13);
            pat::CompositeCandidate patMu1(recoMu1);
// get second muon
            child = tree_combination->movePointerToTheNextChild();
            RefCountedKinematicParticle fitMu2 = tree_combination->currentParticle();
            if (!child) continue;
            float m2_ma_fit = fitMu2->currentState().mass();
            int   m2_ch_fit = fitMu2->currentState().particleCharge();
            float m2_px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
            float m2_py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
            float m2_pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
            float m2_en_fit = sqrt(m2_ma_fit*m2_ma_fit+m2_px_fit*m2_px_fit+m2_py_fit*m2_py_fit+m2_pz_fit*m2_pz_fit);
            reco::CompositeCandidate recoMu2(m2_ch_fit,math::XYZTLorentzVector(m2_px_fit,m2_py_fit,m2_pz_fit,m2_en_fit),
                                             math::XYZPoint(combination_vx_fit,combination_vy_fit,combination_vz_fit),13);
            pat::CompositeCandidate patMu2(recoMu2);

// Define psi from two muons
	    pat::CompositeCandidate patDimuon;
	    patDimuon.addDaughter(patMu1,"muon1");
            patDimuon.addDaughter(patMu2,"muon2");
            patDimuon.setP4(patMu1.p4()+patMu2.p4());

// get tracks
            RefCountedKinematicParticle fitTrk;
            RefCountedKinematicParticle fitTrk2;
            double ditrack_vp_fit = -9999.;
            double ditrack_x2_fit = 10000.;

            if (make_ditrack_vertex_) {
               RefCountedKinematicVertex ditrack_vertex = ditrack_tree->currentDecayVertex();
               if ( ditrack_vertex ) {
                  ditrack_x2_fit = ditrack_vertex->chiSquared();
                  ditrack_vp_fit = ChiSquaredProbability(ditrack_x2_fit,(double)(ditrack_vertex->degreesOfFreedom()));
               }
               if (!ditrack_tree->movePointerToTheFirstChild()) continue;
               fitTrk = ditrack_tree->currentParticle();
               if (!ditrack_tree->movePointerToTheNextChild()) continue;
               fitTrk2 = ditrack_tree->currentParticle();
            } else {
               if (!tree_combination->movePointerToTheNextChild()) continue;
               fitTrk = tree_combination->currentParticle();
               if (!tree_combination->movePointerToTheNextChild()) continue;
               fitTrk2 = tree_combination->currentParticle();
            }

            float tk_ma_fit = fitTrk->currentState().mass();
            int   tk_ch_fit = fitTrk->currentState().particleCharge();
            float tk_px_fit = fitTrk->currentState().kinematicParameters().momentum().x();
            float tk_py_fit = fitTrk->currentState().kinematicParameters().momentum().y();
            float tk_pz_fit = fitTrk->currentState().kinematicParameters().momentum().z();
            float tk_en_fit = sqrt(tk_ma_fit*tk_ma_fit+tk_px_fit*tk_px_fit+tk_py_fit*tk_py_fit+tk_pz_fit*tk_pz_fit);
            reco::CompositeCandidate recoTk(tk_ch_fit,math::XYZTLorentzVector(tk_px_fit,tk_py_fit,tk_pz_fit,tk_en_fit),
                                             math::XYZPoint(combination_vx_fit,combination_vy_fit,combination_vz_fit),321);
            pat::CompositeCandidate patTk(recoTk);

            float tk2_ma_fit = fitTrk2->currentState().mass();
            int   tk2_ch_fit = fitTrk2->currentState().particleCharge();
            float tk2_px_fit = fitTrk2->currentState().kinematicParameters().momentum().x();
            float tk2_py_fit = fitTrk2->currentState().kinematicParameters().momentum().y();
            float tk2_pz_fit = fitTrk2->currentState().kinematicParameters().momentum().z();
            float tk2_en_fit = sqrt(tk2_ma_fit*tk2_ma_fit+tk2_px_fit*tk2_px_fit+tk2_py_fit*tk2_py_fit+tk2_pz_fit*tk2_pz_fit);
            reco::CompositeCandidate recoTk2(tk2_ch_fit,math::XYZTLorentzVector(tk2_px_fit,tk2_py_fit,tk2_pz_fit,tk2_en_fit),
                                             math::XYZPoint(combination_vx_fit,combination_vy_fit,combination_vz_fit),321);
            pat::CompositeCandidate patTk2(recoTk2);

// Define psi from two muons
            pat::CompositeCandidate patDitrack;
            patDitrack.addDaughter(patTk,"trak1");
            patDitrack.addDaughter(patTk2,"trak2");
            patDitrack.setP4(patTk.p4()+patTk2.p4());
            patDitrack.addUserFloat("vProb",ditrack_vp_fit);
            patDitrack.addUserFloat("vChi2",ditrack_x2_fit);

	    patCombo.addDaughter(patDimuon,"dimuon");
	    patCombo.addDaughter(patDitrack,"ditrack");

            CombinationCollection->push_back(patCombo);
          }
	}
  }
// End kinematic fit

  }
// now sort by vProb
  DiMuonDiTrackKF::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  std::sort(CombinationCollection->begin(),CombinationCollection->end(),vPComparator);
  iEvent.put(std::move(CombinationCollection),product_name_);

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonDiTrackKF::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrackKF);
