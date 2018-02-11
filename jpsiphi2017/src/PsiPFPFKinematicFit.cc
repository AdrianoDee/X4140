// -*- C++ -*-
//
// Package:    PsiPFPFKinematicFit
// Class:      PsiPFPFKinematicFit
//
/**\class PsiPFPFKinematicFit Ponia/OniaTrak/src/PsiPFPFKinematicFit.cc

 Description: performs vertex kinematical fit for J/psi Track

**/
//
// Original Author:  Alberto Sanchez-Hernandez
//         Created:  Abril 2014

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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


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

class PsiPFPFKinematicFit : public edm::EDProducer {
   public:
      explicit PsiPFPFKinematicFit(const edm::ParameterSet&);
      ~PsiPFPFKinematicFit() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginJob() override ;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endJob() override ;

      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

// ----------member data ---------------------------
  edm::EDGetTokenT<pat::CompositeCandidateCollection> oniat_cand_;
  double mass_;
  std::vector<double> OniaTrakTrakMassCuts_;
  std::vector<double> MassTraks_;
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
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PsiPFPFKinematicFit::PsiPFPFKinematicFit(const edm::ParameterSet& iConfig) {
  oniat_cand_   = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("PsiPFPF"));
  mass_         = iConfig.getParameter<double>("mass_constraint");
  OniaTrakTrakMassCuts_ = iConfig.getParameter<std::vector<double>>("OniaTrakTrakMassCuts");
  MassTraks_    = iConfig.getParameter<std::vector<double>>("MassTraks");
  product_name_ = iConfig.getParameter<std::string>("product_name");

// kinematic refit collections
  produces<pat::CompositeCandidateCollection>(product_name_);

// now do what ever other initialization is needed
}

PsiPFPFKinematicFit::~PsiPFPFKinematicFit() {
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void PsiPFPFKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  int debug = 0;
  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> PsiTCandHandle;
  iEvent.getByToken(oniat_cand_, PsiTCandHandle);

// Kinematic refit collection
  std::unique_ptr< pat::CompositeCandidateCollection > PsiTCandRefitColl(new pat::CompositeCandidateCollection);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  std::cout<<"-Debug "<<++debug<<std::endl;
  int indexPsiT=-1;
  for (pat::CompositeCandidateCollection::const_iterator oniat=PsiTCandHandle->begin(); oniat!=PsiTCandHandle->end(); ++oniat) {

    const pat::CompositeCandidate *dimuonC = dynamic_cast<const pat::CompositeCandidate *>(oniat->daughter("onia"));
    const pat::CompositeCandidate *ditrakC = dynamic_cast<const pat::CompositeCandidate*>(oniat->daughter("ditrak"));

    if(dimuonC->userFloat("vProb")<0.0)
      continue;

    indexPsiT++;

    std::vector <reco::TrackRef> JpsiTk;

    const pat::PackedCandidate *trak1 = dynamic_cast<const pat::PackedCandidate*>(ditrakC->daughter("trak1"));
    const pat::PackedCandidate *trak2 = dynamic_cast<const pat::PackedCandidate*>(ditrakC->daughter("trak2"));

    JpsiTk.push_back(( dynamic_cast<const pat::Muon*>(oniat->daughter("onia")->daughter("muon1") ) )->innerTrack());
    JpsiTk.push_back(( dynamic_cast<const pat::Muon*>(oniat->daughter("onia")->daughter("muon2") ) )->innerTrack());
    std::cout<<"Debug "<<++debug<<std::endl;
    if(!trak1->hasTrackDetails())
      continue;
    if(!trak2->hasTrackDetails())
      continue;

    const reco::Vertex thePrimaryV = *dimuonC->userData<reco::Vertex>("PVwithmuons");
    std::cout<<"Debug 0"<<++debug<<std::endl;
    std::vector<reco::TransientTrack> MuMuTT;
    MuMuTT.push_back((*theB).build(&JpsiTk[0]));
    MuMuTT.push_back((*theB).build(&JpsiTk[1]));
    MuMuTT.push_back((*theB).build(&(trak1->pseudoTrack()))); // K+
    MuMuTT.push_back((*theB).build(&(trak2->pseudoTrack()))); // K+

    KinematicParticleFactoryFromTransientTrack pFactory;

    std::cout<<"Debug 1"<<++debug<<std::endl;
    const ParticleMass muonMass(0.1056583);
    float muonSigma = muonMass*1E-6;
    const ParticleMass kaonMass1(MassTraks_[0]);
    float kaonSigma1 = kaonMass1*1E-6;
    const ParticleMass kaonMass2(MassTraks_[1]);
    float kaonSigma2 = kaonMass2*1E-6;

    std::vector<RefCountedKinematicParticle> allPsiTDaughters;
    allPsiTDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
    allPsiTDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
    allPsiTDaughters.push_back(pFactory.particle (MuMuTT[2], kaonMass1, float(0), float(0), kaonSigma1));
    allPsiTDaughters.push_back(pFactory.particle (MuMuTT[3], kaonMass2, float(0), float(0), kaonSigma2));

    KinematicConstrainedVertexFitter constVertexFitter;
    MultiTrackKinematicConstraint *onia_mtc = new  TwoTrackMassKinematicConstraint(mass_);
    RefCountedKinematicTree PsiTTree = constVertexFitter.fit(allPsiTDaughters,onia_mtc);

    if (!PsiTTree->isEmpty()) {
       PsiTTree->movePointerToTheTop();
       RefCountedKinematicParticle fitPsiT = PsiTTree->currentParticle();
       RefCountedKinematicVertex PsiTDecayVertex = PsiTTree->currentDecayVertex();
// Get PsiT reffited
       double oniat_ma_fit = 14000.;
       double oniat_vp_fit = -9999.;
       double oniat_x2_fit = 10000.;
       std::cout<<"Debug 2"<<++debug<<std::endl;
       if (fitPsiT->currentState().isValid()) {
         oniat_ma_fit = fitPsiT->currentState().mass();
         oniat_x2_fit = PsiTDecayVertex->chiSquared();
         oniat_vp_fit = ChiSquaredProbability(oniat_x2_fit,
                                              (double)(PsiTDecayVertex->degreesOfFreedom()));
       }
       if ( oniat_ma_fit > OniaTrakTrakMassCuts_[0] && oniat_ma_fit < OniaTrakTrakMassCuts_[1] && oniat_vp_fit > 0.0 ) {
            TVector3 vtx;
            TVector3 pvtx;
            VertexDistanceXY vdistXY;
            int   oniat_ch_fit = oniat->charge();
            double oniat_px_fit = fitPsiT->currentState().kinematicParameters().momentum().x();
            double oniat_py_fit = fitPsiT->currentState().kinematicParameters().momentum().y();
            double oniat_pz_fit = fitPsiT->currentState().kinematicParameters().momentum().z();
            double oniat_en_fit = sqrt(oniat_ma_fit*oniat_ma_fit+oniat_px_fit*oniat_px_fit+
                                      oniat_py_fit*oniat_py_fit+oniat_pz_fit*oniat_pz_fit);
            double oniat_vx_fit = PsiTDecayVertex->position().x();
	          double oniat_vy_fit = PsiTDecayVertex->position().y();
            double oniat_vz_fit = PsiTDecayVertex->position().z();

            vtx.SetXYZ(oniat_vx_fit,oniat_vy_fit,0);
            TVector3 pperp(oniat_px_fit, oniat_py_fit, 0);
            AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
            pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
            TVector3 vdiff = vtx - pvtx;
            double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
            Measurement1D distXY = vdistXY.distance(reco::Vertex(*PsiTDecayVertex), thePrimaryV);
            double ctauPV = distXY.value()*cosAlpha * oniat_ma_fit/pperp.Perp();
            GlobalError v1e = (reco::Vertex(*PsiTDecayVertex)).error();
            GlobalError v2e = thePrimaryV.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
            double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*oniat_ma_fit/(pperp.Perp2());
            std::cout<<"Debug 3"<<++debug<<std::endl;
	    reco::CompositeCandidate recoPsiT(oniat_ch_fit,math::XYZTLorentzVector(oniat_px_fit,oniat_py_fit,oniat_pz_fit,oniat_en_fit),
                                               math::XYZPoint(oniat_vx_fit,oniat_vy_fit,oniat_vz_fit),531);
	    pat::CompositeCandidate patPsiT(recoPsiT);
            patPsiT.addUserFloat("vProb",oniat_vp_fit);
            patPsiT.addUserFloat("vChi2",oniat_x2_fit);
            patPsiT.addUserFloat("cosAlpha",cosAlpha);
            patPsiT.addUserFloat("ctauPV",ctauPV);
            patPsiT.addUserFloat("ctauErrPV",ctauErrPV);

            patPsiT.addUserInt("bIndex",indexPsiT);

// get first muon
            bool child = PsiTTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitMu1 = PsiTTree->currentParticle();
            if (!child) break;
            float m1_ma_fit = fitMu1->currentState().mass();
            int   m1_ch_fit = fitMu1->currentState().particleCharge();
            float m1_px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
            float m1_py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
            float m1_pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
            float m1_en_fit = sqrt(m1_ma_fit*m1_ma_fit+m1_px_fit*m1_px_fit+m1_py_fit*m1_py_fit+m1_pz_fit*m1_pz_fit);
            reco::CompositeCandidate recoMu1(m1_ch_fit,math::XYZTLorentzVector(m1_px_fit,m1_py_fit,m1_pz_fit,m1_en_fit),
                                             math::XYZPoint(oniat_vx_fit,oniat_vy_fit,oniat_vz_fit),13);
            pat::CompositeCandidate patMu1(recoMu1);
// get second muon
            child = PsiTTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitMu2 = PsiTTree->currentParticle();
            if (!child) break;
            float m2_ma_fit = fitMu2->currentState().mass();
            int   m2_ch_fit = fitMu2->currentState().particleCharge();
            float m2_px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
            float m2_py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
            float m2_pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
            float m2_en_fit = sqrt(m2_ma_fit*m2_ma_fit+m2_px_fit*m2_px_fit+m2_py_fit*m2_py_fit+m2_pz_fit*m2_pz_fit);
            reco::CompositeCandidate recoMu2(m2_ch_fit,math::XYZTLorentzVector(m2_px_fit,m2_py_fit,m2_pz_fit,m2_en_fit),
                                             math::XYZPoint(oniat_vx_fit,oniat_vy_fit,oniat_vz_fit),13);
            pat::CompositeCandidate patMu2(recoMu2);
            std::cout<<"Debug 4"<<++debug<<std::endl;
// Define psi from two muons
	          pat::CompositeCandidate psi;
	          psi.addDaughter(patMu1,"muon1");
            psi.addDaughter(patMu2,"muon2");
            psi.setP4(patMu1.p4()+patMu2.p4());
// get kaon
            child = PsiTTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitTrk = PsiTTree->currentParticle();
            if (!child) break;
            float tk_ma_fit = fitTrk->currentState().mass();
            int   tk_ch_fit = fitTrk->currentState().particleCharge();
            float tk_px_fit = fitTrk->currentState().kinematicParameters().momentum().x();
            float tk_py_fit = fitTrk->currentState().kinematicParameters().momentum().y();
            float tk_pz_fit = fitTrk->currentState().kinematicParameters().momentum().z();
            float tk_en_fit = sqrt(tk_ma_fit*tk_ma_fit+tk_px_fit*tk_px_fit+tk_py_fit*tk_py_fit+tk_pz_fit*tk_pz_fit);
            reco::CompositeCandidate recoTk(tk_ch_fit,math::XYZTLorentzVector(tk_px_fit,tk_py_fit,tk_pz_fit,tk_en_fit),
                                             math::XYZPoint(oniat_vx_fit,oniat_vy_fit,oniat_vz_fit),321);
            pat::CompositeCandidate patTk(recoTk);
            std::cout<<"Debug 5"<<++debug<<std::endl;
// get kaon2
            child = PsiTTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitTrk2 = PsiTTree->currentParticle();
            if (!child) break;
            float tk2_ma_fit = fitTrk2->currentState().mass();
            int   tk2_ch_fit = fitTrk2->currentState().particleCharge();
            float tk2_px_fit = fitTrk2->currentState().kinematicParameters().momentum().x();
            float tk2_py_fit = fitTrk2->currentState().kinematicParameters().momentum().y();
            float tk2_pz_fit = fitTrk2->currentState().kinematicParameters().momentum().z();
            float tk2_en_fit = sqrt(tk2_ma_fit*tk2_ma_fit+tk2_px_fit*tk2_px_fit+tk2_py_fit*tk2_py_fit+tk2_pz_fit*tk2_pz_fit);
            reco::CompositeCandidate recoTk2(tk2_ch_fit,math::XYZTLorentzVector(tk2_px_fit,tk2_py_fit,tk2_pz_fit,tk2_en_fit),
                                             math::XYZPoint(oniat_vx_fit,oniat_vy_fit,oniat_vz_fit),321);
            pat::CompositeCandidate patTk2(recoTk2);

// Define psi from two muons
            pat::CompositeCandidate phi;
            phi.addDaughter(patTk,"trak1");
            phi.addDaughter(patTk2,"trak2");
            phi.setP4(patTk.p4()+patTk2.p4());

	          patPsiT.addDaughter(psi,"onia");
	          patPsiT.addDaughter(phi,"ditrak");
            std::cout<<"Debug 6"<<++debug<<std::endl;
            PsiTCandRefitColl->push_back(patPsiT);
          }
	}
  }
// End kinematic fit

// now sort by vProb
  PsiPFPFKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  std::sort(PsiTCandRefitColl->begin(),PsiTCandRefitColl->end(),vPComparator);
  iEvent.put(std::move(PsiTCandRefitColl),product_name_);

  std::cout << "PsiTCandRefitColl size: "<< PsiTCandRefitColl->size()<<std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void PsiPFPFKinematicFit::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PsiPFPFKinematicFit::endJob() {}

// ------------ method called when starting to processes a run  ------------
void PsiPFPFKinematicFit::beginRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void PsiPFPFKinematicFit::endRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void PsiPFPFKinematicFit::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void PsiPFPFKinematicFit::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PsiPFPFKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PsiPFPFKinematicFit);
