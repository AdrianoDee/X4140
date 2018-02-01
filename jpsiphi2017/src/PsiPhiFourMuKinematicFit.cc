// -*- C++ -*-
//
// Package:    PsiPhiFourMuKinematicFit
// Class:      PsiPhiFourMuKinematicFit
//
/**\class PsiPhiFourMuKinematicFit Ponia/OniaTrak/src/PsiPhiFourMuKinematicFit.cc

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

class PsiPhiFourMuKinematicFit : public edm::EDProducer {
   public:
      explicit PsiPhiFourMuKinematicFit(const edm::ParameterSet&);
      ~PsiPhiFourMuKinematicFit() override;

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
  edm::EDGetTokenT<pat::CompositeCandidateCollection> oniafourmu_cand_;
  double mass_phi,mass_jpsi;
  std::vector<double> FourOniaMassCuts_;
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
PsiPhiFourMuKinematicFit::PsiPhiFourMuKinematicFit(const edm::ParameterSet& iConfig) {
  oniafourmu_cand_   = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("PsiTrakTrak"));
  mass_phi         = iConfig.getParameter<double>("phi_constraint");
  mass_jpsi         = iConfig.getParameter<double>("jpsi_constraint");
  FourOniaMassCuts_ = iConfig.getParameter<std::vector<double>>("OniaTrakTrakMassCuts");
  product_name_ = iConfig.getParameter<std::string>("product_name");

// kinematic refit collections
  produces<pat::CompositeCandidateCollection>(product_name_);

// now do what ever other initialization is needed
}

PsiPhiFourMuKinematicFit::~PsiPhiFourMuKinematicFit() {
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void PsiPhiFourMuKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> PsiPhiCandHandle;
  iEvent.getByToken(oniafourmu_cand_, PsiPhiCandHandle);

// Kinematic refit collection
  std::unique_ptr< pat::CompositeCandidateCollection > PsiPhiCandRefitColl(new pat::CompositeCandidateCollection);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  int indexB0s=-1;
  for (pat::CompositeCandidateCollection::const_iterator oniat=PsiPhiCandHandle->begin(); oniat!=PsiPhiCandHandle->end(); ++oniat) {

    indexB0s++;
    reco::TrackRef JpsiPhiTk[4]={
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("jpsi")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("jpsi")->daughter("muon2") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("phi")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("phi")->daughter("muon2") ) )->innerTrack()
    };

    reco::TrackRef JpsiTk[2]={
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("jpsi")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("jpsi")->daughter("muon2") ) )->innerTrack()
    };

    reco::TrackRef PhiTk[2]={
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("phi")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("phi")->daughter("muon2") ) )->innerTrack()
    };

    std::vector<reco::TransientTrack> JpsiTT;
    JpsiTT.push_back((*theB).build(&JpsiTk[0]));
    JpsiTT.push_back((*theB).build(&JpsiTk[1]));

    std::vector<reco::TransientTrack> PhiTT;
    PhiTT.push_back((*theB).build(PhiTk[0]));
    PhiTT.push_back((*theB).build(PhiTk[1]));

    const pat::CompositeCandidate *jpsiCand = dynamic_cast<const pat::CompositeCandidate *>(oniat->daughter("jpsi"));
    const reco::Vertex thePrimaryV = *jpsiCand->userData<reco::Vertex>("PVwithmuons");

    std::vector<reco::TransientTrack> fourMuTT;
    fourMuTT.push_back((*theB).build(&JpsiPhiTk[0]));
    fourMuTT.push_back((*theB).build(&JpsiPhiTk[1]));
    fourMuTT.push_back((*theB).build(&JpsiPhiTk[2]));
    fourMuTT.push_back((*theB).build(&JpsiPhiTk[3]));

    KinematicParticleFactoryFromTransientTrack pFactory;

    const ParticleMass muonMass(0.1056583);
    float muonSigma = muonMass*1E-6;

    std::vector<RefCountedKinematicParticle> allFourMuDaughters;
    allFourMuDaughters.push_back(pFactory.particle (fourMuTT[0], muonMass, float(0), float(0), muonSigma));
    allFourMuDaughters.push_back(pFactory.particle (fourMuTT[1], muonMass, float(0), float(0), muonSigma));
    allFourMuDaughters.push_back(pFactory.particle (fourMuTT[2], muonMass, float(0), float(0), muonSigma));
    allFourMuDaughters.push_back(pFactory.particle (fourMuTT[3], muonMass, float(0), float(0), muonSigma));

    std::vector<RefCountedKinematicParticle> jpsiParticles;
    jpsiParticles.push_back(pFactory.particle(JpsiTT[0],muonMass,float(0),float(0),muonSigma));
    jpsiParticles.push_back(pFactory.particle(JpsiTT[1],muonMass,float(0),float(0),muonSigma));

    std::vector<RefCountedKinematicParticle> phiParticles;
    phiParticles.push_back(pFactory.particle(PhiTT[0],muonMass,float(0),float(0),muonSigma));
    phiParticles.push_back(pFactory.particle(PhiTT[1],muonMass,float(0),float(0),muonSigma));

    KinematicParticleVertexFitter fitter;

    RefCountedKinematicTree jpsiVertexFitTree;
    jpsiVertexFitTree = fitter.fit(jpsiParticles);

    if(jpsiVertexFitTree->isValid())
    {

      const ParticleMass jpsi_mass(mass_jpsi);
      float jpsi_sigma = 1E-6 * jpsi_mass;

      KinematicParticleFitter csFitterJpsi;
      KinematicConstraint * jpsi_c = new MassKinematicConstraint(jpsi_mass,jpsi_sigma);

      jpsiVertexFitTree->movePointerToTheTop();
      jpsiVertexFitTree = csFitterJpsi.fit(jpsi_c,jpsiVertexFitTree);

      if (jpsiVertexFitTree->isValid()) {

        RefCountedKinematicTree phiVertexFitTree;
        phiVertexFitTree = fitter.fit(phiParticles);

        jpsiVertexFitTree->movePointerToTheTop();
      	RefCountedKinematicParticle fitJpsi = jpsiVertexFitTree->currentParticle();

        std::vector<RefCountedKinematicParticle> allB0sDaughters;
        allB0sDaughters.push_back(pFactory.particle (PhiTT[0], muonMass, float(0), float(0), muonSigma));
        allB0sDaughters.push_back(pFactory.particle (PhiTT[1], muonMass, float(0), float(0), muonSigma));
      	allB0sDaughters.push_back(fitJpsi);

        KinematicConstrainedVertexFitter constVertexFitter;

      	MultiTrackKinematicConstraint *phi_mtc = new  TwoTrackMassKinematicConstraint(mass_phi);
      	RefCountedKinematicTree B0sTree = constVertexFitter.fit(allB0sDaughters,phi_mtc);

        if (!B0sTree->isEmpty())
        {
          B0sTree->movePointerToTheTop();
          RefCountedKinematicParticle fitB0s = B0sTree->currentParticle();
          RefCountedKinematicVertex B0sDecayVertex = B0sTree->currentDecayVertex();

          if (fitB0s->currentState().isValid())
          {

            float B0sM_fit  = fitB0s->currentState().mass();
	          float B0sPx_fit = fitB0s->currentState().kinematicParameters().momentum().x();
            float B0sPy_fit = fitB0s->currentState().kinematicParameters().momentum().y();
	          float B0sPz_fit = fitB0s->currentState().kinematicParameters().momentum().z();
            float B0sVtxX_fit = B0sDecayVertex->position().x();
            float B0sVtxY_fit = B0sDecayVertex->position().y();
            float B0sVtxZ_fit = B0sDecayVertex->position().z();
            float B0s_x2_fit = B0sDecayVertex->chiSquared();
            float B0s_nDof_fit = B0sDecayVertex->degreesOfFreedom();
            float B0sVtxP_fit = ChiSquaredProbability((double)(B0s_x2_fit),
                                                       (double)(B0s_nDof_fit));
            float B0s_en_fit = sqrt(B0sM_fit*B0sM_fit+B0sPx_fit*B0sPx_fit+
                                                       B0sPy_fit*B0sPy_fit+B0sPz_fit*B0sPz_fit);

            if ( B0sM_fit < FourOniaMassCuts_[0] || B0sM_fit > FourOniaMassCuts_[1])
            continue;

            TVector3 vtx;
            TVector3 pvtx;
            VertexDistanceXY vdistXY;

            vtx.SetXYZ(B0sVtxX_fit,B0sVtxY_fit,0);

            TVector3 pperp(B0sPx_fit, B0sPy_fit, 0);
            AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
            pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
            TVector3 vdiff = vtx - pvtx;
            double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
            Measurement1D distXY = vdistXY.distance(reco::Vertex(*B0sDecayVertex), thePrimaryV);
            double ctauPV = distXY.value()*cosAlpha * B0sM_fit/pperp.Perp();
            GlobalError v1e = (reco::Vertex(*B0sDecayVertex)).error();
            GlobalError v2e = thePrimaryV.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
            double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*B0sM_fit/(pperp.Perp2());

            reco::CompositeCandidate recoB0s(0, math::XYZTLorentzVector(B0sPx_fit, B0sPy_fit, B0sPz_fitB0s_en_fit), math::XYZPoint(B0sVtxX_fit,
                                              B0sVtxY_fit, B0sVtxZ_fit), 531);



	          pat::CompositeCandidate patB0s(recoB0s);
            patB0s.addUserFloat("vProb",B0sVtxP_fit);
            patB0s.addUserFloat("vChi2",B0s_x2_fit);
            patB0s.addUserFloat("cosAlpha",cosAlpha);
            patB0s.addUserFloat("ctauPV",ctauPV);
            patB0s.addUserFloat("ctauErrPV",ctauErrPV);

            patB0s.addUserInt("bIndex",indexB0s);
            //get first muon
            bool child = B0sTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitPhiMu1 = B0sTree->currentParticle();
            if (!child) break;

            float phiMu1M_fit  = fitPhiMu1->currentState().mass();
            float phiMu1Q_fit  = fitPhiMu1->currentState().particleCharge();
            float phiMu1Px_fit = fitPhiMu1->currentState().kinematicParameters().momentum().x();
	          float phiMu1Py_fit = fitPhiMu1->currentState().kinematicParameters().momentum().y();
            float phiMu1Pz_fit = fitPhiMu1->currentState().kinematicParameters().momentum().z();
	          reco::CompositeCandidate recoPhiMu1(phiMu1Q_fit, math::XYZTLorentzVector(phiMu1Px_fit, phiMu1Py_fit, phiMu1Pz_fit,
                                             sqrt(phiMu1M_fit*phiMu1M_fit + phiMu1Px_fit*phiMu1Px_fit + phiMu1Py_fit*phiMu1Py_fit +
                                             phiMu1Pz_fit*phiMu1Pz_fit)), math::XYZPoint(B0sVtxX_fit, B0sVtxY_fit, B0sVtxZ_fit), 13);
	          pat::CompositeCandidate patPhiMu1(recoPhiMu1);

            //get second muon
            child = B0sTree->movePointerToTheNextChild();
	          RefCountedKinematicParticle fitPhiMu2 = B0sTree->currentParticle();
            if (!child) break;

	          float phiMu2M_fit  = fitPhiMu2->currentState().mass();
            float phiMu2Q_fit  = fitPhiMu2->currentState().particleCharge();
            float phiMu2Px_fit = fitPhiMu2->currentState().kinematicParameters().momentum().x();
	          float phiMu2Py_fit = fitPhiMu2->currentState().kinematicParameters().momentum().y();
	          float phiMu2Pz_fit = fitPhiMu2->currentState().kinematicParameters().momentum().z();
            reco::CompositeCandidate recoPhiMu2(phiMu2Q_fit, math::XYZTLorentzVector(phiMu2Px_fit, phiMu2Py_fit, phiMu2Pz_fit,
                                             sqrt(phiMu2M_fit*phiMu2M_fit + phiMu2Px_fit*phiMu2Px_fit + phiMu2Py_fit*phiMu2Py_fit +
                                             phiMu2Pz_fit*phiMu2Pz_fit)), math::XYZPoint(B0sVtxX_fit, B0sVtxY_fit, B0sVtxZ_fit), 13);
            pat::CompositeCandidate patPhiMu2(recoPhiMu2);

            pat::CompositeCandidate phiRefit;
            phiRefit.addDaughter(patPhiMu1,"muon1");
            phiRefit.addDaughter(patPhiMu2,"muon2");
            phiRefit.setP4(patPhiMu1.p4()+patPhiMu2.p4());

            child = B0sTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitJpsi = B0sTree->currentParticle();
	          if (!child) break;

            float jpsiM_fit  = fitJpsi->currentState().mass();
            float jpsiPx_fit = fitJpsi->currentState().kinematicParameters().momentum().x();
	          float jpsiPy_fit = fitJpsi->currentState().kinematicParameters().momentum().y();
            float jpsiPz_fit = fitJpsi->currentState().kinematicParameters().momentum().z();

            reco::CompositeCandidate recoJpsi(0, math::XYZTLorentzVector(jpsiPx_fit, jpsiPy_fit, jpsiPz_fit,
                                               sqrt(jpsiM_fit*jpsiM_fit + jpsiPx_fit*jpsiPx_fit + jpsiPy_fit*jpsiPy_fit +
                                               jpsiPz_fit*jpsiPz_fit)), math::XYZPoint(B0sVtxX_fit, B0sVtxY_fit, B0sVtxZ_fit), 22);
            pat::CompositeCandidate patJpsi(recoJpsi);

            patB0s.addDaughter(phiRefit,"phi");
            patB0s.addDaughter(patJpsi,"jpsi");

            PsiPhiCandRefitColl->push_back(patB0s);
          }

        }

      }

    }
  }


// now sort by vProb
  PsiPhiFourMuKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  std::sort(PsiPhiCandRefitColl->begin(),PsiPhiCandRefitColl->end(),vPComparator);
  iEvent.put(std::move(PsiPhiCandRefitColl),product_name_);
}

// ------------ method called once each job just before starting event loop  ------------
void PsiPhiFourMuKinematicFit::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PsiPhiFourMuKinematicFit::endJob() {}

// ------------ method called when starting to processes a run  ------------
void PsiPhiFourMuKinematicFit::beginRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void PsiPhiFourMuKinematicFit::endRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void PsiPhiFourMuKinematicFit::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void PsiPhiFourMuKinematicFit::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PsiPhiFourMuKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PsiPhiFourMuKinematicFit);
