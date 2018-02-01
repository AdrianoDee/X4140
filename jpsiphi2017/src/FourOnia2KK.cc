#include "../interface/FourOnia2KK.h"

//Headers for the data items
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

//Headers for services and tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "../interface/FourOniaVtxReProducer.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

FourOnia2KKPAT::FourOnia2KKPAT(const edm::ParameterSet& iConfig):
trakCollection_(consumes<edm::View<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("tracks"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
higherPuritySelection_(iConfig.getParameter<std::string>("higherPuritySelection")),
lowerPuritySelection_(iConfig.getParameter<std::string>("lowerPuritySelection")),
dimuonSelection_(iConfig.existsAs<std::string>("dimuonSelection") ? iConfig.getParameter<std::string>("dimuonSelection") : ""),
addCommonVertex_(iConfig.getParameter<bool>("addCommonVertex")),
addMuonlessPrimaryVertex_(iConfig.getParameter<bool>("addMuonlessPrimaryVertex")),
resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity")),
addMCTruth_(iConfig.getParameter<bool>("addMCTruth")),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("HLTFilters"))
{
  revtxtrks_ = consumes<reco::TrackCollection>((edm::InputTag)"generalTracks"); //if that is not true, we will raise an exception
  revtxbs_ = consumes<reco::BeamSpot>((edm::InputTag)"offlineBeamSpot");
  produces<pat::CompositeCandidateCollection>();
}


FourOnia2KKPAT::~FourOnia2KKPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

UInt_t FourOnia2KKPAT::isTriggerMatched(pat::CompositeCandidate *diMuon_cand) {
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
    // std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
  }

  return matched;
}

// ------------ method called to produce the data  ------------

void
FourOnia2KKPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  vector<double> kMasses;
  kMasses.push_back( 0.493677 );
  kMasses.push_back( 0.493677 );

  const ParticleMass kaon_mass(0.493677);
  float kaon_sigma = 1E-6;

  const ParticleMass phi_mass(1.019461);
  float phi_sigma = 1E-6;

  std::unique_ptr<pat::CompositeCandidateCollection> phiOutput(new pat::CompositeCandidateCollection);
  //std::cout<<"Four muonia producer"<<std::endl;
  Vertex thePrimaryV;
  Vertex theBeamSpotV;


  Handle<BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  Handle<VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
  }

  edm::Handle< View<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_,thePATTrackHandle);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  // JPsi candidates only from muons

  double TrMaxNormChi2 = 10.0;
  double TrMinPt = 0.0;
  bool do_phimass_constrain = false;

  std::cout<<"M"<< std::endl;

  for(View<pat::GenericParticle>::const_iterator kTrack1 = thePATTrackHandle->begin();kTrack1 != thePATTrackHandle->end(); ++kTrack1 )
  {
    if(kTrack1->charge()==0) continue;

    // if ((kTrack1->track() ==  nullptr)) continue;
    if ((kTrack1->track()->chi2() / kTrack1->track()->ndof() > TrMaxNormChi2)  ||  kTrack1->pt() < TrMinPt) continue;

    std::cout << kTrack1->track()->pt() << " - " << kTrack1->pt() << std::endl;

    for(View<pat::GenericParticle>::const_iterator kTrack2 = kTrack1+1; kTrack2 != thePATTrackHandle->end(); ++kTrack2 )
    {
      // if ((kTrack2->track() ==  nullptr)) continue;
      if(kTrack1==kTrack2) continue;
      if(kTrack2->charge()==0) continue;

      if ((kTrack2->track()->chi2() / kTrack2->track()->ndof() > TrMaxNormChi2)  ||  kTrack2->pt() < TrMinPt) continue;

      if (kTrack1->charge() * kTrack2->charge() > 0) continue;//TODO CHECK IF phi->K0K0 ... ?

      std::vector<reco::TransientTrack> phiTracks;
      phiTracks.push_back((*theTTBuilder).build(kTrack1->track()));
      phiTracks.push_back((*theTTBuilder).build(kTrack2->track()));

      KinematicParticleFactoryFromTransientTrack pFactory;
      std::vector<RefCountedKinematicParticle> PhiParticles;
      PhiParticles.push_back(pFactory.particle(phiTracks[0],kaon_mass,float(0),float(0),kaon_sigma));
      PhiParticles.push_back(pFactory.particle(phiTracks[1],kaon_mass,float(0),float(0),kaon_sigma));

      KinematicParticleVertexFitter fitter;
      RefCountedKinematicTree phiVertexFitTree;
      phiVertexFitTree = fitter.fit(PhiParticles);

      if (!phiVertexFitTree->isValid())
      {
        edm::ParameterSet pSet;
        pSet.addParameter<double>("maxDistance", 3);
        pSet.addParameter<int>("maxNbrOfIterations", 10000);
        KinematicParticleVertexFitter fitter2(pSet);
        phiVertexFitTree = fitter2.fit(PhiParticles);
      }


      if (phiVertexFitTree->isValid())
      {
        KinematicParticleFitter fitterPhi;
        KinematicConstraint * phi_const = new MassKinematicConstraint(phi_mass,phi_sigma);

        phiVertexFitTree->movePointerToTheTop();

        if(do_phimass_constrain)
        phiVertexFitTree = fitterPhi.fit(phi_const,phiVertexFitTree);

        if (phiVertexFitTree->isValid())
        {

          phiVertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle kkCandFitted = phiVertexFitTree->currentParticle();
          RefCountedKinematicVertex kkDecayVertex = phiVertexFitTree->currentDecayVertex();

          float kkM_fit  = kkCandFitted->currentState().mass();
          float kkPx_fit = kkCandFitted->currentState().kinematicParameters().momentum().x();
          float kkPy_fit = kkCandFitted->currentState().kinematicParameters().momentum().y();
          float kkPz_fit = kkCandFitted->currentState().kinematicParameters().momentum().z();
          float kkVtxX_fit = kkDecayVertex->position().x();
          float kkVtxY_fit = kkDecayVertex->position().y();
          float kkVtxZ_fit = kkDecayVertex->position().z();
          float kkVtxP_fit = ChiSquaredProbability((double)(kkDecayVertex->chiSquared()),
          (double)(kkDecayVertex->degreesOfFreedom()));

          reco::CompositeCandidate recoKKCand(0, math::XYZTLorentzVector(kkPx_fit, kkPy_fit, kkPz_fit,
            sqrt(kkM_fit*kkM_fit + kkPx_fit*kkPx_fit + kkPy_fit*kkPy_fit +
            kkPz_fit*kkPz_fit)), math::XYZPoint(kkVtxX_fit,
            kkVtxY_fit, kkVtxZ_fit), 50551);

          pat::CompositeCandidate patKKCand(recoKKCand);

          patKKCand.addUserFloat("vProb",kkVtxP_fit);

          phiOutput->push_back(patKKCand);

        }

      }



    }
  }

  iEvent.put(std::move(phiOutput));

}

      bool
      FourOnia2KKPAT::isAbHadron(int pdgID) {

        if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
        return false;

      }

      bool
      FourOnia2KKPAT::isAMixedbHadron(int pdgID, int momPdgID) {

        if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
        (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
        return true;
        return false;

      }

      std::pair<int, float>
      FourOnia2KKPAT::findJpsiMCInfo(reco::GenParticleRef genJpsi) {

        int momJpsiID = 0;
        float trueLife = -99.;

        if (genJpsi->numberOfMothers()>0) {

          TVector3 trueVtx(0.0,0.0,0.0);
          TVector3 trueP(0.0,0.0,0.0);
          TVector3 trueVtxMom(0.0,0.0,0.0);

          trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
          trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

          bool aBhadron = false;
          reco::GenParticleRef Jpsimom = genJpsi->motherRef();       // find mothers
          if (Jpsimom.isNull()) {
            std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
            return result;
          } else {
            reco::GenParticleRef Jpsigrandmom = Jpsimom->motherRef();
            if (isAbHadron(Jpsimom->pdgId())) {
              if (Jpsigrandmom.isNonnull() && isAMixedbHadron(Jpsimom->pdgId(),Jpsigrandmom->pdgId())) {
                momJpsiID = Jpsigrandmom->pdgId();
                trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
              } else {
                momJpsiID = Jpsimom->pdgId();
                trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
              }
              aBhadron = true;
            } else {
              if (Jpsigrandmom.isNonnull() && isAbHadron(Jpsigrandmom->pdgId())) {
                reco::GenParticleRef JpsiGrandgrandmom = Jpsigrandmom->motherRef();
                if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(Jpsigrandmom->pdgId(),JpsiGrandgrandmom->pdgId())) {
                  momJpsiID = JpsiGrandgrandmom->pdgId();
                  trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
                } else {
                  momJpsiID = Jpsigrandmom->pdgId();
                  trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
                }
                aBhadron = true;
              }
            }
            if (!aBhadron) {
              momJpsiID = Jpsimom->pdgId();
              trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
            }
          }

          TVector3 vdiff = trueVtx - trueVtxMom;
          //trueLife = vdiff.Perp()*3.09688/trueP.Perp();
          trueLife = vdiff.Perp()*genJpsi->mass()/trueP.Perp();
        }
        std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
        return result;

      }

      // ------------ method called once each job just before starting event loop  ------------
      void
      FourOnia2KKPAT::beginJob()
      {
      }

      // ------------ method called once each job just after ending the event loop  ------------
      void
      FourOnia2KKPAT::endJob() {
      }

      //define this as a plug-in
      DEFINE_FWK_MODULE(FourOnia2KKPAT);
