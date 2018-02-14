#include "../interface/Phi2KKPAT.h"

#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/PatCandidates/interface/UserData.h>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

Phi2KKPAT::Phi2KKPAT(const edm::ParameterSet& iConfig):
  kaons_(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("kaons"))),
  thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
  theOnias_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("OniaTag"))),
  higherPuritySelection_(iConfig.getParameter<std::string>("higherPuritySelection")),
  lowerPuritySelection_(iConfig.getParameter<std::string>("lowerPuritySelection")),
  dikaonSelection_(iConfig.existsAs<std::string>("dikaonSelection") ? iConfig.getParameter<std::string>("dikaonSelection") : ""),
  addCommonVertex_(iConfig.getParameter<bool>("addCommonVertex")),
  resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity"))
{
    produces<pat::CompositeCandidateCollection>();
}

int Phi2KKPAT::GetPVFromOnia(edm::Event& iEvent, edm::Handle<reco::VertexCollection> priVtxs) {
  int ivertex = -1;
  edm::Handle<pat::CompositeCandidateCollection> onias_;
  iEvent.getByToken(theOnias_, onias_);
  //std::cout << "Phi2KKPAT::GetPVFromOnia: onia size " << onias_->size() << " pv size " << priVtxs->size() << std::endl;
  if (onias_.isValid() && !onias_->empty() && !priVtxs->empty()) {
    const pat::CompositeCandidate *ionia =  &(onias_->at(0));
    if (ionia) {
      const reco::Vertex *ipv = ionia->userData<reco::Vertex>("PVwithmuons");
      ivertex = 0;               // if not vertex available in onia the defakt one is used
      if (ipv) {
        //std::cout << "Phi2KKPAT::GetPVFromOnia: finding pv, from PVwithmuons " << ipv->x() << " " << ipv->y() << " " << ipv->z() << std::endl;
        int index_v = -1;
        for (std::vector<reco::Vertex>::const_iterator iv = priVtxs->begin(), ivend = priVtxs->end(); iv != ivend; ++iv) {
          index_v++;
          if (vertexComparator_(*iv,*ipv)) {
            //std::cout << index_v << ": " << iv->x() << " " << iv->y() << " " << iv->z() << std::endl;
	    ivertex = index_v;
	    break;
          }
        }
	if (index_v < 0) std::cout << "Phi2KKPAT::GetPVFromOnia: *** non matching PV to Onia PV found in PVCollection, using iPV=0" << std::endl;
      } else std::cout << "Phi2KKPAT::GetPVFromOnia: *** no PV associate to Onia object found in composite candidate, using iPV=0" << std::endl;
    } //else std::cout << "Phi2KKPAT::GetPVFromOnia: *** no Onia combination pass selection" << std::endl;
  } else { if (!priVtxs->empty()) ivertex = 0; }  // in case no onia is present we will use the PV==0, so a filter has to be applied befor.
  return ivertex;
}

// ------------ method called to produce the data  ------------
void Phi2KKPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  std::vector<double> kMasses;
  kMasses.push_back( 0.493677 );
  kMasses.push_back( 0.493677 );

  std::unique_ptr<pat::CompositeCandidateCollection> phiOutput(new pat::CompositeCandidateCollection);

  reco::Vertex thePrimaryV;
  reco::Vertex theBeamSpotV;

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  edm::Handle<BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  edm::Handle<VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);

  if ( priVtxs->begin() != priVtxs->end() ) thePrimaryV = Vertex(*(priVtxs->begin()));
  else thePrimaryV = Vertex(bs.position(), bs.covariance3D());

  edm::Handle<pat::GenericParticleCollection> kaons;
  iEvent.getByToken(kaons_,kaons);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  int npv = priVtxs->size();
  int which_vertex = GetPVFromOnia(iEvent,priVtxs);
  if (which_vertex < 0)   std::cout << "Phi2KKPAT::GetPVFromOnia: *** Event does not have enough data to continue == -1" << endl;
  if (which_vertex > npv) std::cout << "Phi2KKPAT::GetPVFromOnia: *** Onia vertex index above bondaries " << which_vertex << " " << npv << endl;

  // filter the tracks ..

  int ntracks_pv = kaons->size();
  int npvtracks = 0;

  // phi candidates only from kaons
  for (size_t ii=0; ii<(size_t)ntracks_pv; ii++) {
    const pat::GenericParticle *it = &(kaons->at(ii));
    // both must pass low quality
    if(!lowerPuritySelection_(*it)) continue;
    for (size_t jj=ii+1; jj<(size_t)ntracks_pv; jj++) {
      const pat::GenericParticle *it2 = &(kaons->at(jj));

      int pv_index = -1;
      // both must pass low quality
      if(!lowerPuritySelection_(*it2)) continue;
      // one must pass tight quality
      if (!(higherPuritySelection_(*it) || higherPuritySelection_(*it2))) continue;

      float deltaphi  = it->p4().Phi() - it2->p4().Phi();
      while (deltaphi < -M_PI) deltaphi += 2*M_PI;
      while (deltaphi >  M_PI) deltaphi -= 2*M_PI;
      float deltaeta  = it->p4().Eta() - it2->p4().Eta();
      float deltar    = sqrt(pow(deltaphi,2) + pow(deltaeta,2));

      pat::CompositeCandidate myPhi;
      vector<TransientVertex> pvs;

      // ---- no explicit order defined ----
      myPhi.addDaughter(*it, "kaon1");
      myPhi.addDaughter(*it2,"kaon2");

      // ---- define and set candidate's 4momentum  ----
      LorentzVector phi = it->p4() + it2->p4();
      myPhi.setP4(phi);
      myPhi.setCharge(it->charge()+it2->charge());
      // ----
      myPhi.addUserFloat("deltar",deltar);

      // ---- apply the dikaon cut ----
      if (!dikaonSelection_(myPhi)) continue;

      if (!(it->track().isNonnull() && it2->track().isNonnull())) continue;

      // ---- fit vertex using Tracker tracks (if they have tracks) ----
      if (it->track().isNonnull() && it2->track().isNonnull()) {

	//build the dikaon secondary vertex
	vector<TransientTrack> t_tks;
	t_tks.push_back(theTTBuilder->build(it->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
	t_tks.push_back(theTTBuilder->build(it2->track())); // otherwise the vertex will have transient refs inside.
	TransientVertex myVertex = vtxFitter.vertex(t_tks);

	CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );

        Measurement1D MassWErr(phi.M(),-9999.);
        if ( field->nominalValue() > 0 ) MassWErr = massCalculator.invariantMass( VtxForInvMass, kMasses );
        else myVertex = TransientVertex();                   // this is an invalid vertex by definition

	myPhi.addUserFloat("MassErr",MassWErr.error());

	if (myVertex.isValid()) {
	  float vChi2 = myVertex.totalChiSquared();
	  float vNDF  = myVertex.degreesOfFreedom();
	  float vProb(TMath::Prob(vChi2,(int)vNDF));

	  myPhi.addUserFloat("vNChi2",vChi2/vNDF);
	  myPhi.addUserFloat("vProb",vProb);

	  TVector3 vtx;
          TVector3 pvtx;
          VertexDistanceXY vdistXY;

	  vtx.SetXYZ(myVertex.position().x(),myVertex.position().y(),0);
	  TVector3 pperp(phi.px(), phi.py(), 0);
	  AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

	  float minDz = 99999.;
	  float extrapZ=-9E20;

	  if (resolveAmbiguity_) {
	    TwoTrackMinimumDistance ttmd;
	    bool status = ttmd.calculate( GlobalTrajectoryParameters(
                                                                     GlobalPoint(myVertex.position().x(), myVertex.position().y(), myVertex.position().z()),
                                                                     GlobalVector(myPhi.px(),myPhi.py(),myPhi.pz()),TrackCharge(0),&(*magneticField)),
					  GlobalTrajectoryParameters(
								     GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
								     GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));
	    if (status) extrapZ=ttmd.points().first.z();

	      int ii_pv = -1;
	      for (VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv) {
		ii_pv++;
		float deltaZ = fabs(extrapZ - itv->position().z()) ;
		if ( deltaZ < minDz ) {
		  minDz = deltaZ;
		  thePrimaryV = Vertex(*itv);
		  pv_index = ii_pv;
		}
	      }
	  } else {
            minDz = -1;
            pv_index = which_vertex;
            thePrimaryV = (*priVtxs)[which_vertex];
            extrapZ = thePrimaryV.position().z();
          }

          myPhi.addUserInt("oniaPV",which_vertex);
	  myPhi.addUserInt("iPV",pv_index);
	  myPhi.addUserFloat("dzPV",minDz);
	  myPhi.addUserFloat("extrapZPV",extrapZ);

	  // count the number of high Purity tracks with pT > 500 MeV attached to the chosen vertex
	  double vertexWeight = 0., sumPTPV = 0.;
	  int countTksOfPV = 0;
          for (size_t kk=1; kk<(size_t)ntracks_pv; kk++) {
              const pat::GenericParticle *it3 = &(kaons->at(kk));
              if (!it3->track().isNonnull())                  continue;
              reco::Track track = *it3->track();
              if(track.pt() < 0.5)                            continue;
              if(!track.quality(reco::TrackBase::highPurity)) continue;
              TransientTrack tt = theTTBuilder->build(track);
              pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,thePrimaryV);
              if (!tkPVdist.first)                  continue;
              if (tkPVdist.second.significance()>3) continue;
              if (track.ptError()/track.pt()>0.1)   continue;
              if (it3 == it2 || it3 == it)          continue;
              countTksOfPV++;
              sumPTPV += track.pt();
              vertexWeight += thePrimaryV.trackWeight(it3->track());
          }

	  myPhi.addUserInt("countTksOfPV", countTksOfPV);
	  myPhi.addUserFloat("vertexWeight", (float) vertexWeight);
	  myPhi.addUserFloat("sumPTPV", (float) sumPTPV);

	  ///DCA
	  TrajectoryStateClosestToPoint k1TS = t_tks[0].impactPointTSCP();
	  TrajectoryStateClosestToPoint k2TS = t_tks[1].impactPointTSCP();
	  float dca = 1E20;
	  if (k1TS.isValid() && k2TS.isValid()) {
	    ClosestApproachInRPhi cApp;
	    cApp.calculate(k1TS.theState(), k2TS.theState());
	    if (cApp.status() ) dca = cApp.distance();
	  }
	  myPhi.addUserFloat("DCA", dca );
	  ///end DCA

	  myPhi.addUserData("PVwithkaons",thePrimaryV);
	  npvtracks = thePrimaryV.nTracks();

	  // lifetime using PV
          pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
	  TVector3 vdiff = vtx - pvtx;
	  double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
	  Measurement1D distXY = vdistXY.distance(Vertex(myVertex), thePrimaryV);
	  double ctauPV = distXY.value()*cosAlpha * myPhi.mass()/pperp.Perp();
	  GlobalError v1e = (Vertex(myVertex)).error();
	  GlobalError v2e = thePrimaryV.error();
          AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
	  double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*myPhi.mass()/(pperp.Perp2());

	  myPhi.addUserFloat("ppdlPV",ctauPV);
          myPhi.addUserFloat("ppdlErrPV",ctauErrPV);
	  myPhi.addUserFloat("cosAlpha",cosAlpha);

	  // lifetime using BS
          pvtx.SetXYZ(theBeamSpotV.position().x(),theBeamSpotV.position().y(),0);
	  vdiff = vtx - pvtx;
	  cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
	  distXY = vdistXY.distance(Vertex(myVertex), theBeamSpotV);
	  double ctauBS = distXY.value()*cosAlpha*myPhi.mass()/pperp.Perp();
	  GlobalError v1eB = (Vertex(myVertex)).error();
	  GlobalError v2eB = theBeamSpotV.error();
          AlgebraicSymMatrix33 vXYeB = v1eB.matrix()+ v2eB.matrix();
	  double ctauErrBS = sqrt(ROOT::Math::Similarity(vpperp,vXYeB))*myPhi.mass()/(pperp.Perp2());

	  myPhi.addUserFloat("ppdlBS",ctauBS);
          myPhi.addUserFloat("ppdlErrBS",ctauErrBS);

	  if (addCommonVertex_) myPhi.addUserData("commonVertex",Vertex(myVertex));
          myPhi.addUserInt("npvtracks", npvtracks);
	  myPhi.addUserInt("ntracks_pv",ntracks_pv );

          // ---- If here push back to output ----
          phiOutput->push_back(myPhi);

	} // vtxfit of two tracks is valid
      }   // both trks have track
    }     // loop trk2
  }       // loop trk1

  if (phiOutput->size() > 1) std::sort(phiOutput->begin(),phiOutput->end(),vPComparator_);
  iEvent.put(std::move(phiOutput));

}

//define this as a plug-in
DEFINE_FWK_MODULE(Phi2KKPAT);
