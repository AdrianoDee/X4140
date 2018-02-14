/**
   \file
   Declaration of DiMuonTracksProducer
   \author
   Alberto Sanchez-Hernandez
   September 2014

*/

#ifndef __DiMuonTracksProducer_h_
#define __DiMuonTracksProducer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <TLorentzVector.h>
#include <vector>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>

#include "FourOniaVtxReProducer.h"

template<typename T>
struct GreaterByVProb {
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};

class DiMuonTracksProducer : public edm::EDProducer {

 public:
  explicit DiMuonTracksProducer(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;
  void endJob() override;



  edm::EDGetTokenT<pat::CompositeCandidateCollection> phi_tracks_Label;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> jpsi_dimuon_Label;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  StringCutObjectSelector<reco::Candidate, true> quadmuonSelection_;

  double dzMax_;
  bool addMuonlessPrimaryVertex_;
  bool resolveAmbiguity_;
  bool addMCTruth_;
  bool triggerMatch_;


  const pat::CompositeCandidate makeCandidate(const pat::CompositeCandidate&,
						 const pat::CompositeCandidate&);

  float Getdz(const pat::CompositeCandidate&, const reco::Candidate::Point &);
  // check if the mass difference is in desired range
  bool cutDeltaMass(const pat::CompositeCandidate&,const pat::CompositeCandidate&);

  bool cutdz(float dz){return dz<dzMax_; }

  bool isOverlappedMuons(const pat::CompositeCandidate *phi,const pat::CompositeCandidate *jpsi);

  edm::EDGetTokenT<reco::TrackCollection> revtxtrks_;
  edm::EDGetTokenT<reco::BeamSpot> revtxbs_;

  int candidates;
  int delta_mass_fail;
  int dz_phi_cut_fail;
  int dz_jps_cut_fail;

  GreaterByVProb<pat::CompositeCandidate> vPComparator_;

  InvariantMassFromVertex massCalculator;

};

#endif // __DiMuonTracksProducer_h_
