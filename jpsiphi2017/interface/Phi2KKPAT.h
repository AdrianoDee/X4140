#ifndef HeavyFlavorAnalysis_Onia2MuMu_Phi2KKPAT_h
#define HeavyFlavorAnalysis_Onia2MuMu_Phi2KKPAT_h

// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
//#include <DataFormats/PatCandidates/interface/Muon.h>

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

template<typename T>
struct GreaterByVProb {
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};

template<typename T>
struct EqualByPosition {
  bool operator() ( T  v1, T  v2 ) {
    return (v1.x() == v2.x() && v1.y() == v2.y() && v1.z() == v2.z());
  }
};


//
// class decleration
//

class Phi2KKPAT : public edm::EDProducer {
 public:
  explicit Phi2KKPAT(const edm::ParameterSet&);
  ~Phi2KKPAT() override {};
 private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  int GetPVFromOnia(edm::Event&, edm::Handle<reco::VertexCollection>);

  // ----------member data ---------------------------
 private:
  edm::EDGetTokenT<std::vector<pat::GenericParticle>> kaons_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> theOnias_;
  StringCutObjectSelector<pat::GenericParticle> higherPuritySelection_;
  StringCutObjectSelector<pat::GenericParticle> lowerPuritySelection_; 
  StringCutObjectSelector<reco::Candidate, true> dikaonSelection_;
  bool addCommonVertex_;
  bool resolveAmbiguity_;

  GreaterByVProb<pat::CompositeCandidate> vPComparator_;
  InvariantMassFromVertex massCalculator;
  EqualByPosition<reco::Vertex> vertexComparator_;
};

#endif
