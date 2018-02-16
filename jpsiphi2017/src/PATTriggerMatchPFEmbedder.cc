

// system include files
#include <memory>

// FW include files
#include "mmkk/mmkk/PATTriggerMatchEmbedder.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "PhysicsTools/PatAlgos/plugins/PATTriggerMatchEmbedder.cc"

namespace pat
{
  typedef PATTriggerMatchEmbedder< CompositeCandidate > PATTriggerMatchPFEmbedder;
}

DEFINE_FWK_MODULE( PATTriggerMatchPFEmbedder );
