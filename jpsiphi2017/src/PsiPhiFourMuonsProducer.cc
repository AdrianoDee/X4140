#include "../interface/PsiPhiFourMuonsProducer.h"

PsiPhiFourMuonsProducer::PsiPhiFourMuonsProducer(const edm::ParameterSet& ps):
  PsiCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("JPsiCollection"))),
  PhiCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("PhiCollection"))),
  JPsiMassCuts_(ps.getParameter<std::vector<double>>("JPsiMassCuts")),
  PhiMassCuts_(ps.getParameter<std::vector<double>>("PhiMassCuts")),
  FourOniaMassCuts_(ps.getParameter<std::vector<double>>("FourOniaMassCuts"))
{
  produces<pat::CompositeCandidateCollection>("PsiPhiFourMuonsCandidates");
  candidates = 0;
  nevents = 0;
  nPhi = 0;
  nJps = 0;
}

void PsiPhiFourMuonsProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::unique_ptr<pat::CompositeCandidateCollection> FourMuCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> psiOnia;
  event.getByToken(PsiCollection_,psiOnia);

  edm::Handle<pat::CompositeCandidateCollection> phiOnia;
  event.getByToken(PhiCollection_,phiOnia);

  float JPsiMassMax_ = JPsiMassCuts_[1];
  float JPsiMassMin_ = JPsiMassCuts_[0];
  float PhiMassMax_ = PhiMassCuts_[1];
  float PhiMassMin_ = PhiMassCuts_[0];

  float FourOniaMassMax_ = FourOniaMassCuts_[1];
  float FourOniaMassMin_ = FourOniaMassCuts_[0];

// Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
  //Looking for J/Psi
  for (pat::CompositeCandidateCollection::const_iterator jPsiCand = psiOnia->begin(); jPsiCand != psiOnia->end(); ++jPsiCand){
     if ( jPsiCand->mass() < JPsiMassMax_  && jPsiCand->mass() > JPsiMassMin_ ) {
       const pat::Muon *jPsiMu1 = dynamic_cast<const pat::Muon*>(jPsiCand->daughter("muon1"));
       const pat::Muon *jPsiMu2 = dynamic_cast<const pat::Muon*>(jPsiCand->daughter("muon2"));


       for (pat::CompositeCandidateCollection::const_iterator phiCand = jPsiCand + 1; phiCand != phiOnia->end(); ++phiCand){
          if ( phiCand->mass() < PhiMassMax_  && phiCand->mass() > PhiMassMin_ ) {
            const pat::Muon *phiMu1 = dynamic_cast<const pat::Muon*>(phiCand->daughter("muon1"));
            const pat::Muon *phiMu2 = dynamic_cast<const pat::Muon*>(phiCand->daughter("muon2"));

            if( phiMu1 == phiMu2 || phiMu1 == jPsiMu1 || phiMu1 == jPsiMu2 ) continue;
            if( phiMu2 == jPsiMu1 || phiMu2 == jPsiMu2 ) continue;
            if( jPsiMu1 == jPsiMu2 ) continue;

            pat::CompositeCandidate fourOniaCandidate = makeCandidate(*phiCand, *jPsiCand);

            if(fourOniaCandidate.charge() != 0.0) continue;

            if ( fourOniaCandidate.mass() < FourOniaMassMax_ && fourOniaCandidate.mass() > FourOniaMassMin_)
              {
                candidates++;
                FourMuCandColl->push_back(fourOniaCandidate);
              }
            }
          }
        }
      }
     // if (OnlyBest_) break;
     if ( !(psiOnia->empty()) )  nPhi++;
     if ( !(phiOnia->empty()) )  nJps++;

     event.put(std::move(FourMuCandColl),"FourOniaCandidates");
     nevents++;
  }


void PsiPhiFourMuonsProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "OniaTrakTrak Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with Phi  candidates " << nPhi << std::endl;
  std::cout << "Events with JPsi candidates " << nJps << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " OniaTrakTrak candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}


const pat::CompositeCandidate PsiPhiFourMuonsProducer::makeCandidate(const pat::CompositeCandidate& phi,
  const pat::CompositeCandidate& jpsi){
    pat::CompositeCandidate xCand;
    xCand.addDaughter(phi,"phi");
    xCand.addDaughter(jpsi,"jpsi");
    reco::Candidate::LorentzVector vX = phi.p4() + jpsi.p4();
    xCand.setP4(vX);
    return xCand;
  }

reco::Candidate::LorentzVector PsiPhiFourMuonsProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(PsiPhiFourMuonsProducer);
