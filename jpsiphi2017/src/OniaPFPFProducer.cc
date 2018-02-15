#include "../interface/OniaPFPFProducer.h"

OniaPFPFProducer::OniaPFPFProducer(const edm::ParameterSet& ps):
  OniaCollection_(consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("Onia"))),
  PFCandCollection_(consumes<std::vector<pat::PackedCandidate>>(ps.getParameter<edm::InputTag>("PFCandidates"))),
  OniaMassCuts_(ps.getParameter<std::vector<double>>("OniaMassCuts")),
  TrakTrakMassCuts_(ps.getParameter<std::vector<double>>("TrakTrakMassCuts")),
  OniaPFPFMassCuts_(ps.getParameter<std::vector<double>>("OniaPFPFMassCuts")),
  MassTraks_(ps.getParameter<std::vector<double>>("MassTraks")),
  OnlyBest_(ps.getParameter<bool>("OnlyBest"))
{
  produces<pat::CompositeCandidateCollection>("OniaPFPFCandidates");
  candidates = 0;
  nevents = 0;
  nonia = 0;
  nreco = 0;
}

void OniaPFPFProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::unique_ptr<pat::CompositeCandidateCollection> OniaTTCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> onia;
  event.getByToken(OniaCollection_,onia);

  edm::Handle<std::vector<pat::PackedCandidate> > trak;
  event.getByToken(PFCandCollection_,trak);

  uint ncombo = 0;
  float OniaMassMax_ = OniaMassCuts_[1];
  float OniaMassMin_ = OniaMassCuts_[0];
  float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
  float TrakTrakMassMin_ = TrakTrakMassCuts_[0];
  float OniaPFPFMassMax_ = OniaPFPFMassCuts_[1];
  float OniaPFPFMassMin_ = OniaPFPFMassCuts_[0];

// Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
  for (pat::CompositeCandidateCollection::const_iterator oniaCand = onia->begin(); oniaCand != onia->end(); ++oniaCand){
     if ( oniaCand->mass() < OniaMassMax_  && oniaCand->mass() > OniaMassMin_ ) {
       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon1"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon2"));

// loop on track candidates, make OniaT candidate, positive charge
       for (std::vector<pat::PackedCandidate>::const_iterator posTrack = trak->begin(), trakend=trak->end(); posTrack!= trakend; ++posTrack){

         if(posTrack->charge()==0) continue;
         if(posTrack->pt()<0.5) continue;
	       if(fabs(posTrack->pdgId())!=211) continue;
	       if(!(posTrack->trackHighPurity())) continue;

         if ( IsTheSame(*posTrack,*pmu1) || IsTheSame(*posTrack,*pmu2) || posTrack->charge() < 0 ) continue;

// loop over second track candidate, negative charge
         for (std::vector<pat::PackedCandidate>::const_iterator negTrack = trak->begin(); negTrack!= trakend; ++negTrack){

           if(negTrack->charge()==0) continue;
           if(negTrack->pt()<0.5) continue;
  	       if(fabs(negTrack->pdgId())!=211) continue;
  	       if(!(negTrack->trackHighPurity())) continue;

           if (negTrack == posTrack) continue;
           if ( IsTheSame(*negTrack,*pmu1) || IsTheSame(*negTrack,*pmu2) || negTrack->charge() > 0 ) continue;

           pat::CompositeCandidate TTCand = makeTTCandidate(*posTrack, *negTrack);

           if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {

           pat::CompositeCandidate OniaTTCand = makeOniaTTCandidate(*oniaCand, *&TTCand);

           if ( OniaTTCand.mass() < OniaPFPFMassMax_ && OniaTTCand.mass() > OniaPFPFMassMin_) {

             OniaTTCandColl->push_back(OniaTTCand);
             candidates++;
             ncombo++;
           }
        }

         }
         } // loop over second track
       }   // loop on track candidates
       if (OnlyBest_) break;
     }

  if ( ncombo != OniaTTCandColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != OniaTT ("<<OniaTTCandColl->size()<<")"<< std::endl;
  if ( !onia->empty() )  nonia++;
  if ( ncombo > 0 ) nreco++;
  event.put(std::move(OniaTTCandColl),"OniaPFPFCandidates");
  nevents++;
}

void OniaPFPFProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "OniaPFPF Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with Onia candidates " << nonia << std::endl;
  std::cout << "Events with OniaPFPF candidates " << nreco << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " OniaPFPF candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

bool OniaPFPFProducer::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

const pat::CompositeCandidate OniaPFPFProducer::makeOniaTTCandidate(
                                          const pat::CompositeCandidate& onia,
				          const pat::CompositeCandidate& tt
                                         ){

  pat::CompositeCandidate OniaTCand;
  OniaTCand.addDaughter(onia,"onia");
  OniaTCand.addDaughter(tt,"ditrak");
  OniaTCand.setVertex(onia.vertex());
  OniaTCand.setCharge(tt.charge());

  reco::Candidate::LorentzVector vOniaT = onia.p4() + tt.p4();
  OniaTCand.setP4(vOniaT);

  return OniaTCand;

}

const pat::CompositeCandidate OniaPFPFProducer::makeTTCandidate(
                                          const pat::PackedCandidate& trakP,
                                          const pat::PackedCandidate& trakN
                                         ){

  pat::CompositeCandidate TTCand;
  TTCand.addDaughter(trakP,"trakP");
  TTCand.addDaughter(trakN,"trakN");
  TTCand.setCharge(trakP.charge()+trakN.charge());

  double m_kaon1 = MassTraks_[0];
  math::XYZVector mom_kaon1 = trakP.momentum();
  double e_kaon1 = sqrt(m_kaon1*m_kaon1 + mom_kaon1.Mag2());
  math::XYZTLorentzVector p4_kaon1 = math::XYZTLorentzVector(mom_kaon1.X(),mom_kaon1.Y(),mom_kaon1.Z(),e_kaon1);
  double m_kaon2 = MassTraks_[1];
  math::XYZVector mom_kaon2 = trakN.momentum();
  double e_kaon2 = sqrt(m_kaon2*m_kaon2 + mom_kaon2.Mag2());
  math::XYZTLorentzVector p4_kaon2 = math::XYZTLorentzVector(mom_kaon2.X(),mom_kaon2.Y(),mom_kaon2.Z(),e_kaon2);
  reco::Candidate::LorentzVector vTT = p4_kaon1 + p4_kaon2;
  TTCand.setP4(vTT);

  return TTCand;
}


reco::Candidate::LorentzVector OniaPFPFProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(OniaPFPFProducer);
