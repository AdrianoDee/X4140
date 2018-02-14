#include "../interface/FourOniaProducer.h"
#include "../interface/FourOniaVtxReProducer.h"


FourOniaProducer::FourOniaProducer(const edm::ParameterSet& iConfig):
phi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("phidimuons"))),
jpsi_dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("jpsidimuons"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
quadmuonSelection_(iConfig.existsAs<std::string>("quadmuonSelection") ? iConfig.getParameter<std::string>("quadmuonSelection") : ""),
dzMax_(iConfig.getParameter<double>("dzmax")),
addMuonlessPrimaryVertex_(iConfig.getParameter<bool>("addMuonlessPrimaryVertex")),
resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity")),
// addMCTruth_(iConfig.getParameter<bool>("addMCTruth")),
triggerMatch_(iConfig.getParameter<bool>("triggerMatch"))
//doCombinatorial(iConfig.existsAs<bool>("doCombinatorial") ? iConfig.getParameter<bool>("doCombinatorial") : false)
{
  revtxtrks_ = consumes<reco::TrackCollection>((edm::InputTag)"generalTracks"); //if that is not true, we will raise an exception
  revtxbs_ = consumes<reco::BeamSpot>((edm::InputTag)"offlineBeamSpot");
  produces<pat::CompositeCandidateCollection>();
  candidates = 0;
  delta_mass_fail = 0;
  dz_phi_cut_fail = 0;
  dz_jps_cut_fail = 0;
}


void FourOniaProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  using namespace edm;
  using namespace std;
  using namespace reco;

  Vertex thePrimaryV;
  Vertex theBeamSpotV;

  vector<double> muMasses;
  muMasses.push_back( 0.1056583715 );
  muMasses.push_back( 0.1056583715 );
  muMasses.push_back( 0.1056583715 );
  muMasses.push_back( 0.1056583715 );

  //std::cout<<"FourOniaProducer inside"<<std::endl;

  std::unique_ptr<pat::CompositeCandidateCollection> xCandColl(new pat::CompositeCandidateCollection);

  ESHandle<MagneticField> magneticField;
  esetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  edm::Handle<pat::CompositeCandidateCollection> dimuonsPhi;
  event.getByToken(phi_dimuon_Label,dimuonsPhi);

  edm::Handle<pat::CompositeCandidateCollection> dimuonsJPsi;
  event.getByToken(jpsi_dimuon_Label,dimuonsJPsi);

  Handle<BeamSpot> theBeamSpot;
  event.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  esetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);
  TrackCollection muonLess;

  Handle<VertexCollection> priVtxs;
  event.getByToken(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
  }

  //std::cout << "JPsis size: " << dimuonsJPsi.product()->size() << std::endl;
  //std::cout << "Phis size: " << dimuonsPhi.product()->size() << std::endl;
  // Note: since Dimuon cand are sorted by decreasing vertex probability then the first chi cand is the one associated with the "best" dimuon
  for (pat::CompositeCandidateCollection::const_iterator  phiCand = dimuonsPhi->begin(); phiCand!= dimuonsPhi->end(); ++phiCand){

    // use only trigger-matched Jpsi or Upsilon if so requested

    if ( ( triggerMatch_ ) && ( !phiCand->userInt("isTriggerMatched") ) )
    continue;
    //std::cout << "Phi muons trigger matched" << std::endl;
    for (pat::CompositeCandidateCollection::const_iterator  jpsiCand = dimuonsJPsi->begin(); jpsiCand!= dimuonsJPsi->end(); ++jpsiCand){

      if (( triggerMatch_ ) && (  !jpsiCand->userInt("isTriggerMatched") ) )
      continue;
      //std::cout << "Jps muons trigger matched" << std::endl;
      if(areOverlappedMuons(&(*phiCand),&(*jpsiCand)))
      continue;

      int pMatch = 0, jMatch = 0;
      float pDeltaR = -1.0, jDeltaR = -1.0;

      pMatch = phiCand->userInt("isTriggerMatched");
      jMatch = jpsiCand->userInt("isTriggerMatched");

      // std::cout << phiCand->userInt("isTriggerMatched") << std::endl;

      pat::CompositeCandidate xCand;

      xCand = makeCandidate(*phiCand, *jpsiCand);

      if(!quadmuonSelection_(xCand)) continue;

      if (((dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon1") )->track()).isNonnull())
      && ((dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon2") )->track()).isNonnull())
      && ((dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon1") )->track()).isNonnull())
      && ((dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon2") )->track()).isNonnull())){

        xCand.addUserInt("phi_isTriggerMatched",pMatch);
        xCand.addUserInt("jpsi_isTriggerMatched",jMatch);

        jDeltaR = reco::deltaR2((dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon1") )->track())->eta(),
            (dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon1") )->track())->phi(),
            (dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon2") )->track())->eta(),
            (dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon2") )->track())->phi());

        pDeltaR = reco::deltaR2((dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon1") )->track())->eta(),
            (dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon1") )->track())->phi(),
            (dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon2") )->track())->eta(),
            (dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon2") )->track())->phi());

        xCand.addUserFloat("phi_deltaR",pDeltaR);
        xCand.addUserFloat("jpsi_deltaR",jDeltaR);

        vector<TransientVertex> pvs;

        vector<TransientTrack> t_tks;
        t_tks.push_back(theTTBuilder->build(dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon1") )->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
        t_tks.push_back(theTTBuilder->build(dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon2") )->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
        t_tks.push_back(theTTBuilder->build(dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon1") )->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
        t_tks.push_back(theTTBuilder->build(dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon2") )->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
        TransientVertex myVertex = vtxFitter.vertex(t_tks);

        CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );

        Measurement1D MassWErr(xCand.mass(),-9999.);
        if ( field->nominalValue() > 0 ) {
          MassWErr = massCalculator.invariantMass( VtxForInvMass, muMasses );
        } else {
          myVertex = TransientVertex();                      // with no arguments it is invalid
        }

        xCand.addUserFloat("MassErr",MassWErr.error());

        if (myVertex.isValid())
        {
          float vChi2 = myVertex.totalChiSquared();
          float vNDF  = myVertex.degreesOfFreedom();
          float vProb(TMath::Prob(vChi2,(int)vNDF));

          xCand.addUserFloat("vNChi2",vChi2/vNDF);
          xCand.addUserFloat("vProb",vProb);

          TVector3 vtx, vtx3D;
          TVector3 pvtx, pvtx3D;
          VertexDistanceXY vdistXY;

          vtx.SetXYZ(myVertex.position().x(),myVertex.position().y(),0);
          TVector3 pperp(xCand.px(), xCand.py(), 0);
          TVector3 pperp3D(xCand.px(), xCand.py(), xCand.pz());
          AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
          AlgebraicVector3 vpperp3D(pperp.x(),pperp.y(),pperp.z());

          if (resolveAmbiguity_) {

            float minDz = 999999.;
            TwoTrackMinimumDistance ttmd;
            bool status = ttmd.calculate(
              GlobalTrajectoryParameters(GlobalPoint(myVertex.position().x(), myVertex.position().y(), myVertex.position().z()),
              GlobalVector(xCand.px(),xCand.py(),xCand.pz()),TrackCharge(0),&(*magneticField)),
              GlobalTrajectoryParameters(GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
              GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));

              float extrapZ=-9E20;
              if (status) extrapZ=ttmd.points().first.z();

              for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv)
              {
                float deltaZ = fabs(extrapZ - itv->position().z()) ;
                if ( deltaZ < minDz ) {
                  minDz = deltaZ;
                  thePrimaryV = Vertex(*itv);
                }
              }
            }


            Vertex theOriginalPV = thePrimaryV;

            muonLess.clear();
            muonLess.reserve(thePrimaryV.tracksSize());
            if( addMuonlessPrimaryVertex_  && thePrimaryV.tracksSize()>2) {
              // Primary vertex matched to the dimuon, now refit it removing the two muons
              FourOniaVtxReProducer revertex(priVtxs, event);
              edm::Handle<reco::TrackCollection> pvtracks;
              event.getByToken(revtxtrks_,   pvtracks);
              if( !pvtracks.isValid()) { std::cout << "pvtracks NOT valid " << std::endl; }
              else {
                edm::Handle<reco::BeamSpot> pvbeamspot;
                event.getByToken(revtxbs_, pvbeamspot);
                if (pvbeamspot.id() != theBeamSpot.id()) edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
                // I need to go back to the reco::Muon object, as the TrackRef in the pat::Muon can be an embedded ref.
                const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon1") )->originalObject());
                const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon2") )->originalObject());
                const reco::Muon *rmu3 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon1") )->originalObject());
                const reco::Muon *rmu4 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon2") )->originalObject());
                // check that muons are truly from reco::Muons (and not, e.g., from PF Muons)
                // also check that the tracks really come from the track collection used for the BS
                if (rmu1 != nullptr && rmu2 != nullptr && rmu1->track().id() == pvtracks.id() && rmu2->track().id() == pvtracks.id() &&
                rmu3 != nullptr && rmu4 != nullptr && rmu3->track().id() == pvtracks.id() && rmu4->track().id() == pvtracks.id() ) {
                  // Save the keys of the tracks in the primary vertex
                  // std::vector<size_t> vertexTracksKeys;
                  // vertexTracksKeys.reserve(thePrimaryV.tracksSize());
                  if( thePrimaryV.hasRefittedTracks() ) {
                    // Need to go back to the original tracks before taking the key
                    std::vector<reco::Track>::const_iterator itRefittedTrack = thePrimaryV.refittedTracks().begin();
                    std::vector<reco::Track>::const_iterator refittedTracksEnd = thePrimaryV.refittedTracks().end();
                    for( ; itRefittedTrack != refittedTracksEnd; ++itRefittedTrack ) {
                      if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu1->track().key() ) continue;
                      if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu2->track().key() ) continue;
                      if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu3->track().key() ) continue;
                      if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu4->track().key() ) continue;
                      // vertexTracksKeys.push_back(thePrimaryV.originalTrack(*itRefittedTrack).key());
                      muonLess.push_back(*(thePrimaryV.originalTrack(*itRefittedTrack)));
                    }
                  }
                  else {
                    std::vector<reco::TrackBaseRef>::const_iterator itPVtrack = thePrimaryV.tracks_begin();
                    for( ; itPVtrack != thePrimaryV.tracks_end(); ++itPVtrack ) if (itPVtrack->isNonnull()) {
                      if( itPVtrack->key() == rmu1->track().key() ) continue;
                      if( itPVtrack->key() == rmu2->track().key() ) continue;
                      if( itPVtrack->key() == rmu3->track().key() ) continue;
                      if( itPVtrack->key() == rmu4->track().key() ) continue;
                      // vertexTracksKeys.push_back(itPVtrack->key());
                      muonLess.push_back(**itPVtrack);
                    }
                  }
                  if (muonLess.size()>1 && muonLess.size() < thePrimaryV.tracksSize()){
                    pvs = revertex.makeVertices(muonLess, *pvbeamspot, esetup) ;
                    if (!pvs.empty()) {
                      Vertex muonLessPV = Vertex(pvs.front());
                      thePrimaryV = muonLessPV;
                    }
                  }
                }
              }
            }

            // count the number of high Purity tracks with pT > 400 MeV attached to the chosen vertex
            double vertexWeight = -1., sumPTPV = -1.;
            int countTksOfPV = -1;
            const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon1"))->originalObject());
            const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*jpsiCand).daughter("muon2"))->originalObject());
            const reco::Muon *rmu3 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon1"))->originalObject());
            const reco::Muon *rmu4 = dynamic_cast<const reco::Muon *>(dynamic_cast<const pat::Muon*>((*phiCand).daughter("muon2"))->originalObject());
            try{
              for(reco::Vertex::trackRef_iterator itVtx = theOriginalPV.tracks_begin(); itVtx != theOriginalPV.tracks_end(); itVtx++) if(itVtx->isNonnull()){
                const reco::Track& track = **itVtx;
                if(!track.quality(reco::TrackBase::highPurity)) continue;
                if(track.pt() < 0.4) continue; //reject all rejects from counting if less than 400 MeV
                TransientTrack tt = theTTBuilder->build(track);
                pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,thePrimaryV);
                if (!tkPVdist.first) continue;
                if (tkPVdist.second.significance()>3) continue;
                if (track.ptError()/track.pt()>0.1) continue;
                // do not count the two muons
                if (rmu1 != nullptr && rmu1->innerTrack().key() == itVtx->key())
                continue;
                if (rmu2 != nullptr && rmu2->innerTrack().key() == itVtx->key())
                continue;

                if (rmu3 != nullptr && rmu3->innerTrack().key() == itVtx->key())
                continue;
                if (rmu4 != nullptr && rmu4->innerTrack().key() == itVtx->key())
                continue;

                vertexWeight += theOriginalPV.trackWeight(*itVtx);
                if(theOriginalPV.trackWeight(*itVtx) > 0.5){
                  countTksOfPV++;
                  sumPTPV += track.pt();
                }
              }
            } catch (std::exception & err) {std::cout << " muon Selection%G�%@failed " << std::endl; return ; }

            xCand.addUserInt("countTksOfPV", countTksOfPV);
            xCand.addUserFloat("vertexWeight", (float) vertexWeight);
            xCand.addUserFloat("sumPTPV", (float) sumPTPV);

            if (addMuonlessPrimaryVertex_)
            xCand.addUserData("muonlessPV",Vertex(thePrimaryV));

            xCand.addUserData("PVwithmuons",Vertex(theOriginalPV));

            // lifetime using PV

            float cosAlpha, cosAlpha3D, ctauPV, ctauErrPV, l_xyz, l_xy, lErr_xyz, lErr_xy;

            pvtx.SetXYZ(theOriginalPV.position().x(),theOriginalPV.position().y(),0);
            TVector3 vdiff = vtx - pvtx;
            cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());

            Measurement1D distXY = vdistXY.distance(Vertex(myVertex), theOriginalPV);
            //double ctauPV = distXY.value()*cosAlpha*3.09688/pperp.Perp();
            ctauPV = distXY.value()*cosAlpha * xCand.mass()/pperp.Perp();

            GlobalError v1e = (Vertex(myVertex)).error();
            GlobalError v2e = theOriginalPV.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
            //double ctauErrPV = sqrt(vXYe.similarity(vpperp))*3.09688/(pperp.Perp2());
            ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*xCand.mass()/(pperp.Perp2());

            AlgebraicVector3 vDiff;
            vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0 ;
            l_xy = vdiff.Perp();
            lErr_xy = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();

            /// 3D
            pvtx3D.SetXYZ(thePrimaryV.position().x(), thePrimaryV.position().y(), thePrimaryV.position().z());
            TVector3 vdiff3D = vtx3D - pvtx3D;
            cosAlpha3D = vdiff3D.Dot(pperp3D)/(vdiff3D.Mag()*vdiff3D.Mag());
            l_xyz = vdiff3D.Mag();

            AlgebraicVector3 vDiff3D;
            vDiff3D[0] = vdiff3D.x(); vDiff3D[1] = vdiff3D.y(); vDiff3D[2] = vdiff3D.z() ;
            lErr_xyz = sqrt(ROOT::Math::Similarity(vDiff3D,vXYe)) / vdiff3D.Mag();

            xCand.addUserFloat("ctauPV",ctauPV);
            xCand.addUserFloat("ctauErrPV",ctauErrPV);
            xCand.addUserFloat("cosAlpha",cosAlpha);
            xCand.addUserFloat("cosAlpha3D",cosAlpha3D);

            xCand.addUserFloat("l_xy",l_xy);
            xCand.addUserFloat("lErr_xy",lErr_xy);

            xCand.addUserFloat("l_xyz",l_xyz);
            xCand.addUserFloat("lErr_xyz",lErr_xyz);

            if (addMuonlessPrimaryVertex_){
              pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
              TVector3 vdiff = vtx - pvtx;
              double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
              Measurement1D distXY = vdistXY.distance(Vertex(myVertex), thePrimaryV);
              //double ctauPV = distXY.value()*cosAlpha*3.09688/pperp.Perp();
              double ctauPV = distXY.value()*cosAlpha * xCand.mass()/pperp.Perp();
              GlobalError v1e = (Vertex(myVertex)).error();
              GlobalError v2e = theOriginalPV.error();
              AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
              //double ctauErrPV = sqrt(vXYe.similarity(vpperp))*3.09688/(pperp.Perp2());
              double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*xCand.mass()/(pperp.Perp2());

              xCand.addUserFloat("ctauPVMuLess",ctauPV);
              xCand.addUserFloat("ctauErrPVMuLess",ctauErrPV);
              xCand.addUserFloat("cosAlphaMuLess",cosAlpha);

            }
            else
            {
              xCand.addUserFloat("ctauPVMuLess",-100.0);
              xCand.addUserFloat("ctauErrPVMuLess",ctauErrPV);
              xCand.addUserFloat("cosAlphaMuLess",cosAlpha);
            }

            // lifetime using BS
            float cosAlphaBS, cosAlphaBS3D, ctauBS, ctauErrBS, l_xyBS, lErr_xyBS, l_xyzBS, lErr_xyzBS;

            pvtx.SetXYZ(theBeamSpotV.position().x(),theBeamSpotV.position().y(),0);
            vdiff = vtx - pvtx;
            cosAlphaBS = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
            distXY = vdistXY.distance(Vertex(myVertex), theBeamSpotV);
            //double ctauBS = distXY.value()*cosAlpha*3.09688/pperp.Perp();

            ctauBS = distXY.value()*cosAlpha*xCand.mass()/pperp.Perp();

            GlobalError v1eB = (Vertex(myVertex)).error();
            GlobalError v2eB = theBeamSpotV.error();
            AlgebraicSymMatrix33 vXYeB = v1eB.matrix()+ v2eB.matrix();
            //double ctauErrBS = sqrt(vXYeB.similarity(vpperp))*3.09688/(pperp.Perp2());
            ctauErrBS = sqrt(ROOT::Math::Similarity(vpperp,vXYeB))*xCand.mass()/(pperp.Perp2());

            vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0 ;
            l_xyBS = vdiff.Perp();
            lErr_xyBS = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();

            /// 3D
            pvtx3D.SetXYZ(theBeamSpotV.position().x(), theBeamSpotV.position().y(), theBeamSpotV.position().z());
            vdiff3D = vtx3D - pvtx3D;
            cosAlphaBS3D = vdiff3D.Dot(pperp3D)/(vdiff3D.Mag()*vdiff3D.Mag());
            l_xyzBS = vdiff3D.Mag();
            vDiff3D[0] = vdiff3D.x(); vDiff3D[1] = vdiff3D.y(); vDiff3D[2] = vdiff3D.z() ;
            lErr_xyzBS = sqrt(ROOT::Math::Similarity(vDiff3D,vXYe)) / vdiff3D.Mag();

            xCand.addUserFloat("ctauBS",ctauBS);
            xCand.addUserFloat("ctauErrBS",ctauErrBS);

            xCand.addUserFloat("cosAlphaBS",cosAlphaBS);
            xCand.addUserFloat("cosAlphaBS3D",cosAlphaBS3D);

            xCand.addUserFloat("l_xyBS",l_xyBS);
            xCand.addUserFloat("lErr_xyBS",lErr_xyBS);

            xCand.addUserFloat("l_xyzBS",l_xyzBS);
            xCand.addUserFloat("lErr_xyzBS",lErr_xyzBS);

            xCand.addUserData("commonVertex",Vertex(myVertex));

          } else {
            xCand.addUserFloat("vNChi2",-1);
            xCand.addUserFloat("vProb", -1);
            xCand.addUserFloat("ctauPV",-100);
            xCand.addUserFloat("ctauErrPV",-100);

            xCand.addUserFloat("ctauBS",-100);
            xCand.addUserFloat("ctauErrBS",-100);

            xCand.addUserInt("countTksOfPV", -1);
            xCand.addUserFloat("vertexWeight", -100.);
            xCand.addUserFloat("sumPTPV", -100.);

            xCand.addUserData("commonVertex",Vertex());
            if (addMuonlessPrimaryVertex_)
            {
              xCand.addUserData("muonlessPV",Vertex());
              xCand.addUserFloat("ctauPVMuLess",-100);
              xCand.addUserFloat("ctauErrPVMuLess",-100);
              xCand.addUserFloat("cosAlphaMuLess",-100);
            }

            xCand.addUserFloat("cosAlpha",-100);
            xCand.addUserFloat("cosAlpha3D",-100);
            xCand.addUserFloat("cosAlphaBS",-100);
            xCand.addUserFloat("cosAlphaBS3D",-100);

            xCand.addUserFloat("l_xy",-100);
            xCand.addUserFloat("lErr_xy",-100);

            xCand.addUserFloat("l_xyz",-100);
            xCand.addUserFloat("lErr_xyz",-100);

            xCand.addUserFloat("l_xyBS",-100);
            xCand.addUserFloat("lErr_xyBS",-100);

            xCand.addUserFloat("l_xyzBS",-100);
            xCand.addUserFloat("lErr_xyzBS",-100);
            xCand.addUserData("PVwithmuons",Vertex());
          }

        }
        else
        continue;


        const reco::Vertex *ipv = phiCand->userData<reco::Vertex>("commonVertex");
        float dzPhi = fabs(Getdz(*phiCand,ipv->position()));              // onia2mumu stores vertex as userData
        xCand.addUserFloat("dzPhi",dzPhi);

        if (!cutdz(dzPhi)){
          dz_phi_cut_fail++;
          continue;
        }
        //std::cout << "Phi dz cut passed" << std::endl;
        ipv = jpsiCand->userData<reco::Vertex>("commonVertex");
        float dzJpsi = fabs(Getdz(*jpsiCand,ipv->position()));              // onia2mumu stores vertex as userData
        xCand.addUserFloat("dzJpsi",dzJpsi);
        //std::cout << "Jps dz cut passed" << std::endl;
        if (!cutdz(dzJpsi)){
          dz_jps_cut_fail++;
          continue;
        }

        float dz = fabs(Getdz(*phiCand,ipv->position()));
        xCand.addUserFloat("dzFourMuons",dz);


        xCandColl->push_back(xCand);
        candidates++;
      }
    }

    std::sort(xCandColl->begin(),xCandColl->end(),vPComparator_);

    event.put(std::move(xCandColl));
  }

  float FourOniaProducer::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {

    const reco::Candidate::LorentzVector& mom = c.p4();
    const reco::Candidate::Point& vtx = c.vertex();

    double dz = (vtx.Z()-p.Z()) - ((vtx.X()-p.X())*mom.X()+(vtx.Y()-p.Y())*mom.Y())/mom.Rho() * mom.Z()/mom.Rho();
    return (float) dz;

  }

  void FourOniaProducer::endJob(){
    std::cout << "###########################" << std::endl;
    std::cout << "X Candidate producer report:" << std::endl;
    std::cout << "###########################" << std::endl;
    // std::cout << "Delta mass fail: " << delta_mass_fail << std::endl;
    std::cout << "Dz phi fail:      " << dz_phi_cut_fail << std::endl;
    std::cout << "Dz jps fail:      " << dz_jps_cut_fail << std::endl;
    std::cout << "###########################" << std::endl;
    std::cout << "Found " << candidates << " X candidates." << std::endl;
    std::cout << "###########################" << std::endl;
  }
  //
  // UInt_t FourOniaProducer::isTriggerMatched(const pat::CompositeCandidate *diMuon_cand) {
  //   const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  //   const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  //   UInt_t matched = 0;  // if no list is given, is not matched
  //
  //   // if matched a given trigger, set the bit, in the same order as listed
  //   for (unsigned int iTr = 0; iTr<HLTFilters_.size(); iTr++ ) {
  //     const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
  //     const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
  //     if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
  //   }
  //   return matched;
  // }

  bool FourOniaProducer::areOverlappedMuons(const pat::CompositeCandidate *phi,const pat::CompositeCandidate *jpsi) {

    bool same = false;

    std::vector<const pat::Muon*> fourMuons;

    fourMuons.push_back(dynamic_cast<const pat::Muon*>(phi->daughter("muon1")) );
    fourMuons.push_back(dynamic_cast<const pat::Muon*>(phi->daughter("muon2")) );

    fourMuons.push_back(dynamic_cast<const pat::Muon*>(jpsi->daughter("muon1")) );
    fourMuons.push_back(dynamic_cast<const pat::Muon*>(jpsi->daughter("muon2")) );

    for (size_t i = 0; i < fourMuons.size(); i++)
    for (size_t j = i+1; j < fourMuons.size(); j++)
    same = same || (fourMuons[i]->track() == fourMuons[j]->track());

    // const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(phi->daughter("muon2"));
    // const pat::Muon* muon3 = dynamic_cast<const pat::Muon*>(jpsi->daughter("muon1"));
    // const pat::Muon* muon4 = dynamic_cast<const pat::Muon*>(jpsi->daughter("muon2"));

    return same;

  }

  const pat::CompositeCandidate FourOniaProducer::makeCandidate(const pat::CompositeCandidate& phi,
    const pat::CompositeCandidate& jpsi){
      pat::CompositeCandidate xCand;
      xCand.addDaughter(phi,"phi");
      xCand.addDaughter(jpsi,"jpsi");
      reco::Candidate::LorentzVector vX = phi.p4() + jpsi.p4();
      xCand.setP4(vX);
      return xCand;
    }


    //define this as a plug-in
    DEFINE_FWK_MODULE(FourOniaProducer);
