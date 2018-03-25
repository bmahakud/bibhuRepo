// -*- C++ -*-
//
// Package:    MuonHLT/MuonHLTDebugger
// Class:      MuonHLTDebugger
// 
/**\class MuonHLTDebugger MuonHLTDebugger.cc MuonHLT/MuonHLTDebugger/plugins/MuonHLTDebugger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Santiago Folgueras
//         Created:  Thu, 22 Sep 2016 13:30:13 GMT
//
//

//#define USEGENINFO 1

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"


//
// class declaration
//
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"

#include <cassert>


#include "RecoTracker/FinalTrackSelectors/interface/TrackMVAClassifier.h"
#include "TLorentzVector.h"
#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"
#include "TPRegexp.h"
#include "TString.h"
//class TrajectoryStateOnSurface;
//class TrajectorySeed;
class MuonServiceProxy;

const double NOMATCH = 999.;
const double NOMATCHITS =  0.;
const std::string EFFICIENCY_SUFFIXES[2] = {"den", "num"};

class MuonHLTDebugger : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit MuonHLTDebugger(const edm::ParameterSet&);
  ~MuonHLTDebugger();
  
  template <class T1, class T2> std::vector<size_t> matchByDeltaR(const std::vector<T1> &, const std::vector<T2> &, const double maxDeltaR = NOMATCH);
  std::vector<size_t> matchBySharedHits(const reco::MuonCollection& Muons, trigger::TriggerObjectCollection& hltl3muons, 
					const edm::Handle<reco::MuonTrackLinksCollection>& links, const double minSharedFrac=NOMATCHITS) ; 

  std::vector<size_t> matchBySharedHits(const reco::TrackExtraCollection& Muons, trigger::TriggerObjectCollection& hltl3muons, 
					const edm::Handle<reco::MuonTrackLinksCollection>& links, const double minSharedFrac=NOMATCHITS) ; 
  
  int sharedHits(const reco::Track& track1, const reco::Track& track2) const;
  int sharedHits(const reco::TrackExtra& track1, const reco::Track& track2) const;
  
  reco::MuonCollection selectedMuons(const reco::MuonCollection &);//, const StringCutObjectSelector<reco::Muon> &);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  int num3DLayers(reco::Track const & tk, bool isHLT);


private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) ;
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup) ;
  virtual void endJob() ;

  // ----------member data ---------------------------
  
  // Extra Methods
  
  // Member Variables
  HLTConfigProvider hltConfig_;

  //  edm::InputTag vertexTag_;
  //  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::InputTag muonTag_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::TrackExtraCollection>  TrackCollectionToken_;
  
  edm::InputTag genTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

  edm::InputTag l3candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l3candToken_; 
  edm::InputTag l2candTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l2candToken_; 
  edm::InputTag l1candTag_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1candToken_; 

  // Trigger process
  edm::InputTag triggerResultsTag_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultsToken_;
  edm::InputTag triggerSummaryTag_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummaryToken_;
  
  // Trigger process
  std::string triggerProcess_;
  std::string triggerName_;
  std::string l1filterLabel_;
  std::string l2filterLabel_;
  std::string l3filterLabel_;

  unsigned int debuglevel_;
  bool isMC_;
  bool runSharedHits_;
  bool UseGenInfo_;
  // Trigger indexes
  int tagTriggerIndex_;
  int triggerIndex_;
  
  edm::Service<TFileService> outfile_;
  std::map<std::string, TH1*> hists_;     

  //  edm::EDGetTokenT<reco::TrackCollection>           theIOTracksToken_;
  //  edm::EDGetTokenT<reco::TrackCollection>           theIOFromL1TracksToken_;
  //  edm::EDGetTokenT<reco::TrackCollection>           theOITracksToken_;
  edm::InputTag  theTracksOICandTag_;
  edm::EDGetTokenT<TrackCandidateCollection>    theTracksOICandToken_;
  edm::InputTag  theTracksOINoHPTag_;
  edm::EDGetTokenT<reco::TrackCollection>           theTracksOINoHPToken_;
  edm::InputTag  theTracksOITag_;
  edm::EDGetTokenT<reco::TrackCollection>           theTracksOIToken_;

  edm::InputTag  recoTrackTag_;
  edm::EDGetTokenT<reco::TrackCollection>           recoTrackToken_;


  edm::EDGetTokenT<reco::MuonTrackLinksCollection>  theLinksToken_;
  edm::EDGetTokenT<reco::TrackExtraCollection>  theMuonsWithHitsToken_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > >  puSummaryInfo_;
  edm::EDGetTokenT<reco::BeamSpot>                     theBeamSpotToken_;
  edm::EDGetTokenT<reco::VertexCollection>  theVertexToken_;
  
  // SEEDS... 
  edm::InputTag  theSeedsOITag_;
  edm::EDGetTokenT<TrajectorySeedCollection> theSeedsOIToken_;  
  edm::EDGetTokenT<L3MuonTrajectorySeedCollection> theSeedsOISToken_;
  edm::EDGetTokenT<L3MuonTrajectorySeedCollection> theSeedsOIHToken_;
  
  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> theB;
  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  edm::ESHandle<GlobalTrackingGeometry> trackingGeometry_;





  
  std::string ttrhbuilder_ ;
  
  //  StringCutObjectSelector<reco::Muon> targetMuonSelector_;
  MuonServiceProxy *theService;
};

//
// constants, enums and typedefs
//
//double pt_bins[12]  = { 20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double pt_bins[19]  = { 3, 5, 7, 10, 12, 15, 17, 20,  24,  27,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };

double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double phi_bins[14] = {-M_PI, -M_PI*(11/12.), -M_PI*(9/12.), -M_PI*(7/12.), -M_PI*(5/12.), -M_PI*(3/12.), -M_PI*(1/12.), M_PI*(1/12.), M_PI*(3/12.), M_PI*(5/12.), M_PI*(7/12.), M_PI*(9/12.), M_PI*(11/12.), M_PI};
using namespace edm;
using namespace std;

//
// static data member definitions
//

//
// constructors and destructor
//
MuonHLTDebugger::MuonHLTDebugger(const edm::ParameterSet& cfg):
  muonTag_                (cfg.getUntrackedParameter<edm::InputTag>("muonTag")),
  muonToken_              (consumes<std::vector<reco::Muon> >(muonTag_)),
  TrackCollectionToken_   (consumes<reco::TrackExtraCollection>(edm::InputTag("globalMuons"))),
  //recoTrackToken_(consumes<reco::TrackCollection>(edm::InputTag("globalMuons"))),



  genTag_                 (cfg.getUntrackedParameter<edm::InputTag>("genParticlesTag")),
  genToken_               (consumes<reco::GenParticleCollection>(genTag_)),
  l3candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L3Candidates")),
  l3candToken_            (consumes<reco::RecoChargedCandidateCollection>(l3candTag_)),
  l2candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L2Candidates")),
  l2candToken_            (consumes<reco::RecoChargedCandidateCollection>(l2candTag_)),
  l1candTag_              (cfg.getUntrackedParameter<edm::InputTag>("L1Candidates")),
  l1candToken_            (consumes<l1t::MuonBxCollection>(l1candTag_)),
  triggerResultsTag_      (cfg.getUntrackedParameter<edm::InputTag>("triggerResults")), 
  triggerResultsToken_    (consumes<edm::TriggerResults>(triggerResultsTag_)),
  triggerSummaryTag_      (cfg.getUntrackedParameter<edm::InputTag>("triggerSummary")), 
  triggerSummaryToken_    (consumes<trigger::TriggerEvent>(triggerSummaryTag_)),
  triggerProcess_         (cfg.getParameter<std::string>("triggerProcess")), 
  triggerName_            (cfg.getParameter<std::string>("triggerName")), 
  l1filterLabel_          (cfg.getParameter<std::string>("l1filterLabel")), 
  l2filterLabel_          (cfg.getParameter<std::string>("l2filterLabel")), 
  l3filterLabel_          (cfg.getParameter<std::string>("l3filterLabel")),
  debuglevel_             (cfg.getUntrackedParameter<unsigned int>("debuglevel")),
  isMC_                   (cfg.getUntrackedParameter<bool>("isMC")),
  runSharedHits_          (cfg.getUntrackedParameter<bool>("runSharedHits")),
  //UseGenInfo_             (cfg.getUntrackedParameter<bool>("UseGenInfo")),
  //  theOITracksToken_       (consumes<reco::TrackCollection>(edm::InputTag("hltIterL3OIMuCtfWithMaterialTracks","","TEST"))),
  //  thePixelTracksIter0Token_ (mayConsume<reco::TrackCollection>(edm::InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks","","TEST"))),
  //  theSeedsIter0Token_     (mayConsume<TrajectorySeedCollection>(edm::InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks","","TEST"))),
  //  theSeedsIter2Token_     (mayConsume<TrajectorySeedCollection>(edm::InputTag("hltIter2IterL3MuonPixelSeeds","","TEST"))),
  //  theTracksCandIter0Token_(mayConsume<TrackCandidateCollection>(edm::InputTag("hltIter0IterL3MuonCkfTrackCandidates","","TEST"))),
  //  theTracksCandIter2Token_(mayConsume<TrackCandidateCollection>(edm::InputTag("hltIter2IterL3MuonCkfTrackCandidates","","TEST"))),
  //  theTracksCandOIToken_   (mayConsume<TrackCandidateCollection>(edm::InputTag("hltIterL3OITrackCandidates","","TEST"))),
  //  theTracksCandOISToken_   (mayConsume<TrackCandidateCollection>(edm::InputTag("hltL3TrackCandidateFromL2OIState","","TEST"))),
  //  theTracksCandOIHToken_   (mayConsume<TrackCandidateCollection>(edm::InputTag("hltL3TrackCandidateFromL2OIHit","","TEST"))),
  //  theTracksNoHPIter0Token_(mayConsume<reco::TrackCollection>(edm::InputTag("hltIter0IterL3MuonCtfWithMaterialTracks","","TEST"))),
  //  theTracksNoHPIter2Token_(mayConsume<reco::TrackCollection>(edm::InputTag("hltIter2IterL3MuonCtfWithMaterialTracks","","TEST"))),
  //  theTracksIter0Token_    (mayConsume<reco::TrackCollection>(edm::InputTag("hltIter0IterL3MuonTrackSelectionHighPurity","","TEST"))),
  //  theTracksIter2Token_    (mayConsume<reco::TrackCollection>(edm::InputTag("hltIter2IterL3MuonTrackSelectionHighPurity","","TEST"))),
  theTracksOICandTag_     (cfg.getUntrackedParameter<edm::InputTag>("theTracksOICand")),
  theTracksOICandToken_   (mayConsume<TrackCandidateCollection>(theTracksOICandTag_)),
  theTracksOINoHPTag_     (cfg.getUntrackedParameter<edm::InputTag>("theTracksOINoHP")),
  theTracksOINoHPToken_   (mayConsume<reco::TrackCollection>(theTracksOINoHPTag_)),
  theTracksOITag_         (cfg.getUntrackedParameter<edm::InputTag>("theTracksOI")),
  theTracksOIToken_       (mayConsume<reco::TrackCollection>(theTracksOITag_)),

  recoTrackTag_        (cfg.getUntrackedParameter<edm::InputTag>("recoTrackTag")),
   recoTrackToken_        (mayConsume<reco::TrackCollection>(recoTrackTag_)),


  
  theLinksToken_          (consumes<reco::MuonTrackLinksCollection>(cfg.getUntrackedParameter<edm::InputTag>("MuonLinksTag"))),
  theMuonsWithHitsToken_  (consumes<reco::TrackExtraCollection>(edm::InputTag("globalMuons"))),
  puSummaryInfo_          (consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("addPileupInfo"))),
  theBeamSpotToken_       (consumes<reco::BeamSpot>(edm::InputTag("hltOnlineBeamSpot"))),
  theVertexToken_         (consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"))),

  theSeedsOITag_          (cfg.getUntrackedParameter<edm::InputTag>("theTracksOISeeds")),
  theSeedsOIToken_        (mayConsume<TrajectorySeedCollection>(theSeedsOITag_)),
  theSeedsOISToken_       (mayConsume<L3MuonTrajectorySeedCollection>(edm::InputTag("hltL3TrajSeedOIState","","TEST"))),
  theSeedsOIHToken_       (mayConsume<L3MuonTrajectorySeedCollection>(edm::InputTag("hltL3TrajSeedOIHit","","TEST")))
  //  targetMuonSelector_     ("isGlobalMuon && abs(eta) < 2.4 && pt > 10")
{
  theService = new MuonServiceProxy(cfg.getParameter<ParameterSet>("ServiceParameters"));
  
  usesResource("TFileService");
}


MuonHLTDebugger::~MuonHLTDebugger()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//
// ------------ method called for each event  ------------
void
MuonHLTDebugger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  theService->update(iSetup);
  
  /// READING THE OBJECTS
  edm::Handle<reco::GenParticleCollection> genParticles;  
 // edm::Handle<reco::MuonCollection> muons;
 // iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(genToken_, genParticles);

    
  int nGoodVtx=0;
  edm::Handle<reco::VertexCollection> pvHandle; 
  iEvent.getByToken(theVertexToken_, pvHandle);
  const reco::VertexCollection & vertices = *pvHandle.product();
  for(reco::VertexCollection::const_iterator it=vertices.begin(); it!=vertices.end(); ++it) {
    if( it->ndof()>4                     && 
	(fabs(it->z())<=24.)             && 
	(fabs(it->position().rho())<=2.)   ) 
      nGoodVtx++;
  }
  if( nGoodVtx==0 ) return;
  //  const reco::Vertex & pv = vertices[0];

  const reco::Vertex           & pv      = pvHandle->at(0);


  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(theBeamSpotToken_,recoBeamSpotHandle);
  const reco::BeamSpot& beamSpot = *recoBeamSpotHandle;

  edm::Handle <reco::RecoChargedCandidateCollection> l2Muons;
  iEvent.getByToken(l2candToken_,l2Muons);
  edm::Handle <reco::RecoChargedCandidateCollection> l3Muons;
  iEvent.getByToken(l3candToken_, l3Muons);

  // Loop over muons and fill histograms: 
  int numGenPerEvent=0;
  std::vector<const reco::GenParticle*> targetMuons;
    //    for (auto& gen: *(genParticles.product())) {
    for (unsigned int g(0); g < genParticles->size(); ++g){
      const reco::GenParticle* gen = &genParticles->at(g);
      if (fabs(gen->pdgId())!=13) continue;
      if (gen->pt()<10)           continue;
      if (gen->status()!=1)       continue;
      if (fabs(gen->eta())>2.4)   continue;
      ++numGenPerEvent;
      targetMuons.push_back(gen);
      
      if (debuglevel_ > 1) 
	std::cout << "gen muon found: pt: " << gen->pt() 
		  << " eta: " << gen->eta() 
		  << " phi: " << gen->phi() << std::endl;
     }

  if ( numGenPerEvent==0) return; //if no st1 muon skip the event.

    edm::Handle<std::vector<reco::Muon> > muons;
    iEvent.getByToken(muonToken_, muons);



  
    edm::Handle<reco::TrackExtraCollection> trackMuons;
    iEvent.getByToken(theMuonsWithHitsToken_, trackMuons);
    
    reco::TrackExtraCollection trackerMuons;
  
  
  /// DEBUGGING THE OUTSIDE-IN COMPONENT   (ITER-L3) 
      edm::Handle<TrajectorySeedCollection> hltL3TrajSeedOI;
      iEvent.getByToken(theSeedsOIToken_, hltL3TrajSeedOI);
      

     edm::Handle<reco::TrackExtraCollection> TrackCollectionGM;
     iEvent.getByToken(TrackCollectionToken_, TrackCollectionGM);  

     edm::Handle<reco::TrackCollection> recoTrack;
     iEvent.getByToken(recoTrackToken_, recoTrack);

 
    /////////////////////////
  

   edm::Handle<reco::TrackCollection> hltL3TkTracksOINoHP;
      iEvent.getByToken(theTracksOINoHPToken_, hltL3TkTracksOINoHP);

       /*
      if (debuglevel_ > 1) {
        cout << "# of hltL3TkTracksOINoHP(Transient Tracks)  = " << hltL3TkTracksOINoHP->size() << endl;
        for (reco::TrackCollection::const_iterator t=hltL3TkTracksOINoHP->begin(); t!=hltL3TkTracksOINoHP->end(); t++) {
          cout<<"hltL3TkTracksOINoHP(Transient Tracks) -->      pt: " << t->pt() << " eta: " << t->eta() << " phi: " << t->phi() <<std::endl;
          cout<<"Pixel Hits(hitPattern().numberOfValidPixelHits()):    "<<t->hitPattern().numberOfValidPixelHits()<<endl;
          cout<<"dz Wrt beam spot position:        "<<t->dz(beamSpot.position())<<endl;
          cout<<"d0:  "<<t->d0()<<endl;
          cout<<"d0 Error:    "<<t->d0Error()<<endl;
          cout<<"num of Degrees of freedom: :    "<<t->ndof()<<std::endl;


          cout<<"Chi2 value(normalizedChi2()):                "<<t->normalizedChi2()<<endl;
          cout<<"Hits:(numberOfValidHits()):      "<<t->numberOfValidHits()<<endl;
          cout<<"TrkLayers with measurement:      "<<t->hitPattern().trackerLayersWithMeasurement()<<endl;
          uint32_t nlayers3D   = t->hitPattern().pixelLayersWithMeasurement();

          size_t count3D = 0;
      for ( auto it = t->recHitsBegin(), et = t->recHitsEnd(); it!=et; ++it) {
        const TrackingRecHit* hit = (*it);
        if ( trackerHitRTTI::isUndef(*hit) ) continue;

        if ( hit->dimension()==2 ) {
          auto const & thit = static_cast<BaseTrackerRecHit const&>(*hit);
          if (thit.isMatched()) count3D++;
        }
      }
      nlayers3D += count3D;

    cout<<"num3D Layers: "<<nlayers3D<<endl;


          cout<<"Trk Layers without measurement:  "<<t->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS)<<endl;
          cout<<"dxy:       "<<t->dxy(beamSpot.position())<<endl;
          cout<<"Number of lost hits: "<<t->numberOfLostHits()<<endl;





          }
       }
   */




   ///////////////////////////









  for(reco::TrackCollection::const_iterator trk=recoTrack->begin();trk!=recoTrack->end(); ++trk){
      reco::Track  track=(*trk);


      

          double DelRRecoMuon=999.0;
          TLorentzVector trkMuon;
          trkMuon.SetPtEtaPhiM(track.pt(),track.eta(),track.phi(),0);


          for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1)
          { //reco muon loop 2

          if(mu1->pt()<1.0)continue;

            double chargedDep_dR04=mu1->pfIsolationR04().sumChargedHadronPt;
            double photonDep_dR04=mu1->pfIsolationR04().sumPhotonEt;
            double neutralDep_dR04=mu1->pfIsolationR04().sumNeutralHadronEt;
            double puPt_dR04=mu1->pfIsolationR04().sumPUPt;



            double isolation= chargedDep_dR04 + std::max(0., photonDep_dR04 + neutralDep_dR04 - 0.5*puPt_dR04);
            double RelIso=isolation/mu1 -> pt();

             if(mu1 -> isGlobalMuon  () && muon::isTightMuon ( (*mu1), pv ) && RelIso<0.15){//

                   TLorentzVector lz;
                   lz.SetPtEtaPhiM(mu1 -> pt(),mu1 -> eta(),mu1 -> phi(),0);
                   
                   if(DelRRecoMuon>lz.DeltaR(trkMuon)){DelRRecoMuon=lz.DeltaR(trkMuon);}

                              }//
           } //reco muon loop 2 

        if(DelRRecoMuon>0.1)continue;

        cout<<"rec::Track muon   (recoMatched, pt,eta, phi)"<<track.pt()<<", "<<track.eta()<<", "<<track.phi()<<endl;
        //cout<<"number of valid hits: "<<track.numberOfValidHits()<<endl;
        int recHitnum=0;
        int hitsTracker=0;
        int hitsMuS=0;
        int hitsTOBTEC=0;
        int hitsTOB=0;
        int hitsTEC=0;
        for (trackingRecHit_iterator hitIt = track.recHitsBegin(); hitIt != track.recHitsEnd(); ++hitIt) {//track rechit
               recHitnum++;

                const TrackingRecHit* rechit = (*hitIt)->hit();
                DetId detid = rechit->geographicalId();
                int subDet = detid.subdetId();
                 if(detid.det()==1 ){hitsTracker++;
                  if(detid.subdetId()==5 || detid.subdetId()==6){
                   hitsTOBTEC++;

                   cout<<"printing raw ID:   "<<detid.rawId()<<endl;
                   //cout<<"recHit global position:   "<<rechit->globalPosition().x()<<endl;

                    }

                  if(detid.subdetId()==5)hitsTOB++;
                  if(detid.subdetId()==6)hitsTEC++;


                  }
                 if(detid.det()==2)hitsMuS++ ;        



                //cout<<"detID:  "<<detid.det()<<endl;
                //cout<<"subdetID:  "<<detid.subdetId()<<endl;
                //cout<<"raw id: "<<detid.rawId()<<endl;

          }
       cout<<"number of  recHits : "<<recHitnum<<endl;
       cout<<"number of tracker hits:  "<<hitsTracker<<endl;
       cout<<"number of muon system hits: "<<hitsMuS<<endl;
       cout<<"number of TOBTEC hits:  "<<hitsTOBTEC<<endl;
       cout<<"number of TOB hits:  "<<hitsTOB<<endl;
       cout<<"number of TEC hits:  "<<hitsTEC<<endl;



       ////////////////////////////////decide good or bad

      bool trackMuonGOOD=false;
      int NumMatch=0;

      for( unsigned int il3 = 0; il3 < l3Muons->size(); ++il3) {//l3 loop


        reco::RecoChargedCandidateRef candref(l3Muons, il3);

        TLorentzVector  OIMuon;
        OIMuon.SetPtEtaPhiM(candref -> pt(),candref -> eta(),candref -> phi(),0);

       double delR=trkMuon.DeltaR(OIMuon);

       if(delR<0.1){trackMuonGOOD=true;}




     }//l3 loop




    //////////////////////////////////////decide good or bad      






          for(TrajectorySeedCollection::const_iterator seed = hltL3TrajSeedOI->begin(); seed != hltL3TrajSeedOI->end(); ++seed){//loop seeds
          // cout<<"New Seed ......................"<<endl;

	    TrajectorySeed::range seedHits = seed->recHits();

            PTrajectoryStateOnDet pTSOD = seed->startingState();

	    DetId seedDetId(pTSOD.detId());

	    const GeomDet* gdet = theService->trackingGeometry()->idToDet( seedDetId );


	    TrajectoryStateOnSurface seedTSOS = trajectoryStateTransform::transientState(pTSOD, &(gdet->surface()), &*(theService)->magneticField());


            for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){//seed hit
	      if (!(*iseed).isValid()) continue; 
             float deltax(9999.);
             float deltay(9999.);
             float deltaz(9999.);

             DetId detidSeed = (*iseed).geographicalId();
	     int subDetSeed = detidSeed.subdetId();

           for (trackingRecHit_iterator hitIt = track.recHitsBegin(); hitIt != track.recHitsEnd(); ++hitIt) {//track rechit


		if (not (*hitIt)->isValid()) continue;

		const TrackingRecHit* rechit = (*hitIt)->hit();
		DetId detid = rechit->geographicalId();
		int subDet = detid.subdetId();

                if (detid.det() != detidSeed.det()) continue;

		if (subDetSeed != subDet) continue;

               //if (detid.rawId() != detidSeed.rawId()) continue;


               //cout<<"DetID:inside seedloop-- "<<detid.det()<<endl;
               //cout<<"subdet iD: inside seedloop-- "<<subDet<<endl;
               //cout<<"raw ID: insideloop--- "<<detid.rawId()<<endl;

                //cout<<"GP: "<<rechit->globalPosition().x();


                deltax = (*iseed).localPosition().x() - rechit->localPosition().x();
                deltay = (*iseed).localPosition().y() - rechit->localPosition().y();
                deltaz = (*iseed).localPosition().z() - rechit->localPosition().z();

               double dslocal=sqrt((deltax)*(deltax)+(deltay)*(deltay)+(deltaz)*(deltaz));

              // cout<<"local distance: "<<dslocal<<endl;

               //double delxglobal=(*iseed).globalPosition().x() - rechit->globalPosition().x();
              // double delyglobal=(*iseed).globalPosition().y() - rechit->globalPosition().y();
              // double delzglobal=(*iseed).globalPosition().z() - rechit->globalPosition().z();

             // double dsglobal=sqrt((delxglobal)*(delxglobal)+(delyglobal)*(delyglobal)+(delzglobal)*(delzglobal));

              if(trackMuonGOOD==false){
              hists_["distSeedhitRechitBADLocal"]->Fill(dslocal);
             // hists_["distSeedhitRechitBADGlobal"]->Fill(dsglobal);
                 }


             if(trackMuonGOOD==true){
              hists_["distSeedhitRechitGOODLocal"]->Fill(dslocal);
              //hists_["distSeedhitRechitGOODGlobal"]->Fill(dsglobal);
                 }




             // SAME layer, rod, module (TOB) or wheel, petal, ring, module (TEC) detector.
                NumMatch++;




  








                if(trackMuonGOOD==false){
              hists_["distSeedhitRechitBADLocalwRaw"]->Fill(dslocal);
              //hists_["distSeedhitRechitBADGlobalwRaw"]->Fill(dsglobal);
                 }


             if(trackMuonGOOD==true){
              hists_["distSeedhitRechitGOODLocalwRaw"]->Fill(dslocal);
             // hists_["distSeedhitRechitGOODGlobalwRaw"]->Fill(dsglobal);
                 }












                }///track rechit loop




         }//seedhits loop
 

      }//loop seeds


       if(trackMuonGOOD==false){
      hists_["nSeedhitRechitPairBAD"]->Fill(NumMatch);
         }

      if(trackMuonGOOD==true){
      hists_["nSeedhitRechitPairGOOD"]->Fill(NumMatch);
        }







     }//lopping over recoTracks










    for(std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1)
  {

   if(mu1->pt()<1.0)continue;

   double chargedDep_dR04=mu1->pfIsolationR04().sumChargedHadronPt;
   double photonDep_dR04=mu1->pfIsolationR04().sumPhotonEt;
   double neutralDep_dR04=mu1->pfIsolationR04().sumNeutralHadronEt;
   double puPt_dR04=mu1->pfIsolationR04().sumPUPt;



    double isolation= chargedDep_dR04 + std::max(0., photonDep_dR04 + neutralDep_dR04 - 0.5*puPt_dR04);
    double RelIso=isolation/mu1 -> pt();

   if(mu1 -> isGlobalMuon  () && muon::isTightMuon ( (*mu1), pv ) && RelIso<0.15){
   std::cout<<"Tight id and Isolated reco muon, pt,eta,phi:  "<<mu1 -> pt()<<", "<<mu1 -> eta()<<", "<<mu1->phi()<<std::endl;
    std::cout<<"Tight id and Isolated reco muon, px,py,pz:  "<<mu1 -> px()<<", "<<mu1 ->py()<<", "<<mu1->pz()<<std::endl;  


    }



   }



    try{
   std::cout<<"Num of L3OI muon cands: "<<l3Muons->size()<<std::endl;
     for( unsigned int il3 = 0; il3 < l3Muons->size(); ++il3) {//l3 loop


    reco::RecoChargedCandidateRef candref(l3Muons, il3);

    std::cout<<"L3 OI cand pt, eta, phi:  "<<candref -> pt()<<"  :  "<<candref -> eta()<<"  : "<<candref -> phi()<<std::endl;


     }//l3 loop




      }
     catch(...){}


















      
}

reco::MuonCollection MuonHLTDebugger::selectedMuons(const reco::MuonCollection & allMuons) { //,  const StringCutObjectSelector<reco::Muon> &selector){ 
  
  reco::MuonCollection reducedMuons;
  for (auto const& mu : allMuons){
    const reco::Track * track = 0;
    if (mu.isTrackerMuon())         track = & * mu.innerTrack();
    else if (mu.isStandAloneMuon()) track = & * mu.outerTrack();

    if (!track) continue;
    if (!mu.isGlobalMuon()) continue;
    if (std::abs(mu.eta()) > 2.4) continue;
    reducedMuons.push_back(mu);
  }
  return reducedMuons;
}

int MuonHLTDebugger::sharedHits(const reco::Track& track1, const reco::Track& track2) const {
  
  int match = 0;
  cout << "Number of RecHits Track1: " << track1.recHitsSize() << endl;
  cout << "Number of RecHits Track2: " << track2.recHitsSize() << endl;
  
  for (trackingRecHit_iterator hit1 = track1.recHitsBegin(); hit1 != track1.recHitsEnd(); ++hit1) {
    if ( !(*hit1)->isValid() ) continue;
    DetId id1 = (*hit1)->geographicalId();
    for (trackingRecHit_iterator hit2 = track2.recHitsBegin(); hit2 != track2.recHitsEnd(); ++hit2) {
      if ( !(*hit2)->isValid() ) continue;
      DetId id2 = (*hit2)->geographicalId();
      if (id2.subdetId() != id1.subdetId()) continue;
      if (id2.rawId() != id1.rawId() ) continue;

      GlobalPoint pos1 = theService->trackingGeometry()->idToDet(id1)->surface().toGlobal((*hit1)->localPosition());
      GlobalPoint pos2 = theService->trackingGeometry()->idToDet(id2)->surface().toGlobal((*hit2)->localPosition());
      double diff = ( pos1 - pos2 ).mag();
      match++;
      cout << "GLOBAL: Pos1 - Pos2 = " << pos1 << " - " << pos2 << " = " << diff << endl;
      cout << "LOCAL:  Pos1 - Pos2 = " << (*hit1)->localPosition() <<" - "<< (*hit2)->localPosition() << " = " << ((*hit1)->localPosition()-(*hit2)->localPosition()).mag()<< endl;
      
//      for(vector<TrackingRecHit*>::iterator hit=track1hits.begin();hit!=track1hits.end();++hit){
// 	for(vector<TrackingRecHit*>::iterator otherhit=track2hits.begin();otherhit!=track2hits.end();++otherhit){
// 	  cout << "Shares Input ? " << (*hit)->sharesInput(*otherhit,TrackingRecHit::some) << endl;
// 	  if((*hit)->sharesInput(*otherhit,TrackingRecHit::some)) match++;
// 	}
//      }
    }
  }
  cout << "N Shared Hits: " << match << " / " << track1.numberOfValidHits() << " --> " << track1.pt() << " , " << track1.eta() << " , " << track1.phi() << endl;
  cout << "N Shared Hits: " << match << " / " << track2.numberOfValidHits() << " --> " << track2.pt() << " , " << track2.eta() << " , " << track2.phi() << endl;
  return match;
}

int MuonHLTDebugger::sharedHits(const reco::TrackExtra& track1, const reco::Track& track2) const {
  
  int match = 0;
  cout << "Number of RecHits Track1: " << track1.recHitsSize() << endl;
  cout << "Number of RecHits Track2: " << track2.recHitsSize() << endl;
  
  for (trackingRecHit_iterator hit1 = track1.recHitsBegin(); hit1 != track1.recHitsEnd(); ++hit1) {
    if ( !(*hit1)->isValid() ) continue;
    DetId id1 = (*hit1)->geographicalId();
    for (trackingRecHit_iterator hit2 = track2.recHitsBegin(); hit2 != track2.recHitsEnd(); ++hit2) {
      if ( !(*hit2)->isValid() ) continue;
      DetId id2 = (*hit2)->geographicalId();
      if (id2.subdetId() != id1.subdetId()) continue;
      if (id2.rawId() != id1.rawId() ) continue;

      GlobalPoint pos1 = theService->trackingGeometry()->idToDet(id1)->surface().toGlobal((*hit1)->localPosition());
      GlobalPoint pos2 = theService->trackingGeometry()->idToDet(id2)->surface().toGlobal((*hit2)->localPosition());
      double diff = ( pos1 - pos2 ).mag();
      match++;
      cout << "GLOBAL: Pos1 - Pos2 = " << pos1 << " - " << pos2 << " = " << diff << endl;
      cout << "LOCAL:  Pos1 - Pos2 = " << (*hit1)->localPosition() <<" - "<< (*hit2)->localPosition() << " = " << ((*hit1)->localPosition()-(*hit2)->localPosition()).mag()<< endl;
      
//      for(vector<TrackingRecHit*>::iterator hit=track1hits.begin();hit!=track1hits.end();++hit){
// 	for(vector<TrackingRecHit*>::iterator otherhit=track2hits.begin();otherhit!=track2hits.end();++otherhit){
// 	  cout << "Shares Input ? " << (*hit)->sharesInput(*otherhit,TrackingRecHit::some) << endl;
// 	  if((*hit)->sharesInput(*otherhit,TrackingRecHit::some)) match++;
// 	}
//      }
    }
  }
  //  cout << "N Shared Hits: " << match << " / " << track1.numberOfValidHits() << " --> " << track1.pt() << " , " << track1.eta() << " , " << track1.phi() << endl;
  //  cout << "N Shared Hits: " << match << " / " << track2.numberOfValidHits() << " --> " << track2.pt() << " , " << track2.eta() << " , " << track2.phi() << endl;
  return match;
}
/*
  int MuonHLTDebugger::sharedHits(const reco::Track& track1, const reco::Track& track2) const {
  int match = 0;
  for (trackingRecHit_iterator hit1 = track1.recHitsBegin(); hit1 != track1.recHitsEnd(); ++hit1) {
  if ( !(*hit1)->isValid() ) continue;
  DetId id1 = (*hit1)->geographicalId();
  GlobalPoint pos1 = theService->trackingGeometry()->idToDet(id1)->surface().toGlobal((*hit1)->localPosition());
  for (trackingRecHit_iterator hit2 = track2.recHitsBegin(); hit2 != track2.recHitsEnd(); ++hit2) {
  if ( !(*hit2)->isValid() ) continue;
  DetId id2 = (*hit2)->geographicalId();
  if (id2.subdetId() != id1.subdetId()) continue;
  if (id2.rawId() != id1.rawId() ) continue;
  ///      if (id2.det() == DetId::Muon ) continue;
  GlobalPoint pos2 = theService->trackingGeometry()->idToDet(id2)->surface().toGlobal((*hit2)->localPosition());
  double diff = ( pos1 - pos2 ).mag();
  cout << "GLOBAL: Pos1 - Pos2 = " << pos1 << " - " << pos2 << " = " << diff << endl;
  cout << "LOCAL:  Pos1 - Pos2 = " << (*hit1)->localPosition() <<" - "<< (*hit2)->localPosition() << " = " << ((*hit1)->localPosition()-(*hit2)->localPosition()).mag()<< endl;
  cout << "Shares input? " << (*hit1)->sharesInput(*hit2,TrackingRecHit::some) << endl;
  //hists_["hltL3mu_DiffHits"] ->Fill(diff);
  if ( diff < 1e-3 ) match++;
  }
  }
  cout << "N Shared Hits: " << match << " / " << track1.numberOfValidHits() << " --> " << track1.pt() << " , " << track1.eta() << " , " << track1.phi() << endl;
  cout << "N Shared Hits: " << match << " / " << track2.numberOfValidHits() << " --> " << track2.pt() << " , " << track2.eta() << " , " << track2.phi() << endl;
  return match;
  }
*/

template <class T1, class T2> vector<size_t> MuonHLTDebugger::matchByDeltaR(const vector<T1> & collection1, const vector<T2> & collection2, const double maxDeltaR)  {
  
  const size_t n1 = collection1.size();
  const size_t n2 = collection2.size();
  
  vector<size_t> result(n1, -1);
  vector<vector<double> > deltaRMatrix(n1, vector<double>(n2, NOMATCH));
  
  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++) {
#ifndef USEGENINFO
      deltaRMatrix[i][j] = deltaR(collection1[i], collection2[j]);
#endif 
#ifdef USEGENINFO
      deltaRMatrix[i][j] = deltaR(*(collection1.at(i)), collection2[j]);
#endif
    }
  
  // Run through the matrix n1 times to make sure we've found all matches.
  for (size_t k = 0; k < n1; k++) {
    size_t i_min = -1;
    size_t j_min = -1;
    double minDeltaR = maxDeltaR;
    // find the smallest deltaR
    for (size_t i = 0; i < n1; i++)
      for (size_t j = 0; j < n2; j++)
	if (deltaRMatrix[i][j] < minDeltaR) {
	  i_min = i;
	  j_min = j;
	  minDeltaR = deltaRMatrix[i][j];
	}
    
    if (minDeltaR < maxDeltaR) {
      result[i_min] = j_min;
      deltaRMatrix[i_min] = vector<double>(n2, NOMATCH);
      for (size_t i = 0; i < n1; i++)
	deltaRMatrix[i][j_min] = NOMATCH;
    }
  }
  return result;
}

vector<size_t> MuonHLTDebugger::matchBySharedHits(const reco::MuonCollection& muons, 
						  trigger::TriggerObjectCollection& hltl3muons,
						  const edm::Handle<reco::MuonTrackLinksCollection>& links, const double minSharedFrac)  
{
  
  const size_t n1 = muons.size();
  const size_t n2 = hltl3muons.size();
  
  vector<size_t> result(n1, -1);
  vector<vector<double> > SharedHitMatrix(n1, vector<double>(n2, NOMATCHITS));
  
  for (size_t i = 0; i < n1; i++) {
    const reco::Muon  recoMu = muons[i];
    reco::TrackRef mutk = recoMu.innerTrack();
    
    for (size_t j = 0; j < n2; j++) {
      trigger::TriggerObject & l3mu = hltl3muons[j];
      
      for(unsigned int l(0); l <links->size(); ++l){
	const reco::MuonTrackLinks* link = &links->at(l);
	const reco::Track& globalTrack = *link->globalTrack();
	float dR2 = deltaR2(l3mu.eta(),l3mu.phi(),globalTrack.eta(),globalTrack.phi());
	float dPt = std::abs(l3mu.pt() - globalTrack.pt())/l3mu.pt();
	
	if (dR2 < 0.02*0.02 and dPt < 0.001) {
	  reco::TrackRef tk = link->trackerTrack();
	  SharedHitMatrix[i][j] = sharedHits(*mutk, *tk)/recoMu.innerTrack()->numberOfValidHits();
	  hists_["hltL3mu_numValidHits"]->Fill(tk->numberOfValidHits());
	  hists_["hltL3mu_normalizedChi2"]->Fill(tk->normalizedChi2());
	  hists_["hltL3mu_numOfLostHits"]->Fill(tk->numberOfLostHits());
	  hists_["hltL3mu_numValidMuonHits"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
	  
	  if (std::abs(l3mu.eta())<0.9) {
	    hists_["hltL3mu_numValidHits_barrel"]->Fill(tk->numberOfValidHits());
	    hists_["hltL3mu_normalizedChi2_barrel"]->Fill(tk->normalizedChi2());
	    hists_["hltL3mu_numOfLostHits_barrel"]->Fill(tk->numberOfLostHits());
	    hists_["hltL3mu_numValidMuonHits_barrel"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
	  }
	  else {
	    hists_["hltL3mu_numValidHits_endcap"]->Fill(tk->numberOfValidHits());
	    hists_["hltL3mu_normalizedChi2_endcap"]->Fill(tk->normalizedChi2());
	    hists_["hltL3mu_numOfLostHits_endcap"]->Fill(tk->numberOfLostHits());
	    hists_["hltL3mu_numValidMuonHits_endcap"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
	  }
	  hists_["hltL3mu_NumSharedHits"]->Fill(SharedHitMatrix[i][j]);
	  hists_["hltL3mu_FracSharedHits"]->Fill(SharedHitMatrix[i][j]/recoMu.innerTrack()->numberOfValidHits());
	}
      }
    }
    // RECO MUO param
    hists_["RecoMu_numValidHits"]->Fill(recoMu.innerTrack()->numberOfValidHits());
    hists_["RecoMu_normalizedChi2"]->Fill(recoMu.innerTrack()->normalizedChi2());
    hists_["RecoMu_numOfLostHits"]->Fill(recoMu.innerTrack()->numberOfLostHits());
    hists_["RecoMu_numValidMuonHits"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
    if (std::abs(recoMu.eta())<0.9) {
      hists_["RecoMu_numValidHits_barrel"]->Fill(recoMu.innerTrack()->numberOfValidHits());
      hists_["RecoMu_normalizedChi2_barrel"]->Fill(recoMu.innerTrack()->normalizedChi2());
      hists_["RecoMu_numOfLostHits_barrel"]->Fill(recoMu.innerTrack()->numberOfLostHits());
      hists_["RecoMu_numValidMuonHits_barrel"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
    }
    else{ 
      hists_["RecoMu_numValidHits_endcap"]->Fill(recoMu.innerTrack()->numberOfValidHits());
      hists_["RecoMu_normalizedChi2_endcap"]->Fill(recoMu.innerTrack()->normalizedChi2());
      hists_["RecoMu_numOfLostHits_endcap"]->Fill(recoMu.innerTrack()->numberOfLostHits());
      hists_["RecoMu_numValidMuonHits_endcap"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
    }
  }
  
  
  // Run through the matrix n1 times to make sure we've found all matches.
  for (size_t k = 0; k < n1; k++) {
    size_t i_min = -1;
    size_t j_min = -1;
    double maxSharedFrac = minSharedFrac;
    // find the smallest deltaR
    for (size_t i = 0; i < n1; i++)
      for (size_t j = 0; j < n2; j++)
	if (SharedHitMatrix[i][j] > maxSharedFrac) {
	  i_min = i;
	  j_min = j;
	  maxSharedFrac = SharedHitMatrix[i][j];
	}

    // If a match has been made, save it and make those candidates unavailable.
    if (maxSharedFrac > minSharedFrac) {
      result[i_min] = j_min;
      SharedHitMatrix[i_min] = vector<double>(n2, NOMATCHITS);
      for (size_t i = 0; i < n1; i++)
	SharedHitMatrix[i][j_min] = NOMATCHITS;
    }
  }
  return result;
}


int MuonHLTDebugger::num3DLayers(reco::Track const & tk, bool isHLT){
 uint32_t nlayers3D   = tk.hitPattern().pixelLayersWithMeasurement();    
    if (!isHLT)
      nlayers3D += tk.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
    else {
      size_t count3D = 0;
      for ( auto it = tk.recHitsBegin(), et = tk.recHitsEnd(); it!=et; ++it) {
	const TrackingRecHit* hit = (*it);
	if ( trackerHitRTTI::isUndef(*hit) ) continue;
	
	if ( hit->dimension()==2 ) {
	  auto const & thit = static_cast<BaseTrackerRecHit const&>(*hit);
	  if (thit.isMatched()) count3D++;
	}
      }
      nlayers3D += count3D;
    }
return nlayers3D;



}




vector<size_t> MuonHLTDebugger::matchBySharedHits(const reco::TrackExtraCollection& muons, 
						  trigger::TriggerObjectCollection& hltl3muons,
						  const edm::Handle<reco::MuonTrackLinksCollection>& links, const double minSharedFrac)  
{
  
  const size_t n1 = muons.size();
  const size_t n2 = hltl3muons.size();
  
  vector<size_t> result(n1, -1);
  vector<vector<double> > SharedHitMatrix(n1, vector<double>(n2, NOMATCHITS));
  
  for (size_t i = 0; i < n1; i++) {
    const reco::TrackExtra  recoMu = muons[i];
    //    reco::TrackExtraRef mutk = recoMu;
    
    for (size_t j = 0; j < n2; j++) {
      trigger::TriggerObject & l3mu = hltl3muons[j];
      
      for(unsigned int l(0); l <links->size(); ++l){
	const reco::MuonTrackLinks* link = &links->at(l);
	const reco::Track& globalTrack = *link->globalTrack();
	float dR2 = deltaR2(l3mu.eta(),l3mu.phi(),globalTrack.eta(),globalTrack.phi());
	float dPt = std::abs(l3mu.pt() - globalTrack.pt())/l3mu.pt();
	
	if (dR2 < 0.02*0.02 and dPt < 0.001) {
	  reco::TrackRef tk = link->trackerTrack();
	  SharedHitMatrix[i][j] = sharedHits(recoMu, *tk); //recoMu.innerTrack()->numberOfValidHits();
//	  hists_["hltL3mu_numValidHits"]->Fill(tk->numberOfValidHits());
//	  hists_["hltL3mu_normalizedChi2"]->Fill(tk->normalizedChi2());
//	  hists_["hltL3mu_numOfLostHits"]->Fill(tk->numberOfLostHits());
//	  hists_["hltL3mu_numValidMuonHits"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
//	  
//	  if (std::abs(l3mu.eta())<0.9) {
//	    hists_["hltL3mu_numValidHits_barrel"]->Fill(tk->numberOfValidHits());
//	    hists_["hltL3mu_normalizedChi2_barrel"]->Fill(tk->normalizedChi2());
//	    hists_["hltL3mu_numOfLostHits_barrel"]->Fill(tk->numberOfLostHits());
//	    hists_["hltL3mu_numValidMuonHits_barrel"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
//	  }
//	  else {
//	    hists_["hltL3mu_numValidHits_endcap"]->Fill(tk->numberOfValidHits());
//	    hists_["hltL3mu_normalizedChi2_endcap"]->Fill(tk->normalizedChi2());
//	    hists_["hltL3mu_numOfLostHits_endcap"]->Fill(tk->numberOfLostHits());
//	    hists_["hltL3mu_numValidMuonHits_endcap"]->Fill(tk->hitPattern().numberOfValidMuonHits());	
//	  }
//	  hists_["hltL3mu_NumSharedHits"]->Fill(SharedHitMatrix[i][j]);
//	  hists_["hltL3mu_FracSharedHits"]->Fill(SharedHitMatrix[i][j]/recoMu.innerTrack()->numberOfValidHits());
	}
      }
    }
//    // RECO MUO param
//    hists_["RecoMu_numValidHits"]->Fill(recoMu.innerTrack()->numberOfValidHits());
//    hists_["RecoMu_normalizedChi2"]->Fill(recoMu.innerTrack()->normalizedChi2());
//    hists_["RecoMu_numOfLostHits"]->Fill(recoMu.innerTrack()->numberOfLostHits());
//    hists_["RecoMu_numValidMuonHits"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
//    if (std::abs(recoMu.eta())<0.9) {
//      hists_["RecoMu_numValidHits_barrel"]->Fill(recoMu.innerTrack()->numberOfValidHits());
//      hists_["RecoMu_normalizedChi2_barrel"]->Fill(recoMu.innerTrack()->normalizedChi2());
//      hists_["RecoMu_numOfLostHits_barrel"]->Fill(recoMu.innerTrack()->numberOfLostHits());
//      hists_["RecoMu_numValidMuonHits_barrel"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
//    }
//    else{ 
//      hists_["RecoMu_numValidHits_endcap"]->Fill(recoMu.innerTrack()->numberOfValidHits());
//      hists_["RecoMu_normalizedChi2_endcap"]->Fill(recoMu.innerTrack()->normalizedChi2());
//      hists_["RecoMu_numOfLostHits_endcap"]->Fill(recoMu.innerTrack()->numberOfLostHits());
//      hists_["RecoMu_numValidMuonHits_endcap"]->Fill(recoMu.innerTrack()->hitPattern().numberOfValidMuonHits());
//    }
  }

  
  // Run through the matrix n1 times to make sure we've found all matches.
  for (size_t k = 0; k < n1; k++) {
    size_t i_min = -1;
    size_t j_min = -1;
    double maxSharedFrac = minSharedFrac;
    // find the smallest deltaR
    for (size_t i = 0; i < n1; i++)
      for (size_t j = 0; j < n2; j++)
	if (SharedHitMatrix[i][j] > maxSharedFrac) {
	  i_min = i;
	  j_min = j;
	  maxSharedFrac = SharedHitMatrix[i][j];
	}

    // If a match has been made, save it and make those candidates unavailable.
    if (maxSharedFrac > minSharedFrac) {
      result[i_min] = j_min;
      SharedHitMatrix[i_min] = vector<double>(n2, NOMATCHITS);
      for (size_t i = 0; i < n1; i++)
	SharedHitMatrix[i][j_min] = NOMATCHITS;
    }
  }
  return result;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonHLTDebugger::beginJob()
{

  hists_["distSeedhitRechitBADLocal"]  = outfile_->make<TH1F>("distSeedhitRechitBADLocal",  "distSeedhitRechitBADLocal; distSeedhitRechitBADLocal", 100,  0, 50);
  hists_["distSeedhitRechitBADGlobal"]  = outfile_->make<TH1F>("distSeedhitRechitBADGlobal",  "distSeedhitRechitBADGlobal; distSeedhitRechitBADGlobal", 100,  0, 50);

  hists_["distSeedhitRechitGOODLocal"]  = outfile_->make<TH1F>("distSeedhitRechitGOODLocal",  "distSeedhitRechitGOODLocal; distSeedhitRechitGOODLocal", 100,  0, 50);
  hists_["distSeedhitRechitGOODGlobal"]  = outfile_->make<TH1F>("distSeedhitRechitGOODGlobal",  "distSeedhitRechitGOODGlobal; distSeedhitRechitGOODGlobal", 100,  0, 50);



  hists_["distSeedhitRechitBADLocalwRaw"]  = outfile_->make<TH1F>("distSeedhitRechitBADLocalwRaw",  "distSeedhitRechitBADLocalwRaw; distSeedhitRechitBADLocalwRaw", 100,  0, 50);
  hists_["distSeedhitRechitBADGlobalwRaw"]  = outfile_->make<TH1F>("distSeedhitRechitBADGlobalwRaw",  "distSeedhitRechitBADGlobalwRaw; distSeedhitRechitBADGlobalwRaw", 100,  0, 50);

  hists_["distSeedhitRechitGOODLocalwRaw"]  = outfile_->make<TH1F>("distSeedhitRechitGOODLocalwRaw",  "distSeedhitRechitGOODLocalwRaw; distSeedhitRechitGOODLocalwRaw", 100,  0, 50);
  hists_["distSeedhitRechitGOODGlobalwRaw"]  = outfile_->make<TH1F>("distSeedhitRechitGOODGlobalwRaw",  "distSeedhitRechitGOODGlobalwRaw; distSeedhitRechitGOODGlobalwRaw", 100,  0, 50);






  
  hists_["nSeedhitRechitPairBAD"]  = outfile_->make<TH1F>("nSeedhitRechitPairBAD",  "nSeedhitRechitPairBAD; nSeedhitRechitPairBAD", 50,  0, 50);
  hists_["nSeedhitRechitPairGOOD"]  = outfile_->make<TH1F>("nSeedhitRechitPairGOOD",  "nSeedhitRechitPairGOOD; nSeedhitRechitPairGOOD", 50,  0, 50);








  hists_["gen_pt"]  = outfile_->make<TH1F>("gen_pt",  "Gen Muon p_{T}; gen #mu p_{T}", 30,  pt_bins[0], pt_bins[11]);
  hists_["gen_eta"] = outfile_->make<TH1F>("gen_eta", "Gen Muon #eta; gen #mu #eta", 30, eta_bins[0], eta_bins[15]);
  hists_["gen_phi"] = outfile_->make<TH1F>("gen_phi", "Gen Muon #phi; gen #mu #phi", 30, -3.3, 3.3);
 
  // L1,L2,L3 values and efficiencies: 
  hists_["hltL1_pt"]     = outfile_->make<TH1F>("hltL1_pt",  "HLT (L1) p_{T}; p_{T} of L1 object", 18, pt_bins );
  hists_["hltL1_eta"]    = outfile_->make<TH1F>("hltL1_eta", "HLT (L1) #eta; #eta of L1 object", 15, eta_bins );
  hists_["hltL1_phi"]    = outfile_->make<TH1F>("hltL1_phi", "HLT (L1) #phi;#phi of L1 object", 13, phi_bins);
  hists_["hltL1_DeltaR"] = outfile_->make<TH1F>("hltL1_DeltaR", "HLT (L1) #Delta R; #Delta wrt L1 object", 15, 0., 1.);
  hists_["hltL1p_resEta"] = outfile_->make<TH1F>("hltL1p_resEta", "L1 Resolution (+);#eta^{reco}-#eta^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1p_resPhi"] = outfile_->make<TH1F>("hltL1p_resPhi", "L1 Resolution (+);#phi^{reco}-#phi^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1m_resEta"] = outfile_->make<TH1F>("hltL1m_resEta", "L1 Resolution (-);#eta^{reco}-#eta^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1m_resPhi"] = outfile_->make<TH1F>("hltL1m_resPhi", "L1 Resolution (-);#phi^{reco}-#phi^{HLT}",  100,  -0.1,   0.1);
  hists_["hltL1p_resEta_barrel"] = outfile_->make<TH1F>("hltL1p_resEta_barrel", "L1 Resolution (+);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resEta_endcap"] = outfile_->make<TH1F>("hltL1p_resEta_endcap", "L1 Resolution (+);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resPhi_barrel"] = outfile_->make<TH1F>("hltL1p_resPhi_barrel", "L1 Resolution (+);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resPhi_endcap"] = outfile_->make<TH1F>("hltL1p_resPhi_endcap", "L1 Resolution (+);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resEta_barrel"] = outfile_->make<TH1F>("hltL1m_resEta_barrel", "L1 Resolution (-);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resEta_endcap"] = outfile_->make<TH1F>("hltL1m_resEta_endcap", "L1 Resolution (-);#eta^{HLT}/#eta^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resPhi_barrel"] = outfile_->make<TH1F>("hltL1m_resPhi_barrel", "L1 Resolution (-);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1m_resPhi_endcap"] = outfile_->make<TH1F>("hltL1m_resPhi_endcap", "L1 Resolution (-);#phi^{HLT}/#phi^{GEN}-1",  100,  -0.25,   0.25);
  hists_["hltL1p_resPt"]  = outfile_->make<TH1F>("hltL1p_resPt",  "L1 Resolution (+);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);
  hists_["hltL1m_resPt"]  = outfile_->make<TH1F>("hltL1p_resPt",  "L1 Resolution (+);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);

  hists_["hltL2p_resEta"]        = outfile_->make<TH1F>("hltL2p_resEta", "L2 Resolution (+);#eta^{reco}-#eta^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2p_resPhi"]        = outfile_->make<TH1F>("hltL2p_resPhi", "L2 Resolution (+);#phi^{reco}-#phi^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2m_resEta"]        = outfile_->make<TH1F>("hltL2m_resEta", "L2 Resolution (-);#eta^{reco}-#eta^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2m_resPhi"]        = outfile_->make<TH1F>("hltL2m_resPhi", "L2 Resolution (-);#phi^{reco}-#phi^{HLT}",  50,  -0.05,   0.05);
  hists_["hltL2p_resEta_barrel"] = outfile_->make<TH1F>("hltL2p_resEta_barrel", "L2 Resolution (+);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resEta_endcap"] = outfile_->make<TH1F>("hltL2p_resEta_endcap", "L2 Resolution (+);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resPhi_barrel"] = outfile_->make<TH1F>("hltL2p_resPhi_barrel", "L2 Resolution (+);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resPhi_endcap"] = outfile_->make<TH1F>("hltL2p_resPhi_endcap", "L2 Resolution (+);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resEta_barrel"] = outfile_->make<TH1F>("hltL2m_resEta_barrel", "L2 Resolution (-);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resEta_endcap"] = outfile_->make<TH1F>("hltL2m_resEta_endcap", "L2 Resolution (-);#eta^{HLT}/#eta^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resPhi_barrel"] = outfile_->make<TH1F>("hltL2m_resPhi_barrel", "L2 Resolution (-);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2m_resPhi_endcap"] = outfile_->make<TH1F>("hltL2m_resPhi_endcap", "L2 Resolution (-);#phi^{HLT}/#phi^{gen}-1",  100,  -0.25,   0.25);
  hists_["hltL2p_resPt"]         = outfile_->make<TH1F>("hltL2p_resPt",         "L2 Resolution (+);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);
  hists_["hltL2m_resPt"]         = outfile_->make<TH1F>("hltL2m_resPt",         "L2 Resolution (-);p_{T}^{reco}-p_{T}^{HLT}", 60,  -0.30,   0.30);

  hists_["hltL2_pt"]            = outfile_->make<TH1F>("hltL2_pt",  "HLT (L2) p_{T}; p_{T} of L2 object", 18, pt_bins );
  hists_["hltL2_eta"]           = outfile_->make<TH1F>("hltL2_eta", "HLT (L2) #eta; #eta of L2 object", 15, eta_bins );
  hists_["hltL2_phi"]           = outfile_->make<TH1F>("hltL2_phi", "HLT (L2) #phi;#phi of L2 object", 13, phi_bins);
  hists_["hltL2_DeltaR"] = outfile_->make<TH1F>("hltL2_DeltaR", "HLT (L2) #Delta R; #Delta wrt L2 object", 15, 0., 1.);
  hists_["hltL2_resEta"] = outfile_->make<TH1F>("hltL2_resEta", "L2 Resolution;#eta^{reco}-#eta^{HLT}",  20,  -0.05,   0.05);
  hists_["hltL2_resPhi"] = outfile_->make<TH1F>("hltL2_resPhi", "L2 Resolution;#phi^{reco}-#phi^{HLT}",  20,  -0.05,   0.05);
  hists_["hltL2_resPt"]  = outfile_->make<TH1F>("hltL2_resPt",  "L2 Resolution;p_{T}^{reco}-p_{T}^{HLT}", 40,  -0.5,   0.5);

  hists_["hltL2L1_resEta"] = outfile_->make<TH1F>("hltL2L1_resEta", "L2/L1 Resolution;#eta^{reco}-#eta^{HLT}",   100,  -0.25,   0.25);
  hists_["hltL2L1_resPhi"] = outfile_->make<TH1F>("hltL2L1_resPhi", "L2/L1 Resolution;#phi^{reco}-#phi^{HLT}",   100,  -0.25,   0.25);
  hists_["hltL2L1_resPt"]  = outfile_->make<TH1F>("hltL2L1_resPt",  "L2/L1 Resolution;p_{T}^{reco}-p_{T}^{HLT}", 100,  -0.5,   0.5);

  hists_["hltL3_pt"]     = outfile_->make<TH1F>("hltL3_pt",  "HLT (L3) p_{T}; p_{T} of L3 object", 18, pt_bins );
  hists_["hltL3_eta"]    = outfile_->make<TH1F>("hltL3_eta", "HLT (L3) #eta; #eta of L3 object", 15, eta_bins );
  hists_["hltL3_phi"]    = outfile_->make<TH1F>("hltL3_phi", "HLT (L3) #phi;#phi of L3 object", 13, phi_bins);
  hists_["hltL3_DeltaR"] = outfile_->make<TH1F>("hltL3_DeltaR", "HLT (L3) #Delta R; #Delta wrt L3 object", 15, 0., 1.);
  hists_["hltL3_resEta"] = outfile_->make<TH1F>("hltL3_resEta", "L3 Resolution;#eta^{reco}-#eta^{HLT}",  20,  -0.01,   0.01);
  hists_["hltL3_resPhi"] = outfile_->make<TH1F>("hltL3_resPhi", "L3 Resolution;#phi^{reco}-#phi^{HLT}",  20,  -0.01,   0.01);
  hists_["hltL3_resPt"]  = outfile_->make<TH1F>("hltL3_resPt",  "L3 Resolution;p_{T}^{reco}-p_{T}^{HLT}", 40,  -0.30,   0.30);
  
  //// COUNTERS
  hists_["hlt_FracL2Match" ] = outfile_->make<TH1F>("hlt_FracL2Match","Fracber of L2 Matched", 500, 0., 1.0);
  hists_["hlt_FracL3Match" ] = outfile_->make<TH1F>("hlt_FracL3Match","Fracber of L3 Matched", 500, 0., 1.0);

  hists_["hlt_NumL1Match" ] = outfile_->make<TH1F>("hlt_NumL1Match","Number of L1 Matched", 5, -0.5, 4.5);
  hists_["hlt_NumL2Match" ] = outfile_->make<TH1F>("hlt_NumL2Match","Number of L2 Matched", 5, -0.5, 4.5);
  hists_["hlt_NumL3Match" ] = outfile_->make<TH1F>("hlt_NumL3Match","Number of L3 Matched", 5, -0.5, 4.5);

  hists_["hlt_NumL1" ] = outfile_->make<TH1F>("hlt_NumL1","Number of L1 Found", 5, -0.5, 4.5);
  hists_["hlt_NumL2" ] = outfile_->make<TH1F>("hlt_NumL2","Number of L2 Found", 5, -0.5, 4.5);
  hists_["hlt_NumL3" ] = outfile_->make<TH1F>("hlt_NumL3","Number of L3 Found", 5, -0.5, 4.5);

  /// OTHER CHECKS:   
  hists_["hlt_numSeeds"] = outfile_->make<TH1F>("hlt_numSeeds","Number of Seeds (Iter0+Iter2)", 50, -0.5, 49.5);  
  hists_["hlt_numSeedsIter0"] = outfile_->make<TH1F>("hlt_numSeedsIter0","Number of Seeds (Iter0)", 50, -0.5, 49.5);
  hists_["hlt_numSeedsIter2"] = outfile_->make<TH1F>("hlt_numSeedsIter2","Number of Seeds (Iter2)", 50, -0.5, 49.5);

  hists_["hlt_numSeeds_PU"] = outfile_->make<TH2F>("hlt_numSeeds_PU","Number of Seeds (Iter0+Iter2) vs NPU"    , 75,  0,   75., 50, -0.5, 49.5);  
  hists_["hlt_numSeedsIter0_PU"] = outfile_->make<TH2F>("hlt_numSeedsIter0_PU","Number of Seeds (Iter0) vs NPU", 75,  0,   75., 50, -0.5, 49.5);
  hists_["hlt_numSeedsIter2_PU"] = outfile_->make<TH2F>("hlt_numSeedsIter2_PU","Number of Seeds (Iter2) vs NPU", 75,  0,   75., 50, -0.5, 49.5);

  hists_["hlt_L3OI_numTracks"] = outfile_->make<TH1F>("hlt_L3OI_numTracks","Number of Tracks (outside-in)", 15, -0.5, 14.5);
  hists_["hlt_L3IO_numTracks"] = outfile_->make<TH1F>("hlt_L3IO_numTracks","Number of Tracks (inside-out)", 15, -0.5, 14.5);
  hists_["hlt_L3IOFromL1_numTracks"] = outfile_->make<TH1F>("hlt_L3IOFromL1_numTracks","Number of Tracks (inside-out fromL1)", 15, -0.5, 14.5);
    
  hists_["hlt_numPixelHits"]      = outfile_->make<TH1F>("hlt_numPixelHits"     ,"Number of PixelHits (Iter0+Iter2)", 50, -0.5, 49.5);
  hists_["hlt_numPixelHitsIter0"] = outfile_->make<TH1F>("hlt_numPixelHitsIter0","Number of PixelHits (Iter0)", 50, -0.5, 49.5);
  hists_["hlt_numPixelHitsIter2"] = outfile_->make<TH1F>("hlt_numPixelHitsIter2","Number of PixelHits (Iter2)", 50, -0.5, 49.5);
  hists_["hlt_numValidHits"]       = outfile_->make<TH1F>("hlt_numValidHits", "N Valid Hits Iter0+Iter2 ", 70,  -0.5, 69.5 ); 
  hists_["hlt_normalizedChi2"]     = outfile_->make<TH1F>("hlt_normalizedChi2", "Normalised Chi2 of Iter0+Iter2 ", 100, 0., 20); 
  hists_["hlt_numValidMuonHits"]   = outfile_->make<TH1F>("hlt_numValidMuonHits", "N Valid Muon Hits of Iter0+Iter2", 70,  -0.5, 69.5 );

  /// for failing L3 
  hists_["hlt_noL3_numPixelTracksIter0"] = outfile_->make<TH1F>("hlt_noL3_numPixelTracksIter0","Number of Pixel Tracks (Iter0)",           50, -0.5, 49.5);
  hists_["hlt_noL3_numSeedsIter0"]     = outfile_->make<TH1F>("hlt_noL3_numSeedsIter0"    ,"Number of Seeds of Failed L3 (Iter0)",           50, -0.5, 49.5);
  hists_["hlt_noL3_numSeedsIter2"]     = outfile_->make<TH1F>("hlt_noL3_numSeedsIter2"    ,"Number of Seeds of Failed L3 (Iter2)",           50, -0.5, 49.5);

  hists_["hlt_noL3_numTracks"]         = outfile_->make<TH1F>("hlt_noL3_numTracks"        ,"Number of Tracks (Iter0+Iter2)",   15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksIter0"]    = outfile_->make<TH1F>("hlt_noL3_numTracksIter0"   ,"Number of Tracks (Iter0)",          15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksIter2"]    = outfile_->make<TH1F>("hlt_noL3_numTracksIter2"   ,"Number of Tracks (Iter2)",          15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksCandIter0"] = outfile_->make<TH1F>("hlt_noL3_numTracksCandIter0","Number of Tracks Candidates (Iter0)", 15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksCandIter2"] = outfile_->make<TH1F>("hlt_noL3_numTracksCandIter2","Number of Tracks Candidates (Iter2)", 15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksNoHPIter0"] = outfile_->make<TH1F>("hlt_noL3_numTracksNoHPIter0","Number of Tracks No HP (Iter0)", 15, -0.5, 14.5);
  hists_["hlt_noL3_numTracksNoHPIter2"] = outfile_->make<TH1F>("hlt_noL3_numTracksNoHPIter2","Number of Tracks No HP (Iter2)", 15, -0.5, 14.5);

  hists_["hlt_noL3_numPixelHits"]      = outfile_->make<TH1F>("hlt_noL3_numPixelHits"     ,"Number of PixelHits of Failed L3 (Iter0+Iter2)", 50, -0.5, 49.5);
  hists_["hlt_noL3_numPixelHitsIter0"] = outfile_->make<TH1F>("hlt_noL3_numPixelHitsIter0","Number of PixelHits of Failed L3 (Iter0)",       50, -0.5, 49.5);
  hists_["hlt_noL3_numPixelHitsIter2"] = outfile_->make<TH1F>("hlt_noL3_numPixelHitsIter2","Number of PixelHits of Failed L3 (Iter2)",       50, -0.5, 49.5);
  hists_["hlt_noL3_numValidHits"]      = outfile_->make<TH1F>("hlt_noL3_numValidHits"     ,"Number of Valid Hits of Failed L3 ",             70, -0.5, 69.5); 
  hists_["hlt_noL3_normalizedChi2"]    = outfile_->make<TH1F>("hlt_noL3_normalizedChi2"   ,"Normalised Chi2 of Failed L3 Muon",             100,  0. , 10. ); 
  hists_["hlt_noL3_numValidMuonHits"]  = outfile_->make<TH1F>("hlt_noL3_numValidMuonHits" ,"Number of Valid Muon Hits of Failed L3 Muon",    70, -0.5, 69.5);

  // other track properties...
  hists_["hlt_noL3_trackpt"]       = outfile_->make<TH1F>("hlt_noL3_trackpt"  ,"Track (L3) p_{T}; p_{T} of track",  18, pt_bins );
  hists_["hlt_noL3_tracketa"]      = outfile_->make<TH1F>("hlt_noL3_tracketa" ,"Track (L3) #eta; #eta of L3 track",  15, eta_bins );
  hists_["hlt_noL3_trackphi"]      = outfile_->make<TH1F>("hlt_noL3_trackphi" ,"Track (L3) #phi; #phi of L3 track",  13, phi_bins);
  hists_["hlt_noL3_DeltaR"]        = outfile_->make<TH1F>("hlt_noL3_DeltaR",   "Iter L3 #Delta R; #Delta R (tk,L2)", 15, 0., 0.2);
  hists_["hlt_noL3_DeltaEta"]      = outfile_->make<TH1F>("hlt_noL3_DeltaEta", "Iter L3 #Delta #eta;#eta^{tk}-#eta^{L2}",  50,  -0.2,  0.2);
  hists_["hlt_noL3_DeltaPhi"]      = outfile_->make<TH1F>("hlt_noL3_DeltaPhi", "Iter L3 #Delta #phi;#phi^{tk}-#phi^{L2}",  50,  -0.2,  0.2);
  hists_["hlt_noL3_DeltaPt"]       = outfile_->make<TH1F>("hlt_noL3_DeltaPt",  "Iter L3 #Delta p_{T};p_{T}^{tk}-p_{T}^{L2}", 60,  -0.30,   0.30);
  
  hists_["hltL3mu_dR2withLink"]      = outfile_->make<TH1F>("hltL3mu_dR2withLink"     , "#Delta R2 with Link", 100, 0, 0.02);
  hists_["hltL3mu_dPtwithLink"]      = outfile_->make<TH1F>("hltL3mu_dPtwithLink"     , "#Delta pt with Link", 100, 0, 0.05);
  hists_["hltL3mu_pt"]               = outfile_->make<TH1F>("hltL3mu_pt"              , "#p_{T} of the L3 Muon; p_{T} L3 muon", 18, pt_bins );
  hists_["hltL3mu_eta"]		     = outfile_->make<TH1F>("hltL3mu_eta"             , "#eta of the L3 Muon; #eta L3 muons", 15, eta_bins ); 
  hists_["hltL3mu_dr"]		     = outfile_->make<TH1F>("hltL3mu_dr"              , "Dr w.r.t. BS", 50, 0., 0.2 ); 
  hists_["hltL3mu_dz"]		     = outfile_->make<TH1F>("hltL3mu_dz"              , "Dz w.r.t. BS", 50, 0., 0.5); 
  hists_["hltL3mu_dxySig"]	     = outfile_->make<TH1F>("hltL3mu_dxySig"          , "Dxy Significance ", 50, 0., 15.); 
  hists_["hltL3mu_dxy"]		     = outfile_->make<TH1F>("hltL3mu_dxy"             , "Dxy w.r.t. BS", 50, 0., 0.2); 
  hists_["hltL3mu_numValidHits"]     = outfile_->make<TH1F>("hltL3mu_numValidHits"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_numOfLostHits"]    = outfile_->make<TH1F>("hltL3mu_numOfLostHits"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_normalizedChi2"]   = outfile_->make<TH1F>("hltL3mu_normalizedChi2"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["hltL3mu_numValidMuonHits"] = outfile_->make<TH1F>("hltL3mu_numValidMuonHits", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["RecoMu_numValidHits"]     = outfile_->make<TH1F>("RecoMu_numValidHits"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_numOfLostHits"]    = outfile_->make<TH1F>("RecoMu_numOfLostHits"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_normalizedChi2"]   = outfile_->make<TH1F>("RecoMu_normalizedChi2"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["RecoMu_numValidMuonHits"] = outfile_->make<TH1F>("RecoMu_numValidMuonHits", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["hltL3mu_numValidHits_barrel"]     = outfile_->make<TH1F>("hltL3mu_numValidHits_barrel"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_numOfLostHits_barrel"]    = outfile_->make<TH1F>("hltL3mu_numOfLostHits_barrel"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_normalizedChi2_barrel"]   = outfile_->make<TH1F>("hltL3mu_normalizedChi2_barrel"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["hltL3mu_numValidMuonHits_barrel"] = outfile_->make<TH1F>("hltL3mu_numValidMuonHits_barrel", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["RecoMu_numValidHits_barrel"]     = outfile_->make<TH1F>("RecoMu_numValidHits_barrel"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_numOfLostHits_barrel"]    = outfile_->make<TH1F>("RecoMu_numOfLostHits_barrel"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_normalizedChi2_barrel"]   = outfile_->make<TH1F>("RecoMu_normalizedChi2_barrel"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["RecoMu_numValidMuonHits_barrel"] = outfile_->make<TH1F>("RecoMu_numValidMuonHits_barrel", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["hltL3mu_numValidHits_endcap"]     = outfile_->make<TH1F>("hltL3mu_numValidHits_endcap"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_numOfLostHits_endcap"]    = outfile_->make<TH1F>("hltL3mu_numOfLostHits_endcap"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_normalizedChi2_endcap"]   = outfile_->make<TH1F>("hltL3mu_normalizedChi2_endcap"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["hltL3mu_numValidMuonHits_endcap"] = outfile_->make<TH1F>("hltL3mu_numValidMuonHits_endcap", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["RecoMu_numValidHits_endcap"]     = outfile_->make<TH1F>("RecoMu_numValidHits_endcap"    , "N Valid Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_numOfLostHits_endcap"]    = outfile_->make<TH1F>("RecoMu_numOfLostHits_endcap"   , "N Lost Hits L3 Muon", 70,  -0.5, 69.5 ); 
  hists_["RecoMu_normalizedChi2_endcap"]   = outfile_->make<TH1F>("RecoMu_normalizedChi2_endcap"  , "Normalised Chi2 of L3 Muon", 100, 0., 20); 
  hists_["RecoMu_numValidMuonHits_endcap"] = outfile_->make<TH1F>("RecoMu_numValidMuonHits_endcap", "N Valid Muon Hits L3 Muon", 70,  -0.5, 69.5 );

  hists_["hltL3mu_DiffHits"]               = outfile_->make<TH1F>("hltL3mu_DiffHits"    , "Diff Shared Hits L3/RECO Muon", 100,  0., 5. ); 
  hists_["hltL3mu_NumSharedHits"]          = outfile_->make<TH1F>("hltL3mu_NumSharedHits"    , "N Shared Hits L3/RECO Muon", 70,  -0.5, 69.5 ); 
  hists_["hltL3mu_FracSharedHits"]         = outfile_->make<TH1F>("hltL3mu_FracValidHits"    , "Fraction Shared Hits L3/RECO Muon", 100,  0., 1.05 ); 
 
  hists_["hlt_OI_NumSeeds"]               = outfile_->make<TH1F>("hlt_OI_NumSeeds","Number of Seeds (OI)", 50, -0.5, 49.5);
  hists_["hlt_OI_NumSeedsVsEta"]          = outfile_->make<TH2F>("hlt_OI_NumSeedsVsEta","Number of Seeds (OI) vs #eta", 15, eta_bins, 50, -0.5, 49.5);
  hists_["hlt_OI_NumberOfRecHitsPerSeed"] = outfile_->make<TH1F>("hlt_OI_NumberOfRecHitsPerSeed","Number Hits per Seed (OI)", 50, -0.5, 49.5);
  hists_["hlt_OI_seedEta"]                = outfile_->make<TH1F>("hlt_OI_seedEta", "Seed eta  (OI)", 15, eta_bins);
  hists_["hlt_OI_seedEtaVsL2Eta"]         = outfile_->make<TH2F>("hlt_OI_seedEtaVsL2Eta", "Seed eta - L2 eta (OI)", 15, eta_bins,100, -1., 1.);
  hists_["hlt_OI_seedPhiVsL2Phi"]         = outfile_->make<TH2F>("hlt_OI_seedPhiVsL2Phi", "Seed phi - L2 phi (OI)", 15, eta_bins,100, -1., 1.);
  hists_["hlt_OI_seedPtErr"]              = outfile_->make<TH1F>("hlt_OI_seedPtErr", "Seed pt Error (OI)", 50, 0., 50. );
  hists_["hlt_OI_seedPtErr_barrel"]       = outfile_->make<TH1F>("hlt_OI_seedPtErr_barrel", "Seed pt Error (OI)", 50, 0., 50. );
  hists_["hlt_OI_seedPtErr_endcap"]       = outfile_->make<TH1F>("hlt_OI_seedPtErr_endcap", "Seed pt Error (OI)", 50, 0., 50. );
  hists_["hlt_OI_seedPtErrVsEta"]         = outfile_->make<TH2F>("hlt_OI_seedPrErrVsEta", "Seed Pt Error vs Eta (OI)", 15, eta_bins, 50, 0., 50. );
  hists_["hlt_OI_seedPhiErr"]             = outfile_->make<TH1F>("hlt_OI_seedPhiErr","Seed Phi error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedPhiErr_barrel"]      = outfile_->make<TH1F>("hlt_OI_seedPhiErr_barrel","Seed Phi error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedPhiErr_endcap"]      = outfile_->make<TH1F>("hlt_OI_seedPhiErr_endcap","Seed Phi error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedEtaErr"]             = outfile_->make<TH1F>("hlt_OI_seedEtaErr","Seed Eta error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedEtaErr_barrel"]      = outfile_->make<TH1F>("hlt_OI_seedEtaErr_barrel","Seed Eta error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_seedEtaErr_endcap"]      = outfile_->make<TH1F>("hlt_OI_seedEtaErr_endcap","Seed Eta error (OI)", 100, 0., 0.1);
  hists_["hlt_OI_hitPt"]                  = outfile_->make<TH1F>("hlt_OI_hitPt","Hit p_{T} (OI);p_{T} hit", 18, pt_bins );
  hists_["hlt_OI_hitPhi"]                 = outfile_->make<TH1F>("hlt_OI_hitPhi","Hit #phi (OI);#phi hit", 13, phi_bins);
  hists_["hlt_OI_hitEta"]                 = outfile_->make<TH1F>("hlt_OI_hitEta","Hit #eta (OI);#eta hit", 15,  eta_bins ) ;
  hists_["hlt_OI_hitx"]                   = outfile_->make<TH1F>("hlt_OI_hitx","Hit x (OI);x hit", 200, -100, 100);
  hists_["hlt_OI_hity"]                   = outfile_->make<TH1F>("hlt_OI_hity","Hit y (OI);y hit", 200, -100, 100);
  hists_["hlt_OI_hitz"]                   = outfile_->make<TH1F>("hlt_OI_hitz","Hit z (OI);z hit", 600, -300, 300);
  hists_["hlt_OI_hitdx"]                  = outfile_->make<TH1F>("hlt_OI_hitdx","Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_hitdy"]                  = outfile_->make<TH1F>("hlt_OI_hitdy","Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);
   
  hists_["hlt_OI_trackPt"]                = outfile_->make<TH1F>("hlt_OI_trackPt","Track (L3) p_{T}; p_{T} of track",  18, pt_bins);
  hists_["hlt_OI_trackEta"]               = outfile_->make<TH1F>("hlt_OI_trackEta","Track (L3) #eta; #eta of L3 track",  15, eta_bins);
  hists_["hlt_OI_trackPhi"]               = outfile_->make<TH1F>("hlt_OI_trackPhi","Track (L3) #phi; #phi of L3 track",  13, phi_bins);
  hists_["hlt_OI_trackChi2"]              = outfile_->make<TH1F>("hlt_OI_trackChi2"     ,"Normalised Chi2 of the Track", 100, 0., 20);
  hists_["hlt_OI_trackDxy"]               = outfile_->make<TH1F>("hlt_OI_trackDxy"     ,"d_xy of the Track w.r.t BS", 100, 0., 1);
  hists_["hlt_OI_trackDz"]                = outfile_->make<TH1F>("hlt_OI_trackDz"     ,"d_z of the Track w.r.t. BS", 100, 0., 1);
  hists_["hlt_OI_trackValidPixelHits"]    = outfile_->make<TH1F>("hlt_OI_trackPixelHits","N Valid Pixel Hits of the Track", 10,  -0.5,  9.5 );
  hists_["hlt_OI_trackValidHits"]         = outfile_->make<TH1F>("hlt_OI_trackValidHits","N Valid Hits of the Track", 30,  -0.5, 29.5 );
  hists_["hlt_OI_trackLayers"]            = outfile_->make<TH1F>("hlt_OI_trackLayers"   ,"N Layers of the Track", 30,  -0.5, 29.5 );
  
  // vs eta
  hists_["hlt_OI_trackChi2VsEta"]           = outfile_->make<TH2F>("hlt_OI_trackChi2VsEta"     ,"Normalised Chi2 of the Track", 15, eta_bins, 100, 0., 20);
  hists_["hlt_OI_trackDxyVsEta"]            = outfile_->make<TH2F>("hlt_OI_trackDxyVsEta"      ,"d_xy of the Track w.r.t BS", 15, eta_bins, 100, 0., 1.);
  hists_["hlt_OI_trackDzVsEta"]             = outfile_->make<TH2F>("hlt_OI_trackDzVsEta"       ,"d_z of the Track w.r.t. BS", 15, eta_bins, 100, 0., 1.);
  hists_["hlt_OI_trackValidPixelHitsVsEta"] = outfile_->make<TH2F>("hlt_OI_trackPixelHitsVsEta","N Valid Pixel Hits of the Track", 15,eta_bins,10,-0.5, 9.5);
  hists_["hlt_OI_trackValidHitsVsEta"]      = outfile_->make<TH2F>("hlt_OI_trackValidHitsVsEta","N Valid Hits of the Track", 15, eta_bins, 30,  -0.5, 29.5 );
  hists_["hlt_OI_trackLayersVsEta"]         = outfile_->make<TH2F>("hlt_OI_trackLayersVsEta"   ,"N Layers of the Track", 15, eta_bins, 30,  -0.5, 29.5 );

  hists_["hlt_OI_trackResChi2VsEta"]           = outfile_->make<TH2F>("hlt_OI_trackResChi2VsEta"     ,"Normalised Chi2 of the Track", 15, eta_bins,  100, -1., 1.);
  hists_["hlt_OI_trackResDxyVsEta"]            = outfile_->make<TH2F>("hlt_OI_trackResDxyVsEta"      ,"d_xy of the Track w.r.t BS", 15, eta_bins,    100, -0.1, 0.1);
  hists_["hlt_OI_trackResDzVsEta"]             = outfile_->make<TH2F>("hlt_OI_trackResDzVsEta"       ,"d_z of the Track w.r.t. BS", 15, eta_bins,    100, -0.1, 0.1);
  hists_["hlt_OI_trackResValidPixelHitsVsEta"] = outfile_->make<TH2F>("hlt_OI_trackResPixelHitsVsEta","N Valid Pixel Hits of the Track", 15,eta_bins,100, -10., 10.);
  hists_["hlt_OI_trackResValidHitsVsEta"]      = outfile_->make<TH2F>("hlt_OI_trackResValidHitsVsEta","N Valid Hits of the Track", 15, eta_bins,     100, -10., 10.);
  hists_["hlt_OI_trackResLayersVsEta"]         = outfile_->make<TH2F>("hlt_OI_trackResLayersVsEta"   ,"N Layers of the Track", 15, eta_bins,         100, -10., 10.);

  /// TOB hits
  hists_["hlt_OI_TOBhitdx"]               = outfile_->make<TH1F>("hlt_OI_TOBhitdx",      "TOB Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TOBhitdy"]               = outfile_->make<TH1F>("hlt_OI_TOBhitdy",      "TOB Hit #Delta y (OI);#Delta y hit", 100, -0.015, 0.015);
  hists_["hlt_OI_TOBmonohitdx"]           = outfile_->make<TH1F>("hlt_OI_TOBmonohitdx",  "TOB Mono Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TOBmonohitdy"]           = outfile_->make<TH1F>("hlt_OI_TOBmonohitdy",  "TOB Mono Hit #Delta y (OI);#Delta y hit", 100, -0.015, 0.015);
  hists_["hlt_OI_TOBstereohitdx"]         = outfile_->make<TH1F>("hlt_OI_TOBstereohitdx","TOB Stereo Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TOBstereohitdy"]         = outfile_->make<TH1F>("hlt_OI_TOBstereohitdy","TOB Stereo Hit #Delta y (OI);#Delta y hit", 100, -0.015, 0.015);

  /// TEC hits
  hists_["hlt_OI_TEChitdx"]               = outfile_->make<TH1F>("hlt_OI_TEChitdx",      "TEC Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TEChitdy"]               = outfile_->make<TH1F>("hlt_OI_TEChitdy",      "TEC Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);
  hists_["hlt_OI_TECmonohitdx"]           = outfile_->make<TH1F>("hlt_OI_TECmonohitdx",  "TEC Mono Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TECmonohitdy"]           = outfile_->make<TH1F>("hlt_OI_TECmonohitdy",  "TEC Mono Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);
  hists_["hlt_OI_TECstereohitdx"]         = outfile_->make<TH1F>("hlt_OI_TECstereohitdx","TEC Stereo Hit #Delta x (OI);#Delta x hit", 100, -5, 5);
  hists_["hlt_OI_TECstereohitdy"]         = outfile_->make<TH1F>("hlt_OI_TECstereohitdy","TEC Stereo Hit #Delta y (OI);#Delta y hit", 100, -0.05, 0.05);

}

void 
MuonHLTDebugger::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) 
{
  // Initialize hltConfig
  bool changed = true;
  if( hltConfig_.init(run, eventSetup, triggerProcess_, changed) ) {
  }
  else {
    std::cout << "Warning, didn't find process " << triggerProcess_.c_str() << std::endl;
    // Now crash
    assert(false);
  }
  
  triggerIndex_ = -1; 
  
  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    std::string tempName = hltConfig_.triggerName(iHltPath);
    if(tempName.find(triggerName_) != std::string::npos) {
      triggerIndex_ = int(iHltPath);
    }
    
    if( triggerIndex_>-1) break; 
  } // end for each path
  
  if( triggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find trigger " <<  triggerName_.c_str() << std::endl;
    // Now crash
    assert(false);    
  }  
}
// ------------ method called once each job just after ending the event loop  ------------
void MuonHLTDebugger::endJob() {}
void MuonHLTDebugger::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonHLTDebugger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonHLTDebugger);

//  LocalWords:  endl
