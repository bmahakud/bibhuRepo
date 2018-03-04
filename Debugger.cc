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
  edm::Handle<reco::MuonCollection> muons;
#ifndef USEGENINFO
  iEvent.getByToken(muonToken_, muons);
#endif
  if (isMC_) {
    iEvent.getByToken(genToken_, genParticles);
  }

 // if (debuglevel_ > 1)
   // cout << "#####################################################################################" << endl;
   // cout << "[EVENT] Run:Event --> " << iEvent.id().run() << " : " << iEvent.id().event() << endl;
  
  //########################### Trigger Info ###########################
  // Get objects from the event.  
  Handle<trigger::TriggerEvent> triggerSummary;
  iEvent.getByToken(triggerSummaryToken_, triggerSummary);

  if(!triggerSummary.isValid()) 
    {
      LogError("MuonHLTDebugger")<<"Missing triggerSummary collection" << endl;
      return;
    }

  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken_, triggerResults);

  if(!triggerResults.isValid()) 
    {
      LogError("MuonHLTDebugger")<<"Missing triggerResults collection" << endl;
      return;
    }
    

//  edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
//  iEvent.getByToken(puSummaryInfo_,puInfo);
//  int nbMCvtx = -1;
//  if (puInfo.isValid()) {
//    std::vector<PileupSummaryInfo>::const_iterator PVI;
//    for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
//      if(PVI->getBunchCrossing()==0){
//	nbMCvtx = PVI->getPU_NumInteractions();
//	break;
//      }
//    }
//  }

  unsigned int nGoodVtx = 0; 
#ifndef USEGENINFO
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

  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(theBeamSpotToken_,recoBeamSpotHandle);
  const reco::BeamSpot& beamSpot = *recoBeamSpotHandle;
#endif

  //Get the Reconstructed object collections:
  edm::Handle<l1t::MuonBxCollection> l1Muons;
  iEvent.getByToken(l1candToken_,l1Muons);
  edm::Handle <reco::RecoChargedCandidateCollection> l2Muons;
  iEvent.getByToken(l2candToken_,l2Muons);
  edm::Handle <reco::RecoChargedCandidateCollection> l3Muons;
  iEvent.getByToken(l3candToken_, l3Muons);

  //Get filter objects, these are the names of the paths for the Mu50 path:
  size_t L1MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag( l1filterLabel_,"",triggerProcess_));//The L1 Filter
  size_t L2MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag( l2filterLabel_,"",triggerProcess_));//The L2 Filter
  size_t L3MuonFilterIndex = (*triggerSummary).filterIndex(edm::InputTag( l3filterLabel_,"",triggerProcess_));//The L3 Filter

  trigger::TriggerObjectCollection L1MuonTrigObjects;
  trigger::TriggerObjectCollection L2MuonTrigObjects;
  trigger::TriggerObjectCollection L3MuonTrigObjects;
  trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();

  if (L1MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    //save the trigger objects corresponding to muon leg
    const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L1MuonFilterIndex);
    for (size_t j = 0; j < keysMuons.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
      L1MuonTrigObjects.push_back(foundObject);
    }
  }
  if (L2MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    //save the trigger objects corresponding to muon leg
    const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L2MuonFilterIndex);
    for (size_t j = 0; j < keysMuons.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
      L2MuonTrigObjects.push_back(foundObject);
    }
  }
  if (L3MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    //save the trigger objects corresponding to muon leg
    const trigger::Keys &keysMuons = (*triggerSummary).filterKeys(L3MuonFilterIndex);
    for (size_t j = 0; j < keysMuons.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysMuons[j]];
      L3MuonTrigObjects.push_back(foundObject);
    }
  }
  
  
  // Loop over muons and fill histograms: 
  int numGenPerEvent=0;
#ifdef USEGENINFO
  std::vector<const reco::GenParticle*> targetMuons;
#endif
  if (isMC_) { 
    //    for (auto& gen: *(genParticles.product())) {
    for (unsigned int g(0); g < genParticles->size(); ++g){
      const reco::GenParticle* gen = &genParticles->at(g);
      if (fabs(gen->pdgId())!=13) continue;
      if (gen->pt()<10)           continue;
      if (gen->status()!=1)       continue;
      if (fabs(gen->eta())>2.4)   continue;
      ++numGenPerEvent;
#ifdef USEGENINFO  
      targetMuons.push_back(gen);
#endif
      hists_["gen_pt"] ->Fill(gen->pt());
      hists_["gen_eta"]->Fill(gen->eta());
      hists_["gen_phi"]->Fill(gen->phi());
      
      if (debuglevel_ > 1) 
	std::cout << "gen muon found: pt: " << gen->pt() 
		  << " eta: " << gen->eta() 
		  << " phi: " << gen->phi() << std::endl;

      /// DRAWING RESOLUTION PLOTS: 
      for (int ibx = l1Muons->getFirstBX(); ibx <= l1Muons->getLastBX(); ++ibx) {
	if (ibx != 0) continue; 
	for (auto it = l1Muons->begin(ibx); it != l1Muons->end(ibx); it++){
	  l1t::MuonRef l1muon(l1Muons, distance(l1Muons->begin(l1Muons->getFirstBX()),it) );
	  if (gen->charge() > 0 && l1muon->charge()>0) {
	    hists_["hltL1p_resEta"]->Fill(l1muon->eta()-gen->eta()); 
	    hists_["hltL1p_resPhi"]->Fill(l1muon->phi()-gen->phi());
	    hists_["hltL1p_resPt"] ->Fill(l1muon->pt() -gen->pt() );
	    
	    if (fabs(gen->eta()) < 0.9) {
	      hists_["hltL1p_resEta_barrel"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1p_resPhi_barrel"]->Fill(l1muon->phi()-gen->phi());
	    }
	    if (fabs(gen->eta()) > 0.9) {
	      hists_["hltL1p_resEta_endcap"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1p_resPhi_endcap"]->Fill(l1muon->phi()-gen->phi());
	    }
	  }
	  else if (gen->charge() < 0 && l1muon->charge()<0)  {
	    hists_["hltL1m_resEta"]->Fill(l1muon->eta()-gen->eta()); 
	    hists_["hltL1m_resPhi"]->Fill(l1muon->phi()-gen->phi());
	    hists_["hltL1m_resPt"] ->Fill(l1muon->pt() -gen->pt() );
	    
	    if (fabs(gen->eta()) < 0.9) {
	      hists_["hltL1m_resEta_barrel"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1m_resPhi_barrel"]->Fill(l1muon->phi()-gen->phi());
	    }
	    if (fabs(gen->eta()) > 0.9) {
	      hists_["hltL1m_resEta_endcap"]->Fill(l1muon->eta()-gen->eta());
	      hists_["hltL1m_resPhi_endcap"]->Fill(l1muon->phi()-gen->phi());
	    }
	  }
	} //l1Muons->begin(ibx)
      } //L1 l1Muons->getFirstBX()
      
      for (unsigned int t(0); t < L2MuonTrigObjects.size(); ++t){
	trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(t);
	
	if (gen->charge() > 0) {
	  hists_["hltL2p_resEta"]->Fill(l2mu->eta()-gen->eta()); 
	  hists_["hltL2p_resPhi"]->Fill(l2mu->phi()-gen->phi());
	  hists_["hltL2p_resPt"] ->Fill(l2mu->pt() -gen->pt() );
	  
	  if (fabs(gen->eta()) < 0.9) {
	    hists_["hltL2p_resEta_barrel"]->Fill(l2mu->eta()-gen->eta());
	    hists_["hltL2p_resPhi_barrel"]->Fill(l2mu->phi()-gen->phi());
	  }
	  if (fabs(gen->eta()) > 0.9) {
	    hists_["hltL2p_resEta_endcap"]->Fill(l2mu->eta()-gen->eta());
	    hists_["hltL2p_resPhi_endcap"]->Fill(l2mu->phi()-gen->phi());
	  }
	}
	else if (gen->charge() < 0)  {
	  hists_["hltL2m_resEta"]->Fill(l2mu->eta()-gen->eta()); 
	  hists_["hltL2m_resPhi"]->Fill(l2mu->phi()-gen->phi());
	  hists_["hltL2m_resPt"] ->Fill(l2mu->pt() -gen->pt() );
	  
	  if (fabs(gen->eta()) < 0.9) {
	    hists_["hltL2m_resEta_barrel"]->Fill(l2mu->eta()-gen->eta());
	    hists_["hltL2m_resPhi_barrel"]->Fill(l2mu->phi()-gen->phi());
	  }
	  if (fabs(gen->eta()) > 0.9) {
	    hists_["hltL2m_resEta_endcap"]->Fill(l2mu->eta()-gen->eta());
	    hists_["hltL2m_resPhi_endcap"]->Fill(l2mu->phi()-gen->phi());
	  }
	}
      } // L2
    } //genParticles
  }//isMC

  if (isMC_ && numGenPerEvent==0) return; //if no st1 muon skip the event.

#ifndef USEGENINFO
  reco::MuonCollection targetMuons = selectedMuons(*muons); //, targetMuonSelector_);
#endif
  edm::Handle<reco::MuonTrackLinksCollection> links;
  iEvent.getByToken(theLinksToken_, links);
  
  vector<size_t> matchesL1 = matchByDeltaR(targetMuons,L1MuonTrigObjects,0.4); 
  vector<size_t> matchesL2 = matchByDeltaR(targetMuons,L2MuonTrigObjects,0.3); 
  vector<size_t> matchesL3 = matchByDeltaR(targetMuons,L3MuonTrigObjects,0.2); 
  
  int NumL2Matched(0), NumL3Matched(0);
  if (runSharedHits_) {
    edm::Handle<reco::TrackExtraCollection> trackMuons;
    iEvent.getByToken(theMuonsWithHitsToken_, trackMuons);
    
    reco::TrackExtraCollection trackerMuons;
    for (auto const& mu : *trackMuons) trackerMuons.push_back(mu);
    
    vector<size_t> matchesL3_hits = matchBySharedHits(trackerMuons,L3MuonTrigObjects,links, 0.1); // fill with shared hits...    
  }

  for (size_t i = 0; i < targetMuons.size(); i++) {
#ifdef USEGENINFO
    const reco::GenParticle & recoMu = *(targetMuons.at(i));
#endif
#ifndef USEGENINFO
    reco::Muon & recoMu = targetMuons[i];    
#endif
    if (debuglevel_ > 1) { 
      cout << "RECO Muon - pT, eta, phi: " << recoMu.pt() << " , " << recoMu.eta() << " , " << recoMu.phi() << endl;
    }
    
    if (matchesL2[i] < targetMuons.size()) { 
      trigger::TriggerObject & l2mu = L2MuonTrigObjects[matchesL2[i]];      
      NumL2Matched++;
      hists_["hltL2_DeltaR"]->Fill(deltaR(recoMu,l2mu));
      hists_["hltL2_pt"]    ->Fill(l2mu.pt());
      hists_["hltL2_eta"]   ->Fill(l2mu.eta());
      hists_["hltL2_phi"]   ->Fill(l2mu.phi());
      hists_["hltL2_resEta"]->Fill( recoMu.eta()-l2mu.eta());
      hists_["hltL2_resPhi"]->Fill( recoMu.phi()-l2mu.phi());
      hists_["hltL2_resPt"] ->Fill((recoMu.pt() -l2mu.pt())/recoMu.pt());
      if (matchesL1[i] < targetMuons.size()) {
	trigger::TriggerObject & l1mu = L1MuonTrigObjects[matchesL1[i]];      
	hists_["hltL2L1_resEta"]->Fill(l1mu.eta()-l2mu.eta());
	hists_["hltL2L1_resPhi"]->Fill(l1mu.phi()-l2mu.phi());
	hists_["hltL2L1_resPt"] ->Fill((l1mu.pt()-l2mu.pt())/l2mu.pt());
      }
      if (debuglevel_ > 1) { 
	cout << "L2 Muon - pT, eta, phi, DR(recomu,L2): " << l2mu.pt() << " , " << l2mu.eta() << " , "
	     << l2mu.phi() << " , " << deltaR(recoMu,l2mu) << endl;
      }
    }
    if (matchesL3[i] < targetMuons.size()) {
      trigger::TriggerObject & l3mu = L3MuonTrigObjects[matchesL3[i]];
      NumL3Matched++;
      hists_["hltL3_DeltaR"]->Fill(deltaR(recoMu,l3mu));
      hists_["hltL3_pt"]    ->Fill(l3mu.pt());
      hists_["hltL3_eta"]   ->Fill(l3mu.eta());
      hists_["hltL3_phi"]   ->Fill(l3mu.phi());
      hists_["hltL3_resEta"]->Fill( recoMu.eta()-l3mu.eta());
      hists_["hltL3_resPhi"]->Fill( recoMu.phi()-l3mu.phi());
      hists_["hltL3_resPt"] ->Fill((recoMu.pt() -l3mu.pt())/recoMu.pt());
      if (debuglevel_ > 1) { 
//B	cout << "L3 Muon - pT, eta, phi, DR(mu,L3): " << l3mu.pt() << " , " << l3mu.eta() << " , "
//B	     << l3mu.phi() << " , " << deltaR(recoMu,l3mu) << endl;
      }
    }    
  }
  
  hists_["hlt_NumL1" ]      ->Fill(L1MuonTrigObjects.size());
  hists_["hlt_NumL2" ]      ->Fill(L2MuonTrigObjects.size());
  hists_["hlt_NumL3" ]      ->Fill(L3MuonTrigObjects.size());
  hists_["hlt_NumL1Match" ] ->Fill(matchesL1.size());
  hists_["hlt_NumL2Match" ] ->Fill(NumL2Matched);
  hists_["hlt_NumL3Match" ] ->Fill(NumL3Matched);
  if (L2MuonTrigObjects.size()>0) hists_["hlt_FracL2Match" ]->Fill((float)NumL2Matched/(float)L2MuonTrigObjects.size());
  if (L3MuonTrigObjects.size()>0) hists_["hlt_FracL3Match" ]->Fill((float)NumL3Matched/(float)L3MuonTrigObjects.size());

  if (debuglevel_ > 1) {
  //B  std::cout << "Number of L1s passing filter = " << L1MuonTrigObjects.size() << " ( " << matchesL1.size() << " ) " << std::endl;
  //B  std::cout << "Number of L2s passing filter = " << L2MuonTrigObjects.size() << " ( " << NumL2Matched << " ) " << std::endl;
  //B  std::cout << "Number of L3s passing filter = " << L3MuonTrigObjects.size() << " ( " << NumL3Matched << " ) " << std::endl;
  }
  
 // if (debuglevel_ > 1) {
 //   if (L2MuonTrigObjects.size() > 0 && L3MuonTrigObjects.size()!=L2MuonTrigObjects.size())   cout << "[FAILING EVENT - IO] Run:Event --> " << iEvent.id().run() << " : " << iEvent.id().event() << endl;
 // }
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ////////          EVENT DEBUGGING     
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
//  edm::Handle<reco::TrackCollection> hltL3OITracks;
//  iEvent.getByToken(theOITracksToken_, hltL3OITracks);    
//  hists_["hlt_L3OI_numTracks"]->Fill(hltL3OITracks->size());
//
//  edm::Handle<reco::TrackCollection> hltL3IOTracks;
//  iEvent.getByToken(theIOTracksToken_, hltL3IOTracks);    
//  hists_["hlt_L3IO_numTracks"]->Fill(hltL3IOTracks->size());
//
//  edm::Handle<reco::TrackCollection> hltL3IOFromL1Tracks;
//  iEvent.getByToken(theIOFromL1TracksToken_, hltL3IOFromL1Tracks);    
//  hists_["hlt_L3IOFromL1_numTracks"]->Fill(hltL3IOFromL1Tracks->size());
  

  /// CHECKING THE INSIDE-OUT SEQUENCE: 
  /*
  try { 
    if (NumL2Matched>0){ 
      if (NumL3Matched>0) {
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter0;
	iEvent.getByToken(theSeedsIter0Token_, hltL3TrajSeedIter0);
	hists_["hlt_numSeedsIter0"]->Fill(hltL3TrajSeedIter0->size());
	hists_["hlt_numSeedsIter0_PU"]->Fill(nbMCvtx, hltL3TrajSeedIter0->size());
      
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter2;
	iEvent.getByToken(theSeedsIter2Token_, hltL3TrajSeedIter2);
	hists_["hlt_numSeedsIter2"]->Fill(hltL3TrajSeedIter2->size());
	hists_["hlt_numSeedsIter2_PU"]->Fill(nbMCvtx, hltL3TrajSeedIter2->size());	
	hists_["hlt_numSeeds"]->Fill(hltL3TrajSeedIter0->size()+hltL3TrajSeedIter2->size());
	hists_["hlt_numSeeds_PU"]->Fill(nbMCvtx, hltL3TrajSeedIter0->size()+hltL3TrajSeedIter2->size());
      
	// now for the number of tracks
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter0;
	iEvent.getByToken(theTracksIter0Token_, hltL3TkTracksIter0);
	hists_["hlt_numTracksIter0"]->Fill(hltL3TkTracksIter0->size());
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter0->begin(); t!=hltL3TkTracksIter0->end(); t++) {
	  hists_["hlt_numPixelHitsIter0"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
      
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter2;
	iEvent.getByToken(theTracksIter2Token_, hltL3TkTracksIter2);
	hists_["hlt_numTracksIter2"]->Fill(hltL3TkTracksIter2->size());
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter2->begin(); t!=hltL3TkTracksIter2->end(); t++) {
	  hists_["hlt_numPixelHitsIter2"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
      
	edm::Handle<reco::TrackCollection> hltL3TkTracks;
	iEvent.getByToken(theTracksToken_, hltL3TkTracks);    
	hists_["hlt_numTracks"]->Fill(hltL3TkTracks->size());
	for (reco::TrackCollection::const_iterator t=hltL3TkTracks->begin(); t!=hltL3TkTracks->end(); t++) {
	  hists_["hlt_numPixelHits"]->Fill(t->hitPattern().numberOfValidPixelHits());
	  hists_["hlt_normalizedChi2"]->Fill(t->normalizedChi2()); 
	  hists_["hlt_numValidHits"]->Fill(t->numberOfValidHits());	
	  hists_["hlt_numValidMuonHits"]->Fill(t->hitPattern().numberOfValidMuonHits());
	}
	edm::Handle<reco::MuonTrackLinksCollection> IOlink;
	iEvent.getByToken(theIOLinksToken_, IOlink);
      
	for(reco::RecoChargedCandidateCollection::const_iterator cand=l3Muons->begin(); cand!=l3Muons->end(); cand++){
	  for(unsigned int l(0); l <IOlink->size(); ++l){
	    const reco::MuonTrackLinks* link = &IOlink->at(l);
	    bool useThisLink=false;
	    reco::TrackRef tk = cand->track();
	    reco::TrackRef trkTrack = link->trackerTrack();
	    // Using the same method that was used to create the links
	    // ToDo: there should be a better way than dR,dPt matching
	    const reco::Track& globalTrack = *link->globalTrack();
	    float dR2 = deltaR2(tk->eta(),tk->phi(),globalTrack.eta(),globalTrack.phi());
	    float dPt = std::abs(tk->pt() - globalTrack.pt())/tk->pt();
	  
	    hists_["hltL3mu_dR2withLink"]->Fill(dR2);
	    hists_["hltL3mu_dPtwithLink"]->Fill(dPt);
	  
	    if (dR2 < 0.02*0.02 and dPt < 0.001) {
	      useThisLink=true;
	    }
	  
	    if (useThisLink) {
	      const reco::TrackRef staTrack = link->standAloneTrack();
	      if (abs(cand->eta()) < 1.2) continue;
	      hists_["hltL3mu_pt"]->Fill(cand->pt());
	      hists_["hltL3mu_eta"]->Fill(cand->eta());
	      hists_["hltL3mu_numValidHits"]->Fill(tk->numberOfValidHits());
	      //	    hists_["hltL3mu_dr"]->Fill(std::abs( (- (cand->vx()-beamSpot.x0()) * cand->py() + (cand->vy()-beamSpot.y0()) * cand->px() ) / cand->pt() ));
	      ///	    hists_["hltL3mu_dz"]->Fill(std::abs((cand->vz()-beamSpot.z0()) - ((cand->vx()-beamSpot.x0())*cand->px()+(cand->vy()-beamSpot.y0())*cand->py())/cand->pt() * cand->pz()/cand->pt()));
	      //	    hists_["hltL3mu_dxySig"]->Fill(std::abs(tk->dxy(beamSpot.position())/tk->dxyError()));
	      hists_["hltL3mu_normalizedChi2"]->Fill(tk->normalizedChi2());
	      //	    hists_["hltL3mu_dxy"]->Fill(std::abs(tk->dxy(beamSpot.position())));
	      hists_["hltL3mu_numValidMuonHits"]->Fill(tk->hitPattern().numberOfValidMuonHits());
	    } //end of useThisLink
	  } //end of muons in links collection
	} //end of RecoCand collection
      }
      else { 
	trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(0);  
	if (debuglevel_ > 1)   cout << "[FAILING EVENT - IO] Run:Event --> " << iEvent.id().run() << " : " << iEvent.id().event() << endl;
	
	/// Store some information of such events failing the L3 reconstruction: 
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter0;
	iEvent.getByToken(theSeedsIter0Token_, hltL3TrajSeedIter0);
	hists_["hlt_noL3_numSeedsIter0"]->Fill(hltL3TrajSeedIter0->size());
	
	edm::Handle<reco::TrackCollection> hltL3PixelTracsIter0;
	iEvent.getByToken(thePixelTracksIter0Token_, hltL3PixelTracsIter0);
	hists_["hlt_noL3_numPixelTracksIter0"]->Fill(hltL3TrajSeedIter0->size());
	
	edm::Handle<TrajectorySeedCollection> hltL3TrajSeedIter2;
	iEvent.getByToken(theSeedsIter2Token_, hltL3TrajSeedIter2);
	hists_["hlt_noL3_numSeedsIter2"]->Fill(hltL3TrajSeedIter2->size());
	
	// now for the number of tracks
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter0;
	iEvent.getByToken(theTracksIter0Token_, hltL3TkTracksIter0);
	hists_["hlt_noL3_numTracksIter0"]->Fill(hltL3TkTracksIter0->size());
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter0->begin(); t!=hltL3TkTracksIter0->end(); t++) {
	  hists_["hlt_noL3_numPixelHitsIter0"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
	
	edm::Handle<reco::TrackCollection> hltL3TkTracksIter2;
	iEvent.getByToken(theTracksIter2Token_, hltL3TkTracksIter2);
	hists_["hlt_noL3_numTracksIter2"]->Fill(hltL3TkTracksIter2->size());
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter2->begin(); t!=hltL3TkTracksIter2->end(); t++) {
	  hists_["hlt_noL3_numPixelHitsIter2"]->Fill(t->hitPattern().numberOfValidPixelHits());
	}
	
	edm::Handle<reco::TrackCollection> hltL3TkTracks;
	iEvent.getByToken(theTracksToken_, hltL3TkTracks);    
	hists_["hlt_noL3_numTracks"]->Fill(hltL3TkTracks->size());
	
	edm::Handle<TrackCandidateCollection> hltL3TkTracksCandIter0;
	iEvent.getByToken(theTracksCandIter0Token_, hltL3TkTracksCandIter0);
	hists_["hlt_noL3_numTracksCandIter0"]->Fill(hltL3TkTracksCandIter0->size());
	
	edm::Handle<reco::TrackCollection> hltL3TkTracksNoHPIter0;
	iEvent.getByToken(theTracksNoHPIter0Token_, hltL3TkTracksNoHPIter0);
	hists_["hlt_noL3_numTracksNoHPIter0"]->Fill(hltL3TkTracksCandIter0->size());
	
	edm::Handle<TrackCandidateCollection> hltL3TkTracksCandIter2;
	iEvent.getByToken(theTracksCandIter2Token_, hltL3TkTracksCandIter2);
	hists_["hlt_noL3_numTracksCandIter2"]->Fill(hltL3TkTracksCandIter2->size());

	edm::Handle<reco::TrackCollection> hltL3TkTracksNoHPIter2;
	iEvent.getByToken(theTracksNoHPIter2Token_, hltL3TkTracksNoHPIter2);
	hists_["hlt_noL3_numTracksNoHPIter2"]->Fill(hltL3TkTracksCandIter2->size());

	if (debuglevel_ > 1) { 
	  cout << "# of hltL3TrajSeedIter0: "  << hltL3TrajSeedIter0->size() << endl;
	  cout << "# of hltL3TkTracksCandIter0 = " << hltL3TkTracksCandIter0->size() << endl; 
	  cout << "# of hltL3TkTracksNoHPIter0 = " << hltL3TkTracksNoHPIter0->size() << endl; 
	}   
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksNoHPIter0->begin(); t!=hltL3TkTracksNoHPIter0->end(); t++) {
	  if (debuglevel_ > 1) cout << "    hltL3TkTracksNoHPIter0 " << t->pt() << " " << t->eta() << " " << t->phi() << " -> " 
				    << t->hitPattern().numberOfValidPixelHits() << " " << t->numberOfValidHits() << " " << t->hitPattern().numberOfValidMuonHits()  << " " 
				    << deltaR(*t,*l2mu) << " " << t->eta()-l2mu->eta() << " " << t->phi()-l2mu->phi() << endl;
	}

	if (debuglevel_ > 1)  cout << "# of hltL3TkTracksIter0 = " << hltL3TkTracksIter0->size() << endl;    
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter0->begin(); t!=hltL3TkTracksIter0->end(); t++) {
	  if (debuglevel_ > 1)  cout << "    hltL3TkTracksIter0 " << t->pt() << " " << t->eta() << " " << t->phi() << " -> " 
				     << t->hitPattern().numberOfValidPixelHits() << " " << t->numberOfValidHits() << " " << t->hitPattern().numberOfValidMuonHits()  << " " 
				     << deltaR(*t,*l2mu) << " " << t->eta()-l2mu->eta() << " " << t->phi()-l2mu->phi() << endl;
	}
	
	if (debuglevel_ > 1) {
	  cout << "# of hltL3TrajSeedIter2: "  << hltL3TrajSeedIter2->size() << endl;
	  cout << "# of hltL3TkTracksCandIter2 = " << hltL3TkTracksCandIter2->size() << endl; 
	  cout << "# of hltL3TkTracksNoHPIter2 = " << hltL3TkTracksNoHPIter2->size() << endl;    
	}
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksNoHPIter2->begin(); t!=hltL3TkTracksNoHPIter2->end(); t++) {
	  if (debuglevel_ > 1)  cout << "    hltL3TkTracksNoHPIter2 " << t->pt() << " " << t->eta() << " " << t->phi() << " -> " 
				     << t->hitPattern().numberOfValidPixelHits() << " " << t->numberOfValidHits() << " " << t->hitPattern().numberOfValidMuonHits()  << " " 
				     << deltaR(*t,*l2mu) << " " << t->eta()-l2mu->eta() << " " << t->phi()-l2mu->phi() << endl;
	}
	if (debuglevel_ > 1)   cout << "# of hltL3TkTracksIter2 = " << hltL3TkTracksIter2->size() << endl;    
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksIter2->begin(); t!=hltL3TkTracksIter2->end(); t++) {
	  if (debuglevel_ > 1)  cout << "    hltL3TkTracksIter2 " << t->pt() << " " << t->eta() << " " << t->phi() << " -> " 
				     << t->hitPattern().numberOfValidPixelHits() << " " << t->numberOfValidHits() << " " << t->hitPattern().numberOfValidMuonHits()  << " " 
				     << deltaR(*t,*l2mu) << " " << t->eta()-l2mu->eta() << " " << t->phi()-l2mu->phi() << endl;
	}	    

	for (reco::TrackCollection::const_iterator t=hltL3TkTracks->begin(); t!=hltL3TkTracks->end(); t++) {
	  if (debuglevel_ > 1)  cout << "    hltL3TkTracks " << t->pt() << " " << t->eta() << " " << t->phi() << " -> " 
				     << t->hitPattern().numberOfValidPixelHits() << " " << t->numberOfValidHits() << " " << t->hitPattern().numberOfValidMuonHits() 
				     << deltaR(*t,*l2mu) << " " << t->eta()-l2mu->eta() << " " << t->phi()-l2mu->phi() << endl;
	  
	    hists_["hlt_noL3_numPixelHits"]    ->Fill(t->hitPattern().numberOfValidPixelHits());
	    hists_["hlt_noL3_normalizedChi2"]  ->Fill(t->normalizedChi2()); 
	    hists_["hlt_noL3_numValidHits"]    ->Fill(t->numberOfValidHits());
	    hists_["hlt_noL3_numValidMuonHits"]->Fill(t->hitPattern().numberOfValidMuonHits());
	    hists_["hlt_noL3_trackpt"]         ->Fill(t->pt());
	    hists_["hlt_noL3_tracketa"]        ->Fill(t->eta());
	    hists_["hlt_noL3_trackphi"]        ->Fill(t->phi());
	    hists_["hlt_noL3_DeltaR"]          ->Fill(deltaR(*t,*l2mu));
	    hists_["hlt_noL3_DeltaEta"]        ->Fill(t->eta()-l2mu->eta());
	    hists_["hlt_noL3_DeltaPhi"]        ->Fill(t->phi()-l2mu->phi());
	    hists_["hlt_noL3_DeltaPt"]         ->Fill((t->pt()-l2mu->pt())/l2mu->pt());
	}
      } // NO L3 Objects
    }//L2Objects
  }
  catch (...) {
  }
  */
  
  /// CASCADE!!! 
  try {	
    if (l2Muons->size()>0) {
      edm::Handle<L3MuonTrajectorySeedCollection> hltL3TrajSeedOIS;
      iEvent.getByToken(theSeedsOISToken_, hltL3TrajSeedOIS);
      
      edm::Handle<L3MuonTrajectorySeedCollection> hltL3TrajSeedOIH;
      iEvent.getByToken(theSeedsOIHToken_, hltL3TrajSeedOIH);
      
      hists_["hlt_OI_NumSeeds"]->Fill((hltL3TrajSeedOIS->size()+hltL3TrajSeedOIH->size())/l2Muons->size());
      for (unsigned int t(0); t < l2Muons->size(); ++t){
       	reco::RecoChargedCandidateRef candref(l2Muons, t);
       	hists_["hlt_OI_NumSeedsVsEta"]->Fill(candref->eta(),(hltL3TrajSeedOIS->size()+hltL3TrajSeedOIH->size())/l2Muons->size());
      }
      
      for (L3MuonTrajectorySeedCollection::const_iterator seed=hltL3TrajSeedOIS->begin(); seed!=hltL3TrajSeedOIS->end();seed++){
	PTrajectoryStateOnDet ptod = seed->startingState();
	DetId id(ptod.detId());
	const GeomDet * g = theService->trackingGeometry()->idToDet( id );
	const Surface * surface=&g->surface();
	
	TrajectoryStateOnSurface seedTSOS = trajectoryStateTransform::transientState(ptod, surface,  &*(theService)->magneticField());
	AlgebraicSymMatrix66 errors = seedTSOS.cartesianError().matrix();
	double partialPterror = errors(3,3)*pow(seedTSOS.globalMomentum().x(),2) + errors(4,4)*pow(seedTSOS.globalMomentum().y(),2);
	float seedeta = seedTSOS.globalMomentum().eta(); 
	float seedphi = seedTSOS.globalMomentum().phi(); 
	/// Match to L2? 
	float l2eta  = 9999.;
	float l2phi  = 9999.;
	float mindR = 9999.;
	for (unsigned int t(0); t < l2Muons->size(); ++t){
	  reco::RecoChargedCandidateRef l2(l2Muons, t);
	  float deta = l2->eta() - seedeta;
	  float dphi = l2->phi() - seedphi;
	  float dist = sqrt(deta*deta+dphi*dphi);
	  if (dist < mindR) {
	    mindR = dist; 
	    l2eta = l2->eta();
	    l2phi = l2->phi();
	  }
	}
	hists_["hlt_OI_seedEtaVsL2Eta"]       ->Fill(seedeta,seedeta-l2eta);
	hists_["hlt_OI_seedPhiVsL2Phi"]       ->Fill(seedeta,seedphi-l2phi);
	//  End match to L2
	
	hists_["hlt_OI_seedEta"]       ->Fill(seedeta);
	hists_["hlt_OI_seedPtErrVsEta"]->Fill(seedeta,sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	hists_["hlt_OI_seedPtErr"]     ->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	hists_["hlt_OI_seedPhiErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	hists_["hlt_OI_seedEtaErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	if (std::abs(seedeta) < 0.9) {
	  hists_["hlt_OI_seedPtErr_barrel"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPhiErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	  hists_["hlt_OI_seedEtaErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	}
	else { 
	  hists_["hlt_OI_seedPtErr_endcap"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPhiErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	  hists_["hlt_OI_seedEtaErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	}	    
//B	if (debuglevel_ > 1) { 
//B	  cout << " TSOS of OIState     = \n"
//B	       <<"x: "<<seedTSOS.globalPosition()<< " --> " << seedTSOS.localPosition() << "\n" 
//B	       <<"p: "<<seedTSOS.globalMomentum()<< "\n"
//B	       <<"pt: " << seedTSOS.globalMomentum().perp() << " +/- " << sqrt(partialPterror)/seedTSOS.globalMomentum().perp() << "\n"
//B	       <<"eta: " << seedTSOS.globalMomentum().eta() << " +/- " << sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())) << "\n"
//B	       <<"phi: " << seedTSOS.globalMomentum().phi() << " +/- " << sqrt(seedTSOS.curvilinearError().matrix()(2,2)) << "\n"
//B	       << id.subdetId() << " " << id.rawId() << endl;
	  
//B	  TrajectorySeed::range seedHits = seed->recHits();
//B	  for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){
//B	    cout << " hits of hltL3TrajSeedOI -->  x: " << (*iseed).globalPosition() << " - " << (*iseed).localPosition() << "\n"
//B		 << (*iseed).geographicalId().subdetId() << " " << (*iseed).geographicalId().rawId() << endl;
//B	  }

//B	}
      }
      for (L3MuonTrajectorySeedCollection::const_iterator seed=hltL3TrajSeedOIH->begin(); seed!=hltL3TrajSeedOIH->end();seed++){
	PTrajectoryStateOnDet ptod = seed->startingState();
	DetId id(ptod.detId());
	const GeomDet * g = theService->trackingGeometry()->idToDet( id );
	const Surface * surface=&g->surface();

	TrajectoryStateOnSurface seedTSOS = trajectoryStateTransform::transientState(ptod, surface,  &*(theService)->magneticField());
	
	TrajectorySeed::range seedHits = seed->recHits();
	hists_["hlt_OI_NumberOfRecHitsPerSeed"]->Fill(seed->nHits());	    
	
	AlgebraicSymMatrix66 errors = seedTSOS.cartesianError().matrix();
	double partialPterror = errors(3,3)*pow(seedTSOS.globalMomentum().x(),2) + errors(4,4)*pow(seedTSOS.globalMomentum().y(),2);
	float seedeta = seedTSOS.globalMomentum().eta(); 
	float seedphi = seedTSOS.globalMomentum().phi(); 
	/// Match to L2? 
	float l2eta  = 9999.;
	float l2phi  = 9999.;
	float mindR = 9999.;
	for (unsigned int t(0); t < l2Muons->size(); ++t){
	  reco::RecoChargedCandidateRef l2(l2Muons, t);
	  float deta = l2->eta() - seedeta;
	  float dphi = l2->phi() - seedphi;
	  float dist = sqrt(deta*deta+dphi*dphi);
	  if (dist < mindR) {
	    mindR = dist; 
	    l2eta = l2->eta();
	    l2phi = l2->phi();
	  }
	}
	hists_["hlt_OI_seedEtaVsL2Eta"]       ->Fill(seedeta,seedeta-l2eta);
	hists_["hlt_OI_seedPhiVsL2Phi"]       ->Fill(seedeta,seedphi-l2phi);
       	//  End match to L2

	hists_["hlt_OI_seedEta"]       ->Fill(seedeta);
	hists_["hlt_OI_seedPtErrVsEta"]->Fill(seedeta,sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	hists_["hlt_OI_seedPtErr"]     ->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	hists_["hlt_OI_seedPhiErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	hists_["hlt_OI_seedEtaErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	if (std::abs(seedeta) < 0.9) {
	  hists_["hlt_OI_seedPtErr_barrel"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPhiErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	  hists_["hlt_OI_seedEtaErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	}
	else { 
	  hists_["hlt_OI_seedPtErr_endcap"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPhiErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	  hists_["hlt_OI_seedEtaErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	}	    
	
	for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){
	  if (!(*iseed).isValid()) continue;
	  hists_["hlt_OI_hitPt"]->Fill((*iseed).globalPosition().perp());
	  hists_["hlt_OI_hitPhi"]->Fill((*iseed).globalPosition().phi());
	  hists_["hlt_OI_hitEta"]->Fill((*iseed).globalPosition().eta());
	  hists_["hlt_OI_hitx"]->Fill((*iseed).globalPosition().x());
	  hists_["hlt_OI_hity"]->Fill((*iseed).globalPosition().y());
	  hists_["hlt_OI_hitz"]->Fill((*iseed).globalPosition().z());
	}

//B	if (debuglevel_ > 1) { 
//B	  cout << " TSOS of OIHit     = \n"
//B	       <<"x: "<<seedTSOS.globalPosition()<< " --> " << seedTSOS.localPosition() << "\n" 
//B	       <<"p: "<<seedTSOS.globalMomentum()<< "\n"
//B	       <<"pt: " << seedTSOS.globalMomentum().perp() << " +/- " << sqrt(partialPterror)/seedTSOS.globalMomentum().perp() << "\n"
//B	       <<"eta: " << seedTSOS.globalMomentum().eta() << " +/- " << sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())) << "\n"
//B	       <<"phi: " << seedTSOS.globalMomentum().phi() << " +/- " << sqrt(seedTSOS.curvilinearError().matrix()(2,2)) << "\n"
//B	       << id.subdetId() << " " << id.rawId() << endl;
	  
//B	  TrajectorySeed::range seedHits = seed->recHits();
//B	  for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){
//B	    cout << " hits of hltL3TrajSeedOI -->  x: " << (*iseed).globalPosition() << " - " << (*iseed).localPosition() << "\n"
//B		 << (*iseed).geographicalId().subdetId() << " " << (*iseed).geographicalId().rawId() << endl;
//B	  }
//B	}
      }
    }
  }
  catch (...) {
  }

  /*
    try {
    if (NumL2Matched>0) {
    trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(0);  
      
    edm::Handle<L3MuonTrajectorySeedCollection> hltL3TrajSeedOIH;
    iEvent.getByToken(theSeedsOIHToken_, hltL3TrajSeedOIH);
    edm::Handle<TrackCandidateCollection> hltL3TkTracksCandOIH;
    iEvent.getByToken(theTracksCandOIHToken_, hltL3TkTracksCandOIH);
    edm::Handle<reco::TrackCollection> hltL3TkTracksOIH;
      iEvent.getByToken(theTracksOIHToken_, hltL3TkTracksOIH);
      
      if (debuglevel_ > 1) {
	cout << "# of hltL3TrajSeedOIH     = " << hltL3TrajSeedOIH->size() << endl;
	for (L3MuonTrajectorySeedCollection::const_iterator seed=hltL3TrajSeedOIH->begin(); seed!=hltL3TrajSeedOIH->end();seed++){
	  PTrajectoryStateOnDet ptod = seed->startingState();
	  DetId id(ptod.detId());
	  const GeomDet * g = theService->trackingGeometry()->idToDet( id );
	  const Surface * surface=&g->surface();
	  TrajectoryStateOnSurface currentState(trajectoryStateTransform::transientState(ptod, surface,  &*(theService)->magneticField()));
	  cout << " TSOS of hltL3TrajSeedOIH     = \n"
	       <<"x: "<<currentState.globalPosition()<< " --> " << currentState.localPosition() << "\n" 
	       <<"p: "<<currentState.globalMomentum()<<"\n" 
	       << id.subdetId() << " " << id.rawId() << endl;
	  TrajectorySeed::range seedHits = seed->recHits();
	  for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){
	    cout << " hits of hltL3TrajSeedOIH -->  x: " << (*iseed).globalPosition() << " - " << (*iseed).localPosition() << "\n"
		 << (*iseed).geographicalId().subdetId() << " " << (*iseed).geographicalId().rawId() << endl;
	  }
	}
	cout << "# of hltL3TkTracksCandOIH = " << hltL3TkTracksCandOIH->size() << endl; 
	cout << "# of hltL3TkTracksOIH     = " << hltL3TkTracksOIH->size() << endl;    
	for (reco::TrackCollection::const_iterator t=hltL3TkTracksOIH->begin(); t!=hltL3TkTracksOIH->end(); t++) {
	  cout << "    hltL3TkTracksOIH " << t->pt() << " " << t->eta() << " " << t->phi() << " -> " 
	       << t->normalizedChi2() << " " << t->hitPattern().numberOfValidPixelHits() << " " << t->numberOfValidHits() << " " << t->hitPattern().numberOfValidMuonHits()  << " " 
	       << deltaR(*t,*l2mu) << " " << t->eta()-l2mu->eta() << " " << t->phi()-l2mu->phi() << endl;
	}
      }
    }
  }
  catch (...) {
  }
  
  if (NumL2Matched>0 && NumL3Matched==0 && debuglevel_ > 1) {
    trigger::TriggerObject* l2mu = &L2MuonTrigObjects.at(0);
    if (abs(l2mu->eta()) > 2.0)     cout << "[FAILING EVENT] Run:Event --> " << iEvent.id().run() << " : " << iEvent.id().event() << endl;
  }
  */
  
  
  /// DEBUGGING THE OUTSIDE-IN COMPONENT   (ITER-L3) 
  try {
    if (l2Muons->size()>0) { 
      edm::Handle<TrajectorySeedCollection> hltL3TrajSeedOI;
      iEvent.getByToken(theSeedsOIToken_, hltL3TrajSeedOI);
      hists_["hlt_OI_NumSeeds"]->Fill(hltL3TrajSeedOI->size()/l2Muons->size());
      
      for (unsigned int t(0); t < l2Muons->size(); ++t){
	reco::RecoChargedCandidateRef candref(l2Muons, t);
	hists_["hlt_OI_NumSeedsVsEta"]->Fill(candref->eta(),hltL3TrajSeedOI->size()/l2Muons->size());
      }
      
      for(TrajectorySeedCollection::const_iterator seed = hltL3TrajSeedOI->begin(); seed != hltL3TrajSeedOI->end(); ++seed){
	PTrajectoryStateOnDet ptod = seed->startingState();
	DetId id(ptod.detId());
	const GeomDet * g = theService->trackingGeometry()->idToDet( id );
	const Surface * surface=&g->surface();
	TrajectoryStateOnSurface seedTSOS = trajectoryStateTransform::transientState(ptod, surface,  &*(theService)->magneticField());
	
	TrajectorySeed::range seedHits = seed->recHits();
	hists_["hlt_OI_NumberOfRecHitsPerSeed"]->Fill(seed->nHits());	    
       //std::cout<<"Seed TSOS position(x,y,z): "<<seedTSOS.globalPosition().x()<<", "<<seedTSOS.globalPosition().y()<<", "<<seedTSOS.globalPosition().z()<<std::endl;
       //std::cout<<"Seed TSOS momentum:(px,py,pz) :   "<<seedTSOS.globalMomentum().x()<<", "<<seedTSOS.globalMomentum().y()<<", "<<seedTSOS.globalMomentum().z()<<std::endl;



	
	AlgebraicSymMatrix66 errors = seedTSOS.cartesianError().matrix();
	double partialPterror = errors(3,3)*pow(seedTSOS.globalMomentum().x(),2) + errors(4,4)*pow(seedTSOS.globalMomentum().y(),2);
	float seedeta = seedTSOS.globalMomentum().eta(); 
	float seedphi = seedTSOS.globalMomentum().phi(); 
	/// Match to L2? 
	float l2eta  = 9999.;
	float l2phi  = 9999.;
	float mindR = 9999.;
	for (unsigned int t(0); t < l2Muons->size(); ++t){
	  reco::RecoChargedCandidateRef l2(l2Muons, t);
	  float deta = l2->eta() - seedeta;
	  float dphi = l2->phi() - seedphi;
	  float dist = sqrt(deta*deta+dphi*dphi);
	  if (dist < mindR) {
	    mindR = dist; 
	    l2eta = l2->eta();
	    l2phi = l2->phi();
	  }
	}
	hists_["hlt_OI_seedEtaVsL2Eta"]       ->Fill(seedeta,seedeta-l2eta);
	hists_["hlt_OI_seedPhiVsL2Phi"]       ->Fill(seedeta,seedphi-l2phi);
       	//  End match to L2

	hists_["hlt_OI_seedEta"]       ->Fill(seedeta);
	hists_["hlt_OI_seedPtErrVsEta"]->Fill(seedeta,sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	hists_["hlt_OI_seedPtErr"]     ->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	hists_["hlt_OI_seedPhiErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	hists_["hlt_OI_seedEtaErr"]    ->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	if (std::abs(seedeta) < 0.9) {
	  hists_["hlt_OI_seedPtErr_barrel"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPhiErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	  hists_["hlt_OI_seedEtaErr_barrel"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	}
	else { 
	  hists_["hlt_OI_seedPtErr_endcap"]->Fill(sqrt(partialPterror)/seedTSOS.globalMomentum().perp());
	  hists_["hlt_OI_seedPhiErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(2,2)));
	  hists_["hlt_OI_seedEtaErr_endcap"]->Fill(sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta())));
	}	    
	
	for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){
	  if (!(*iseed).isValid()) continue;
	  hists_["hlt_OI_hitPt"]->Fill((*iseed).globalPosition().perp());
	  hists_["hlt_OI_hitPhi"]->Fill((*iseed).globalPosition().phi());
	  hists_["hlt_OI_hitEta"]->Fill((*iseed).globalPosition().eta());
	  hists_["hlt_OI_hitx"]->Fill((*iseed).globalPosition().x());
	  hists_["hlt_OI_hity"]->Fill((*iseed).globalPosition().y());
	  hists_["hlt_OI_hitz"]->Fill((*iseed).globalPosition().z());
	}
      
//	if (debuglevel_ >1) {
	 // cout << " TSOS of hltL3TrajSeedOI     = \n"
	 //      <<"x: "<<seedTSOS.globalPosition()<< " --> " << seedTSOS.localPosition() << "\n" 
	 //      <<"p: "<<seedTSOS.globalMomentum()<< "\n"
	 //      <<"pt: " << seedTSOS.globalMomentum().perp() << " +/- " << sqrt(partialPterror)/seedTSOS.globalMomentum().perp()  << "\n"
	 //      <<"eta: " << seedTSOS.globalMomentum().eta() << " +/- " << sqrt(seedTSOS.curvilinearError().matrix()(1,1))*abs(sin(seedTSOS.globalMomentum().theta()))  << "\n"
	 //      <<"phi: " << seedTSOS.globalMomentum().phi() << " +/- " << sqrt(seedTSOS.curvilinearError().matrix()(2,2)) << "\n"
	 //      << id.subdetId() << " " << id.rawId() << endl;
	  
//	  for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){

  //          cout<<"looping over the seed hits .."<<endl;

//	    cout << " hits of hltL3TrajSeedOI -->  global position: " << (*iseed).globalPosition() <<endl;// " - " << (*iseed).localPosition() << "\n"
   
  //         cout<<"Subdet ID: "<< (*iseed).geographicalId().subdetId() << ":   rawID:   " << (*iseed).geographicalId().rawId() << endl;
//	  }


//	}


      }

   edm::Handle<reco::TrackExtraCollection> TrackCollectionGM;
   iEvent.getByToken(TrackCollectionToken_, TrackCollectionGM);  

   cout<<"Reco Track Size:  "<<TrackCollectionGM->size()<<endl;


   for (reco::TrackExtraCollection::const_iterator trk = TrackCollectionGM->begin(); trk!=TrackCollectionGM->end(); ++trk){//looping over reco Tracks
	  reco::TrackExtra track = (*trk);

          cout<<"testing reco tracks"<<endl;

          for(TrajectorySeedCollection::const_iterator seed = hltL3TrajSeedOI->begin(); seed != hltL3TrajSeedOI->end(); ++seed){//loop seeds
          cout<<"test2"<<endl;

	    TrajectorySeed::range seedHits = seed->recHits();
           cout<<"test3"<<endl;

            PTrajectoryStateOnDet pTSOD = seed->startingState();
          cout<<"test4"<<endl;

	    DetId seedDetId(pTSOD.detId());
         cout<<"test5"<<endl;

	    const GeomDet* gdet = theService->trackingGeometry()->idToDet( seedDetId );

           cout<<"test6"<<endl;

	    TrajectoryStateOnSurface seedTSOS = trajectoryStateTransform::transientState(pTSOD, &(gdet->surface()), &*(theService)->magneticField());

           cout<<"test7"<<endl;

            for ( TrajectorySeed::const_iterator iseed=seedHits.first; iseed!=seedHits.second; ++iseed){//seed hit
	      if (!(*iseed).isValid()) continue; 


             cout<<"test8"<<endl;
             float deltax(9999.);
             float deltay(9999.);
             float deltaz(9999.);

             DetId detidSeed = (*iseed).geographicalId();
	     int subDetSeed = detidSeed.subdetId();

             
           for (trackingRecHit_iterator hitIt = track.recHitsBegin(); hitIt != track.recHitsEnd(); ++hitIt) {//track rechit
		if (not (*hitIt)->isValid()) continue;

                  cout<<"test9"<<endl;
		const TrackingRecHit* rechit = (*hitIt)->hit();
		DetId detid = rechit->geographicalId();
		int subDet = detid.subdetId();

              cout<<"test10"<<endl;
                if (detid.det() != detidSeed.det()) continue;

                cout<<"test11"<<endl;
		if (subDetSeed != subDet) continue;

                cout<<"test12"<<endl;
		if (detid.rawId() != detidSeed.rawId()) continue; // SAME layer, rod, module (TOB) or wheel, petal, ring, module (TEC) detector.
                 cout<<"test13"<<endl;
                deltax = (*iseed).localPosition().x() - rechit->localPosition().x();
		deltay = (*iseed).localPosition().y() - rechit->localPosition().y();
  		deltaz = (*iseed).localPosition().z() - rechit->localPosition().z();

               cout<<"delX: "<<deltax<<endl;
               cout<<"delY: "<<deltay<<endl;
               cout<<"delZ: "<<deltaz<<endl;



                }///track rechit loop







         }//seedhits loop
 

      }//loop seeds

     }//lopping over recoTracks

    }
  }
  catch (...) {
  }

  /// NOW for TRACKS
  try {
    if (l2Muons->size()>0) { 
      edm::Handle<TrackCandidateCollection> hltL3TkTracksCandOI;
      iEvent.getByToken(theTracksOICandToken_, hltL3TkTracksCandOI);
      if (debuglevel_ > 1) {
	cout << "---------------------------------------------" << endl;
	cout << "# of hltL3TkTracksCandOI(Trackcandidates)  = " << hltL3TkTracksCandOI->size() << endl;
	for( TrackCandidateCollection::const_iterator cand = hltL3TkTracksCandOI->begin(); cand != hltL3TkTracksCandOI->end(); ++cand) {
	  auto const & candSS = cand->trajectoryStateOnDet();
	  DetId id(candSS.detId());
	  const GeomDet * g = theService->trackingGeometry()->idToDet( id );
	  const Surface * surface=&g->surface();
	  TrajectoryStateOnSurface candTSOS = trajectoryStateTransform::transientState(candSS, surface,  &*(theService)->magneticField());
	  
	  /// TSOS errors: 
	  AlgebraicSymMatrix66 errors = candTSOS.cartesianError().matrix();
	  double partialPterror = errors(3,3)*pow(candTSOS.globalMomentum().x(),2) + errors(4,4)*pow(candTSOS.globalMomentum().y(),2);
	  double numberOfHits = cand->recHits().second-cand->recHits().first;
	  unsigned int stopReason = static_cast<unsigned int>(cand->stopReason());
/*B	  if (debuglevel_ >1) {
	    cout << " TSOS of hltL3TkTracksCandOI     = \n"
		 <<"x: "<<candTSOS.globalPosition()<< " --> " << candTSOS.localPosition() << "\n" 
		 <<"p: "<<candTSOS.globalMomentum()<< "\n"
		 <<"pt: " << candTSOS.globalMomentum().perp() << " +/- " << sqrt(partialPterror)/candTSOS.globalMomentum().perp()  << "\n"
		 <<"eta: " << candTSOS.globalMomentum().eta() << " +/- " << sqrt(candTSOS.curvilinearError().matrix()(1,1))*abs(sin(candTSOS.globalMomentum().theta()))  << "\n"
		 <<"phi: " << candTSOS.globalMomentum().phi() << " +/- " << sqrt(candTSOS.curvilinearError().matrix()(2,2)) << "\n"
		 << id.subdetId() << " " << id.rawId() << "\n"
	         <<"nhits: " << numberOfHits << "; stop reason: " << stopReason << endl;
	    
	  }
B*/

	}
      }
      edm::Handle<reco::TrackCollection> hltL3TkTracksOINoHP;
      iEvent.getByToken(theTracksOINoHPToken_, hltL3TkTracksOINoHP);
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
    //if (!isHLT)
     // nlayers3D += t.hitPattern().numberOfValidStripLayersWithMonoAndStereo();
    //else {
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
    //}
   
        cout<<"num3D Layers: "<<nlayers3D<<endl; 


          cout<<"Trk Layers without measurement:  "<<t->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS)<<endl;
          cout<<"dxy:
