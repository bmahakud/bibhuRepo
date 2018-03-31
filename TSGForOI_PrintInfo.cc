/**
  \class    TSGForOI
  \brief    Create L3MuonTrajectorySeeds from L2 Muons updated at vertex in an outside in manner
  \author   Benjamin Radburn-Smith, Santiago Folgueras
 */

#include "RecoMuon/TrackerSeedGenerator/plugins/TSGForOI.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include <memory>

using namespace edm;
using namespace std;

TSGForOI::TSGForOI(const edm::ParameterSet & iConfig) :
  src_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
  numOfMaxSeedsParam_(iConfig.getParameter<uint32_t>("maxSeeds")),
  numOfLayersToTry_(iConfig.getParameter<int32_t>("layersToTry")),
  numOfHitsToTry_(iConfig.getParameter<int32_t>("hitsToTry")),
  fixedErrorRescalingForHits_(iConfig.getParameter<double>("fixedErrorRescaleFactorForHits")),
  fixedErrorRescalingForHitless_(iConfig.getParameter<double>("fixedErrorRescaleFactorForHitless")),
  adjustErrorsDynamicallyForHits_(iConfig.getParameter<bool>("adjustErrorsDynamicallyForHits")),
  adjustErrorsDynamicallyForHitless_(iConfig.getParameter<bool>("adjustErrorsDynamicallyForHitless")),
  estimatorName_(iConfig.getParameter<std::string>("estimator")),
  minEtaForTEC_(iConfig.getParameter<double>("minEtaForTEC")),
  maxEtaForTOB_(iConfig.getParameter<double>("maxEtaForTOB")),
  useHitLessSeeds_(iConfig.getParameter<bool>("UseHitLessSeeds")),
  useStereoLayersInTEC_(iConfig.getParameter<bool>("UseStereoLayersInTEC")),
  updator_(new KFUpdator()),
  measurementTrackerTag_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("MeasurementTrackerEvent"))),
  pT1_(iConfig.getParameter<double>("pT1")),
  pT2_(iConfig.getParameter<double>("pT2")),
  pT3_(iConfig.getParameter<double>("pT3")),
  eta1_(iConfig.getParameter<double>("eta1")),
  eta2_(iConfig.getParameter<double>("eta2")),
  SF1_(iConfig.getParameter<double>("SF1")),
  SF2_(iConfig.getParameter<double>("SF2")),
  SF3_(iConfig.getParameter<double>("SF3")),
  SF4_(iConfig.getParameter<double>("SF4")),
  SF5_(iConfig.getParameter<double>("SF5")),
  tsosDiff_(iConfig.getParameter<double>("tsosDiff")),
  propagatorName_(iConfig.getParameter<std::string>("propagatorName")), 
  theCategory(string("Muon|RecoMuon|TSGForOI"))
{
  produces<std::vector<TrajectorySeed> >();
}


TSGForOI::~TSGForOI(){
}


void TSGForOI::produce(edm::StreamID sid, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  /// Init variables
  using namespace std;
  unsigned int numOfMaxSeeds = numOfMaxSeedsParam_;
  unsigned int numSeedsMade=0;
  bool analysedL2 = false;
  unsigned int layerCount = 0;




  /// Surface used to make a TSOS at the PCA to the beamline
  Plane::PlanePointer dummyPlane = Plane::build(Plane::PositionType(), Plane::RotationType());

  /// Read ESHandles
  edm::Handle<MeasurementTrackerEvent>          measurementTrackerH;
  edm::ESHandle<Chi2MeasurementEstimatorBase>   estimatorH;
  edm::ESHandle<MagneticField>                  magfieldH;
  edm::ESHandle<Propagator>                     propagatorAlongH;
  edm::ESHandle<Propagator>                     propagatorOppositeH;
  edm::ESHandle<TrackerGeometry>                tmpTkGeometryH;  
  edm::ESHandle<GlobalTrackingGeometry>         geometryH;

  iSetup.get<IdealMagneticFieldRecord>().get(magfieldH);
  iSetup.get<TrackingComponentsRecord>().get(propagatorName_, propagatorOppositeH);
  iSetup.get<TrackingComponentsRecord>().get(propagatorName_, propagatorAlongH);
  iSetup.get<GlobalTrackingGeometryRecord>().get(geometryH);
  iSetup.get<TrackerDigiGeometryRecord>().get(tmpTkGeometryH);
  iSetup.get<TrackingComponentsRecord>().get(estimatorName_,estimatorH);
  iEvent.getByToken(measurementTrackerTag_, measurementTrackerH);


  //edm::ESHandle<TrackerGeometry> theG;
   // iSetup.get<TrackerDigiGeometryRecord>().get(theG);



  /// Read L2 track collection
  edm::Handle<reco::TrackCollection> l2TrackCol;
  iEvent.getByToken(src_, l2TrackCol);

  //	The product:
  std::unique_ptr<std::vector<TrajectorySeed> > result(new std::vector<TrajectorySeed>());

  //	Get vector of Detector layers once:
  std::vector<BarrelDetLayer const*> const& tob = measurementTrackerH->geometricSearchTracker()->tobLayers();
  std::vector<ForwardDetLayer const*> const& tecPositive = tmpTkGeometryH->isThere(GeomDetEnumerators::P2OTEC) ? 
                                                                measurementTrackerH->geometricSearchTracker()->posTidLayers() : 
                                                                measurementTrackerH->geometricSearchTracker()->posTecLayers(); 
  std::vector<ForwardDetLayer const*> const& tecNegative = tmpTkGeometryH->isThere(GeomDetEnumerators::P2OTEC) ? 
                                                                measurementTrackerH->geometricSearchTracker()->negTidLayers() : 
                                                                measurementTrackerH->geometricSearchTracker()->negTecLayers();
  edm::ESHandle<TrackerTopology> tTopo_handle;
  iSetup.get<TrackerTopologyRcd>().get(tTopo_handle);
  const TrackerTopology* tTopo = tTopo_handle.product();

  //	Get the suitable propagators: 
  std::unique_ptr<Propagator> propagatorAlong = SetPropagationDirection(*propagatorAlongH,alongMomentum);
  std::unique_ptr<Propagator> propagatorOpposite = SetPropagationDirection(*propagatorOppositeH,oppositeToMomentum);

  edm::ESHandle<Propagator> SmartOpposite;
  edm::ESHandle<Propagator> SHPOpposite;
  iSetup.get<TrackingComponentsRecord>().get("hltESPSmartPropagatorAnyOpposite", SmartOpposite);
  iSetup.get<TrackingComponentsRecord>().get("hltESPSteppingHelixPropagatorOpposite", SHPOpposite);

  int il2=0;
  int totalSeeds=0;

  //	Loop over the L2's and make seeds for all of them:
  LogTrace(theCategory) << "TSGForOI::produce: Number of L2's: " << l2TrackCol->size();
  for (unsigned int l2TrackColIndex(0);l2TrackColIndex!=l2TrackCol->size();++l2TrackColIndex){
     if(tsosDiff_==0.2){//

     il2++;

   unsigned int numHitlessSeedsMade=0;
   unsigned int numHitbasedSeedsMade=0;



    const reco::TrackRef l2(l2TrackCol, l2TrackColIndex);
    std::unique_ptr<std::vector<TrajectorySeed> > out(new std::vector<TrajectorySeed>());
    LogTrace("TSGForOI") << "TSGForOI::produce: L2 muon pT, eta, phi --> " << l2->pt() << " , " << l2->eta() << " , " << l2->phi() << endl;
   
    std::cout<<"L2 Muon pt,eta, phi:   "<<l2->pt() << " , " << l2->eta() << " , " << l2->phi() << endl;
    std::cout<<"L2 Muon px,py, pz:   "<<l2->px() << " , " << l2->py() << " , " << l2->pz() << endl;
      


 
    FreeTrajectoryState fts = trajectoryStateTransform::initialFreeState(*l2, magfieldH.product());


    dummyPlane->move(fts.position() - dummyPlane->position());
    TrajectoryStateOnSurface tsosAtIP = TrajectoryStateOnSurface(fts, *dummyPlane);
    LogTrace("TSGForOI") << "TSGForOI::produce: Created TSOSatIP: " << tsosAtIP << std::endl;
    
    cout<<"TSOS IP:position-> "<<tsosAtIP.globalPosition()<<endl;
    cout<<"TSOS IP:momentum-> "<<tsosAtIP.globalMomentum()<<endl;


    // get the TSOS on the innermost layer of the L2. 
    TrajectoryStateOnSurface tsosAtMuonSystem = trajectoryStateTransform::innerStateOnSurface(*l2, *geometryH, magfieldH.product());
    LogTrace("TSGForOI") << "TSGForOI::produce: Created TSOSatMuonSystem: " << tsosAtMuonSystem <<endl;
    
    LogTrace("TSGForOI") << "TSGForOI::produce: Check the error of the L2 parameter and use hit seeds if big errors" << endl;


    cout<<"TSOS mu Station:position-> "<<tsosAtMuonSystem.globalPosition()<<endl;
    cout<<"TSOS mu station :momentum-> "<<tsosAtMuonSystem.globalMomentum()<<endl;





    StateOnTrackerBound fromInside(propagatorAlong.get());
    TrajectoryStateOnSurface outerTkStateInside = fromInside(fts);
      
    StateOnTrackerBound fromOutside(&*SHPOpposite);//implemented by me

   // StateOnTrackerBound fromOutside(&*SmartOpposite);
    TrajectoryStateOnSurface outerTkStateOutside = fromOutside(tsosAtMuonSystem);


    //cout<<"TSOS at mu station after tracker state found: position "<<tsosAtMuonSystem.globalPosition()<<endl;
    //cout<<"TSOS at mu station after tracker state found: momentum "<<tsosAtMuonSystem.globalMomentum()<<endl;


      // for now only checking if the two positions (using updated and not-updated) agree withing certain extent, 
      // will probably have to design something fancier for the future. 
    auto dist=0.0;
    bool useBoth = false;
    if (outerTkStateInside.isValid() && outerTkStateOutside.isValid()){
	dist = match_Chi2(outerTkStateInside,outerTkStateOutside);
    }
    if (dist>tsosDiff_){
	useBoth = true; 
    }

    numSeedsMade=0;
    analysedL2 = false;

                                
                                unsigned int NumCompatibleHitsTECtsIP=0;
                                unsigned int NumCompatibleHitsTEBtsIP=0;
                                unsigned int NumCompatibleHitsTECtsMuS=0;
                                unsigned int NumCompatibleHitsTEBtsMuS=0;









    // if both TSOSs agree, use only the one at vertex, as it uses more information. If they do not agree, search for seed based on both
     int numTOBLayers=0;
    //		BARREL
    if (std::abs(l2->eta()) < maxEtaForTOB_) {
      layerCount = 0;
      for (auto it=tob.rbegin(); it!=tob.rend(); ++it) {	//This goes from outermost to innermost layer
       numTOBLayers++;
	LogTrace("TSGForOI") << "TSGForOI::produce: looping in TOB layer " << layerCount << endl; 
        int numseed=out->size();
      //  cout<<" TOB Layer number:  "<<numTOBLayers<<endl;
	findSeedsOnLayer(iSetup,"TEBip",NumCompatibleHitsTECtsIP,NumCompatibleHitsTEBtsIP,NumCompatibleHitsTECtsMuS,NumCompatibleHitsTEBtsMuS,numHitlessSeedsMade,numHitbasedSeedsMade,tTopo, **it, tsosAtIP,  *(propagatorAlong.get()), *(propagatorOpposite.get()), l2, 
			 estimatorH, measurementTrackerH, numSeedsMade, numOfMaxSeeds, layerCount, analysedL2, out);

      //cout<<"       num seed madeTOB(usingIP): "<<out->size()-numseed<<endl;
      }

      if (useBoth){

    //  layerCount = 0;

	      numSeedsMade=0;
	      for (auto it=tob.rbegin(); it!=tob.rend(); ++it) {	//This goes from outermost to innermost layer
			LogTrace("TSGForOI") << "TSGForOI::produce: looping in TOB layer " << layerCount << endl; 
      int numseed=out->size();
			findSeedsOnLayer(iSetup,"TEBmus",NumCompatibleHitsTECtsIP,NumCompatibleHitsTEBtsIP,NumCompatibleHitsTECtsMuS,NumCompatibleHitsTEBtsMuS,numHitlessSeedsMade,numHitbasedSeedsMade,tTopo, **it, outerTkStateOutside,  *(propagatorOpposite.get()), *(propagatorOpposite.get()), l2, 
			 estimatorH, measurementTrackerH, numSeedsMade, numOfMaxSeeds, layerCount, analysedL2, out);
      //cout<<"       num seed madeTOB(using mu station): "<<out->size()-numseed<<endl;

      		}
	}
    }
    //		Reset Number of seeds if in overlap region:

    if (std::abs(l2->eta())>minEtaForTEC_ && std::abs(l2->eta())<maxEtaForTOB_){
      numSeedsMade=0;
    }


    int nTECLayers=0;
    //		ENDCAP+
    if (l2->eta() > minEtaForTEC_) {
      layerCount = 0;
      for (auto it=tecPositive.rbegin(); it!=tecPositive.rend(); ++it) {
	LogTrace("TSGForOI") << "TSGForOI::produce: looping in TEC+ layer " << layerCount << endl; 
        nTECLayers++;
        int numseed=out->size();
	findSeedsOnLayer(iSetup,"TECip",NumCompatibleHitsTECtsIP,NumCompatibleHitsTEBtsIP,NumCompatibleHitsTECtsMuS,NumCompatibleHitsTEBtsMuS,numHitlessSeedsMade,numHitbasedSeedsMade,tTopo, **it, tsosAtIP, *(propagatorAlong.get()), *(propagatorOpposite.get()), l2, 
			 estimatorH, measurementTrackerH, numSeedsMade, numOfMaxSeeds, layerCount, analysedL2, out);

      //cout<<"       num seed madeTEC+(using IP): "<<out->size()-numseed<<endl;

	}
   //cout<<"nTEC layers: "<<nTECLayers<<endl;

	if (useBoth){

    //layerCount = 0;
		numSeedsMade=0;
      		for (auto it=tecPositive.rbegin(); it!=tecPositive.rend(); ++it) {
			LogTrace("TSGForOI") << "TSGForOI::produce: looping in TEC+ layer " << layerCount << endl; 

                  int numseed=out->size();
			findSeedsOnLayer(iSetup,"TECmus",NumCompatibleHitsTECtsIP,NumCompatibleHitsTEBtsIP,NumCompatibleHitsTECtsMuS,NumCompatibleHitsTEBtsMuS,numHitlessSeedsMade,numHitbasedSeedsMade,tTopo, **it, outerTkStateOutside,  *(propagatorOpposite.get()), *(propagatorOpposite.get()), l2, 
			 estimatorH, measurementTrackerH, numSeedsMade, numOfMaxSeeds, layerCount, analysedL2, out);
           //cout<<"       num seed madeTEC+(using mu station): "<<out->size()-numseed<<endl;

		}
     	}
    }

    //		ENDCAP-
    if (l2->eta() < -minEtaForTEC_) {
      layerCount = 0;
      for (auto it=tecNegative.rbegin(); it!=tecNegative.rend(); ++it) {
	LogTrace("TSGForOI") << "TSGForOI::produce: looping in TEC- layer " << layerCount << endl; 
        int numseed=out->size();
	findSeedsOnLayer(iSetup,"TECip",NumCompatibleHitsTECtsIP,NumCompatibleHitsTEBtsIP,NumCompatibleHitsTECtsMuS,NumCompatibleHitsTEBtsMuS,numHitlessSeedsMade,numHitbasedSeedsMade,tTopo, **it, tsosAtIP,  *(propagatorAlong.get()), *(propagatorOpposite.get()), l2, 
			 estimatorH, measurementTrackerH, numSeedsMade, numOfMaxSeeds, layerCount, analysedL2, out);

       //cout<<"       num seed madeTEC-(using IP): "<<out->size()-numseed<<endl;
	}
	if (useBoth){

   //layerCount = 0;
		numSeedsMade=0;
	      	for (auto it=tecNegative.rbegin(); it!=tecNegative.rend(); ++it) {
			LogTrace("TSGForOI") << "TSGForOI::produce: looping in TEC- layer " << layerCount << endl; 

             int numseed=out->size();
			findSeedsOnLayer(iSetup,"TECmus",NumCompatibleHitsTECtsIP,NumCompatibleHitsTEBtsIP,NumCompatibleHitsTECtsMuS,NumCompatibleHitsTEBtsMuS,numHitlessSeedsMade,numHitbasedSeedsMade,tTopo, **it, outerTkStateOutside,  *(propagatorOpposite.get()), *(propagatorOpposite.get()), l2, estimatorH, measurementTrackerH, numSeedsMade, numOfMaxSeeds, layerCount, analysedL2, out);

            //cout<<"       num seed madeTEC-(using Mu station): "<<out->size()-numseed<<endl;
		}
      	}
    }

    

    for (std::vector<TrajectorySeed>::iterator it=out->begin(); it!=out->end(); ++it){
      result->push_back(*it);
    }


              cout<<"	Number of Compatible hits TEC(TSOSIP):   "<< NumCompatibleHitsTECtsIP<<endl;
              cout<<"	Number of Compatible hits TEC (TSOSmuS): "<<NumCompatibleHitsTECtsMuS<<endl;
              cout<<"	Number of Compatible hits TEB (TSOSIP): "<<NumCompatibleHitsTEBtsIP<<endl;
              cout<<"	Number of Compatible hits TEB (TSOSMuS): "<<NumCompatibleHitsTEBtsMuS<<endl;




  
  cout<<"      Number of hitless seeds: "<<numHitlessSeedsMade<<endl;
  cout<<"      Number of hitbased seeds: "<<numHitbasedSeedsMade<<endl;

   cout<<"total num of seeds made:  "<<numHitlessSeedsMade+numHitbasedSeedsMade<<endl;  







   }//ending dist==0.2
  } //L2Collection

  cout<<"total number of seeds in this event: "<<result->size()<<endl;
  edm::LogInfo(theCategory) << "TSGForOI::produce: number of seeds made: " << result->size();





  iEvent.put(std::move(result));
}

void TSGForOI::findSeedsOnLayer(const edm::EventSetup& iSetup,
                                std::string SubDet,
                                unsigned int& NumCompatibleHitsTECtsIP,
                                unsigned int& NumCompatibleHitsTEBtsIP,
                                unsigned int& NumCompatibleHitsTECtsMuS,
                                unsigned int& NumCompatibleHitsTEBtsMuS, 
                                unsigned int& numHitlessSeedsMade,
                                unsigned int& numHitbasedSeedsMade,
				const TrackerTopology* tTopo,
				const GeometricSearchDet &layer,
				const TrajectoryStateOnSurface &tsosAtIP,
				const Propagator& propagatorAlong,
				const Propagator& propagatorOpposite,
				const reco::TrackRef l2,
				edm::ESHandle<Chi2MeasurementEstimatorBase>& estimatorH,
				edm::Handle<MeasurementTrackerEvent>& measurementTrackerH,
				unsigned int& numSeedsMade,
				unsigned int& numOfMaxSeeds,
				unsigned int& layerCount,
				bool& analysedL2,
				std::unique_ptr<std::vector<TrajectorySeed> >& out)  const{
  
  if (numSeedsMade>numOfMaxSeeds) return;
  LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: numSeedsMade = " << numSeedsMade << " , layerCount = " <<  layerCount << endl;
  
  double errorSFHits=1.0;
  double errorSFHitless=1.0;
  if (!adjustErrorsDynamicallyForHits_) errorSFHits = fixedErrorRescalingForHits_;
  else                                  errorSFHits = calculateSFFromL2(l2);
  if (!adjustErrorsDynamicallyForHitless_) errorSFHitless = fixedErrorRescalingForHitless_;
  // else if(!analysedL2 && adjustErrorsDynamicallyForHitless_) {
   //       errorSFHitless=calculateSFFromL2(l2);
   //       analysedL2=true;
  // }



  // Hitless:  TO Be discarded from here at some point. 
  if (useHitLessSeeds_) {
    LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: Start hitless" << endl;
    std::vector< GeometricSearchDet::DetWithState > dets;
    //TrajectoryStateOnSurface onLayer(tsosAtIP);
    //onLayer.rescaleError(errorSFHitless);

    layer.compatibleDetsV(tsosAtIP, propagatorAlong, *estimatorH, dets);
    if (!dets.empty()) { 
      //cout<<"     comaptible det found for hitless: size"<<dets.size()<<endl; 
      auto const& detOnLayer = dets.front().first;
      auto const& tsosOnLayer = dets.front().second;
      LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: tsosOnLayer " << tsosOnLayer << endl;
      if (!tsosOnLayer.isValid()){
	edm::LogInfo(theCategory) << "ERROR!: Hitless TSOS is not valid!";
      }
      else{
	// calculate SF from L2 (only once -- if needed)
	if (!analysedL2 && adjustErrorsDynamicallyForHitless_) {
	  errorSFHitless=calculateSFFromL2(l2);
	  analysedL2=true;
	}
	
        dets.front().second.rescaleError(errorSFHitless);
        //tsosOnLayer.rescaleError(errorSFHitless);

        //cout<<"      hitlessseed tsos position: "<<tsosOnLayer.globalPosition()<<endl;
        //cout<<"      hitlessseed tsos momentum:  "<<tsosOnLayer.globalMomentum()<<endl;
	PTrajectoryStateOnDet const& ptsod = trajectoryStateTransform::persistentState(tsosOnLayer,detOnLayer->geographicalId().rawId());
	TrajectorySeed::recHitContainer rHC;
	out->push_back(TrajectorySeed(ptsod,rHC,oppositeToMomentum));
        numHitlessSeedsMade++;
	LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: TSOD (Hitless) done " << endl;
        numSeedsMade++;
      }
    }
  }




  // Hits:
  if (layerCount>numOfLayersToTry_) {
  //cout<<"     not going to find hitbased seed in this layer"<<endl;
   return;}
  LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: Start Hits" <<endl;  
  if (makeSeedsFromHits(iSetup,SubDet,NumCompatibleHitsTECtsIP,NumCompatibleHitsTEBtsIP,NumCompatibleHitsTECtsMuS,NumCompatibleHitsTEBtsMuS,numHitbasedSeedsMade,tTopo, layer, tsosAtIP, *out, propagatorAlong, *measurementTrackerH, estimatorH, numSeedsMade, errorSFHits, l2->eta()))  ++layerCount; 
}

double TSGForOI::calculateSFFromL2(const reco::TrackRef track) const{

  double theSF=1.0;
  //	L2 direction vs pT blowup - as was previously done:
  //	Split into 4 pT ranges: <pT1_, pT1_<pT2_, pT2_<pT3_, <pT4_: 13,30,70
  //	Split into 2 eta ranges for the middle two pT ranges: 1.0,1.4
  double abseta = std::abs(track->eta());
  if (track->pt()<=pT1_) theSF=SF1_;
  if (track->pt()>pT1_ && track->pt()<=pT2_){
    if (abseta<=eta1_) theSF=SF3_;
    if (abseta>eta1_ && abseta<=eta2_) theSF=SF2_;
    if (abseta>eta2_) theSF=SF3_;
  }
  if (track->pt()>pT2_ && track->pt()<=pT3_){
    if (abseta<=eta1_) theSF=SF5_;
    if (abseta>eta1_ && abseta<=eta2_) theSF=SF4_;
    if (abseta>eta2_) theSF=SF5_;
  }
  if (track->pt()>pT3_) theSF=SF5_;

  LogTrace(theCategory) << "TSGForOI::calculateSFFromL2: SF has been calculated as: " << theSF;
  return theSF;
}


int TSGForOI::makeSeedsFromHits(const edm::EventSetup& iSetup,
                                std::string SubDet,
                                unsigned int& NumCompatibleHitsTECtsIP,
                                unsigned int& NumCompatibleHitsTEBtsIP,
                                unsigned int& NumCompatibleHitsTECtsMuS,
                                unsigned int& NumCompatibleHitsTEBtsMuS,
                                unsigned int& numHitbasedSeedsMade,
				const TrackerTopology* tTopo,
				const GeometricSearchDet &layer,
				const TrajectoryStateOnSurface &tsosAtIP,
				std::vector<TrajectorySeed> &out,
				const Propagator& propagatorAlong,
				const MeasurementTrackerEvent &measurementTracker,
				edm::ESHandle<Chi2MeasurementEstimatorBase>& estimatorH,
				unsigned int& numSeedsMade,
				const double errorSF,
				const double l2Eta)  const{

  //		Error Rescaling:
 cout<<"starting findhitbased seed ............................................................................with "<<SubDet<<endl;
 
 
  TrajectoryStateOnSurface onLayer(tsosAtIP);
  onLayer.rescaleError(errorSF);    

  std::vector< GeometricSearchDet::DetWithState > dets;
  layer.compatibleDetsV(onLayer, propagatorAlong, *estimatorH, dets);
  //if(dets.size() !=0){ 
  // cout<<"   compatible det for hitbased : size= "<<dets.size()<<endl;
  // } 

    edm::ESHandle<TrackerGeometry> theG;
    iSetup.get<TrackerDigiGeometryRecord>().get(theG);
 

  //	Find Measurements on each DetWithState:
  LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: find measurements on each detWithState  " << dets.size() << endl;

  int countm=0;
  std::vector<TrajectoryMeasurement> meas;

  cout<<"detsize   on layer : "<<dets.size()<<endl;
  for (std::vector<GeometricSearchDet::DetWithState>::iterator it=dets.begin(); it!=dets.end(); ++it) {
    MeasurementDetWithData det = measurementTracker.idToDet(it->first->geographicalId());
    if (det.isNull()) {
      continue;
    }
    if (!it->second.isValid()) continue;	//Skip if TSOS is not valid


     cout<<"Raw ID : "<<it->first->geographicalId().rawId()<<endl;

     if(it->first->geographicalId().det()==1 && it->first->geographicalId().subdetId()==5 ){
     TOBDetId tobdetID(it->first->geographicalId().rawId());
     cout<<"LayerNumber: "<<tobdetID.layerNumber()<<endl; 

      
    
     //cout <<"\t\t\tRecHit in GP "<<theG->idToDet((*hitIt)->geographicalId())->surface().toGlobal((*hitIt)->localPosition()) <<endl; 


     }


     cout<<"TSOS global position:-----------------------    "<<it->second.globalPosition()<<endl;
    //cout<<"TSOS at IP:--------------------------------- "<<onLayer.globalPosition()<<endl;

     std::vector < TrajectoryMeasurement > mymeas = det.fastMeasurements(it->second, onLayer, propagatorAlong, *estimatorH);	//Second TSOS is not used
 
    //cout<<"size my meas:  "<<mymeas.size()<<endl;
   for (std::vector<TrajectoryMeasurement>::const_iterator it2 = mymeas.begin(), ed2 = mymeas.end(); it2 != ed2; ++it2) {
     if (it2->recHit()->isValid()){ meas.push_back(*it2);	//Only save those which are valid

      TrajectoryMeasurement::ConstRecHitPointer recP=it2->recHit();

      cout <<"\t\t\tRecHit globalposition: "<<theG->idToDet(recP->geographicalId())->surface().toGlobal(recP->localPosition()) <<endl;

       cout<<"local position: "<<recP->localPosition()<<endl;
      //const TrackingRecHit* rechit=
      cout<<"          Estimate: "<<it2->estimate()<<endl;

   //  if(SubDet=="TECip")NumCompatibleHitsTECtsIP++;
   //  if(SubDet=="TECmus")NumCompatibleHitsTECtsMuS++;
   //  if(SubDet=="TEBip")NumCompatibleHitsTEBtsIP++;
   //  if(SubDet=="TEBmus")NumCompatibleHitsTEBtsMuS++;
   countm++;
    }
    }
  }
  


  //	Update TSOS using TMs after sorting, then create Trajectory Seed and put into vector:
  LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: Update TSOS using TMs after sorting, then create Trajectory Seed, number of TM = " << meas.size() << endl;
  unsigned int found = 0;
  std::sort(meas.begin(), meas.end(), TrajMeasLessEstim());
  
  //cout<<"measue size:   "<<meas.size()<<endl;

  int countlater=0;
  for (std::vector<TrajectoryMeasurement>::const_iterator it=meas.begin(); it!=meas.end(); ++it) {
  
     countlater++;
     if(SubDet=="TECip")NumCompatibleHitsTECtsIP++;
     if(SubDet=="TECmus")NumCompatibleHitsTECtsMuS++;
     if(SubDet=="TEBip")NumCompatibleHitsTEBtsIP++;
     if(SubDet=="TEBmus")NumCompatibleHitsTEBtsMuS++;

  //  cout<<"TSOS global position before update:   "<<it->forwardPredictedState().globalPosition()<<endl;
    

    TrajectoryStateOnSurface updatedTSOS = updator_->update(it->forwardPredictedState(), *it->recHit());
    LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: TSOS for TM " << found << endl;


   // cout<<"Updated TSOS from forward predited state:    "<<updatedTSOS.globalPosition()<<endl;

    // if(SubDet=="TECip")NumCompatibleHitsTECtsIP++;
    // if(SubDet=="TECmus")NumCompatibleHitsTECtsMuS++;
    // if(SubDet=="TEBip")NumCompatibleHitsTEBtsIP++;
    // if(SubDet=="TEBmus")NumCompatibleHitsTEBtsMuS++;



    if (not updatedTSOS.isValid()) {
     cout<<"TSOS not valid: "<<endl; 
     continue;

      }

    // CHECK if is StereoLayer: 
    if (useStereoLayersInTEC_ && (fabs(l2Eta) > 0.8 && fabs(l2Eta) < 1.6)) { 
     DetId detid = ((*it).recHit()->hit())->geographicalId();
     if (detid.subdetId() == StripSubdetector::TEC) {
    		if (!tTopo->tecIsStereo(detid.rawId())) break;  // try another layer
      }
   }
    
    edm::OwnVector<TrackingRecHit> seedHits;
    seedHits.push_back(*it->recHit()->hit());
     //if(SubDet=="TECip")NumCompatibleHitsTECtsIP++;
     //if(SubDet=="TECmus")NumCompatibleHitsTECtsMuS++;
     //if(SubDet=="TEBip")NumCompatibleHitsTEBtsIP++;
     //if(SubDet=="TEBmus")NumCompatibleHitsTEBtsMuS++;




   //cout<<"   hitbased Seed Updated TSOS position: "<<updatedTSOS.globalPosition()<<endl;
   //cout<<"   hitbased seed Updated TSOS momentum:  "<<updatedTSOS.globalMomentum()<<endl;

    PTrajectoryStateOnDet const& pstate = trajectoryStateTransform::persistentState(updatedTSOS, it->recHit()->geographicalId().rawId());
    TrajectorySeed seed(pstate, std::move(seedHits), oppositeToMomentum);
    LogTrace("TSGForOI") << "TSGForOI::findSeedsOnLayer: number of seedHits: " << seedHits.size() << endl;
    out.push_back(seed);
    numHitbasedSeedsMade++;
    found++;
    numSeedsMade++;
    if (found == numOfHitsToTry_) break;
  }

  //cout<<"countlater: "<<countlater<<endl;

  return found;
}
//
//// calculate Chi^2 of two trajectory states
////

double TSGForOI::match_Chi2(const TrajectoryStateOnSurface& tsos1,
                  const TrajectoryStateOnSurface& tsos2) const{

  if ( !tsos1.isValid() || !tsos2.isValid() ) return -1.;

  AlgebraicVector5 v(tsos1.localParameters().vector() - tsos2.localParameters().vector());
  AlgebraicSymMatrix55 m(tsos1.localError().matrix() + tsos2.localError().matrix());

  bool ierr = !m.Invert();

  if ( ierr ) {
    edm::LogInfo("TSGForOI") << "Error inverting covariance matrix";
    return -1;
  }

  double est = ROOT::Math::Similarity(v,m);

  return est;

}


void TSGForOI::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src",edm::InputTag("hltL2Muons","UpdatedAtVtx"));
  desc.add<int>("layersToTry",1);
  desc.add<double>("fixedErrorRescaleFactorForHitless",2.0);
  desc.add<int>("hitsToTry",1);
  desc.add<bool>("adjustErrorsDynamicallyForHits",false);
  desc.add<bool>("adjustErrorsDynamicallyForHitless",false);
  desc.add<edm::InputTag>("MeasurementTrackerEvent",edm::InputTag("hltSiStripClusters"));
  desc.add<bool>("UseHitLessSeeds",true);
  desc.add<bool>("UseStereoLayersInTEC",false);
  desc.add<std::string>("estimator","hltESPChi2MeasurementEstimator100");
  desc.add<double>("maxEtaForTOB",1.2);
  desc.add<double>("minEtaForTEC",0.8);
  desc.addUntracked<bool>("debug",true);
  desc.add<double>("fixedErrorRescaleFactorForHits",2.0);
  desc.add<unsigned int>("maxSeeds",1);
  desc.add<double>("pT1",13.0);
  desc.add<double>("pT2",30.0);
  desc.add<double>("pT3",70.0);
  desc.add<double>("eta1",1.0);
  desc.add<double>("eta2",1.4);
  desc.add<double>("SF1",3.0);
  desc.add<double>("SF2",4.0);
  desc.add<double>("SF3",5.0);
  desc.add<double>("SF4",7.0);
  desc.add<double>("SF5",10.0);
  desc.add<double>("tsosDiff",0.03);
  desc.add<std::string>("propagatorName","PropagatorWithMaterial");
  descriptions.add("TSGForOI",desc);
}

DEFINE_FWK_MODULE(TSGForOI);

