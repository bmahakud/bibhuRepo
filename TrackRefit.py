import FWCore.ParameterSet.Config as cms

process = cms.Process("Refitting")

### Standard Configurations
#process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# choose!
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data_GRun', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')
#auto:run2_mc_FULL



### Track refitter specific stuff
import RecoTracker.TrackProducer.TrackRefitter_cfi
import CommonTools.RecoAlgos.recoTrackRefSelector_cfi
process.mytkselector = CommonTools.RecoAlgos.recoTrackRefSelector_cfi.recoTrackRefSelector.clone()
process.mytkselector.quality = ['highPurity']
process.mytkselector.min3DLayer = 2
process.mytkselector.ptMin = 0.5
process.mytkselector.tip = 1.0
process.vertexConstraint = RecoTracker.TrackProducer.VertexConstraintProducer_cfi.vertexConstraint.clone()

process.myRefittedTracks = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone()
process.myRefittedTracks.src= 'mytkselector'
process.myRefittedTracks.NavigationSchool = ''
process.myRefittedTracks.Fitter = 'FlexibleKFFittingSmoother'
#process.myRefittedTracks.constraint = 'vertex'
#process.myRefittedTracks.srcConstr = 'vertexConstraint'


### and an analyzer
process.trajCout = cms.EDAnalyzer('TrajectoryAnalyzer',
   trajectoryInput=cms.InputTag('myRefittedTracks')
)



process.vertexConstraint = cms.EDProducer(
    "VertexConstraintProducer",
    srcTrk = cms.InputTag("AlignmentTrackSelector"),
    srcVtx = cms.InputTag("offlinePrimaryVertices")
)


process.source = cms.Source ("PoolSource",
                             fileNames=cms.untracked.vstring('file:step3_RAW2DIGI_L1Reco_RECO_.root',
                              #fileNames=cms.untracked.vstring('/store/mc/RunIISummer16DR80Premix/BsToMuMuPhi_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/004FAC41-C9CF-E611-9812-0242AC130004.root',
                            ),
                             skipEvents=cms.untracked.uint32(0)
                             )


process.TRACKS = cms.OutputModule("PoolOutputModule",
                                 outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                                                        'keep recoTracks_*_*_*',
                                                                        'keep recoTrackExtras_*_*_*'),

                                fileName = cms.untracked.string('refittingMuon.root')
                               )







process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.Path = cms.Path(process.mytkselector+process.myRefittedTracks+process.trajCout)

#process.Path = cms.Path(process.mytkselector+process.vertexConstraint*process.myRefittedTracks+process.trajCout)


process.outpath = cms.EndPath(process.TRACKS)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
