# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --filein /store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+120000+E5F6DFC8-65CF-2A41-B0BE-E82E041CB012.root --fileout file:output.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v20 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename RunIIAutumn18NanoAOD.py --no_exec -n 1000 --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeBxToMuMu --customise=Bmm5/NanoAOD/nano_cff.nanoAOD_customizeV0ForMuonFake --customise_commands=process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))

#import sys
#num = sys.argv[2]


import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NANO',eras.Run2_2018,eras.run2_nanoAOD_102Xv1)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

#InFile = "file:/eos/user/b/bimahaku/BstoMuMu_RECO_miniAOD_Files/BstoMuMuFile-RunIIAutumn18MiniAOD_10_2_18_3.root"#%(str(num))
InFile ="file:/afs/cern.ch/user/b/bimahaku/public/BeamSpotRefit/Test2/New2/CMSSW_10_2_18/src/PrimarySecondary_RecoAndMiniAOD/BPH-RunIIAutumn18DRPremix-00048_MiniAODSIM.root"
SecFile = 'file:/afs/cern.ch/user/b/bimahaku/public/BeamSpotRefit/Test2/New2/CMSSW_10_2_18/src/PrimarySecondary_RecoAndMiniAOD/BPH-RunIIAutumn18DRPremix-00048_RECO.root'
# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('/store/user/dmytro/tmp/store+mc+RunIIAutumn18MiniAOD+QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8+MINIAODSIM+102X_upgrade2018_realistic_v15-v1+120000+E5F6DFC8-65CF-2A41-B0BE-E82E041CB012.root'),
    #fileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18MiniAOD/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/0DF86428-7817-0B40-AF6D-421E65B64647.root'),
   fileNames =cms.untracked.vstring(InFile),

   secondaryFileNames = cms.untracked.vstring(SecFile)

   #secondaryFileNames = cms.untracked.vstring('file:/eos/user/b/bimahaku/BstoMuMu_RECOFiles/step3_RAW2DIGI_L1Reco_RECO_corr_3.root')#'file:/eos/user/b/bimahaku/BstoMuMu_RECOFiles/step3_RAW2DIGI_L1Reco_RECO_1.root')
  
  #secondaryFileNames = cms.untracked.vstring('/store/mc/RunIIAutumn18DRPremix/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/102X_upgrade2018_realistic_v15-v2/00000/A17D7FBA-C7C3-014E-AA6F-1B653BC16087.root')


)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

OutFile = 'file:BstoMuMuFile-RunIIAutumn18_NanoAOD_Tuple.root'#%(str(num))

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(OutFile),#'file:/eos/user/b/bimahaku/BstoMuMu_NanoAOD/BstoMuMuFile-RunIIAutumn18_NanoAOD_1.root'),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v20', '')

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(2)
process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

# Automatic addition of the customisation function from Bmm5.NanoAOD.nano_cff
#from Bmm5.NanoAOD.nano_cff import nanoAOD_customizeBxToMuMu, nanoAOD_customizeV0ForMuonFake 

#call to customisation function nanoAOD_customizeBxToMuMu imported from Bmm5.NanoAOD.nano_cff
#process = nanoAOD_customizeBxToMuMu(process)

#call to customisation function nanoAOD_customizeV0ForMuonFake imported from Bmm5.NanoAOD.nano_cff
#process = nanoAOD_customizeV0ForMuonFake(process)

# End of customisation functions

# Customisation from command line

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
