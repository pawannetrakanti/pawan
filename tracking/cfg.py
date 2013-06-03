# Auto generated configuration file
# using: 
# Revision: 1.341 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: reco --conditions auto:starthi -s RAW2DIGI,L1Reco,RECO --scenario HeavyIons --datatier GEN-SIM-RECO --eventcontent RECODEBUG --filein file:mixed_440_B_01.root --fileout file:test.root --no_exec



import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')

ivars.register ('randomNumber',
                1,
                ivars.multiplicity.singleton,
                ivars.varType.int,
                "Random Seed")

ivars.register ('skipEvents',
                1,
                ivars.multiplicity.singleton,
                ivars.varType.int,
                "Random Seed")

ivars.skipEvents = 0
ivars.randomNumber = 1
ivars.inputFiles = 'file:/mnt/hadoop/cms/store/user/appeltel/HYDJET_b0_61X/hydjet_b0__15.root'
ivars.outputFile = './testReco_612.root'

ivars.parseArguments()


import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# To be uncommented in CMSSW6
#process.Timing = cms.Service("Timing")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(ivars.inputFiles)
)

# drop these, running HLT again
process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_TriggerResults_*_*")

# To be un-commented in CMSSW6
process.options = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.341 $'),
    annotation = cms.untracked.string('step3 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECODEBUGEventContent.outputCommands,
    fileName = cms.untracked.string(ivars.outputFile),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
        )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'STARTHI61_V13::All'


from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
    centralitySrc = cms.InputTag("hiCentrality")
)



process.load("RecoHI.HiTracking.hiIterTracking_cff")
process.heavyIonTracking *= process.hiIterTracking

process.particleFlowClusterPS.thresh_Pt_Seed_Endcap = cms.double(99999.)
process.pfTrack.UseQuality = True
process.pfTrack.TrackQuality = cms.string('highPurity')
process.pfTrack.TkColList = cms.VInputTag("hiGeneralTracks")
#process.pfTrack.GsfTracksInEvents = cms.bool(False)
process.particleFlowBlock.RecMuons = 'muons'
process.particleFlowTmp.muons = 'muons'
process.particleFlowTmp.postMuonCleaning = False

process.globalMuons.TrackerCollectionLabel = "hiGeneralTracks"
process.muons.TrackExtractorPSet.inputTrackCollection = "hiGeneralTracks"
process.muons.inputCollectionLabels = ["hiGeneralTracks", "globalMuons", "standAloneMuons:UpdatedAtVtx", "tevMuons:firstHit", "tevMuons:picky", "tevMuons:dyt"]

process.RECODEBUGoutput.outputCommands.extend(["keep *_hiGeneralTracks_*_*"])
process.RECODEBUGoutput.outputCommands.extend(["keep *_particleTowerProducer_*_*"])

process.load("RecoHI.HiJetAlgos.ParticleTowerProducer_cff")

process.reconstructionHeavyIons_withPF *= process.particleTowerProducer

# load extra trigger stuff
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
#process.load('HLTrigger.Configuration.HLT_HIon_cff')




# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstructionHeavyIons_withPF)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)

from L1Trigger.Configuration.customise_l1EmulatorFromRaw import customise
process = customise(process)

# customize the HLT to use the emulated results
#import HLTrigger.Configuration.customizeHLTforL1Emulator
#process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
#process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimGtDigis( process )
#process.hltTrigReport.HLTriggerResults = cms.InputTag('TriggerResults','','RECO' )

#process.HLTSchedule.remove(process.hltTrigReport)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECODEBUGoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.L1simulation_step,process.L1Reco_step,process.reconstruction_step)
#process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RECODEBUGoutput_step])

