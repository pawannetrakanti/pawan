import FWCore.ParameterSet.Config as cms

import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
hiGeneralTracks = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = (
                      cms.InputTag('hiGlobalPrimTracks'),
                      cms.InputTag('hiLowPtTripletStepTracks'),
                      cms.InputTag('hiSecondPixelTripletGlobalPrimTracks'),
                      cms.InputTag('hiPixelPairGlobalPrimTracks')
                      ),
#    hasSelector=cms.vint32(1,1,1), # Default
    hasSelector=cms.vint32(1,1,1,1), # included the hiLowPtTripletStep
    selectedTrackQuals = cms.VInputTag(
    cms.InputTag("hiInitialStepSelector","hiInitialStep"),
    cms.InputTag("hiLowPtTripletStepSelector","hiLowPtTripletStep"),     #Added by Pawan
    cms.InputTag("hiSecondPixelTripletStepSelector","hiSecondPixelTripletStep"),
    cms.InputTag("hiPixelPairStepSelector","hiPixelPairStep"),
    ),                    
#    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2), pQual=cms.bool(True)),  # should this be False? Default
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3), pQual=cms.bool(True)),  # should this be False? included hiLowPtTripletStep
      ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)
