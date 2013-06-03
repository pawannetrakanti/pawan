import FWCore.ParameterSet.Config as cms

#Added by pawan
hiFirstStepFilter = cms.EDProducer("QualityFilter",
                                 TrackQuality = cms.string('highPurity'),
                                 recTracks = cms.InputTag("hiSelectedTracks")
)


# NEW CLUSTERS (remove previously used clusters)
hiLowPtTripletStepClusters = cms.EDProducer("TrackClusterRemover",
                                            clusterLessSolution= cms.bool(True),
                                            # Heavy-ion
                                            trajectories = cms.InputTag("hiFirstStepFilter"),
                                            TrackQuality = cms.string('highPurity'),
                                            minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
                                            pixelClusters = cms.InputTag("siPixelClusters"),
                                            stripClusters = cms.InputTag("siStripClusters"),
                                            Common = cms.PSet(
    # Added by pawan
    maxSize = cms.uint32(2),
    maxChi2 = cms.double(9.0)
  )
)

# SEEDING LAYERS
import RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi
hiLowPtTripletStepSeedLayers = RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi.pixellayertriplets.clone(
    ComponentName = 'hiLowPtTripletStepSeedLayers'
)
hiLowPtTripletStepSeedLayers.BPix.skipClusters = cms.InputTag('hiLowPtTripletStepClusters')
hiLowPtTripletStepSeedLayers.FPix.skipClusters = cms.InputTag('hiLowPtTripletStepClusters')

# SEEDS
import RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock
hiLowPtTripletStepSeeds = RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff.globalSeedsFromTriplets.clone(
    RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
    ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
    RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(

# Default
    ptMin = cms.double(0.2),
#Playing around
#    ptMin = 0.25,
#    ptMin = 0.30,
#    ptMin = 1.0,

# Default
     originRadius = cms.double(0.02),
#Playing around
#    originRadius = cms.double(0.01),
#    originRadius = cms.double(0.005),
#    originRadius = 0.001,

# Default
    nSigmaZ = cms.double(4.0),
  )
 )
)

hiLowPtTripletStepSeeds.OrderedHitsFactoryPSet.SeedingLayers = 'hiLowPtTripletStepSeedLayers'

#Added by pawan
# to avoid 'too many clusters'
#lowPtTripletStepSeeds.ClusterCheckPSet.doClusterCheck = cms.bool(False)
#lowPtTripletStepSeeds.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = cms.uint32(0)

hiLowPtTripletStepSeeds.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = cms.uint32(5000000)
hiLowPtTripletStepSeeds.ClusterCheckPSet.MaxNumberOfPixelClusters = cms.uint32(5000000)
hiLowPtTripletStepSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = cms.uint32(50000000)

from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
hiLowPtTripletStepSeeds.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'


# QUALITY CUTS DURING TRACK BUILDING
import TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi
hiLowPtTripletStepStandardTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.clone(
    ComponentName = 'hiLowPtTripletStepStandardTrajectoryFilter',
    filterPset = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.filterPset.clone(
    minimumNumberOfHits = cms.int32(3),
    minPt = cms.double(0.075)
 )
)

from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeTrajectoryFilterESProducer_cfi import *
# Composite filter
import TrackingTools.TrajectoryFiltering.CompositeTrajectoryFilterESProducer_cfi
hiLowPtTripletStepTrajectoryFilter = TrackingTools.TrajectoryFiltering.CompositeTrajectoryFilterESProducer_cfi.compositeTrajectoryFilterESProducer.clone(
    ComponentName = cms.string('hiLowPtTripletStepTrajectoryFilter'),
    filterNames   = cms.vstring('lowPtTripletStepStandardTrajectoryFilter',
                                'clusterShapeTrajectoryFilter'
 )
)

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi
hiLowPtTripletStepChi2Est = TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('hiLowPtTripletStepChi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(9.0)
)

# TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi
hiLowPtTripletStepTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi.GroupedCkfTrajectoryBuilder.clone(
    ComponentName = 'hiLowPtTripletStepTrajectoryBuilder',
    MeasurementTrackerName = '',
    trajectoryFilterName = 'hiLowPtTripletStepTrajectoryFilter',
    clustersToSkip = cms.InputTag('hiLowPtTripletStepClusters'),
# Added by Pawan
#    maxLostHits = cms.int32(1),
#    Default    
    maxCand = 4,

# Added by Pawan
#    maxCand = cms.int32(1), # same as in hiSecondPixelTripletStep
    estimator = cms.string('hiLowPtTripletStepChi2Est'),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    # 0.63 GeV is the maximum pT for a charged particle to loop within the 1.1m radius
    # of the outermost Tracker barrel layer (with B=3.8T)
    maxPtForLooperReconstruction = cms.double(0.7) 
)

# MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
hiLowPtTripletStepTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('hiLowPtTripletStepSeeds'),
    ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),
    TrajectoryBuilder = 'hiLowPtTripletStepTrajectoryBuilder',
    doSeedingRegionRebuilding = True,
    useHitsSplitting = True,
# added by pawan
    maxNSeeds = cms.uint32(500000),  # default was 10000
    maxSeedsBeforeCleaning = cms.uint32(10000) # default was 1000
)


# TRACK FITTING
import RecoTracker.TrackProducer.TrackProducer_cfi
hiLowPtTripletStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'hiLowPtTripletStepTrackCandidates',
    AlgorithmName = cms.string('iter1'),
    Fitter = cms.string('FlexibleKFFittingSmoother')
)

from TrackingTools.TrajectoryCleaning.TrajectoryCleanerBySharedHits_cfi import trajectoryCleanerBySharedHits
hiLowPtTripletStepTrajectoryCleanerBySharedHits = trajectoryCleanerBySharedHits.clone(
        ComponentName = cms.string('hiLowPtTripletStepTrajectoryCleanerBySharedHits'),
        fractionShared = cms.double(0.16),
        allowSharedFirstHit = cms.bool(True)
)
hiLowPtTripletStepTrackCandidates.TrajectoryCleaner = 'hiLowPtTripletStepTrajectoryCleanerBySharedHits'

# Final selection
import RecoHI.HiTracking.hiMultiTrackSelector_cfi
hiLowPtTripletStepSelector = RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiMultiTrackSelector.clone(
    src='hiLowPtTripletStepTracks',
    trackSelectors= cms.VPSet(
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiLooseMTS.clone(
    name = 'hiLowPtTripletStepLoose',
    ), #end of pset
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
    name = 'hiLowPtTripletStepTight',
    preFilterName = 'hiLowPtTripletStepLoose',
    ),
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
    name = 'hiLowPtTripletStep',
    preFilterName = 'hiLowPtTripletStepTight',
  ),
 ) #end of vpset
) #end of clone


# Final sequence
hiLowPtTripletStep = cms.Sequence(hiFirstStepFilter
                                  *hiLowPtTripletStepClusters
                                  *hiLowPtTripletStepSeeds
                                  *hiLowPtTripletStepTrackCandidates
                                  *hiLowPtTripletStepTracks
                                  *hiLowPtTripletStepSelector
)
