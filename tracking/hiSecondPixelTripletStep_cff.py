import FWCore.ParameterSet.Config as cms

# Filter on quality tracks
# Commented by Pawan
#hiFirstStepFilter = cms.EDProducer("QualityFilter",
#                                   TrackQuality = cms.string('highPurity'),
#                                   recTracks = cms.InputTag("hiSelectedTracks")
#)

# NEW CLUSTERS (remove previously used clusters)
hiSecondPixelTripletClusters = cms.EDProducer("TrackClusterRemover",
                                              clusterLessSolution= cms.bool(True),
#                                              trajectories = cms.InputTag("hiFirstStepFilter"),
                                              oldClusterRemovalInfo = cms.InputTag("hiLowPtTripletStepClusters"),
                                              trajectories = cms.InputTag("hiLowPtTripletStepTracks"),
                                              overrideTrkQuals = cms.InputTag('hiLowPtTripletStepSelector','hiLowPtTripletStep'),
                                              TrackQuality = cms.string('highPurity'),
                                              pixelClusters = cms.InputTag("siPixelClusters"),
                                              stripClusters = cms.InputTag("siStripClusters"),
                                              Common = cms.PSet(
    maxChi2 = cms.double(9.0)
    ),
    Strip = cms.PSet(
    #Yen-Jie's mod to preserve merged clusters
    maxSize = cms.uint32(2),
    maxChi2 = cms.double(9.0)
    )
)


# SEEDING LAYERS
import RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi
hiSecondPixelTripletSeedLayers = RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi.pixellayertriplets.clone(
        ComponentName = 'hiSecondPixelTripletSeedLayers'
            )
hiSecondPixelTripletSeedLayers.BPix.skipClusters = cms.InputTag('hiSecondPixelTripletClusters')
hiSecondPixelTripletSeedLayers.FPix.skipClusters = cms.InputTag('hiSecondPixelTripletClusters')

# SEEDS
import RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock
hiSecondPixelTripletSeeds = RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff.globalSeedsFromTriplets.clone(
    RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
    ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
    RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
    ptMin = cms.double(4.0),
    originRadius = cms.double(0.005),
    nSigmaZ = cms.double(4.0)
    )
    )
)

hiSecondPixelTripletSeeds.OrderedHitsFactoryPSet.SeedingLayers = 'hiSecondPixelTripletSeedLayers'
hiSecondPixelTripletSeeds.OrderedHitsFactoryPSet.GeneratorPSet.maxElement = cms.uint32(5000000)
hiSecondPixelTripletSeeds.ClusterCheckPSet.MaxNumberOfPixelClusters = cms.uint32(5000000)
hiSecondPixelTripletSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = cms.uint32(50000000)

from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
hiSecondPixelTripletSeeds.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet.ComponentName = 'LowPtClusterShapeSeedComparitor'


# QUALITY CUTS DURING TRACK BUILDING
import TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi
hiSecondPixelTripletTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.clone(
    ComponentName = 'hiSecondPixelTripletTrajectoryFilter',
    filterPset = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.filterPset.clone(
    maxLostHits = cms.int32(1),
    minimumNumberOfHits = cms.int32(6),
    minPt = cms.double(1.0)
  )
)

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi
hiSecondPixelTripletChi2Est = TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('hiSecondPixelTripletChi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(9.0)
)


# TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi
hiSecondPixelTripletTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi.GroupedCkfTrajectoryBuilder.clone(
    ComponentName = 'hiSecondPixelTripletTrajectoryBuilder',
    MeasurementTrackerName = '',
    trajectoryFilterName = 'hiSecondPixelTripletTrajectoryFilter',
    clustersToSkip = cms.InputTag('hiSecondPixelTripletClusters'),
    maxCand = cms.int32(3),
    #estimator = cms.string('hiSecondPixelTripletChi2Est')
)


# MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
hiSecondPixelTripletTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('hiSecondPixelTripletSeeds'),
    TrajectoryBuilder = 'hiSecondPixelTripletTrajectoryBuilder',
    doSeedingRegionRebuilding = True,
    useHitsSplitting = True
)

# TRACK FITTING
import RecoTracker.TrackProducer.TrackProducer_cfi
hiSecondPixelTripletGlobalPrimTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'hiSecondPixelTripletTrackCandidates',
    AlgorithmName = cms.string('iter1')
)


# Final selection
import RecoHI.HiTracking.hiMultiTrackSelector_cfi
hiSecondPixelTripletStepSelector = RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiMultiTrackSelector.clone(
    src='hiSecondPixelTripletGlobalPrimTracks',
    trackSelectors= cms.VPSet(
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiLooseMTS.clone(
    name = 'hiSecondPixelTripletStepLoose',
    ), #end of pset
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiTightMTS.clone(
    name = 'hiSecondPixelTripletStepTight',
    preFilterName = 'hiSecondPixelTripletStepLoose',
    ),
    RecoHI.HiTracking.hiMultiTrackSelector_cfi.hiHighpurityMTS.clone(
    name = 'hiSecondPixelTripletStep',
    preFilterName = 'hiSecondPixelTripletStepTight',
    min_nhits = cms.uint32(14)
    ),
 ) #end of vpset
) #end of clone


import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
hiSecondQual = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('hiSecondPixelTripletGlobalPrimTracks')),
    hasSelector=cms.vint32(1),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("hiSecondPixelTripletStepSelector","hiSecondPixelTripletStep")),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False),
    #writeOnlyTrkQuals = True
)

# Final sequence
# commented by Pawan
#hiSecondPixelTripletStep = cms.Sequence(hiFirstStepFilter
#                                        *hiSecondPixelTripletClusters
hiSecondPixelTripletStep = cms.Sequence(hiSecondPixelTripletClusters
                                        *hiSecondPixelTripletSeeds
                                        *hiSecondPixelTripletTrackCandidates
                                        *hiSecondPixelTripletGlobalPrimTracks
                                        *hiSecondPixelTripletStepSelector
                                        *hiSecondQual
)
