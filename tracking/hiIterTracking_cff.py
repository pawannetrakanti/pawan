from RecoHI.HiTracking.hiLowPtTripletStep_cff import *    #Added by Pawan
from RecoHI.HiTracking.hiSecondPixelTripletStep_cff import *
from RecoHI.HiTracking.hiMixedTripletStep_cff import *
from RecoHI.HiTracking.hiPixelPairStep_cff import *
from RecoHI.HiTracking.MergeTrackCollectionsHI_cff import *

#hiIterTracking = cms.Sequence(hiSecondPixelTripletStep  
hiIterTracking = cms.Sequence(hiLowPtTripletStep
                              *hiSecondPixelTripletStep
                              *hiPixelPairStep
                              *hiGeneralTracks
)
