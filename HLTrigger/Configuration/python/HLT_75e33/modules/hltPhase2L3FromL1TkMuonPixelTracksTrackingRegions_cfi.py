import FWCore.ParameterSet.Config as cms

hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions = cms.EDProducer("CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet(
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        deltaEta = cms.double(0.035),
        deltaPhi = cms.double(0.02),
        input = cms.InputTag("l1tTkMuonsGmt"),
        maxNRegions = cms.int32(10000),
        maxNVertices = cms.int32(1),
        measurementTrackerName = cms.InputTag(""),
        mode = cms.string('BeamSpotSigma'),
        nSigmaZBeamSpot = cms.double(4.0),
        nSigmaZVertex = cms.double(3.0),
        originRadius = cms.double(0.2),
        precise = cms.bool(True),
        ptMin = cms.double(2.0),
        searchOpt = cms.bool(False),
        vertexCollection = cms.InputTag("notUsed"),
        whereToUseMeasurementTracker = cms.string('Never'),
        zErrorBeamSpot = cms.double(24.2),
        zErrorVetex = cms.double(0.2)
    )
)

from Configuration.ProcessModifiers.phase2L3MuonsOIFirst_cff import phase2L3MuonsOIFirst
phase2L3MuonsOIFirst.toModify(
    hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions.RegionPSet,
    input = "hltPhase2L3MuonFilter:L1TkMuToReuse"
)
