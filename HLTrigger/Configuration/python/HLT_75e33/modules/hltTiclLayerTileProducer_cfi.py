import FWCore.ParameterSet.Config as cms

hltTiclLayerTileProducer = cms.EDProducer("TICLLayerTileProducer",
    detector = cms.string('HGCAL'),
    layer_HFNose_clusters = cms.InputTag("hgcalLayerClustersHFNose"),
    layer_clusters = cms.InputTag("hltMergeLayerClusters"),
    mightGet = cms.optional.untracked.vstring
)
