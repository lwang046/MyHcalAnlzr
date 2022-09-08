import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process('MyHcalAnlzr',eras.Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
#process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")

process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/030/00000/38c224b2-b659-4acb-90f9-705bf9cc120a.root",
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/092/00000/12c18f49-2893-412a-9bdc-d5424e89f3f6.root",
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/152/00000/08d8d6df-d609-4b03-9188-6834ae88a976.root",
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/155/00000/0a4e7f59-3130-4a3d-9897-ca3ec7e69121.root",
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/164/00000/0afce2c7-1c39-4e04-9a89-267896298b9f.root",
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/209/00000/6e99adcd-19d3-42ff-b82c-7a6380699c49.root",
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/253/00000/a8000a34-8d76-49ea-8557-3b76cdfedcae.root",
		"file:/eos/cms/tier0/store/data/Run2022D/Cosmics/RAW/v1/000/358/283/00000/74d64b7e-478d-4707-a274-757bbbad24a2.root",
        )
)

process.options = cms.untracked.PSet(
#	SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run3_hlt']

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)


process.MyHcalAnlzr = cms.EDAnalyzer('MyHcalAnlzr',
        #tagRecHit = cms.untracked.InputTag("hbheprereco"),
	#tagQIE11 = cms.InputTag("simHcalDigis", "HBHEQIE11DigiCollection")
	tagQIE11 = cms.untracked.InputTag("hcalDigis"),
        #EEdigiCollection = cms.InputTag("ecalDigis","eeDigis"),
	#EERecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        runtype = cms.untracked.string("Global")
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("output_CosmicsRuns.root"),
      closeFileFast = cms.untracked.bool(True)
)

process.ecalRecoSequence = cms.Sequence((process.ecalMultiFitUncalibRecHit+process.ecalDetIdToBeRecovered+process.ecalRecHit))

process.p = cms.Path(
	process.bunchSpacingProducer*
	#process.ecalDigis*
	#process.ecalRecoSequence*
        process.hcalDigis*
	#process.hcalLocalRecoSequence*
	process.MyHcalAnlzr
)

