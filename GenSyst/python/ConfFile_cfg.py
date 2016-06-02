import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('doExclusive', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Apply Exclusive Selection')
options.register ('jetPtCut', 30, VarParsing.multiplicity.singleton, VarParsing.varType.float, 'Jet Pt Cut')
#options.register ('jetHigherPtCut', 30, VarParsing.multiplicity.singleton, VarParsing.varType.float, 'Jet Pt Cut')
options.parseArguments()
print 'doExclusive =  ', options.doExclusive
print 'Jet Pt Cut = ', options.jetPtCut

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
)

process.demo = cms.EDAnalyzer('GenSystInclusive',
			      jetPt = cms.double(options.jetPtCut),
			      doExclusive = cms.bool(options.doExclusive))
#else:
#        process.demo = cms.EDAnalyzer('GenSystExclusive',
#                                      jetPt = cms.double(options.jetPtCut),
#                                      jetHigherPt = cms.double(options.jetHigherPtCut))
        
process.p = cms.Path(process.demo)
