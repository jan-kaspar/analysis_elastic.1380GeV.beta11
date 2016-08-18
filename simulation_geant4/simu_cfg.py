import FWCore.ParameterSet.Config as cms

process = cms.Process("elasticSimulation")

#----------------------------------------------------------------------------------------------------
# general configuration

seed = 1
events = 1E3
tag = "_1E3_seed1"

#----------------------------------------------------------------------------------------------------

# random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")
process.RandomNumberGeneratorService.g4SimHits.initialSeed = seed
process.RandomNumberGeneratorService.SimG4Object.initialSeed = seed
process.RandomNumberGeneratorService.RPSiDetDigitizer.initialSeed = seed
process.RandomNumberGeneratorService.generator.initialSeed = seed
process.RandomNumberGeneratorService.SmearingGenerator.initialSeed = seed
process.RandomNumberGeneratorService.mix.initialSeed = seed
process.RandomNumberGeneratorService.LHCTransport.initialSeed = seed

# reasonable log
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

#----------------------------------------------------------------------------------------------------
# SIMULATION

# specify the maximum events to simulate
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(int(events))
)

# elastic generator
energy = "1380"
import IOMC.Elegent.ElegentSource_cfi
process.generator = IOMC.Elegent.ElegentSource_cfi.generator
process.generator.fileName = IOMC.Elegent.ElegentSource_cfi.ElegentDefaultFileName(energy)
process.generator.t_min = 0.04
process.generator.t_max = 1.

# smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")

# magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_1380GeV_11_cfi")

# RP geometry 
process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('../user/elastic_analyses/1380GeV,beta11/simulation_geant4/RP_Dist_Beam_Cent.xml')

# G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
process.g4SimHits.Generator.HepMCProductLabel = 'generator'    # The input source for G4 module is connected to "process.source".

##  process.g4SimHits.Watchers = cms.VPSet(
##          cms.PSet(
##              type = cms.string('TotemRP'),
##              TotemRP = cms.PSet(
##                  Names = cms.vstring('TotemHitsRP'),
##                  FileName = cms.string('TotemTestRP_Hits.root'),
##                  FileNameOLD = cms.string('ntuple_sim_hits'+tag+'.root'),
##                  Verbosity = cms.bool(True)
##              )
##          )
##      )


# no pile up for the mixing module
process.load("Configuration.TotemCommon.mixNoPU_cfi")

# RP strip digitization
process.load("SimTotem.RPDigiProducer.RPSiDetConf_cfi")
process.RPSiDetDigitizer.RPVerbosity = 0

#----------------------------------------------------------------------------------------------------
# GEANT4 RECONSTRUCTION

# clusterisation
process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
process.RPClustProd.Verbosity = 0

# reco hit production
process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")
process.RPHecoHitProd.Verbosity = 0

# single-RP pattern recognition
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.verbosity = 0
process.NonParallelTrackFinder.maxHitsPerPlaneToSearch = 5
process.NonParallelTrackFinder.minPlanesPerProjectionToSearch = 3
process.NonParallelTrackFinder.minPlanesPerProjectionToFit = 3
process.NonParallelTrackFinder.threshold = 2.99

# single-RP track fitting
process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 1
process.RPSingleTrackCandCollFit.RPTrackCandCollProducer = 'NonParallelTrackFinder'

#----------------------------------------------------------------------------------------------------
# FAST (IDEAL) SIMULATION AND RECONSTRUCTION

process.load("TotemAlignment.RPFastSimulation.RPFastFullSimulation_cfi")

process.RPSingleTrackCandCollFitIdeal = process.RPSingleTrackCandCollFit.copy()
process.RPSingleTrackCandCollFitIdeal.RPTrackCandCollProducer = 'RPFastFullSimulation'

#----------------------------------------------------------------------------------------------------
# NTUPLIZATION

process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.TotemNtuplizer.outputFileName = "ntuple_full"+tag+".root"
process.TotemNtuplizer.includeDigi = True
process.TotemNtuplizer.includePatterns = True
process.TotemNtuplizer.primaryProtons = cms.bool(True)
process.TotemNtuplizer.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")

process.TotemNtuplizerIdeal = process.TotemNtuplizer.copy()
process.TotemNtuplizerIdeal.outputFileName = "ntuple_ideal"+tag+".root"
process.TotemNtuplizerIdeal.includeDigi = False
process.TotemNtuplizerIdeal.includePatterns = False
process.TotemNtuplizerIdeal.primaryProtons = False
process.TotemNtuplizerIdeal.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFitIdeal")

#----------------------------------------------------------------------------------------------------
# DEBUG TOOLS

process.OptInfo = cms.EDAnalyzer("OpticsInformation")

process.GeomInfo = cms.EDAnalyzer("GeometryInfoModule")

process.EventContentAnalyzer = cms.EDAnalyzer('EventContentAnalyzer')

#----------------------------------------------------------------------------------------------------

process.p1 = cms.Path(
    process.OptInfo
    *process.GeomInfo
    *process.generator
    *process.SmearingGenerator
    *process.g4SimHits
    *process.mix
    *process.RPSiDetDigitizer
    *process.RPClustProd
    *process.RPHecoHitProd
    *process.NonParallelTrackFinder
    *process.RPSingleTrackCandCollFit
    *process.TotemNtuplizer

    *process.RPFastFullSimulation
    *process.RPSingleTrackCandCollFitIdeal
    *process.TotemNtuplizerIdeal

    #*process.EventContentAnalyzer
)

#----------------------------------------------------------------------------------------------------
# OUTPUT

#process.o1 = cms.OutputModule("PoolOutputModule",
#    #outputCommands = cms.untracked.vstring('keep *', 'drop *_*mix*_*_*','drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*'),
#    outputCommands = cms.untracked.vstring('keep *', 'drop *_*mix*_*_*'),
#    fileName = cms.untracked.string('file:elastic_simu.root')
#)
#
#process.outpath = cms.EndPath(process.o1)
