import FWCore.ParameterSet.Config as cms

process = cms.Process("BuildElasticCorrectionsFile")
# minimum of logs
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# empty source
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# geometry
process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/$dataset/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

# include alignments, if any
process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring(
    "TotemAlignment/RPData/LHC/$dataset/sr+hsx/45_220.xml",
    "TotemAlignment/RPData/LHC/$dataset/sr+hsx/56_220.xml",
)

process.builder = cms.EDAnalyzer("BuildElasticCorrectionsFile",
    inputFileName = cms.string("$dir/alignment_fit.out"),
    outputFileName = cms.string("$dir/corrections.xml")
)

process.p = cms.Path(process.builder)
