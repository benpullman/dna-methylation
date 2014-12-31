package org.bp.models

import spray.json._
import DefaultJsonProtocol._

object JsonProtocol extends DefaultJsonProtocol {
  /*
  implicit val nucleobaseFormat = jsonFormat1(Nucleobase.apply)
  implicit val dnaFormat = jsonFormat4(DNA.apply)
  implicit val alignmentFormat = jsonFormat4(AlignmentAnalysis)
  implicit val bisulfiteConversionFormat = jsonFormat3(BisulfiteConversion)
  implicit val analysisFormat = jsonFormat10(Analysis)
  */
  //case class Analysis(barcode: String, sampleName: String, referenceName: String, direction: String, conversion: String, alignment: AlignmentParameters, methylation: Methylation)
  implicit val alignmentParametersFormat = jsonFormat6(AlignmentParameters.apply)
  implicit val methylationFormat = jsonFormat4(Methylation.apply)
  implicit val analysisFormat = jsonFormat8(Analysis.apply)
  implicit val sampleSummaryFormat = jsonFormat1(SampleSummary.apply)
  implicit val sampleListFormat = jsonFormat3(SampleList.apply)
  implicit val regionListFormat = jsonFormat1(RegionList.apply)
  implicit val multiplexedFormat = jsonFormat3(MultiplexedSequence.apply)

}
