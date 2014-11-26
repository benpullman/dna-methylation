/*package org.bp.models

import spray.json._
import DefaultJsonProtocol._

object JsonProtocol extends DefaultJsonProtocol {
  implicit val nucleobaseFormat = jsonFormat1(Nucleobase.apply)
  implicit val dnaFormat = jsonFormat4(DNA.apply)
  implicit val alignmentFormat = jsonFormat4(AlignmentAnalysis)
  implicit val methylationFormat = jsonFormat4(Methylation)
  implicit val bisulfiteConversionFormat = jsonFormat3(BisulfiteConversion)
  implicit val analysisFormat = jsonFormat10(Analysis)
  implicit val sampleListFormat = jsonFormat2(SampleList.apply)
  implicit val regionListFormat = jsonFormat2(RegionList.apply)
  implicit val multiplexedFormat = jsonFormat1(MultiplexedSequence.apply)
}
*/