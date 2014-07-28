package org.bp.models

import spray.json._
import DefaultJsonProtocol._

object JsonProtocol extends DefaultJsonProtocol {
  implicit val alignmentFormat = jsonFormat4(Alignment)
  implicit val methylationFormat = jsonFormat4(Methylation)
  implicit val bisulfiteConversionFormat = jsonFormat3(BisulfiteConversion)
  implicit val analysisFormat = jsonFormat8(Analysis)
}