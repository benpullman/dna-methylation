package org.bp.models

import scala.io.Source

case class MultiplexedSequence(regions: List[RegionList]){
	def addSample(newSequence: DNA, references: List[DNA], sampleMap: Map[String,String]) = {
		val tag = newSequence.nucleobases.take(18).foldLeft("")(_ + _)
		val nucleobases = newSequence.nucleobases.drop(18)
		val barcode = sampleMap.getOrElse(tag,"none")
		//println(tag + " " + barcode)
		val analysis = Sequence.toAnalysis(newSequence,references,5,-3,-20,-2)
		val updatedRegion = regions.filter(_.name == analysis.referenceName).head.addSample(barcode,analysis)
		this.copy(regions = updatedRegion :: this.regions.filter(_.name != analysis.referenceName))
	}
}

case class RegionList(name: String, samples: List[SampleList]){
	def addSample(sampleName: String, sample: Analysis): RegionList = {
		val newSample = samples.filter(_.name == sampleName).head.addSample(sample)
		this.copy(samples =  newSample :: this.samples.filter(_.name != sampleName))
	}
}

case class SampleList(name: String, analyses: List[Analysis]){
	def addSample(analysis: Analysis): SampleList = {
		this.copy(analyses = analysis :: this.analyses)
	}
}

object MultiplexedSequence {
	def create(references: List[DNA], sampleMap: Map[String,String]) = {
		val regions = references.map(reference => RegionList(reference.name,SampleList("none",List()) :: sampleMap.keys.toList.map(barcode => SampleList(barcode,List()))))
		MultiplexedSequence(regions)
	}
	def importMap(filename: String): Map[String, String] = {
		Source.fromFile(filename).getLines.foldLeft[Map[String,String]](Map())((values,toAdd) => {
			val toAddArray = toAdd.split(",")
			values ++ Map(toAddArray(0) -> toAddArray(1).toUpperCase)
		})
	}
}