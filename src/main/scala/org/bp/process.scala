package org.bp

import org.bp.models._
import org.bp.methylation._
import org.bp.fasta._
import org.bp.alignment._
import scala.util.{Try,Success,Failure}


case class Methylation(unmethylated: Int, methylated: Int, sequence: IndexedSeq[String],reference: IndexedSeq[Int]){
	def methylationPercent = Try(this.methylated / (this.unmethylated + this.methylated))
}

case class Analysis(barcode: String, sampleName: String, referenceName: String, direction: String, conversion: String, alignment: AlignmentParameters, methylation: Methylation){
	def length: Int = this.alignment.length
	def percentAligned: Double = this.alignment.percentAligned
	def methylationPercent = this.methylation.methylationPercent
	override def toString(): String = {
		"Analysis for " + sampleName + "--\n" +
		"\tBarcode: " + barcode + "\n" +
		"\tReference: " + referenceName + "\n" +
		"\tPercent Methylated: " + methylation.methylated + "/" + (methylation.methylated + methylation.unmethylated) + "\n" +
		"\tPercent Aligned: " + percentAligned
	}
	
}

/*
case class DNAAlignment(sample: DNA, reference: DNA, params: AlignmentParameters){
	def length: Int = this.params.length
	def percentAligned: Double = this.params.percentAligned
}
case class AlignmentParameters(score: Int, start: Int, end: Int, refStart: Int, gaps: Int, mismatches: Int){
	def length: Int = this.start - this.end
	def percentAligned: Double = this.mismatches.toDou


//read in map and reference files
*/
class Process(name: String, references: List[DNA], map: Map[String,String]){

	val multiplex = MultiplexedSequence.create(name, references, map)

	val ref = references.flatMap(_.generateBisulfiteQuartet).map(_.addMethylation)
	//alignment and methylation

	def pre(input: DNA): Option[(DNA,String)] = {
		val tag = input.nucleobases.take(18).foldLeft("")(_ + _)
    	val nucleobases = input.nucleobases.drop(18)
    	//val barcode = sampleMap.getOrElse(tag,"none")
    	val barcode = Sequence.findBestBarcode(tag,map)
    	//println(tag + " " + barcode)
    	if (barcode == "error"){
      		None
    	}else{
      		Some(input,barcode)
    	}
	}

	def align(input: DNA, barcode: String): Analysis = {
		val (alignment,methylationList) = ref.map(r => (DNAAlignment.run(input,r._1),r._2)).maxBy(_._1.params.score) //Alignment
		println(methylationList)
		println(alignment.sample.nucleobases.length)
		println(alignment.reference.nucleobases.length)
		val methylationPoints = MethylationProcess.checkSites(methylationList,alignment.sample,alignment.params.refStart)
		println(methylationPoints)
		val grouped = methylationPoints.groupBy(t=>t).map(t => (t._1,t._2.length))
		println(grouped)
		val (methylated,unmethylated) = (grouped.getOrElse("M",0),grouped.getOrElse("U",0))
		val methylation = Methylation(unmethylated, methylated, methylationPoints, methylationList.map(_+1))
		val analysis = Analysis(
			barcode = barcode,
			sampleName = alignment.sample.name,
			referenceName = alignment.reference.name,
			direction = alignment.reference.direction.getOrElse(""),
			conversion = alignment.reference.bisulfiteConversion.getOrElse(""),
			alignment = alignment.params,
			methylation = methylation
		)
		println(analysis + "\n")
		analysis
  	}

	//combine back

	def combine(multiplexedSeq: MultiplexedSequence, analysis: Analysis) = {
		multiplexedSeq.addSample(analysis)
	}

	def run(samples: List[DNA]) = samples.map(pre(_)).flatMap(_.map{ case (sample,barcode) => align(sample,barcode)}).filter(_.percentAligned > .9).foldLeft(multiplex)(combine(_,_))

}

object Run{

	def main(args: Array[String]): Unit = {

	    val references = FASTA.read("data/mssm-samples/new_refer.txt")
	    val samples = FASTA.read("data/mssm-samples/new_mini.txt")
	    val map = MultiplexedSequence.importMap("data/mssm-samples/barcode_new.csv")

	    val process = new Process("sample",references,map)
	    println(process.run(samples))

	}

}