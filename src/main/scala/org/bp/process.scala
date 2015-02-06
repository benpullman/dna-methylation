package org.bp

import org.bp.models._
import org.bp.methylation._
import org.bp.fasta._
import org.bp.alignment._
import scala.util.{Try,Success,Failure}
import spray.json._
import JsonProtocol._
import java.io._
import org.apache.spark.rdd.RDD
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf

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

case class Process(name: String, references: List[DNA], map: Map[String,String], barcodeLength: Int) extends Serializable{

//class Process(name: String, references: List[DNA], map: Map[String,String], barcodeLength: Int) extends Serializable{
	def multiplex = MultiplexedSequence.create(name, references, map)

	def ref = references.flatMap(_.generateBisulfiteQuartet).map(_.addMethylation)

	//alignment and methylation

	def pre(input: DNA): Option[(DNA,String)] = {
		val tag = input.nucleobases.take(barcodeLength).foldLeft("")(_ + _)
    	val nucleobases = input.nucleobases.drop(barcodeLength)
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
		println(input.name)
		val (alignment,methylationList) = ref.map(r => (DNAAlignment.run(input,r._1),r._2)).maxBy(_._1.params.score) //Alignment
		//println(methylationList)
		//println(alignment.sample.nucleobases.length)
		//println(alignment.reference.nucleobases.length)
		val methylationPoints = MethylationProcess.checkSites(methylationList,alignment.sample,alignment.params.refStart)
		//println(methylationPoints)
		val grouped = methylationPoints.groupBy(t=>t).map(t => (t._1,t._2.length))
		//println(grouped)
		val (methylated,unmethylated) = (grouped.getOrElse("M",0),grouped.getOrElse("U",0))
		val methylation = Methylation(unmethylated, methylated, methylationPoints, methylationList.map(_+1))
		val analysis = Analysis(
			barcode = barcode,
			referenceLength = alignment.reference.length,
			sampleName = alignment.sample.name,
			referenceName = alignment.reference.name,
			direction = alignment.reference.direction.getOrElse(""),
			conversion = alignment.reference.bisulfiteConversion.getOrElse(""),
			alignment = alignment.params,
			methylation = methylation
		)
		//Sequence.compare(seq1 = alignment.sample.reverse, seq2 = alignment.reference.reverse, lineLength = 60)
		//println(analysis + "\n")
		analysis
  	}

	//combine back

	def combine(multiplexedSeq: MultiplexedSequence, analysis: Analysis) = {
		//println(analysis)
		//println(multiplexedSeq)
		multiplexedSeq.addSample(analysis)
	}

	//def run(samples: RDD[DNA]) = samples.map(pre(_)).flatMap(_.map{ case (sample,barcode) => align(sample,barcode)}).filter(_.percentAligned > .9)//.fold(multiplex)(combine(_,_))

}

object Run{

	//def run(samples: RDD[DNA], process: Process) = samples.map(process.pre(_)).flatMap(_.map{ case (sample,barcode) => process.align(sample,barcode)}).filter(_.percentAligned > .9)//.fold(multiplex)(combine(_,_))

	def main(args: Array[String]): Unit = {

  		val conf = new SparkConf().setMaster("local[20]").setAppName("HITMAP")
    	val sc = new SparkContext(conf)

		val out = new PrintWriter(new File("results.txt"))

		val sampleFile = "data/mssm-samples/new_mini.txt"

		//val sampleData = sc.textFile(sampleFile, 2).cache()


	    val references = FASTA.read("data/mssm-samples/new_refer.txt")
	    val samples = FASTA.read("data/mssm-samples/new_mini.txt")
	    val map = MultiplexedSequence.importMap("data/mssm-samples/barcode_new.csv")

	    //val references = FASTA.read("data/quma-samples/sample_genome_fasta.txt")
	    //val samples = FASTA.read("data/quma-samples/Gm9_J1_seq_fasta.txt")
	    //val map = MultiplexedSequence.importMap("data/quma-samples/easy_example_barcode.csv")
	    //.map(sample => out.write(sample.toString))
	    // out.close()
	    //sampleData.map(samples => println(

	    //val process = new Process("sample",references,map,4)
	    val process = Process("sample",references,map,4)
	    //process.run(sc.parallelize(samples)).foreach(println)
	    sc.parallelize(samples,20).map(sample => process.pre(sample)).flatMap(_.map{ case (sample,barcode) => process.align(sample,barcode)}).foreach(println)
	    //sampleData.map(samples => println(process.run(FASTA.read(samples)).toJson))

	}

}