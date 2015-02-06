package org.bp.models

import scala.io.Source
import spray.json._
import DefaultJsonProtocol._
import org.bp._
import scala.util.{Try,Success,Failure}


case class AlignmentParameters(score: Int, start: Int, end: Int, refStart: Int, gaps: Int, mismatches: Int) extends Serializable{
	def length: Int = this.end - this.start
	def percentAligned: Double = (1 - this.mismatches.toDouble/this.length.toDouble)
}

case class Methylation(unmethylated: Int, methylated: Int, sequence: IndexedSeq[String],reference: IndexedSeq[Int]) extends Serializable{
	def methylationPercent = Try(this.methylated / (this.unmethylated + this.methylated))
}

case class Analysis(barcode: String, referenceLength: Int, sampleName: String, referenceName: String, direction: String, conversion: String, alignment: AlignmentParameters, methylation: Methylation) extends Serializable{
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

case class AlignmentProcess(scoringMatrix: Array[Array[Int]],directionMatrix: Array[Array[String]], sample: String, endX:Int,endY:Int, score: Int) extends Serializable

case class MultiplexedSequence(id: Option[Int], name: String, regions: Map[String,RegionList]) extends Serializable{
	def addSample(analysis: Analysis) = {
		//val tag = newSequence.nucleobases.take(18).foldLeft("")(_ + _)
		//val nucleobases = newSequence.nucleobases.drop(18)
		//val barcode = sampleMap.getOrElse(tag,"none")
		//val barcode = MultiplexedSequence.findBestBarcode(tag,sampleMap)
		//println(tag + " " + barcode)
		//val analysis = Sequence.toAnalysis(newSequence,references,5,-3,-20,-2,barcode)
		val updatedRegion = regions.get(analysis.referenceName).map(_.addSample(analysis)).get
		//println(tag)
		this.copy(regions = this.regions + (analysis.referenceName -> updatedRegion))
	}
}
/*
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
		val regions = references.map(reference => RegionList(reference.name,SampleList("error",List()) :: sampleMap.keys.toList.map(barcode => SampleList(barcode,List()))))
		MultiplexedSequence(regions)
	}
	
	def jsonImport(filename: String, references: List[DNA], sampleMap: Map[String,String]): MultiplexedSequence = {
		import JsonProtocol._
		val rawJson = Source.fromFile(filename).mkString
		val analysisList = rawJson.parseJson.convertTo[List[Analysis]]
		val multi = create(references,sampleMap)
		analysisList.foldLeft(multi)((multiplexed,analysis) => {
			val updatedRegion = multiplexed.regions.filter(_.name == analysis.referenceName).head.addSample(analysis.barcode,analysis)
		//println(tag)
			multiplexed.copy(regions = updatedRegion :: multiplexed.regions.filter(_.name != analysis.referenceName))
		})
	}
	*/

case class RegionList(samples: Map[String,SampleList]) extends Serializable {
	def addSample(sample: Analysis): RegionList = {
		samples.get(sample.barcode).map(_.addSample(sample)).map(newSample => this.copy(samples =  this.samples + (sample.barcode -> newSample))).getOrElse(this)
	}
}

case class SampleSummary(name: String) extends Serializable

case class SampleList(summary: Option[SampleSummary] = None, methylation: Option[Double] = None, analyses: List[Analysis])extends Serializable{
	def addSample(analysis: Analysis): SampleList = {
		this.copy(analyses = analysis :: this.analyses)
	}
}

object MultiplexedSequence {
	def create(name: String, references: List[DNA], sampleMap: Map[String,String]) = {
		//val regions = references.map(reference => RegionList((reference.name -> sampleMap.keys.toList.map(barcode => SampleList((barcode -> List()).toMap)).toMap)))
		val regions = references.foldLeft[Map[String,RegionList]](Map())((racc,rnew) => {
			racc + (rnew.name -> RegionList({
				sampleMap.keys.toList.foldLeft[Map[String,SampleList]](Map())((sacc,snew) => {
					sacc + (snew -> SampleList(None,None,List()))
				})
			}))
		})
		MultiplexedSequence(None, name, regions)
	}

	def importMap(filename: String): Map[String, String] = {
		Source.fromFile(filename).getLines.foldLeft[Map[String,String]](Map())((values,toAdd) => {
			val toAddArray = toAdd.split(",")
			values ++ Map(toAddArray(0) -> toAddArray(1).toUpperCase)
		})
	}
	/*
	def combine(seq1: MultiplexedSequence, seq2: MultiplexedSequence) = {
		MultiplexedSequence(seq1.regions ::: seq2.regions)
	}
	def refWithError(input: String, sampleMap: Map[String,String]) = {
		val samples = sampleMap.keys
		case class SampleCheck(sample: String, score: Int)
		def checkSample(input: String, sample: String): SampleCheck = {
			val score = (0 to 17).foldLeft(0)((acc,next) => if(input(next)==sample(next)) acc+1 else acc)
			SampleCheck(sample,score)
		}
		val bestSample = samples.map(checkSample(input,_)).maxBy(_.score)
		//println(bestSample.score)
		sampleMap.getOrElse(bestSample.sample,"none")
	}
	*/
	def findBestBarcode(input: String, sampleMap: Map[String,String]) = {
		val bestSample = sampleMap.keys.map(barcodeAlignment(input,_,10,-1,-5,-5)).maxBy(_.score)
		def completeAlignment(input: String, alignmentProcess: AlignmentProcess) = {
			var start = 0
			var end = alignmentProcess.endX
			val sample = alignmentProcess.sample
			def read(i: Int,j: Int,directionMatrix:Array[Array[String]]): List[(Char, Char)] = {
				directionMatrix(i)(j) match {
					case "diag" => (input(i-1),sample(j-1)) :: read(i-1,j-1,directionMatrix)
					case "upup" => (input(i-1),'-') :: read(i-1,j,directionMatrix)
					case "left" => ('-',sample(j-1)) :: read(i,j-1,directionMatrix)
					case _ => {
						start = i+1
						List(('s','s'))
					}
				}
			}
			read(alignmentProcess.endX,alignmentProcess.endY,alignmentProcess.directionMatrix).reverse.drop(1)
		}

		def getErrors(pairs: List[(Char,Char)]): Int = {
			val error = pairs.foldLeft(0)((tally,next) => {
				val sample = next._2
				val input = next._1
				if (sample==input){
					tally
				} else {
					tally + 1
				}
			})
			error
		}

		val error = getErrors(completeAlignment(input,bestSample))
		//println(error)
		if (error < 4) {
			sampleMap.getOrElse(bestSample.sample,"error")
		} else {
			"error"
		}
	}
	def barcodeAlignment(input: String, sample: String, hit: Int, miss: Int, gapUp: Int, gapLeft: Int) = {
		var largest = (0,(0,0))
		var start = 0
		var end = 0
		var directionMatrix = Array.ofDim[String](input.length + 1,sample.length+1)
		var scoringMatrix = Array.ofDim[Int](input.length + 1,sample.length+1)
		def initializeScoringMatrix: Unit = {
			for (i <- 1 to input.length; j <- 1 to sample.length){
				scoringMatrix(i)(j) = -1
			}
		}
		def matchScore(i: Int, j: Int): Int = {
			if (input(i-1) == sample(j-1)) hit else miss
		}
		def score(i:Int,j:Int): Int = {
			scoringMatrix(i)(j) match {
				//case -1 => List(score(i-1,j-1) + matchScore(i,j), score(i-1,j) + gap, score(i,j-1) + gap, 0).max
				case -1 => {
					val diag = score(i-1,j-1) + matchScore(i,j)
					val left = score(i,j-1) + gapLeft
					val up = score(i-1,j) + gapUp
					val best = List(diag, left, up, 0).max
					if (best > largest._1) largest = (best,(i,j))
					best match {
						case d if d == diag => {
							directionMatrix(i)(j) = "diag"
							d
						}
						case l if l == left => {
							directionMatrix(i)(j) = "left"
							l
						}
						case u if u == up => {
							directionMatrix(i)(j) = "upup"
							u
						}
						case 0 => {
							directionMatrix(i)(j) = "stop"
							0
						}
						case _ => -100000
					}
				}
				case a => a
			}
		}
		def setScore: Unit = {
			for (i <- 0 to input.length; j <- 0 to sample.length){
				scoringMatrix(i)(j) = score(i,j)
			}
		}

		initializeScoringMatrix
		setScore
		AlignmentProcess(scoringMatrix, directionMatrix, sample, largest._2._1,largest._2._2, largest._1)
		//directionMatrix.foreach{case i => i foreach {j => print(j + " ")}; print('\n')}
		//end = largest._2._1
		//println(largest._1)
		//(read(largest._2._1,largest._2._2).reverse.drop(1),(start,end))
	}
}