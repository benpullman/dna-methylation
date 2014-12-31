package org.bp.alignment

import org.bp.models._

case class DNAAlignment(sample: DNA, reference: DNA, params: AlignmentParameters){
	def length: Int = this.params.length
	def percentAligned: Double = this.params.percentAligned
}

case class Alignment[T](seq1: IndexedSeq[T], seq2: IndexedSeq[T], params: AlignmentParameters){
	def length: Int = this.params.length
	def percentAligned: Double = this.params.percentAligned
}

case class ScoredPoint(score: Int, point: Point)

case class Point(x: Int, y: Int)

object DNAAlignment{

	val sw = new SmithWaterman[Nucleobase](Nucleobase("-"),5,-3,-5,-2,matchBase)

	val swString = new SmithWaterman[String]("-",5,-3,-2,-2,matchBaseString)

	def matchBaseString(i: String, j: String): Boolean = {
      j match {
        case "cT" => "T" == i || "C" == i
        case "cA" => "A" == i || "G" == i
        case nucleobase => (nucleobase == i)
      }
    }

    def matchBase(i: Nucleobase, j: Nucleobase): Boolean = {
      j match {
        case Nucleobase("cT") => (Nucleobase("T") == i || Nucleobase("C") == i)
        case Nucleobase("cA") => (Nucleobase("A") == i || Nucleobase("G") == i)
        case nucleobase => (nucleobase == i)
      }
    }

    def alignmentShiftString(sampleName: String, referenceName: String, bisulfiteConversion: Option[String], direction: Option[String], alignment: Alignment[String]): DNAAlignment = {
    	val sample = DNA(sampleName, alignment.seq1.map(Nucleobase(_)), bisulfiteConversion, direction)
    	val reference = DNA(referenceName, alignment.seq2.map(Nucleobase(_)), bisulfiteConversion, direction)
    	DNAAlignment(sample, reference, alignment.params)
    }

    def alignmentShift(sampleName: String, referenceName: String, bisulfiteConversion: Option[String], direction: Option[String], alignment: Alignment[Nucleobase]): DNAAlignment = {
    	val sample = DNA(sampleName, alignment.seq1, bisulfiteConversion, direction)
    	val reference = DNA(sampleName, alignment.seq2, bisulfiteConversion, direction)
    	DNAAlignment(sample, reference, alignment.params)
    }

   	def run(sample: DNA, reference: DNA): DNAAlignment = {
   		val sampleString = sample.nucleobases.map(_.toString)
   		val referenceString = reference.nucleobases.map(_.toString)
   		val nucleobaseAlignment = swString.run(sampleString, referenceString)
   		alignmentShiftString(sample.name, reference.name, reference.bisulfiteConversion, reference.direction, nucleobaseAlignment)
   	}

   	def runN(sample: DNA, reference: DNA): DNAAlignment = {
   		val nucleobaseAlignment = sw.run(sample.nucleobases, reference.nucleobases)
   		alignmentShift(sample.name, reference.name, reference.bisulfiteConversion, reference.direction, nucleobaseAlignment)
   	}
}

class SmithWaterman[T](gapVal: T, hit: Int, miss: Int, gapUp: Int, gapLeft: Int, comparison: (T,T) => Boolean){

	def run(seq1: IndexedSeq[T], seq2: IndexedSeq[T]): Alignment[T] = {
		val (directionMatrix,score) = calculate(seq1, seq2)
		val ((seq1Aligned,seq2Aligned),(start,end,refStart),(mismatches,gaps)) = readDirectionMatrix(seq1, seq2, directionMatrix, score)
		val params = AlignmentParameters(score.score,start,end,refStart,mismatches,gaps)
		Alignment[T](seq1Aligned,seq2Aligned,params)
	}

	def calculate(seq1: IndexedSeq[T], seq2: IndexedSeq[T]): (Array[Array[String]],ScoredPoint) = {
		var scoringMatrix = createScoringMatrix(seq1.length, seq2.length)
		var directionMatrix = Array.ofDim[String](seq1.length + 1, seq2.length + 1)
		var largest = ScoredPoint(0,Point(0,0))

		for (i <- 0 to seq1.length; j <- 0 to seq2.length){
			scoringMatrix(i)(j) = score(i,j)
		}

		def score(i:Int, j:Int): Int = {
			scoringMatrix(i)(j) match {
				case -1 => {
					val diag = score(i-1,j-1) + matchScore(seq1(i-1),seq2(j-1))
					val left = score(i,j-1) + gapLeft
					val up = score(i-1,j) + gapUp
					val best = List(diag, left, up, 0).max
					if (best > largest.score) largest = ScoredPoint(best,Point(i,j))
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
		(directionMatrix, largest)
	}

	def createScoringMatrix(length1: Int, length2: Int): Array[Array[Int]] = {
		var scoringMatrix = Array.ofDim[Int](length1 + 1 ,length2 + 1)
		for (i <- 1 to length1; j <- 1 to length2){
			scoringMatrix(i)(j) = -1
		}
		scoringMatrix
	}

	def readDirectionMatrix(seq1: IndexedSeq[T], seq2: IndexedSeq[T], directionMatrix: Array[Array[String]], largest: ScoredPoint): ((IndexedSeq[T],IndexedSeq[T]),(Int,Int,Int),(Int,Int)) = {
		var start = 0
		var refStart = 0
		var end = largest.point.x
		var mismatches = 0
		var gaps = 0
		def read(i: Int,j: Int): (IndexedSeq[T],IndexedSeq[T]) = {
			directionMatrix(i)(j) match {
				case "diag" => {
					mismatches += matchError(seq1(i-1),seq2(j-1))
					val prev = read(i-1,j-1)
					(seq1(i-1) +: prev._1, seq2(j-1) +: prev._2)					
				}
				case "upup" => {
					gaps += 1
					mismatches += 1
					val prev = read(i-1,j)
					(seq1(i-1) +: prev._1, gapVal +: prev._2)
					//(prev._1,prev._2)					
				}
				case "left" => {
					gaps += 1
					mismatches += 1
					val prev = read(i,j-1)
					//(gapVal +: prev._1, seq2(j-1) +: prev._2)		
					(prev._1,prev._2)					
			
				}
				case _ => {
					start = i+1
					refStart = j+1
					(IndexedSeq(),IndexedSeq())
				}
			}
		}
		(read(largest.point.x,largest.point.y),(start,end,refStart),(mismatches,gaps))
	}

	def matchScore(i: T, j: T): Int = {
		if (comparison(i,j)) hit else miss
	}

	def matchError(i: T, j: T): Int = {
		if (comparison(i,j)) 0 else 1
	}

}