package org.bp.alignment

import org.bp.models._

case class Alignment(sample: DNA, reference: DNA, params: AlignmentParameters){
	def length: this.params.length
	def percentAligned: this.params.percentAligned
}

case class AlignmentParameters(score: Int, start: Int, end: Int, gaps: Int, mismatches: Int){
	def length: this.start-this.end
	def percentAligned: this.mismatches.toDouble/this.length.toDouble
}

object SmithWaterman{

	def process(sample: DNA, reference: DNA, hit: Int, miss: Int, gapUp: Int, gapLeft: Int) = {

		//begin timing
		val startClock: Double = System.currentTimeMillis / 1000.0
		val finishClock: Double = System.currentTimeMillis / 1000.0
		println("Alignment time + " + (finishClock.toDouble - startClock.toDouble))

		//intialize variables
		val s = sample.nucleobases
		val r = references.nucleobases
		val direction = references.direction
		val bisulfiteConversion = references.bisulfiteConversion
		case class ScoredPoint(score: Int, point: Point)
		case class Point(x: Int, y: Int)
		var largest = ScoredPoint(0,(0,0))
		var directionMatrix = Array.ofDim[String](s.length + 1 ,r.length + 1)
		var scoringMatrix = Array.ofDim[Int](s.length + 1 ,r.length + 1)

		//Operations
		setScoringMatrix
		val aligned = readDirectionMatrix
		val params = AlignmentParameters(largest.score,aligned._2._1,aligned._2._2,aligned._3._1,aligned._3._2)
		val alignedSample = DNA("sample",aligned._1._1,direction,bisulfiteConversion)
		val alignedReference = DNA("reference",aligned._1._2,direction,bisulfiteConversion)
		Alignment(alignedSample,alignedReference,params)

		//return
		//AlignmentProcess(scoringMatrix, directionMatrix, y, largest._2._1,largest._2._2, largest._1)

		def initializeScoringMatrix: Unit = {
			for (i <- 1 to a.length; j <- 1 to b.length){
				scoringMatrix(i)(j) = -1
			}
		}

		def match(i: Int, j: Int): Boolean = {
			r(j-1) match {
				case Nucleobase("cT") => if (Nucleobase("T") == s(i-1) || Nucleobase("C") == s(i-1))
				case Nucleobase("cA") => if (Nucleobase("A") == s(i-1) || Nucleobase("G") == s(i-1))
				case nucleobase => if (nucleobase == s(i-1))
			}
		}

		def matchScore(i: Int, j: Int): Int = {
			if match(i,j) hit else miss
		}

		def matchError(i: Int, j: Int): Int = {
			if match(i,j) 0 else 1
		}

		def setScoringMatrix: Unit = {

			initializeScoringMatrix

			for (i <- 0 to a.length; j <- 0 to b.length){
				scoringMatrix(i)(j) = score(i,j)
			}

			def initializeScoringMatrix: Unit = {
				for (i <- 1 to a.length; j <- 1 to b.length){
					scoringMatrix(i)(j) = -1
				}
			}

			def score(i:Int,j:Int): Int = {
				scoringMatrix(i)(j) match {
					case -1 => {
						val diag = score(i-1,j-1) + matchScore(i,j)
						val left = score(i,j-1) + gapLeft
						val up = score(i-1,j) + gapUp
						val best = List(diag, left, up, 0).max
						if (best > largest._1) largest = ScoredPoint(best,(i,j))
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
		}

		def readDirectionMatrix: ((Vector,Vector),(Int,Int),(Int,Int)) = {
			var start = 0
			var end = largest.point.x
			var mismatches = 0
			var gaps = 0
			def read(i: Int,j: Int): (Vector,Vector) = {
				directionMatrix(i)(j) match {
					case "diag" => {
						mismatches += matchError(i,j)
						val prev = read(i-1,j-1,directionMatrix)
						((a(i-1) +: prev._1, b(j-1) +: prev._2)					
					}
					case "upup" => {
						gaps += 1
						mismatches += 1
						val prev = read(i-1,j,directionMatrix)
						((a(i-1) +: prev._1, Nucleobase("-") +: prev._2)					
					}
					case "left" => {
						gaps += 1
						mismatches += 1
						val prev = read(i,j-1,directionMatrix)
						(Nucleobase("-") +: prev._1, b(j-1) +: prev._2)					
					}
					case _ => {
						start = i+1
						(Vector(),Vector())
					}
				}
			}
			(read(largest.point.x,largest.point.y),(start,end),(mismatches,gaps))
		}

	}

	

	def getErrors(pairs: List[(Nucleobase,Nucleobase)],reference: DNA) = {
		val alignment = pairs.foldLeft((0,0))((tally,next) => {
			val genome = next._2
			val bisulfite = next._1
			bisulfite match {
				case Nucleobase("-") => (tally._1 + 1, tally._2 + 1)
				case Nucleobase("T") => if (Nucleobase("T") == genome || Nucleobase("C") == genome) tally else (tally._1 + 1, tally._2)
				case nucleobase => if (genome != nucleobase) (tally._1 + 1, tally._2) else tally
			}
		})
		Alignment(alignment._1,alignment._2,reference.direction.getOrElse(""),pairs.length)
	}

	def getSequences(pairs: List[(Nucleobase,Nucleobase)], reference: DNA) = {
		val both = pairs.unzip
		//scoringMatrix.foreach{case i => i foreach {j => print(j + " ")}; print('\n')}
		(DNA("bisulfite",both._1.toVector,reference.bisulfiteConversion,reference.direction),DNA(reference.name,both._2.toVector,reference.bisulfiteConversion,reference.direction))
	}
}