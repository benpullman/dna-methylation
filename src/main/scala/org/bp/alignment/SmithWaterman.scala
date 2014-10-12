package org.bp.alignment

import org.bp.models._

case class Alignment[T](seq1: Seq[T], seq2: Seq[T], params: AlignmentParameters){
	def length: Int = this.params.length
	def percentAligned: Double = this.params.percentAligned
}

case class AlignmentParameters(score: Int, start: Int, end: Int, gaps: Int, mismatches: Int){
	def length: Int = this.start - this.end
	def percentAligned: Double = this.mismatches.toDouble/this.length.toDouble
}

case class ScoredPoint(score: Int, point: Point)

case class Point(x: Int, y: Int)

class SmithWaterman[T](gapVal: T, hit: Int, miss: Int, gapUp: Int, gapLeft: Int, comparison: (T,T) => Boolean){

	def run(seq1: Seq[T], seq2: Seq[T]): Alignment {
		val directionScore = calculate(seq1, seq2)
		val aligned = readDirectionMatrix(seq1, seq2, directionScore._1, directionScore._2)
		val params = AlignmentParameters(directionScore._2.score,aligned._2._1,aligned._2._2,aligned._3._1,aligned._3._2)
		Alignment(aligned._1._1,aligned._1._2,params)
	}

	def calculate(seq1: Seq[T], seq2: Seq[T]): (Array[Array[String]],ScoredPoint) {
		var scoringMatrix = createScoringMatrix(seq1.length, seq2.length)
		var directionMatrix = Array.ofDim[String](seq1.length + 1, seq2.length + 1)
		var largest = ScoredPoint(0,Point(0,0))

		for (i <- 0 to seq1.length; j <- 0 to seq2.length){
			scoringMatrix(i)(j) = score(i,j)
		}

		def score(i:Int,j:Int): Int = {
			scoringMatrix(i)(j) match {
				case -1 => {
					val diag = score(i-1,j-1) + matchScore(seq1(i),seq2(j))
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

	def readDirectionMatrix(seq1: Seq[T], seq2: Seq[T], directionMatrix: Array[Array[String]], largest: ScoredPoint): ((Vector[T],Vector[T]),(Int,Int),(Int,Int)) = {
		var start = 0
		var end = largest.point.x
		var mismatches = 0
		var gaps = 0
		def read(i: Int,j: Int): (Vector[Nucleobase],Vector[Nucleobase]) = {
			directionMatrix(i)(j) match {
				case "diag" => {
					mismatches += matchError(i,j)
					val prev = read(i-1,j-1)
					(s(i-1) +: prev._1, r(j-1) +: prev._2)					
				}
				case "upup" => {
					gaps += 1
					mismatches += 1
					val prev = read(i-1,j)
					(s(i-1) +: prev._1, Nucleobase("-") +: prev._2)					
				}
				case "left" => {
					gaps += 1
					mismatches += 1
					val prev = read(i,j-1)
					(Nucleobase("-") +: prev._1, r(j-1) +: prev._2)					
				}
				case _ => {
					start = i+1
					(Vector(),Vector())
				}
			}
		}
		(read(largest.point.x,largest.point.y),(start,end),(mismatches,gaps))
	}

	def matchScore(i: T, j: T): Int = {
		if (comparison(i,j)) hit else miss
	}

	def matchError(i: T, j: T): Int = {
		if (comparison(i,j)) 0 else 1
	}

}