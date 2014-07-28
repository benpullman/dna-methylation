package org.bp.alignment

import org.bp.models._

object SmithWaterman {
	def generateScoringMatrix(x: DNA, y: DNA, hit: Int, miss: Int, gapUp: Int, gapLeft: Int) = {
		val a = x.nucleobases
		val b = y.nucleobases
		var largest = (0,(0,0))
		var start = 0
		var end = 0
		var directionMatrix = Array.ofDim[String](a.length + 1 ,b.length + 1)
		var scoringMatrix = Array.ofDim[Int](a.length + 1 ,b.length + 1)
		def initializeScoringMatrix: Unit = {
			for (i <- 1 to a.length; j <- 1 to b.length){
				scoringMatrix(i)(j) = -1
			}
		}
		def matchScore(i: Int, j: Int) = {
			b(j-1) match {
				case Nucleobase("C") => if (Nucleobase("T") == a(i-1) || Nucleobase("C") == a(i-1)) hit else miss
				case nucleobase => if (nucleobase == a(i-1)) hit else miss
			}
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
			for (i <- 0 to a.length; j <- 0 to b.length){
				scoringMatrix(i)(j) = score(i,j)
			}
		}
		def read(i: Int,j: Int): List[(Nucleobase, Nucleobase)] = {
			directionMatrix(i)(j) match {
				case "diag" => (a(i-1),b(j-1)) :: read(i-1,j-1)
				case "upup" => (a(i-1),Nucleobase("-")) :: read(i-1,j)
				case "left" => (Nucleobase("-"),b(j-1)) :: read(i,j-1)
				case _ => {
					start = i+1
					List((Nucleobase(""),Nucleobase("")))
				}
			}
		}
		initializeScoringMatrix
		setScore
		//directionMatrix.foreach{case i => i foreach {j => print(j + " ")}; print('\n')}
		end = largest._2._1
		(read(largest._2._1,largest._2._2).reverse.drop(1),(start,end))
	}

	def getErrors(pairs: List[(Nucleobase,Nucleobase)]) = {
		val alignment = pairs.foldLeft((0,0))((tally,next) => {
			val genome = next._2
			val bisulfite = next._1
			bisulfite match {
				case Nucleobase("-") => (tally._1 + 1, tally._2 + 1)
				case Nucleobase("T") => if (Nucleobase("T") == genome || Nucleobase("C") == genome) tally else (tally._1 + 1, tally._2)
				case nucleobase => if (genome != nucleobase) (tally._1 + 1, tally._2) else tally
			}
		})
		Alignment(alignment._1,alignment._2,"forward",pairs.length)
	}

	def getSequences(pairs: List[(Nucleobase,Nucleobase)]) = {
		val both = pairs.unzip
		//scoringMatrix.foreach{case i => i foreach {j => print(j + " ")}; print('\n')}
		(DNA("bisulfite",both._1),DNA("genome",both._2))
	}
}