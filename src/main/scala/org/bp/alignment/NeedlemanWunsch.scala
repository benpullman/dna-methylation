package org.bp.alignment

import org.bp.models._

object NeedlemanWunsch {
	def generateScoringMatrix(x: DNA, y: DNA, hit: Int, miss: Int, gap: Int) = {
		val a = x.nucleobases
		val b = y.nucleobases
		var directionMatrix = Array.ofDim[String](a.length + 1 ,b.length + 1)
		var scoringMatrix = Array.ofDim[Int](a.length + 1 ,b.length + 1)
		def initializeScoringMatrix: Unit = {
			for (i <- 1 to a.length; j <- 1 to b.length){
				scoringMatrix(i)(j) = -10000000
			}
			for (i <- 0 to a.length){
				scoringMatrix(i)(0) = gap * i
			}
			for (j <- 0 to b.length){
				scoringMatrix(0)(j) = gap * j
			}
		}
		def matchScore(i: Int, j: Int) = {
			if (a(i-1) == b(j-1)) hit else miss
		}
		def score(i:Int,j:Int): Int = {
			scoringMatrix(i)(j) match {
				//case -1 => List(score(i-1,j-1) + matchScore(i,j), score(i-1,j) + gap, score(i,j-1) + gap, 0).max
				case -10000000 => {
					val diag = score(i-1,j-1) + matchScore(i,j)
					val left = score(i,j-1) + gap
					val up = score(i-1,j) + gap
					List(diag, left, up).max match {
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
				case _ => List((Nucleobase(""),Nucleobase("")))
			}
		}
		initializeScoringMatrix
		setScore
		//directionMatrix.foreach{case i => i foreach {j => print(j + " ")}; print('\n')}
		val both = read(a.length,b.length).reverse.drop(1).unzip
		(DNA("compare",both._1),DNA("reference",both._2))
		//scoringMatrix
	}
}