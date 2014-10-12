/*
package org.bp.alignment

import org.bp.models._

case class Largest(x: Int, y: Int, score: Int)

case class AlignmentObject(score: Int, scoringMatrix: Vector[Vector[Option[Int]]], directionMatrix: Vector[Vector[Option[String]]], largest: Largest, reference: DNA)

case class ReadObject(pairs: List[(Nucleobase,Nucleobase)], start: Int, end: Int)

object SW {
	def setMatrix[T](x: Int, y: Int, newValue: T, matrix: Vector[Vector[Option[T]]]): Vector[Vector[Option[T]]] = {
		val inner = matrix(x)
		val newInner = inner updated (y, Some(newValue))
		matrix updated (x, newInner)
	}

	def initDirectionMatrix(x: Int, y: Int) = {
		Vector.fill(x + 1,y + 1)(None)
	}

	def initScoringMatrix(x: Int, y: Int) = {
		val score = Vector.fill(x + 1,y + 1)(None)
		val zeros = Vector.fill(y + 1)(Some(0))
		val padded = score updated (0, zeros) map (row => row updated (0, Some(0)))
		padded
	}

	def matchForward(sequence: DNA, reference: DNA, hit: Int, miss: Int)(i: Int, j: Int): Int = {
		//need to rewrite
		val seq = sequence.nucleobases
		val ref = reference.nucleobases
		ref(j-1) match {
			case Nucleobase("T") => if (Nucleobase("T") == seq(i-1) || Nucleobase("C") == seq(i-1)) hit else miss
			case nucleobase => if (nucleobase == seq(i-1)) hit else miss
		}
	}

	def run(sequence: DNA, references: List[DNA], hit: Int, miss: Int, gapUp: Int, gapLeft: Int) = {
		val go = runOnce(sequence = sequence, hit = hit, miss = miss, gapUp = gapUp, gapLeft = gapLeft)_
		//val matches = List[matchForward(sequence = sequence, reference = reference, hit = hit, miss = miss)]
		//val scoring = matches.map(score(sequence = sequence, reference = reference, matchScore = _))
		references.map(go(_))
	}

	def runOnce(sequence: DNA, hit: Int, miss: Int, gapUp: Int, gapLeft: Int)(reference: DNA) = {
		val seq = sequence.nucleobases
		val ref = reference.nucleobases
		val directionMatrix = initDirectionMatrix(seq.length,ref.length)
		val scoringMatrix = initScoringMatrix(seq.length,ref.length)
		getMatrix(sequence = sequence, reference = reference, scoringMatrix = scoringMatrix, directionMatrix = directionMatrix, matchScore = matchForward(sequence,reference,hit,miss)_, gapUp = gapUp, gapLeft = gapLeft)
	}

	def score(sequence: DNA, reference: DNA, matchScore: (Int,Int) => Int, gapUp: Int, gapLeft: Int)(i: Int, j: Int)(alignment: AlignmentObject): AlignmentObject = {
		alignment.scoringMatrix(i)(j) match {
			//case -1 => List(score(i-1,j-1) + matchScore(i,j), score(i-1,j) + gap, score(i,j-1) + gap, 0).max
			case None => {
				val setScore = score(sequence,reference,matchScore,gapUp,gapLeft)_
				val diag = setScore(i-1, j-1)(alignment)
				val diagScore = diag.score + matchScore(i,j)
				val left = setScore(i, j-1)(diag)
				val leftScore = left.score + gapLeft
				val up = setScore(i-1, j)(left)
				val upScore = up.score + gapUp
				val bestList = List(diagScore, leftScore, upScore, 0)
				val best = bestList.max
				val newLargest = if (best > alignment.largest.score) Largest(i,j,best) else alignment.largest
				val update = best match {
					case d if d == diagScore => {
						(d,"diag")
					}
					case l if l == leftScore => {
						(l,"left")
					}
					case u if u == upScore => {
						(u,"upup")
					}
					case 0 => {
						(0,"stop")
					}
					case a => {
						println(a + " " + bestList)
						(-100000,"error")
					}
				}
				AlignmentObject(score = update._1, scoringMatrix = setMatrix(i,j,update._1,up.scoringMatrix), directionMatrix = setMatrix(i,j,update._2,up.directionMatrix), largest = newLargest, reference = reference)
			}
			case Some(score) => alignment
		}
	}

	def getMatrix(sequence: DNA, reference: DNA, scoringMatrix: Vector[Vector[Option[Int]]], directionMatrix: Vector[Vector[Option[String]]],matchScore: (Int,Int) => Int,gapUp:Int,gapLeft:Int): AlignmentObject = {
		val setScore = score(sequence,reference,matchScore,gapUp,gapLeft)_
		val xRange = 0 to scoringMatrix.length-1
		val yRange = 0 to scoringMatrix(0).length-1
		val initAlignment = AlignmentObject(score = 0, scoringMatrix = scoringMatrix, directionMatrix = directionMatrix, largest = Largest(0,0,0), reference = reference)
		xRange.foldLeft(initAlignment)((alignment,x) => yRange.foldLeft(alignment)((alignment,y) => setScore(x,y)(alignment)))
	}

	def readAlignment(sequence: DNA, reference: DNA, alignment: AlignmentObject) = {
		val ref = reference.nucleobases
		val seq = sequence.nucleobases
		val go = read(seq,ref)_
		//println(alignment)
		val result = go(alignment.largest.x,alignment.largest.y,alignment)
		(result.pairs,(result.start,result.end))
	}

	def read(seq: Vector[Nucleobase],ref: Vector[Nucleobase])(i: Int, j: Int, alignment: AlignmentObject): ReadObject = {
		alignment.directionMatrix(i)(j) match {
			case Some("diag") => ReadObject((seq(i-1),ref(j-1)) :: read(seq,ref)(i-1,j-1,alignment).pairs,0,0)
			case Some("upup") => ReadObject((seq(i-1),Nucleobase("-")) :: read(seq,ref)(i-1,j,alignment).pairs,0,0)
			case Some("left") => ReadObject((Nucleobase("-"),ref(j-1)) :: read(seq,ref)(i,j-1,alignment).pairs,0,0)
			case _ => {
				ReadObject(List((Nucleobase(""),Nucleobase(""))),i+1,alignment.largest.x)
			}
		}
	}

	def readErrors(reference: DNA, read: List[(Nucleobase,Nucleobase)]) = {
		val alignment = read.foldLeft((0,0))((tally,next) => {
			val sequence = next._1
			val reference = next._2
			sequence match {
				case Nucleobase("-") => (tally._1 + 1, tally._2 + 1)
				case Nucleobase("T") => if (Nucleobase("T") == reference || Nucleobase("C") == reference) tally else (tally._1 + 1, tally._2)
				case nucleobase => if (reference != nucleobase) (tally._1 + 1, tally._2) else tally
			}
		})
		Alignment(alignment._1,alignment._2,reference.direction.getOrElse(""),read.length)
	}

	def getSequences(pairs: List[(Nucleobase,Nucleobase)]) = {
		val both = pairs.unzip
		//scoringMatrix.foreach{case i => i foreach {j => print(j + " ")}; print('\n')}
		(DNA("bisulfite",both._1.reverse.toVector),DNA("genome",both._2.reverse.toVector))
	}
}
*/
