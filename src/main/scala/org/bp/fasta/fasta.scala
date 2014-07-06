package org.bp.fasta

import scala.io.Source
import org.bp.models._

object FASTA {
	def read(filename: String) = {
		//reading by char
		// two char sequences - 
		val output = Source.fromFile(filename).foldLeft[(Boolean,String,List[Nucleobase],List[DNA])]((true, "", List(), List()))((out, add) => {
			add match {
				case '>' => {
					println(DNA(out._2, out._3))
					(true, "", List(),  DNA(out._2, out._3.reverse) :: out._4)
				}
				case '\n' => if (out._1) (false, out._2, out._3, out._4) else out
				case '\r' => if (out._1) (false, out._2, out._3, out._4) else out
				case other => if (out._1) (out._1, out._2 + add, out._3, out._4) else (out._1, out._2, Nucleobase.charToNucleobase(add) :: out._3, out._4)
				}
			}
		)
		"Done"
		//DNA(output._2,output._3) :: output._4.drop(1)
	}
	def readStream(filename: String) = {
		//reading by char
		// two char sequences - 
		val output = Source.fromFile(filename).foldLeft[(Boolean,String,List[Nucleobase],Stream[DNA],Int)]((true, "", List(), Stream(),0))((out, add) => {
			add match {
				case '>' => {
					println(out._5 + 1)
					(true, "", List(),  DNA(out._2, out._3.reverse) #:: out._4,out._5 + 1)
				}
				case '\n' => if (out._1) (false, out._2, out._3, out._4, out._5) else out
				case '\r' => if (out._1) (false, out._2, out._3, out._4, out._5) else out
				case other => if (out._1) (out._1, out._2 + add, out._3, out._4, out._5) else (out._1, out._2, Nucleobase.charToNucleobase(add) :: out._3, out._4, out._5)
				}
			}
		)
		val outstream = DNA(output._2,output._3) #:: output._4.drop(1)
		outstream(10000) + "\n" + outstream(10001)
	}
	def readNoSave(filename: String) = {
		//reading by char
		// two char sequences - 
		Source.fromFile(filename).toStream.foldLeft[(Boolean,String,List[Nucleobase],Int)]((true, "", List(),0))((out, add) => {
			add match {
				case '>' => {
					println(out._4 + 1)
					//println(DNA(out._2, out._3))
					(true, "", List(), out._4 + 1)
				}
				case '\n' => if (out._1) (false, out._2, out._3, out._4) else out
				case '\r' => if (out._1) (false, out._2, out._3, out._4) else out
				case other => if (out._1) (out._1, out._2 + add, out._3, out._4) else (out._1, out._2, Nucleobase.charToNucleobase(add) :: out._3, out._4)
				}
			}
		)
		"Done"
		//DNA(output._2,output._3) :: output._4.drop(1)
	}
}