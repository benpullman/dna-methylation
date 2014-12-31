package org.bp.methylation

import org.bp.fasta._
import org.bp.models._

object MethylationProcess {
	def checkSites(cpg: IndexedSeq[Int], dna: DNA, start: Int) = {
		val reversedDNA = dna.nucleobases.reverse
		//val start = 0
		//println(dna)
		val end = start + dna.nucleobases.length
		//println(start + " " + end)
		val (pre,rest) = cpg.span(_ <= start)
		val (read,post) = rest.span(_ < end)
		//println(pre.length + " " + read.length + " " + post.length)
		//println(read)
		val bisulfiteConversion = dna.bisulfiteConversion.getOrElse("")
		val direction = dna.direction.getOrElse("")
		val case1 = ((bisulfiteConversion == "C => T") && (direction == "Forward")) || ((bisulfiteConversion == "G => A") && (direction == "Reverse Complement"))
        val case2 = ((bisulfiteConversion == "G => A") && (direction == "Forward")) || ((bisulfiteConversion =="C => T") && (direction == "Reverse Complement"))
		pre.map(m => "X") ++ read.map(site => {
			if (site-start+1<=reversedDNA.length-1){
				if (case1) {
					//println("case1")
					//println(site + " " + reversedDNA(site-start+1))
					reversedDNA(site-start+1) match {
						case Nucleobase("C") => "M"
						case Nucleobase("T") => "U"
						case base: Nucleobase => base.displayString
					}
				}
				else if (case2) {
					reversedDNA(site+1-start) match {
						case Nucleobase("G") => "M"
						case Nucleobase("A") => "U"
						case base: Nucleobase => base.displayString
					}
				} else {
					"error"
				}
			}else{
				"error"
			}
		}) ++ post.map(m => "X")
	}
	/*
	def checkSites(CpgSites: IndexedSeq[Int],dna: dna, start: Int) = {
		println(dna.nucleobases)
		val bisulfiteConversion = dna.bisulfiteConversion.getOrElse("")
		val direction = dna.direction.getOrElse("")
		val case1 = ((bisulfiteConversion == "C => T") && (direction == "Forward")) || ((bisulfiteConversion == "G => A") && (direction == "Reverse Complement"))
        val case2 = ((bisulfiteConversion == "G => A") && (direction == "Forward")) || ((bisulfiteConversion =="C => T") && (direction == "Reverse Complement"))
		CpgSites.map(site => {
			if (case1) {
				dna.nucleobases(site) match {
					case Nucleobase("C") => "M"
					case Nucleobase("T") => "U"
					case base: Nucleobase => base.toString
				}
			}
			else if (case2) {
				dna.nucleobases(site+1) match {
					case Nucleobase("G") => "M"
					case Nucleobase("A") => "U"
					case base: Nucleobase => base.toString
				}
			} else {
				"error"
			}
		})
	}
	*/
}