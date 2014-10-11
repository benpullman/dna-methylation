package org.bp.methylation

import org.bp.fasta._
import org.bp.models._

object MethylationProcess {
	def checkSites(CpgSites: Vector[Int],DNA: DNA) = {
		val bisulfiteConversion = DNA.bisulfiteConversion.getOrElse("")
		val direction = DNA.direction.getOrElse("")
		val case1 = ((bisulfiteConversion == "C => T") && (direction == "Forward")) || ((bisulfiteConversion == "G => A") && (direction == "Reverse Complement"))
        val case2 = ((bisulfiteConversion == "G => A") && (direction == "Forward")) || ((bisulfiteConversion =="C => T") && (direction == "Reverse Complement"))
		CpgSites.map(site => {
			if (case1) {
				DNA.nucleobases(site) match {
					case Nucleobase("C") => "M"
					case Nucleobase("T") => "U"
					case base: Nucleobase => base.toString
				}
			}
			else if (case2) {
				DNA.nucleobases(site+1) match {
					case Nucleobase("G") => "M"
					case Nucleobase("A") => "U"
					case base: Nucleobase => base.toString
				}
			} else {
				"error"
			}
		})
	}
}