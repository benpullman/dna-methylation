package org.bp.methylation

import org.bp.fasta._
import org.bp.models._

object MethylationProcess {
	def checkSites(CpgSites: List[Int],DNA: DNA) = {
		CpgSites.map(site => {
			DNA.nucleobases(site) match {
				case Nucleobase("C") => "M"
				case Nucleobase("T") => "U"
				case base: Nucleobase => base.toString
			}
		})
	}
}