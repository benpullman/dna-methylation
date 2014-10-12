package org.bp.models

case class AlignmentAnalysis(mismatches: Int, gaps: Int, direction: String, length: Int)
case class Methylation(cpGSites: Int, methylatedCpGSites: Int, sequence: Vector[String],reference: Vector[Int])
case class BisulfiteConversion(convertedCpH: Int, unconvertedCpH: Int, direction: String)
case class Analysis(
	sequenceName: String,
	sequenceLength: Int,
	referenceLength: Int,
	referenceName: String,
	alignment: AlignmentAnalysis,
	bisulfite: BisulfiteConversion,
	methylation: Methylation,
	seqStart: Int,
	seqEnd: Int,
	barcode: String
	){
	def toString2(): String = {
		"Sequence name:\t\t\t" + sequenceName + "\n" +
		"Sequence length:\t\t" + sequenceLength + "\n" +
		"Number of mismatches:\t\t" + alignment.mismatches + "\n" +
		"Number of gaps:\t\t\t" + alignment.gaps + "\n" +
		"Alignment identity:\t\t" + ((sequenceLength-alignment.mismatches.toFloat)/sequenceLength)*100 + "%\n" + 
		"Number of CpGs:\t\t\t" + (methylation.cpGSites + methylation.methylatedCpGSites) + "\n" + 
		"Number of methylated CpGs:\t" + methylation.methylatedCpGSites + "\n" + 
		"Percent methylated:\t\t" + (methylation.methylatedCpGSites.toFloat/(methylation.cpGSites + methylation.methylatedCpGSites)) + "%\n" + 
		"Number of unconverted CpHs (CpA/CpT/CpC):\t" + bisulfite.unconvertedCpH + "\n" +
		"Number of CpHs:\t" + bisulfite.convertedCpH + "\n" +
		"Percent converted CpHs:\t\t" + (bisulfite.convertedCpH.toFloat/(bisulfite.convertedCpH + bisulfite.unconvertedCpH))+ "%\n" +
		"Methylation pattern:\t\t" + methylation.sequence.foldLeft("")(_+_) + "\n"
	}
	override def toString(): String = {
		"--------------------------\n" +
		"Summary of information\n" +
		"--------------------------\n" +
		"Bisulfite sequence name\t\t\t\t" + 					sequenceName + "\n" +
		"Conversion\t\t\t\t\t" + 								bisulfite.direction +  " conversion\n" +
		"(conversion of forward strand of the genomic sequence)\n" +
		"Length of bisulfite sequence\t\t\t" + 					sequenceLength + "\n" +
		"Length of target genome sequence\t\t"+					referenceLength + "\n" +
		"Aligned region of bisulfite sequence\t\t" +			seqStart + " - " + seqEnd +"\n" +
		"Alignment direction\t\t\t\t" +							alignment.direction + "\n" +
		"Number of CpGs:\t\t\t\t\t" + 							(methylation.cpGSites + methylation.methylatedCpGSites) + "\n" + 
		"Number of methylated CpGs\t\t\t" + 					methylation.methylatedCpGSites + "\n" + 
		"Percent methylated\t\t\t\t" + 							(methylation.methylatedCpGSites.toFloat/(methylation.cpGSites + methylation.methylatedCpGSites))*100 + "%\n" + 
		"Number of unconverted CpHs (CpA/CpT/CpC)\t" + 			bisulfite.unconvertedCpH + "\n" +
		"Number of CpHs\t\t\t\t\t" + 							(bisulfite.unconvertedCpH + bisulfite.convertedCpH) + "\n" +
		"Percent converted CpHs (CpA/CpT/CpC)\t\t" + 			(bisulfite.convertedCpH.toFloat/(bisulfite.convertedCpH + bisulfite.unconvertedCpH))*100+ "%\n" +
		"Number of mismatches (include gaps)\t\t" + 			alignment.mismatches + "\n" +
		"Number of gaps:\t\t\t\t\t" + 							alignment.gaps + "\n" +
		"Alignment length:\t\t\t\t" + 							alignment.length + "\n" +
		"Percent identity:\t\t\t\t" + 							((sequenceLength-alignment.mismatches.toFloat)/sequenceLength)*100 + "%\n" + 
		"Methylation pattern:\t\t\t\t" + 						methylation.sequence.foldLeft("")(_+_) + "\n" +
		"(U: unmethylated, M: methylated, A,C,G,T,N: mismatch, -: gap)"
	}
}