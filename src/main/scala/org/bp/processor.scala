package org.bp

import org.bp.models._
import org.bp.fasta._
import org.bp.methylation._
import org.bp.alignment._
import spray.json._
import DefaultJsonProtocol._

object Processor {

  def main(args: Array[String]): Unit = {
    import JsonProtocol._
    //val b = DNA("second",List(Nucleobase("A"),Nucleobase("C"),Nucleobase("T"),Nucleobase("G"),Nucleobase("A"),Nucleobase("T"),Nucleobase("T"),Nucleobase("C"),Nucleobase("A")))
    //val a = DNA("first",List(Nucleobase("A"),Nucleobase("C"),Nucleobase("G"),Nucleobase("C"),Nucleobase("A"),Nucleobase("T"),Nucleobase("C"),Nucleobase("A")))
    //val matrix = SmithWaterman.generateScoringMatrix(a,b,2,-3,-2)
    //matrix.foreach{case a => a foreach {b => print(b.toString + "\t")}; print('\n')}
    val reference = FASTA.read("data/quma-samples/sample_genome_fasta.txt").head
    val toCheck = FASTA.read("data/quma-samples/Gm9_J1_seq_fasta.txt")
    val results = toCheck.map(sequence => {
      /*val check = SmithWaterman.getSequences(SmithWaterman.generateScoringMatrix(sequence,reference,5,-3,-20,-2))
      val bisulfite = check._1
      val genome = check._2
      println(sequence.name)
      println("Total length: " + genome.nucleobases.length + " Compared length: " + bisulfite.nucleobases.length)
      val comparison = Sequence.generateComparison(bisulfite,genome)
      Sequence.compare(bisulfite,genome,comparison,60)
      Sequence.methylation(comparison)*/
      Sequence.toAnalysis(sequence,reference,5,-3,-20,-2)
      //println(bisulfite)
      //println(genome)
      })
    print(results.toJson.prettyPrint)
    //println("Length of sequence to check " + toCheck.length)
    //toCheck.map(dna => println(dna.name + "\n" + Methylation.checkSites(reference)(dna)))
    //println("done")
    //println(FASTA.readStream("data/mssm-samples/-home-sbsuser-SMRT-smrtanalysis-current-common-jobs-018-018823-data-reads_of_insert.fasta.txt"))
  }

}