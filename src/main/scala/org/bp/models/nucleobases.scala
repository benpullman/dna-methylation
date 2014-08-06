package org.bp.models

import org.bp.alignment._
import org.bp.methylation._

case class Nucleobase(name: String){
  override def toString = this.name
  def toRNA = this match {
    case Nucleobase("T") => Nucleobase("U")
    case other => other
  }
  def complement = this match {
    case Nucleobase("A") => Nucleobase("T")
    case Nucleobase("C") => Nucleobase("G")
    case Nucleobase("G") => Nucleobase("C")
    case Nucleobase("T") => Nucleobase("A")
  }
}

class Sequence(name: String, nucleobases: List[Nucleobase], bisulfiteConversion: Option[String] = None, direction: Option[String] = None){
    override def toString = "Name: " + name + bisulfiteConversion.map(" " + _ + " ").getOrElse("") + direction.map(_ + " ").getOrElse("") + "\nSequence: " + this.nucleobases.foldLeft("")((nucleotide,nucleobase) => nucleotide + nucleobase)
    def findCpgSites = this.nucleobases.foldLeft[(Boolean,Int,List[Int])]((false,0,List()))((acc,next) => {
      if (acc._1 && next == Nucleobase("G")) {
        (false, acc._2 + 1, (acc._2 - 1) :: acc._3)
      } else if (next == Nucleobase("C")){
        (true, acc._2 + 1, acc._3)
      } else {
        (false, acc._2 + 1, acc._3)
      }
    })._3.reverse
}

case class DNA(name: String, nucleobases: List[Nucleobase], bisulfiteConversion: Option[String] = None, direction: Option[String] = None) extends Sequence(name: String, nucleobases: List[Nucleobase],bisulfiteConversion: Option[String], direction: Option[String]){
  def toRNA = RNA(name = name, nucleobases = this.nucleobases.map(nucleobase => nucleobase.toRNA))
  def reverseComplement = this.nucleobases.reverse.map(nucleobase => nucleobase.complement)
  def reverse = this.copy(nucleobases = this.nucleobases.reverse)
  def add(nucleobase: Nucleobase) = this.copy(nucleobases = nucleobase::nucleobases)
  def generateBisulfiteQuartet = {
      val CtoT = Map(Nucleobase("C") -> Nucleobase("T"))
      val GtoA = Map(Nucleobase("G") -> Nucleobase("A"))
      def bisulfiteConversionForward(input: DNA, conversion: Map[Nucleobase,Nucleobase], conversionString: String) = {
        val sequence = input.nucleobases.foldLeft[(List[Nucleobase],Option[Nucleobase])]((List(),None))((bases,next) => {
          next match {
            case Nucleobase("C") => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(Nucleobase("C"), Nucleobase("C")) :: bases._1, Some(Nucleobase("C")))
              case _ => (bases._1, Some(Nucleobase("C")))
            }
            case Nucleobase("G") => bases._2 match {
              case Some(Nucleobase("C")) => (Nucleobase("G") :: Nucleobase("C") :: bases._1,None)
              case _ => (conversion.getOrElse(Nucleobase("G"), Nucleobase("G")) :: bases._1,None)
            }
            case base => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(base,base) :: conversion.getOrElse(Nucleobase("C"),Nucleobase("C")) :: bases._1,None)
              case _ => (conversion.getOrElse(base,base) :: bases._1,None)
            }
          }
        })._1.reverse
        DNA(input.name,sequence,Some(conversionString),Some("Forward"))
      }
      def bisulfiteConversionReverseComplement(input: DNA, conversion: Map[Nucleobase,Nucleobase], conversionString: String) = {
        val sequence = input.nucleobases.foldLeft[(List[Nucleobase],Option[Nucleobase])]((List(),None))((bases,next) => {
          next match {
            case Nucleobase("C") => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(Nucleobase("C"), Nucleobase("C")).complement :: bases._1, Some(Nucleobase("C")))
              case _ => (bases._1, Some(Nucleobase("C")))
            }
            case Nucleobase("G") => bases._2 match {
              case Some(Nucleobase("C")) => (Nucleobase("C") :: Nucleobase("G") :: bases._1,None)
              case _ => (conversion.getOrElse(Nucleobase("G"), Nucleobase("G")).complement :: bases._1,None)
            }
            case base => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(base,base).complement :: conversion.getOrElse(Nucleobase("C"),Nucleobase("C")).complement :: bases._1,None)
              case _ => (conversion.getOrElse(base,base).complement :: bases._1,None)
            }
          }
        })._1
        DNA(input.name,sequence,Some(conversionString),Some("Reverse Complement"))
      }
      List(
        bisulfiteConversionForward(this,CtoT,"C => T"),
        bisulfiteConversionForward(this,GtoA,"G => A"),
        bisulfiteConversionReverseComplement(this,CtoT,"C => T"),
        bisulfiteConversionReverseComplement(this,GtoA,"G => A")
      )
    }
}

case class RNA(name: String, nucleobases: List[Nucleobase], bisulfiteConversion: Option[String] = None, direction: Option[String] = None) extends Sequence(name: String, nucleobases: List[Nucleobase],bisulfiteConversion: Option[String], direction: Option[String])

object Nucleobase{
	def byteToNucleobase(byte: Byte): Nucleobase = {
    byte match {
        case 65 => Nucleobase("A")
        case 67 => Nucleobase("C")
        case 71 => Nucleobase("G")
        case 84 => Nucleobase("T")
      }
  	}
  def charToNucleobase(char: Char): Nucleobase = {
    char match {
        case 'A' => Nucleobase("A")
        case 'C' => Nucleobase("C")
        case 'G' => Nucleobase("G")
        case 'T' => Nucleobase("T")
        case 'a' => Nucleobase("A")
        case 'c' => Nucleobase("C")
        case 'g' => Nucleobase("G")
        case 't' => Nucleobase("T")
      }
    }
}

object Sequence{
  def makeSequence(added: Nucleobase, before: DNA): DNA = {
    before.add(added)
  }
  def fromString(name:String, input:String) = {
    val sequence = input.foldLeft[List[Nucleobase]](List())((bases,base) => Nucleobase.charToNucleobase(base) :: bases)
    DNA(name,sequence.reverse)
  }
  def GCCount(added: Nucleobase, before: (Int,Int)): (Int,Int) = {
    added match {
      case Nucleobase("C") => (before._1 + 1,before._2 + 1)
      case Nucleobase("G") => (before._1 + 1,before._2 + 1)
      case _ => (before._1,before._2 + 1)
    }
  }
  def compare(seq1: DNA, seq2: DNA, comparison: List[String], lineLength: Int): Unit = {
    def getNameLength(name1: String, name2: String): (String,String) = {
      val tabList = List("\t","\t","\t","\t")
      val tab1 = name1.length%5
      val tab2 = name2.length%5
      val difference = tab1 - tab2
      if (difference > 0) {
        (name1 + "\t",name2 + tabList.take(difference-1).foldLeft("")(_ + _))
      }else{
        (name1  + tabList.take(difference-1).foldLeft("")(_ + _),name2 + "\t")
      }
    }
    def printBases(bases1: List[Nucleobase], bases2: List[Nucleobase], comparison: List[String], name1: String, name2: String): Unit = {
      if (bases1.length != 0 && bases2.length != 0) {
        println(name2 + bases2.take(lineLength).foldLeft("")(_ + _))
        println(name1 + comparison.take(lineLength).foldLeft("")(_ + _))
        println(name1 + bases1.take(lineLength).foldLeft("")(_ + _))
        println("")
        printBases(bases1.drop(lineLength),bases2.drop(lineLength),comparison.drop(lineLength),name1,name2)
      }
    }
    if (seq1.nucleobases.length == seq2.nucleobases.length) {
      val bases1 = seq1.nucleobases
      val bases2 = seq2.nucleobases
      val names = getNameLength(seq1.name,seq2.name)
      printBases(bases1,bases2,comparison,names._1,names._2)
    }else{
      println("Sequences aren\'t the same length")
    }
  }
  def generateComparison(seq1: DNA, seq2: DNA): List[String] = {
    val bases1 = seq1.nucleobases
    val bases2 = seq2.nucleobases
    def generate(bases1: List[Nucleobase], bases2: List[Nucleobase]): List[String] = {
      if (bases1.length > 0 && bases2.length > 0) {
        bases1(0) match {
          case Nucleobase("C") => {
            bases1(1) match {
              case Nucleobase("G") => {
                bases2(0) match {
                  case Nucleobase("C") => "*" :: generate(bases1.drop(1),bases2.drop(1))
                  case Nucleobase("T") => "#" :: generate(bases1.drop(1),bases2.drop(1))
                  case _ => "X" :: generate(bases1.drop(1),bases2.drop(1))
                }
              }
              case notG => {
                bases2(0) match {
                  case Nucleobase("T") => ":" :: generate(bases1.drop(1),bases2.drop(1))
                  case Nucleobase("C") => "!" :: generate(bases1.drop(1),bases2.drop(1))
                  case _ => "X" :: generate(bases1.drop(1),bases2.drop(1))
                }
              }
            }
          }
          case other => {
            if (other == bases2(0)) {
              "|" :: generate(bases1.drop(1),bases2.drop(1))
            } else {
              "X" :: generate(bases1.drop(1),bases2.drop(1))
            }
          }
        }
      }else{
        List("")
      }
    }
    generate(bases2,bases1)
  }
  def methylation(list: List[String]): List[String] = {
    val methylationList = list.filter(comp => comp == "*" || comp == "#").map(comp => if (comp == "*") "M" else "U")
    methylationList
  }
  def toAnalysis(s: DNA, rlist: List[DNA], hit: Int, miss: Int, gapUp: Int, gapLeft: Int): Analysis = {
    //val best = rlist.flatMap(_.generateBisulfiteQuartet).map(ref => SmithWaterman.generateScoringMatrix(s,ref,hit,miss,gapUp,gapLeft)).maxBy(_.score)
    val best = rlist.map(ref => SmithWaterman.generateScoringMatrix(s,ref,hit,miss,gapUp,gapLeft)).maxBy(_.score)
    val r = best.reference
    val sw = SmithWaterman.completeAlignment(s,r,best)
    val input = sw._1
    val sequences = SmithWaterman.getSequences(input)
    val alignmentErrors = SmithWaterman.getErrors(input,r)
    val sequence = sequences._1
    val reference = sequences._2
    val methylationList = sequences._2.findCpgSites
    val methylationPoints = MethylationProcess.checkSites(sequences._2.findCpgSites,sequence)
    val comparison = generateComparison(sequence,reference)
    val sequenceLength = s.nucleobases.length
    val referenceLength = r.nucleobases.length
    val name = s.name
    val bisulfiteMap = comparison.groupBy(t => t).map(t => (t._1,t._2.length))
    val methylation = Methylation(bisulfiteMap.getOrElse("#",0),bisulfiteMap.getOrElse("*",0),methylationPoints, methylationList.map(_+1))
    val bisulfite = BisulfiteConversion(bisulfiteMap.getOrElse(":",0),bisulfiteMap.getOrElse("!",0),r.bisulfiteConversion.getOrElse(""))
    Analysis(sequenceName = name, sequenceLength = sequenceLength, referenceLength = referenceLength, referenceName = r.name, alignment = alignmentErrors, bisulfite = bisulfite, methylation = methylation, seqStart=sw._2._1, seqEnd=sw._2._2)
  }

}