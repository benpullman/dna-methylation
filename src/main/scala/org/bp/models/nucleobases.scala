package org.bp.models

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

class Sequence(name: Option[String] = None, nucleobases: List[Nucleobase]){
    override def toString = name.getOrElse("") + " " + this.nucleobases.foldLeft("")((nucleotide,nucleobase) => nucleotide + nucleobase)
}

case class DNA(name: Option[String] = None, nucleobases: List[Nucleobase]) extends Sequence(name: Option[String], nucleobases: List[Nucleobase]){
  def toRNA = RNA(nucleobases = this.nucleobases.map(nucleobase => nucleobase.toRNA))
  def reverseComplement = this.nucleobases.reverse.map(nucleobase => nucleobase.complement)
  def reverse = this.copy(nucleobases = this.nucleobases.reverse)
  def add(nucleobase: Nucleobase) = this.copy(nucleobases = nucleobase::nucleobases)
}

case class RNA(name: Option[String] = None, nucleobases: List[Nucleobase]) extends Sequence(name: Option[String], nucleobases: List[Nucleobase])

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
      }
    }
}

object Sequence{
  def FASTA(input: Array[Byte]) = {
    val output = input.foldLeft[(Boolean,List[Char],List[Nucleobase],List[DNA])]((true,List(),List(),List()))((current,newByte) =>
      newByte.toChar match {
        case '>' => current._2 match {
          case Nil => (true,List(),List(),current._4)
          case _ => (true,List(),List(),DNA(Some(current._2.reverse.mkString),current._3) :: current._4)
          }
        case '\n' => current._1 match {
          case true => (false,current._2,List(),current._4)
          case false => (false,current._2,current._3,current._4)
          }
        case b => current._1 match {
          case true => (true,b :: current._2,List(),current._4)
          case false => (false,current._2,Nucleobase.charToNucleobase(b) :: current._3,current._4)
          }
        }
    )
    DNA(Some(output._2.reverse.mkString),output._3) :: output._4
  }
  def makeSequence(added: Nucleobase,before: DNA): DNA = {
    before.add(added)
  }
  def GCCount(added: Nucleobase, before: (Int,Int)): (Int,Int) = {
    added match {
      case Nucleobase("C") => (before._1 + 1,before._2 + 1)
      case Nucleobase("G") => (before._1 + 1,before._2 + 1)
      case _ => (before._1,before._2 + 1)
    }
  }
}