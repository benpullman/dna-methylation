package org.bp.models

import org.bp.alignment._
import org.bp.methylation._
import scala.io.Source

case class Nucleobase(name: String){
  /*override def toString = this.name match {
    case "cT" => "T"
    case "cA" => "A"
    case other => other
  }*/
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
    case Nucleobase("cT") => Nucleobase("cA")
    case Nucleobase("cA") => Nucleobase("cT")
  }
}

class Sequence(name: String, nucleobases: IndexedSeq[Nucleobase], bisulfiteConversion: Option[String] = None, direction: Option[String] = None){
    override def toString = "Name: " + name + bisulfiteConversion.map(" " + _ + " ").getOrElse("") + direction.map(_ + " ").getOrElse("") + "\nSequence: " + this.nucleobases.foldLeft("")((nucleotide,nucleobase) => nucleotide + nucleobase)
    def findCpgSites = this.nucleobases.foldLeft[(Boolean,Int,IndexedSeq[Int])]((false,0,IndexedSeq()))((acc,next) => {
      if (acc._1 && next == Nucleobase("G")) {
        (false, acc._2 + 1, (acc._2 - 1) +: acc._3)
      } else if (next == Nucleobase("C")){
        (true, acc._2 + 1, acc._3)
      } else {
        (false, acc._2 + 1, acc._3)
      }
    })._3.reverse
    def length = this.nucleobases.length
}

case class DNA(name: String, nucleobases: IndexedSeq[Nucleobase], bisulfiteConversion: Option[String] = None, direction: Option[String] = None) extends Sequence(name: String, nucleobases: IndexedSeq[Nucleobase],bisulfiteConversion: Option[String], direction: Option[String]){
  //def toRNA = RNA(name = name, nucleobases = this.nucleobases.map(nucleobase => nucleobase.toRNA))
  def reverseComplement = this.nucleobases.reverse.map(nucleobase => nucleobase.complement)
  def reverse = this.copy(nucleobases = this.nucleobases.reverse)
  def add(nucleobase: Nucleobase) = this.copy(nucleobases = nucleobase+:nucleobases)
  def bisulfite: List[DNA] = {
    val forward = Map(IndexedSeq(Nucleobase("C"),Nucleobase("G")) -> IndexedSeq(Nucleobase("G"), Nucleobase("C")))
    val reverse = Map(IndexedSeq(Nucleobase("C"),Nucleobase("G")) -> IndexedSeq(Nucleobase("C"), Nucleobase("G")))
    val CtoT = Map(Nucleobase("C") -> Nucleobase("cT"))
    val GtoA = Map(Nucleobase("G") -> Nucleobase("cA"))
    def forwardCtoT = Sequence.convert(forward,false,CtoT)_
    def reverseComplementCtoT = Sequence.convert(reverse,true,CtoT)_
    def forwardGtoA = Sequence.convert(forward,false,GtoA)_
    def reverseComplementGtoA = Sequence.convert(reverse,true,GtoA)_
    List(
      DNA(this.name,forwardCtoT(this).reverse,Some("C => T"),Some("Forward")),
      DNA(this.name,reverseComplementCtoT(this),Some("C => T"),Some("Reverse Complement")),
      DNA(this.name,forwardGtoA(this).reverse,Some("G => A"),Some("Forward")),
      DNA(this.name,reverseComplementGtoA(this),Some("G => A"),Some("Reverse Complement"))
    )
  }
  def addMethylation = {
    (this,this.findCpgSites)
  }
  
  def generateBisulfiteQuartet = {
      val CtoT = Map(Nucleobase("C") -> Nucleobase("T"))
      val GtoA = Map(Nucleobase("G") -> Nucleobase("A"))
      def bisulfiteConversionForward(input: DNA, conversion: Map[Nucleobase,Nucleobase], conversionString: String) = {
        val sequence = input.nucleobases.foldLeft[(IndexedSeq[Nucleobase],Option[Nucleobase])]((IndexedSeq(),None))((bases,next) => {
          next match {
            case Nucleobase("C") => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(Nucleobase("C"), Nucleobase("C")) +: bases._1, Some(Nucleobase("C")))
              case _ => (bases._1, Some(Nucleobase("C")))
            }
            case Nucleobase("G") => bases._2 match {
              case Some(Nucleobase("C")) => (Nucleobase("G") +: Nucleobase("C") +: bases._1,None)
              case _ => (conversion.getOrElse(Nucleobase("G"), Nucleobase("G")) +: bases._1,None)
            }
            case base => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(base,base) +: conversion.getOrElse(Nucleobase("C"),Nucleobase("C")) +: bases._1,None)
              case _ => (conversion.getOrElse(base,base) +: bases._1,None)
            }
          }
        })._1.reverse
        DNA(input.name,sequence,Some(conversionString),Some("Forward"))
      }
      def bisulfiteConversionReverseComplement(input: DNA, conversion: Map[Nucleobase,Nucleobase], conversionString: String) = {
        val sequence = input.nucleobases.foldLeft[(IndexedSeq[Nucleobase],Option[Nucleobase])]((IndexedSeq(),None))((bases,next) => {
          next match {
            case Nucleobase("C") => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(Nucleobase("C"), Nucleobase("C")).complement +: bases._1, Some(Nucleobase("C")))
              case _ => (bases._1, Some(Nucleobase("C")))
            }
            case Nucleobase("G") => bases._2 match {
              case Some(Nucleobase("C")) => (Nucleobase("C") +: Nucleobase("G") +: bases._1,None)
              case _ => (conversion.getOrElse(Nucleobase("G"), Nucleobase("G")).complement +: bases._1,None)
            }
            case base => bases._2 match {
              case Some(Nucleobase("C")) => (conversion.getOrElse(base,base).complement +: conversion.getOrElse(Nucleobase("C"),Nucleobase("C")).complement +: bases._1,None)
              case _ => (conversion.getOrElse(base,base).complement +: bases._1,None)
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

  object DNA{
    def read(input: String) = {
    //reading by char
    // two char sequences - 
      val output = input.foldLeft[(Boolean,String,IndexedSeq[Nucleobase],List[DNA])]((true, "", IndexedSeq(), List()))((out, add) => {
      add match {
        case '>' => {
          (true, "", IndexedSeq(),  DNA(out._2, out._3.reverse) :: out._4)
        }
        case '\n' => if (out._1) (false, out._2, out._3, out._4) else out
        case '\r' => if (out._1) (false, out._2, out._3, out._4) else out
        case other => if (out._1) (out._1, out._2 + add, out._3, out._4) else (out._1, out._2, Nucleobase.charToNucleobase(add) +: out._3, out._4)
        }
      }
    )
    (DNA(output._2,output._3.reverse) :: output._4).reverse.drop(1)

    }
    def readFromFile(filename: String) = {
    //reading by char
    // two char sequences - 
    val output = Source.fromFile(filename).foldLeft[(Boolean,String,IndexedSeq[Nucleobase],List[DNA])]((true, "", IndexedSeq(), List()))((out, add) => {
      add match {
        case '>' => {
          (true, "", IndexedSeq(),  DNA(out._2, out._3.reverse) :: out._4)
        }
        case '\n' => if (out._1) (false, out._2, out._3, out._4) else out
        case '\r' => if (out._1) (false, out._2, out._3, out._4) else out
        case other => if (out._1) (out._1, out._2 + add, out._3, out._4) else (out._1, out._2, Nucleobase.charToNucleobase(add) +: out._3, out._4)
        }
      }
    )
    (DNA(output._2,output._3.reverse) :: output._4).reverse.drop(1)
  }
    def CpGSiteMap(sequences: List[DNA]): Map[(String,String),IndexedSeq[Int]] = {
      sequences.map(sequence => ((sequence.bisulfiteConversion.get,sequence.direction.get),sequence.findCpgSites)).toMap
    }
  }

//case class RNA(name: String, nucleobases: List[Nucleobase], bisulfiteConversion: Option[String] = None, direction: Option[String] = None) extends Sequence(name: String, nucleobases: List[Nucleobase],bisulfiteConversion: Option[String], direction: Option[String])

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
    val sequence = input.foldLeft[IndexedSeq[Nucleobase]](IndexedSeq())((bases,base) => Nucleobase.charToNucleobase(base) +: bases)
    DNA(name,sequence.reverse)
  }
  def GCCount(added: Nucleobase, before: (Int,Int)): (Int,Int) = {
    added match {
      case Nucleobase("C") => (before._1 + 1,before._2 + 1)
      case Nucleobase("G") => (before._1 + 1,before._2 + 1)
      case _ => (before._1,before._2 + 1)
    }
  }
  
  def convert(direction: Map[IndexedSeq[Nucleobase],IndexedSeq[Nucleobase]], complement: Boolean, bisulfite: Map[Nucleobase,Nucleobase])(sequence: DNA) = {
    def applyComplement(base: Nucleobase): Nucleobase = {
      if(complement){
        base.complement
      }else{
        base
      }
    }
    sequence.nucleobases.foldLeft[(Option[Nucleobase],IndexedSeq[Nucleobase])]((None,IndexedSeq()))((acc,next) => {
      acc._1 match {
        case Some(nucleobase) => {
          if (next == Nucleobase("C")){
            (Some(Nucleobase("C")),applyComplement(bisulfite.getOrElse(nucleobase,nucleobase)) +: acc._2)
          }else{
            val next2 = direction.getOrElse(IndexedSeq(nucleobase,next),IndexedSeq(bisulfite.getOrElse(next,next),bisulfite.getOrElse(nucleobase,nucleobase)).map(applyComplement(_)))
            (None, next2 ++ acc._2)
          }
        }
        case None => {
          next match {
            case Nucleobase("C") => (Some(Nucleobase("C")),acc._2)
            case nucleobase => (None,applyComplement(bisulfite.getOrElse(nucleobase,nucleobase)) +: acc._2)
          }
        }
      }
    })._2
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
    def printBases(bases1: IndexedSeq[Nucleobase], bases2: IndexedSeq[Nucleobase], comparison: List[String], name1: String, name2: String): Unit = {
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
      val names = getNameLength("genome\t\t","bisulfite")
      printBases(bases1,bases2,comparison,names._1,names._2)
    }else{
      println("Sequences aren\'t the same length")
    }
  }
  /*
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
                  case Nucleobase("C") => "*" +: generate(bases1.drop(1),bases2.drop(1))
                  case Nucleobase("T") => "#" +: generate(bases1.drop(1),bases2.drop(1))
                  case _ => "X" +: generate(bases1.drop(1),bases2.drop(1))
                }
              }
              case notG => {
                bases2(0) match {
                  case Nucleobase("T") => ":" +: generate(bases1.drop(1),bases2.drop(1))
                  case Nucleobase("C") => "!" +: generate(bases1.drop(1),bases2.drop(1))
                  case _ => "X" +: generate(bases1.drop(1),bases2.drop(1))
                }
              }
            }
          }
          case other => {
            if (other == bases2(0)) {
              "|" +: generate(bases1.drop(1),bases2.drop(1))
            } else {
              "X" +: generate(bases1.drop(1),bases2.drop(1))
            }
          }
        }
      }else{
        List("")
      }
    }
    generate(bases2,bases1)
  }
  */
  def generateComparison(sequence: DNA, reference: DNA): List[String] = {
    val ref = reference.nucleobases
    val seq = sequence.nucleobases
    def checkConversion(check:Nucleobase,reference:Nucleobase, back: Map[Nucleobase,Nucleobase]): String =  {
      if (check == reference){
        ":"
      }else if (check == back.getOrElse(reference,reference)){
        "!"
      }else{
        "X"
      }
    }
    def generate(bisulfiteConversion: String, direction: String)(seq: IndexedSeq[Nucleobase], ref: IndexedSeq[Nucleobase]): List[String] = {
      val case1 = ((bisulfiteConversion == "C => T") && (direction == "Forward")) || ((bisulfiteConversion == "G => A") && (direction == "Reverse Complement"))
      val case2 = ((bisulfiteConversion == "G => A") && (direction == "Forward")) || ((bisulfiteConversion =="C => T") && (direction == "Reverse Complement"))
      val back = if (case1) {
        Map(Nucleobase("T") -> Nucleobase("C"))
      }else{
        Map(Nucleobase("A") -> Nucleobase("G"))
      }
      val gen = generate(bisulfiteConversion,direction)_
      if (seq.length > 1 && ref.length > 1) {
        (ref(0),ref(1)) match {
          case (Nucleobase("C"),Nucleobase("G")) => {
            (seq(0),seq(1)) match {
              case (Nucleobase("C"),_) if case1 => "*" +: gen(seq.drop(1),ref.drop(1))
              case (_,Nucleobase("G")) if case2 => "*" +: gen(seq.drop(1),ref.drop(1))
              case (Nucleobase("T"),Nucleobase("G")) if case1 => "#" +: gen(seq.drop(1),ref.drop(1))
              case (Nucleobase("C"),Nucleobase("A")) if case2 => "#" +: gen(seq.drop(1),ref.drop(1))
              case (_,_) => "X" +: gen(seq.drop(1),ref.drop(1))
            }
          }
          case (Nucleobase("cT"),_) if case1 => {
            checkConversion(seq(0),Nucleobase("T"),back) +: gen(seq.drop(1),ref.drop(1))
          }
          case (Nucleobase("cA"),_) if case2 => {
            checkConversion(seq(0),Nucleobase("A"),back) +: gen(seq.drop(1),ref.drop(1))
          }
          case (notC,notG) => {
            if (notC == seq(0)) {
              "|" +: gen(seq.drop(1),ref.drop(1))
            } else {
              "X" +: gen(seq.drop(1),ref.drop(1))
            }
          }
        }
      }else if (seq.length==1 && ref.length==1){
        if(seq(0)==ref(0)) List("|") else List("X")
      }else{
        List("")
      }
    }
    generate(reference.bisulfiteConversion.getOrElse(""),reference.direction.getOrElse(""))(seq,ref)
  }
  /*
  def errorsFromComparison(comparison:List[String],sequence:DNA) = {
    val direction = sequence.direction.getOrElse("")
    val errors = comparison.foldLeft(0)((acc,next) => if (next=="X") acc + 1 else acc)
    val gaps = sequence.nucleobases.foldLeft(0)((acc,next) => if (next==Nucleobase("-")) acc + 1 else acc)
    Alignment(mismatches = errors, gaps = gaps, direction = direction, length = sequence.nucleobases.length)
  }
  */
  def methylation(list: List[String]): List[String] = {
    val methylationList = list.filter(comp => comp == "*" || comp == "#").map(comp => if (comp == "*") "M" else "U")
    methylationList
  }
  /*
  def toAnalysis(s: DNA, rlist: List[DNA], hit: Int, miss: Int, gapUp: Int, gapLeft: Int, barcode: String): Analysis = {
    val start: Double = System.currentTimeMillis / 1000.0
    val options = SW.run(s,rlist,hit,miss,gapUp,gapLeft)
    val best = options.maxBy(_.score)
    val endSequences: Double = System.currentTimeMillis / 1000.0
    val r = best.reference
    val sw = SW.readAlignment(s,r,best)
    val input = sw._1
    val sequences = SW.getSequences(input)
    val alignmentErrors = SW.readErrors(r,input)
    val sequence = sequences._1
    val reference = sequences._2
    println(sequence)
    println(reference)
    val methylationList = sequences._2.findCpgSites
    val methylationPoints = MethylationProcess.checkSites(sequences._2.findCpgSites,sequence)
    val comparison = generateComparison(sequence,reference)
    val sequenceLength = s.nucleobases.length
    val referenceLength = r.nucleobases.length
    val name = s.name
    val bisulfiteMap = comparison.groupBy(t => t).map(t => (t._1,t._2.length))
    val methylation = Methylation(bisulfiteMap.getOrElse("#",0),bisulfiteMap.getOrElse("*",0),methylationPoints, methylationList.map(_+1))
    val bisulfite = BisulfiteConversion(bisulfiteMap.getOrElse(":",0),bisulfiteMap.getOrElse("!",0),r.bisulfiteConversion.getOrElse(""))
    val finish: Double = System.currentTimeMillis / 1000.0
    println("Sequence process time: " + (endSequences.toDouble-start.toDouble) + " Other process time: " + (finish.toDouble - endSequences.toDouble))
    Analysis(sequenceName = name, sequenceLength = sequenceLength, referenceLength = referenceLength, referenceName = r.name, alignment = alignmentErrors, bisulfite = bisulfite, methylation = methylation, seqStart=sw._2._1, seqEnd=sw._2._2,barcode=barcode)
  }
  */
  /*def toAnalysis(s: DNA, rlist: List[DNA], hit: Int, miss: Int, gapUp: Int, gapLeft: Int, barcode: String): Analysis = {
    val start: Double = System.currentTimeMillis / 1000.0
    //val best = rlist.flatMap(_.generateBisulfiteQuartet).map(ref => SmithWaterman.generateScoringMatrix(s,ref,hit,miss,gapUp,gapLeft)).maxBy(_.score)
    //val all = rlist.par.map(ref => SmithWaterman.generateScoringMatrix(s,ref,hit,miss,gapUp,gapLeft))
    //all.map(a => println(a.score + " " + a.reference.bisulfiteConversion))
    val best = rlist.par.map(ref => SmithWaterman.generateScoringMatrix(s,ref,hit,miss,gapUp,gapLeft)).maxBy(_.score)
    val endSequences: Double = System.currentTimeMillis / 1000.0
    val r = best.reference
    val sw = SmithWaterman.completeAlignment(s,r,best)
    val input = sw._1
    val sequences = SmithWaterman.getSequences(input,r)
    val sequence = sequences._1
    val reference = sequences._2
    val methylationList = sequences._2.findCpgSites
    val methylationPoints = MethylationProcess.checkSites(sequences._2.findCpgSites,sequence)
    val comparison = generateComparison(sequence,reference)
    val alignmentErrors = errorsFromComparison(comparison,sequence)
    //compare(reference,sequence,comparison,60)
    val sequenceLength = s.nucleobases.length
    val referenceLength = r.nucleobases.length
    val name = s.name
    val bisulfiteMap = comparison.groupBy(t => t).map(t => (t._1,t._2.length))
    val methylation = Methylation(bisulfiteMap.getOrElse("#",0),bisulfiteMap.getOrElse("*",0),methylationPoints, methylationList.map(_+1))
    val bisulfite = BisulfiteConversion(bisulfiteMap.getOrElse(":",0),bisulfiteMap.getOrElse("!",0),r.bisulfiteConversion.getOrElse(""))
    val finish: Double = System.currentTimeMillis / 1000.0
    println("Sequence process time: " + (endSequences.toDouble-start.toDouble) + " Other process time: " + (finish.toDouble - endSequences.toDouble))
    Analysis(sequenceName = name, sequenceLength = sequenceLength, referenceLength = referenceLength, referenceName = r.name, alignment = alignmentErrors, bisulfite = bisulfite, methylation = methylation, seqStart=sw._2._1, seqEnd=sw._2._2,barcode=barcode)
  }
  
  def bestAlignment(s: DNA, rlist: List[DNA], hit: Int, miss: Int, gapUp: Int, gapLeft: Int) = {
    val best = rlist.map(r => DNAAlignment.run(s,r)).maxBy(_.params.score)
    best
  }
  def methylationProcess(s: DNA, r: DNA, bisulfiteMap: Map[String,Int]) = {
    val methylationList = r.findCpgSites
    val methylationPoints = MethylationProcess.checkSites(methylationList,s)
    Methylation(bisulfiteMap.getOrElse("#",0),bisulfiteMap.getOrElse("*",0),methylationPoints, methylationList.map(_+1))
  }
  def bisulfiteProcess(s: DNA, r: DNA, bisulfiteMap: Map[String,Int]) = {
    BisulfiteConversion(bisulfiteMap.getOrElse(":",0),bisulfiteMap.getOrElse("!",0),r.bisulfiteConversion.getOrElse(""))
  }
  def toAnalysis(s: DNA, rlist: List[DNA], hit: Int, miss: Int, gapUp: Int, gapLeft: Int, barcode: String): Analysis = {
    val start: Double = System.currentTimeMillis / 1000.0
    val best = bestAlignment(s,rlist,hit,miss,gapUp,gapLeft)
    val endSequences: Double = System.currentTimeMillis / 1000.0
    val bisulfiteMap = generateComparison(best.sample, best.reference).groupBy(t => t).map(t => (t._1,t._2.length))
    val methylation = methylationProcess(best.sample, best.reference, bisulfiteMap)
    val bisulfite = bisulfiteProcess(best.sample, best.reference, bisulfiteMap)
    val finish: Double = System.currentTimeMillis / 1000.0
    println("Sequence process time: " + (endSequences.toDouble-start.toDouble) + " Other process time: " + (finish.toDouble - endSequences.toDouble))
    Analysis(
      sequenceName = best.sample.name,
      sequenceLength = best.sample.nucleobases.length,
      referenceLength = best.reference.nucleobases.length,
      referenceName = best.reference.name,
      alignment = AlignmentAnalysis(best.params.mismatches, best.params.gaps, best.sample.direction.getOrElse("None"), best.length),
      bisulfite = bisulfite,
      methylation = methylation,
      seqStart= best.params.start,
      seqEnd= best.params.end,
      barcode=barcode
    )
  }
  def process(newSequence: DNA, references: List[DNA], sampleMap: Map[String,String]) = {
    println(newSequence.name)
    val tag = newSequence.nucleobases.take(18).foldLeft("")(_ + _)
    val nucleobases = newSequence.nucleobases.drop(18)
    //val barcode = sampleMap.getOrElse(tag,"none")
    val barcode = Sequence.findBestBarcode(tag,sampleMap)
    //println(tag + " " + barcode)
    if (barcode == "error"){
      None
    }else{
      Some(Sequence.toAnalysis(newSequence,references,5,-3,-20,-2,barcode))
    }
    //println(tag)
  }
  */
  def findBestBarcode(input: String, sampleMap: Map[String,String]) = {
    val bestSample = sampleMap.keys.map(barcodeAlignment(input,_,10,-1,-5,-5)).maxBy(_.score)
    def completeAlignment(input: String, alignmentProcess: AlignmentProcess) = {
      var start = 0
      var end = alignmentProcess.endX
      val sample = alignmentProcess.sample
      def read(i: Int,j: Int,directionMatrix:Array[Array[String]]): List[(Char, Char)] = {
        directionMatrix(i)(j) match {
          case "diag" => (input(i-1),sample(j-1)) +: read(i-1,j-1,directionMatrix)
          case "upup" => (input(i-1),'-') +: read(i-1,j,directionMatrix)
          case "left" => ('-',sample(j-1)) +: read(i,j-1,directionMatrix)
          case _ => {
            start = i+1
            List(('s','s'))
          }
        }
      }
      read(alignmentProcess.endX,alignmentProcess.endY,alignmentProcess.directionMatrix).reverse.drop(1)
    }

    def getErrors(pairs: List[(Char,Char)]): Int = {
      val error = pairs.foldLeft(0)((tally,next) => {
        val sample = next._2
        val input = next._1
        if (sample==input){
          tally
        } else {
          tally + 1
        }
      })
      error
    }

    val error = getErrors(completeAlignment(input,bestSample))
    //println(error)
    if (error < 4) {
      sampleMap.getOrElse(bestSample.sample,"error")
    } else {
      "error"
    }
  }

  def barcodeAlignment(input: String, sample: String, hit: Int, miss: Int, gapUp: Int, gapLeft: Int) = {
    var largest = (0,(0,0))
    var start = 0
    var end = 0
    var directionMatrix = Array.ofDim[String](input.length + 1,sample.length+1)
    var scoringMatrix = Array.ofDim[Int](input.length + 1,sample.length+1)
    def initializeScoringMatrix: Unit = {
      for (i <- 1 to input.length; j <- 1 to sample.length){
        scoringMatrix(i)(j) = -1
      }
    }
    def matchScore(i: Int, j: Int): Int = {
      if (input(i-1) == sample(j-1)) hit else miss
    }
    def score(i:Int,j:Int): Int = {
      scoringMatrix(i)(j) match {
        //case -1 => List(score(i-1,j-1) + matchScore(i,j), score(i-1,j) + gap, score(i,j-1) + gap, 0).max
        case -1 => {
          val diag = score(i-1,j-1) + matchScore(i,j)
          val left = score(i,j-1) + gapLeft
          val up = score(i-1,j) + gapUp
          val best = List(diag, left, up, 0).max
          if (best > largest._1) largest = (best,(i,j))
          best match {
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
            case 0 => {
              directionMatrix(i)(j) = "stop"
              0
            }
            case _ => -100000
          }
        }
        case a => a
      }
    }
    def setScore: Unit = {
      for (i <- 0 to input.length; j <- 0 to sample.length){
        scoringMatrix(i)(j) = score(i,j)
      }
    }

    initializeScoringMatrix
    setScore
    AlignmentProcess(scoringMatrix, directionMatrix, sample, largest._2._1,largest._2._2, largest._1)
    //directionMatrix.foreach{case i => i foreach {j => print(j + " ")}; print('\n')}
    //end = largest._2._1
    //println(largest._1)
    //(read(largest._2._1,largest._2._2).reverse.drop(1),(start,end))
  }

}