name := "dna-methylation-analyze"

version := "1.0"

scalaVersion := "2.11.1"

resolvers += "Scalaz Bintray Repo" at "http://dl.bintray.com/scalaz/releases"

libraryDependencies ++= {
  Seq(
		"org.spire-math" %% "spire" % "0.7.5"
  )
}