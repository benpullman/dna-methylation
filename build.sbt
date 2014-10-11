name := "dna-methylation-analyze"

version := "1.0"

scalaVersion := "2.11.2"

resolvers += "Scalaz Bintray Repo" at "http://dl.bintray.com/scalaz/releases"

resolvers += "spray" at "http://repo.spray.io/"

libraryDependencies ++= {
  Seq(
		"org.spire-math" %% "spire" % "0.7.5",
		"io.spray" %%  "spray-json" % "1.2.6"
  )
}