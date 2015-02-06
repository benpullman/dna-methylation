name := "dna-methylation-analyze"

version := "1.0"

scalaVersion := "2.11.2"

resolvers += "Scalaz Bintray Repo" at "http://dl.bintray.com/scalaz/releases"

resolvers += "spray" at "http://repo.spray.io/"

libraryDependencies ++= {
  Seq(
		"io.spray" %%  "spray-json" % "1.2.6",
		"org.apache.spark" %% "spark-core" % "1.2.0"
  )
}