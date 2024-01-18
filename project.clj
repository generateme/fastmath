(defproject generateme/fastmath "2.3.0"
  :description "Fast and primitive math library"
  :url "https://github.com/generateme/fastmath"
  :license {:name "The MIT Licence"
            :url "https://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.11.1"]
                 [net.jafama/jafama "2.3.2"]
                 [org.apache.commons/commons-math3 "3.6.1"]
                 [com.github.haifengl/smile-base "3.0.2"]
                 [com.github.haifengl/smile-core "3.0.2"]
                 [com.github.haifengl/smile-mkl "3.0.2"]
                 
                 [org.bytedeco/arpack-ng "3.9.0-1.5.9"]
                 [org.bytedeco/arpack-ng-platform "3.9.0-1.5.9"]
                 [org.bytedeco/openblas "0.3.23-1.5.9"]
                 [org.bytedeco/openblas-platform "0.3.23-1.5.9"]
                 [org.bytedeco/javacpp "1.5.9"]

                 [de.sciss/jwave "1.0.3"]
                 [ca.umontreal.iro.simul/ssj "3.3.1"]

                 [org.clojure/data.int-map "1.2.1"]
                 [org.clojure/math.combinatorics "0.2.0"]]
  :pedantic? false
  :resource-path "resources/"
  :java-source-paths ["src" "LBFGSBJava/src"]
  :javac-options ["--release" "8"  "-Xlint:unchecked"]
  :scm {:name "git"
        :url "https://github.com/generateme/fastmath/"}  
  :profiles {:1.10 {:dependencies [[org.clojure/clojure "1.10.3"]]}
             :eastwood {:plugins [[jonase/eastwood "1.4.0"]]
                        :eastwood {:add-linters [:performance :boxed-math]
                                   :exclude-namespaces [stats random utils bootstrap core complex-quaternion
                                                        calculus]}}
             :dev {:dependencies [[io.github.nextjournal/clerk "0.15.957"]
                                  [clojure2d "1.4.6-SNAPSHOT" :exclusions [generateme/fastmath]]
                                  [org.clojure/data.csv "1.0.1"]]
                   :source-paths ["notebooks"]}
             :dev-codox {:codox {:source-uri "https://github.com/generateme/fastmath/blob/master/{filepath}#L{line}"
                                 :namespaces [#"^fastmath\.(?!fields\.[a-z])"]}}})
