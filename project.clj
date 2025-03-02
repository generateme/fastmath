(defproject generateme/fastmath "3.0.0-alpha4-SNAPSHOT"
  :description "Fast and primitive math library"
  :url "https://github.com/generateme/fastmath"
  :license {:name "The MIT Licence"
            :url "https://opensource.org/licenses/MIT"}
  :dependencies [[net.jafama/jafama "2.3.2"]
                 [org.apache.commons/commons-math3 "3.6.1"]

                 [de.sciss/jwave "1.0.3"]
                 [com.github.wendykierp/JTransforms "3.1"]
                 
                 [ca.umontreal.iro.simul/ssj "3.3.1"]

                 [org.clojure/data.int-map "1.3.0"]
                 [org.clojure/math.combinatorics "0.3.0"]]
  :pedantic? false
  :resource-path "resources/"
  :java-source-paths ["src" "LBFGSBJava/src"]
  :javac-options ["--release" "8"  "-Xlint:unchecked"]
  :scm {:name "git"
        :url "https://github.com/generateme/fastmath/"}  
  :profiles {:1.10 {:dependencies [[org.clojure/clojure "1.10.3"]]}
             :1.11 {:dependencies [[org.clojure/clojure "1.11.3"]]}
             :1.12 {:dependencies [[org.clojure/clojure "1.12.0"]]}
             :eastwood {:plugins [[jonase/eastwood "1.4.2"]]
                        :dependencies [[org.clojure/data.csv "1.1.0"]]
                        :eastwood {:add-linters [:performance :boxed-math :wrong-tag]
                                   :source-paths ["src"]
                                   :exclude-namespaces [:test-paths]}}
             :dev {:dependencies [[org.clojure/clojure "1.12.0"]
                                  [org.clojure/data.csv "1.1.0"]
                                  [org.scicloj/clay "2-beta30"]
                                  [scicloj/clojisr "1.0.0"]
                                  [com.github.lread/test-doc-blocks "1.1.20"]]
                   :source-paths ["notebooks" "utils"]}
             :dev-codox {:codox {:source-uri "https://github.com/generateme/fastmath/blob/master/{filepath}#L{line}"
                                 :namespaces [#"^fastmath\.(?!fields\.[a-z])"]}}}
  :aliases {"tests-with-md" ["with-profile" "dev" "do"
                             ["run" "-m" "lread.test-doc-blocks" "gen-tests"
                              "--platform" "clj"
                              "test/docs/*.md"]
                             ["test"]]})
