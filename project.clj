(defproject generateme/fastmath "3.0.0-alpha8-SNAPSHOT"
  :description "Fast and primitive math library"
  :url "https://github.com/generateme/fastmath"
  :license {:name "The MIT Licence"
            :url "https://opensource.org/licenses/MIT"}
  :dependencies [;; all general stuff
                 [net.jafama/jafama "2.3.2"]
                 [org.apache.commons/commons-math3 "3.6.1"]

                 ;; Wavelets and FFT
                 [de.sciss/jwave "1.0.3"]
                 [com.github.wendykierp/JTransforms "3.2"]

                 ;; IIR filters
                 [uk.me.berndporr/iirj "1.7"]

                 ;; Some distributions, interpolations, Eigendecomposition
                 [ca.umontreal.iro.simul/ssj "3.3.2"
                  :exclusions [org.jfree/jfreechart
                               org.jfree/jcommon]]

                 ;; discrete distribution
                 [org.clojure/data.int-map "1.3.1"]

                 ;; integration (permutations)
                 [org.clojure/math.combinatorics "0.3.2"]]
  :pedantic? false
  :resource-path "resources/"
  :java-source-paths ["src" "LBFGSBJava/src"]
  :javac-options ["--release" "8"  "-Xlint:unchecked"]
  :jvm-opts ["-Djdk.attach.allowAttachSelf=true"]
  :scm {:name "git"
        :url "https://github.com/generateme/fastmath/"}  
  :profiles {:1.10 {:dependencies [[org.clojure/clojure "1.10.3"]]}
             :1.11 {:dependencies [[org.clojure/clojure "1.11.4"]]}
             :1.12 {:dependencies [[org.clojure/clojure "1.12.5"]]}
             :1.13 {:dependencies [[org.clojure/clojure "1.13.0-alpha4"]]}
             :eastwood {:plugins [[jonase/eastwood "1.4.3"]]
                        :dependencies [[org.clojure/data.csv "1.1.0"]]
                        :eastwood {:add-linters [:performance :boxed-math :wrong-tag]
                                   :source-paths ["src"]
                                   :exclude-namespaces [:test-paths]}}
             :dev {:dependencies [;;[org.clojure/clojure "1.13.0-alpha4"]
                                  [org.clojure/clojure "1.12.5"]
                                  [org.clojure/data.csv "1.1.0"]
                                  [org.clojure/data.json "2.5.2"]
                                  [org.scicloj/clay "2.0.13"]
                                  [zprint "1.3.0"]
                                  [scicloj/clojisr "1.0.0"]
                                  ;; [scicloj/tablecloth "7.059"]
                                  [virgil "0.5.0"]
                                  [org.ow2.asm/asm "9.9"]
                                  [clojure2d/clojure2d "1.5.0-alpha2-SNAPSHOT"]
                                  [org.scicloj/plotje "0.8.1" :exclusions [clojure2d/clojure2d
                                                                           generateme/fastmath]]]
                   :source-paths ["notebooks" "utils"]}
             :dev-codox {:codox {:source-uri "https://github.com/generateme/fastmath/blob/master/{filepath}#L{line}"
                                 :namespaces [#"^fastmath\.(?!fields\.[a-z])"]}}}
  :aliases {"tests-with-md" ["with-profile" "dev" "do"
                             ["run" "-m" "lread.test-doc-blocks" "gen-tests"
                              "--platform" "clj"
                              "test/docs/*.md"]
                             ["test"]]})
