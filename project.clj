(defproject generateme/fastmath "2.1.2"
  :description "Fast and primitive math library"
  :url "https://github.com/generateme/fastmath"
  :license {:name "The MIT Licence"
            :url "https://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.10.2"]
                 [net.jafama/jafama "2.3.1"]
                 [org.apache.commons/commons-math3 "3.6.1"]
                 [com.github.haifengl/smile-interpolation "2.6.0"]
                 [com.github.haifengl/smile-core "2.6.0"]
                 [com.github.haifengl/smile-mkl "2.6.0"]
                 
                 [org.bytedeco/arpack-ng "3.7.0-1.5.4"]
                 [org.bytedeco/arpack-ng-platform "3.7.0-1.5.4"]
                 [org.bytedeco/openblas "0.3.10-1.5.4"]
                 [org.bytedeco/openblas-platform "0.3.10-1.5.4"]
                 [org.bytedeco/javacpp "1.5.4"]

                 [de.sciss/jwave "1.0.3"]
                 [ca.umontreal.iro.simul/ssj "3.3.1"]]
  :pedantic? false
  :resource-path "resources/"
  :java-source-paths ["src"]
  :javac-options ["-target" "1.8" "-source" "1.8"]
  :scm {:name "git"
        :url "https://github.com/generateme/fastmath/"}  
  :profiles {:dev-codox {:codox {:source-uri "https://github.com/generateme/fastmath/blob/master/{filepath}#L{line}"
                                 :namespaces [#"^fastmath\."]}}})
