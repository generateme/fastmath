(defproject generateme/fastmath "1.1.2"
  :description "Fast and primitive math library"
  :url "https://github.com/generateme/fastmath"
  :license {:name "The Unlicence"
            :url "http://unlicense.org"}
  :dependencies [[org.clojure/clojure "1.9.0"]
                 [net.jafama/jafama "2.3.1"]
                 [org.apache.commons/commons-math3 "3.6.1"]
                 [com.github.haifengl/smile-interpolation "1.5.1"]
                 [com.github.haifengl/smile-core "1.5.1"]
                 [org.slf4j/slf4j-simple "1.7.25"]
                 [de.sciss/jwave "1.0.3"]]
  :resource-path "resources/"
  :java-source-paths ["src"]
  :javac-options ["-target" "1.8" "-source" "1.8"]
  :scm {:name "git"
        :url "https://github.com/generateme/fastmath/"}
  :profiles {:uberjar {:aot :all}
             :dev {:plugins [[refactor-nrepl "2.4.0-SNAPSHOT"]
                             [cider/cider-nrepl "0.18.0-SNAPSHOT"]
                             [com.jakemccrary/lein-test-refresh "0.22.0"]
                             [lein-codox "0.10.3"]]
                   :source-paths ["example"]
                   :dependencies [[codox-theme-rdash "0.1.2"]
                                  [clojure2d "1.0.0"]
                                  [criterium "0.4.4"]
                                  [incanter/incanter-charts "1.9.2"]
                                  [incanter/incanter-core "1.9.2"]
                                  [metadoc "0.2.3"]]
                   :codox {:themes [:rdash]
                           :metadata {:doc/format :markdown}
                           :output-path "docs/"
                           :source-paths ["src"]
                           :source-uri "https://github.com/generateme/fastmath/blob/master/{filepath}#L{line}"
                           :doc-paths ["docs/tutorials/"]
                           :writer metadoc.writers.codox/write-docs
                           :html {:transforms [[:head] [:append [:script {:type "text/javascript",
                                                                          :async ""
                                                                          :src "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML"}]]]}}}})
