(defproject generateme/fastmath "1.3.0-SNAPSHOT"
  :description "Fast and primitive math library"
  :url "https://github.com/generateme/fastmath"
  :license {:name "The Unlicence"
            :url "http://unlicense.org"}
  :dependencies [[org.clojure/clojure "1.10.0"]
                 [net.jafama/jafama "2.3.1"]
                 [org.apache.commons/commons-math3 "3.6.1"]
                 [com.github.haifengl/smile-interpolation "1.5.2" ]
                 [com.github.haifengl/smile-core "1.5.2"]
                 [com.github.haifengl/smile-netlib "1.5.2"]
                 [org.slf4j/slf4j-nop "1.7.25"]
                 [de.sciss/jwave "1.0.3"]
                 [net.mikera/core.matrix "0.62.0"]
                 [clj-boost "0.0.4"]
                 [de.bwaldvogel/liblinear "2.21"]]
  :exclusions [[org.slf4j/slf4j-simple]]
  :resource-path "resources/"
  :java-source-paths ["src"]
  :javac-options ["-target" "1.8" "-source" "1.8"]
  :scm {:name "git"
        :url "https://github.com/generateme/fastmath/"}
  :profiles {:uberjar {:aot :all}
             :dev {:plugins [[refactor-nrepl "2.4.0"]
                             [cider/cider-nrepl "0.20.0"]
                             [com.jakemccrary/lein-test-refresh "0.23.0"]
                             [lein-codox "0.10.5"]]
                   :source-paths ["example"]
                   :dependencies [[codox-theme-rdash "0.1.2"]
                                  [clojure2d "1.0.2"]
                                  [criterium "0.4.4"]
                                  [incanter/incanter-charts "1.9.2"]
                                  [incanter/incanter-core "1.9.2"]
                                  [metadoc "0.2.4"]]
                   :exclusions [[asm]]
                   ;; :exclusions [[org.slf4j/slf4j-simple] [ml.dmlc/xgboost4j] [asm]]
                   ;; :resource-paths ["resources/" "lib/xgboost4j-0.81-criteo-20180821_2.11-win64.jar"]
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
