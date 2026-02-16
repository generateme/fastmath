(ns fastmath.dev.wavelets
  "Parse txt files with wavel"
  (:require [clojure.string :as str]
            [clojure.java.io :as io]
            [fastmath.vector :as v]))

(def names '(beyl.txt     bl10.txt   db17.txt  db31.txt  db5.txt     la10.txt    mb8.3.txt    rbio6.8.txt  sym24.txt  sym39.txt
                        bior1.1.txt  bl7.txt    db18.txt  db32.txt  db6.txt     la12.txt    mb8.4.txt    sym10.txt    sym25.txt  sym4.txt
                        bior1.3.txt  bl9.txt    db19.txt  db33.txt  db7.txt     la14.txt    rbio1.1.txt  sym11.txt    sym26.txt  sym40.txt
                        bior1.5.txt  coif1.txt  db2.txt   db34.txt  db8.txt     la16.txt    rbio1.3.txt  sym12.txt    sym27.txt  sym41.txt
                        bior2.2.txt  coif2.txt  db20.txt  db35.txt  db9.txt     la18.txt    rbio1.5.txt  sym13.txt    sym28.txt  sym42.txt
                        bior2.4.txt  coif3.txt  db21.txt  db36.txt  dmey.txt    la20.txt    rbio2.2.txt  sym14.txt    sym29.txt  sym43.txt
                        bior2.6.txt  coif4.txt  db22.txt  db37.txt  fk14.txt    la8.txt     rbio2.4.txt  sym15.txt    sym3.txt   sym44.txt
                        bior2.8.txt  coif5.txt  db23.txt  db38.txt  fk18.txt    mb10.3.txt  rbio2.6.txt  sym16.txt    sym30.txt  sym45.txt
                        bior3.1.txt  db1.txt    db24.txt  db39.txt  fk22.txt    mb12.3.txt  rbio2.8.txt  sym17.txt    sym31.txt  sym5.txt
                        bior3.3.txt  db10.txt   db25.txt  db4.txt   fk4.txt     mb14.3.txt  rbio3.1.txt  sym18.txt    sym32.txt  sym6.txt
                        bior3.5.txt  db11.txt   db26.txt  db40.txt  fk6.txt     mb16.3.txt  rbio3.3.txt  sym19.txt    sym33.txt  sym7.txt
                        bior3.7.txt  db12.txt   db27.txt  db41.txt  fk8.txt     mb18.3.txt  rbio3.5.txt  sym2.txt     sym34.txt  sym8.txt
                        bior3.9.txt  db13.txt   db28.txt  db42.txt  han2.3.txt  mb24.3.txt  rbio3.7.txt  sym20.txt    sym35.txt  sym9.txt
                        bior4.4.txt  db14.txt   db29.txt  db43.txt  han3.3.txt  mb32.3.txt  rbio3.9.txt  sym21.txt    sym36.txt  vaid.txt
                        bior5.5.txt  db15.txt   db3.txt   db44.txt  han4.5.txt  mb4.2.txt   rbio4.4.txt  sym22.txt    sym37.txt
                        bior6.8.txt  db16.txt   db30.txt  db45.txt  han5.5.txt  mb8.2.txt   rbio5.5.txt  sym23.txt    sym38.txt))

(defn parse-file
  [file-symbol]
  (let [fs (str file-symbol)
        v (subs fs 0 (- (count fs) 4))
        [dl dh rl rh] (->> (str "wavelets/" fs)
                           (io/reader)
                           (line-seq)
                           (map #(mapv parse-double (str/split % #"[,\s]"))))]
    [v (if (or (and (not rl) (not rh)) ;; empty reconstruction
               (and (= dl (reverse rl)) (= dh (reverse rh)))) ;; not biorthogonal
         [dl dh]
         #_[dl (v/sub dh) (reverse rl) (v/sub (reverse rh))]
         [dl dh (reverse rl) (reverse rh)])]))

(spit "resources/wavelets/coeffs.edn" (pr-str (mapv parse-file names)))
