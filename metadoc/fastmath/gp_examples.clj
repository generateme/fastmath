(ns fastmath.gp-examples
  (:require [fastmath.gp :refer :all]
            [metadoc.examples :refer :all]
            [fastmath.kernel :as k]))

(add-examples gaussian-process
  (example-session "Object creation"
    (gaussian-process [1 2 3] [-1 2 1])
    (gaussian-process [1 2 3] [-1 2 1] {:normalize? true
                                        :kernel (k/kernel :gaussian 0.5)
                                        :kscale 2.0
                                        :noise 0.1})))

(add-examples predict
  (example (let [gp (gaussian-process [-5 1 2] [17 10 12])]
             [(gp 1.1)
              (predict gp 1.1)]))
  (example "Predict and return standard deviation"
    (let [gp (gaussian-process [-5 1 2] [17 10 12])]
      [(gp 1.1 true)
       (predict gp 1.1 true)]))
  (example-image "Gaussian process with confidence intervals" "images/gp/ci.jpg")
  (example "Predict and return standard deviation for normalized process"
    (let [gp (gaussian-process [-5 1 2] [17 10 12] {:normalize? true})]
      [(gp 1.1 true)
       (predict gp 1.1 true)]))
  (example-image "Gaussian process with confidence intervals for normalized process" "images/gp/cin.jpg"))

(add-examples prior-samples
  (example (let [gp (gaussian-process [0 1 -2 -2.001] [-2 3 0.5 -0.6])]
             (prior-samples gp (range 0 1 0.1))))
  (example-image "Plot of 10 priors" "images/gp/prior.jpg")
  (example "With added noise" (let [gp (gaussian-process [0 1 -2 -2.001] [-2 3 0.5 -0.6] {:noise 0.1})]
                                (prior-samples gp (range 0 1 0.1))))
  (example-image "Plot of 10 priors (with noise)" "images/gp/priorn.jpg")
  ;; periodic causes LAPACK GETRF error in prior...
  #_(example "With periodic kernel" (let [gp (gaussian-process [0 1 -2 -2.001] [-2 3 0.5 -0.6] {:noise 0.01
                                                                                                :kernel (k/kernel :periodic 0.2)})]
                                      (prior-samples gp (range 0 1 0.1))))
  #_(example-image "Plot of 10 priors (with periodic kernel)" "images/gp/priorp.jpg")
  )

(add-examples posterior-samples
  (example "With gaussian kernel"
    (let [gp (gaussian-process [-5 1 2] [17 10 12] {:kernel (k/kernel :gaussian 0.5)})]
      (posterior-samples gp (range -5 2 0.9))))
  (example-image "Plot of 10 posteriors (with gaussian kernel)" "images/gp/posteriork.jpg")
  (example "With periodic kernel"
    (let [gp (gaussian-process [-5 1 2] [17 10 12] {:kernel (k/kernel :periodic 0.2 6.5)})]
      (posterior-samples gp (range -5 2 0.9))))
  (example-image "Plot of 10 posteriors (with periodic kernel)" "images/gp/posteriorkp.jpg"))

(add-examples predict-all
  (example (let [gp (gaussian-process [-5 1 2] [17 10 12])]
             (predict-all gp [-5 -2 1 2 3])))
  (example "Predict and return standard deviation"
    (let [gp (gaussian-process [-5 1 2] [17 10 12])]
      (predict-all gp [-5 -2 1 2 3] true))))
