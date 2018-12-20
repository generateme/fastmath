package fastmath.java;

import org.apache.commons.math3.random.RandomVectorGenerator;

public class R2 implements RandomVectorGenerator {

    private final static double phi1 = 1.0 / 1.6180339887498948482;
    private final static double phi2 = 1.0 / 1.3247179572447460259609088544780973;
    private final static double phi3 = 1.0 / 1.2207440846057594753616853491088319;
    private final static double phi4 = 1.0 / 1.16730397826141868425604589985484218072;

    private final static double[] phi1_t = { phi1 };
    private final static double[] phi2_t = { phi2, phi2*phi2 };
    private final static double[] phi3_t = { phi3, phi3*phi3, phi3*phi3*phi3 };
    private final static double[] phi4_t = { phi4, phi4*phi4, phi4*phi4*phi4, phi4*phi4*phi4*phi4};
    
    private int dimensions;
    private double[] current, phis;
    
    public R2(int dimensions) {
        this.dimensions = dimensions;
        current = new double[dimensions];

        for(int i=0;i<dimensions;i++) {
            current[i] = 0.5; // seed
        } 
        
        phis = phi1_t;
        switch (dimensions) {
          case 2: phis = phi2_t; break;
          case 3: phis = phi3_t; break;
          case 4: phis = phi4_t; break;
        }
    }

    public double[] nextVector() {
        for(int i=0;i<dimensions;i++) {
            current[i] = (current[i] + phis[i])%1.0;
        }
        return current;
    }
}
