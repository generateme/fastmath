package fastmath.java.noise;

public final class RidgedMulti {
    private static double ridged(double v, double weight) {
        double a = 1.0-Math.abs(v);
        return a * a * weight;
    }
    
    public static double noise(NoiseConfig cfg, double x) {
        double sum = 0.0;
        double weight = 1.0;

        double xx = x;

        for(int i=0; i < cfg.octaves; i++) {
            double signal = ridged(Noise.value(cfg, cfg.perm[i], xx), weight);
            sum += signal * cfg.spectralWeights[i];

            weight = signal / cfg.gain;
            weight = weight > 1.0 ? 1.0 : (weight < 0.0 ? 0.0 : weight);
            xx *= cfg.lacunarity;
        }

        return cfg.normalize ? sum * cfg.fractalBounding : (sum * cfg.fractalBounding * 2.0 - 1.0);
    }

    public static double noise(NoiseConfig cfg, double x, double y) {
        double sum = 0.0;
        double weight = 1.0;
        
        double xx = x;
        double yy = y;

        for(int i=0; i < cfg.octaves; i++) {
            double signal = ridged(Noise.value(cfg, cfg.perm[i], xx, yy), weight);
            sum += signal * cfg.spectralWeights[i];
            
            weight = signal / cfg.gain;
            weight = weight > 1.0 ? 1.0 : (weight < 0.0 ? 0.0 : weight);
            xx *= cfg.lacunarity;
            yy *= cfg.lacunarity;
        }

        return cfg.normalize ? sum * cfg.fractalBounding : (sum * cfg.fractalBounding * 2.0 - 1.0);
    }

    public static double noise(NoiseConfig cfg, double x, double y, double z) {
        double sum = 0.0;
        double weight = 1.0;
        
        double xx = x;
        double yy = y;
        double zz = z;

        for(int i = 0; i < cfg.octaves; i++) {
            double signal = ridged(Noise.value(cfg, cfg.perm[i], xx, yy, zz), weight);
            sum += signal * cfg.spectralWeights[i];

            weight = signal / cfg.gain;
            weight = weight > 1.0 ? 1.0 : (weight < 0.0 ? 0.0 : weight);
            xx *= cfg.lacunarity;
            yy *= cfg.lacunarity;
            zz *= cfg.lacunarity;
        }

        return cfg.normalize ? sum * cfg.fractalBounding : (sum * cfg.fractalBounding * 2.0 - 1.0);
    }
 
}
