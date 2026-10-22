package fastmath.java.noise;

public final class Billow {
    private static double billow(double v) {
        return Math.abs(v) * 2.0 - 1.0;
    }
    
    public static double noise(NoiseConfig cfg, double x) {
        if (!(Math.abs(x) <= cfg.maxCoord)) return Double.NaN;
        double sum = billow(Noise.value(cfg, cfg.perm[0], x));
        double amp = 1.0;

        double xx = x;

        int i=0;
        while(++i < cfg.octaves) {
            xx *= cfg.lacunarity;
            
            amp *= cfg.gain;
            sum += billow(Noise.value(cfg, cfg.perm[i], xx)) * amp;
        }

        return cfg.normalize ? ((sum * cfg.fractalBounding) + 1.0) * 0.5 : sum * cfg.fractalBounding;
    }

    public static double noise(NoiseConfig cfg, double x, double y) {
        if (!(Math.abs(x) <= cfg.maxCoord && Math.abs(y) <= cfg.maxCoord)) return Double.NaN;
        double sum = billow(Noise.value(cfg, cfg.perm[0], x, y));
        double amp = 1.0;

        double xx = x;
        double yy = y;

        int i=0;
        while(++i < cfg.octaves) {
            xx *= cfg.lacunarity;
            yy *= cfg.lacunarity;
            
            amp *= cfg.gain;
            sum += billow(Noise.value(cfg, cfg.perm[i], xx, yy)) * amp;
        }

        return cfg.normalize ? ((sum * cfg.fractalBounding) + 1.0) * 0.5 : sum * cfg.fractalBounding;
    }

    public static double noise(NoiseConfig cfg, double x, double y, double z) {
        if (!(Math.abs(x) <= cfg.maxCoord && Math.abs(y) <= cfg.maxCoord && Math.abs(z) <= cfg.maxCoord)) return Double.NaN;
        double sum = billow(Noise.value(cfg, cfg.perm[0], x, y, z));
        double amp = 1.0;

        double xx = x;
        double yy = y;
        double zz = z;

        int i=0;
        while(++i < cfg.octaves) {
            xx *= cfg.lacunarity;
            yy *= cfg.lacunarity;
            zz *= cfg.lacunarity;
            
            amp *= cfg.gain;
            sum += billow(Noise.value(cfg, cfg.perm[i], xx, yy, zz)) * amp;
        }
 
        return cfg.normalize ? ((sum * cfg.fractalBounding) + 1.0) * 0.5 : sum * cfg.fractalBounding;
    }

}
