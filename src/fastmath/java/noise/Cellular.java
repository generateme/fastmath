package fastmath.java.noise;

public final class Cellular {

    private static int FastRound(double f) {
        return (f >= 0) ? (int)(f + 0.5) : (int)(f - 0.5);
    }

    private static int lut(NoiseConfig cfg, int offset, int x, int y) {
	return cfg.perm[(x & 0xff) + cfg.perm[(y & 0xff) + offset]];
    }

    private static double cellLut(NoiseConfig cfg, int offset, int x, int y) {
        return cfg.cell2dX[NoiseConfig.hash(cfg, offset, x, y) & 0xff];
    }
    
    public static double value(NoiseConfig cfg, int offset, double x, double y) {
        int xr = FastRound(x);
        int yr = FastRound(y);

        double distance = Double.MAX_VALUE;

        int xc = xr;
        int yc = yr;

        for (int xi = xr - 1; xi <= xr + 1; xi++) {
            for (int yi = yr - 1; yi <= yr + 1; yi++) {
                int lutPos = lut(cfg, offset, xi, yi);

                double vecX = xi - x + cfg.cell2dX[lutPos] * 0.45;
                double vecY = yi - y + cfg.cell2dY[lutPos] * 0.45;
                
                double newDistance = vecX * vecX + vecY * vecY;

                if (newDistance < distance) {
                    distance = newDistance;
                    xc = xi;
                    yc = yi;
                }
            }
        }

        return cellLut(cfg, offset, xc, yc);
    }
    
}
