package fastmath.java.noise;

public final class Discrete {

    private static final double AM = 1.0 / 2147483647;

    public static double value (int X, int Y) {
        int n = X + Y * 57;
        n = (n << 13) ^ n;
        return ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) * AM;
    }
}
