package fastmath.java.dist;

import org.apache.commons.math3.distribution.LogisticDistribution;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Logistic distribution with the variance corrected to pi^2 s^2 / 3
 * (Commons Math 3.6.1 returns pi^2 / (3 s^2)). Everything else is inherited.
 */
public class LogisticFixed extends LogisticDistribution {
    private static final long serialVersionUID = 1L;

    public LogisticFixed(double mu, double s) {
        super(mu, s);
    }

    public LogisticFixed(RandomGenerator rng, double mu, double s) {
        super(rng, mu, s);
    }

    @Override
    public double getNumericalVariance() {
        final double s = getScale();
        return Math.PI * Math.PI * s * s / 3.0;
    }
}
