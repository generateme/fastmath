package fastmath.java.dist;

import org.apache.commons.math3.distribution.NakagamiDistribution;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Nakagami distribution with cumulative probability equal to 0 for x &lt;= 0
 * (Commons Math 3.6.1 squares x and returns a positive value for x &lt; 0).
 * Density is 0 for x &lt;= 0 in the parent already.
 */
public class NakagamiFixed extends NakagamiDistribution {
    private static final long serialVersionUID = 1L;

    public NakagamiFixed(double mu, double omega) {
        super(mu, omega);
    }

    public NakagamiFixed(double mu, double omega, double inverseAbsoluteAccuracy) {
        super(mu, omega, inverseAbsoluteAccuracy);
    }

    public NakagamiFixed(RandomGenerator rng, double mu, double omega, double inverseAbsoluteAccuracy) {
        super(rng, mu, omega, inverseAbsoluteAccuracy);
    }

    @Override
    public double cumulativeProbability(double x) {
        return x <= 0.0 ? 0.0 : super.cumulativeProbability(x);
    }
}
