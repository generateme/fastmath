package fastmath.java.dist;

import org.apache.commons.math3.distribution.LevyDistribution;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Levy distribution with density / log-density / cumulative probability equal to 0 / -Inf / 0
 * for x &lt;= mu (Commons Math 3.6.1 returns NaN for x &lt; mu and for the density at x = mu).
 */
public class LevyFixed extends LevyDistribution {
    private static final long serialVersionUID = 1L;

    public LevyFixed(double mu, double c) {
        super(mu, c);
    }

    public LevyFixed(RandomGenerator rng, double mu, double c) {
        super(rng, mu, c);
    }

    @Override
    public double density(double x) {
        return x <= getLocation() ? 0.0 : super.density(x);
    }

    @Override
    public double logDensity(double x) {
        return x <= getLocation() ? Double.NEGATIVE_INFINITY : super.logDensity(x);
    }

    @Override
    public double cumulativeProbability(double x) {
        return x <= getLocation() ? 0.0 : super.cumulativeProbability(x);
    }
}
