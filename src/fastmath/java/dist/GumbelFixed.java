package fastmath.java.dist;

import org.apache.commons.math3.distribution.GumbelDistribution;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Gumbel distribution with the mean corrected to mu + beta * gamma, gamma = 0.5772156649015329
 * (Commons Math 3.6.1 uses 0.5778636748954609). Everything else is inherited.
 */
public class GumbelFixed extends GumbelDistribution {
    private static final long serialVersionUID = 1L;

    private static final double EULER_GAMMA = 0.5772156649015329;

    public GumbelFixed(double mu, double beta) {
        super(mu, beta);
    }

    public GumbelFixed(RandomGenerator rng, double mu, double beta) {
        super(rng, mu, beta);
    }

    @Override
    public double getNumericalMean() {
        return getLocation() + getScale() * EULER_GAMMA;
    }
}
