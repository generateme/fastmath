package fastmath.java.dist;

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * F distribution with density and log-density defined for x &lt;= 0:
 * 0 / -Inf for x &lt; 0 and the limit for x = 0 (+Inf for d1 &lt; 2, 1 for d1 = 2, 0 for d1 &gt; 2).
 * Commons Math 3.6.1 returns NaN there. The cumulative probability is already 0 for x &lt;= 0.
 */
public class FFixed extends FDistribution {
    private static final long serialVersionUID = 1L;

    public FFixed(double numeratorDegreesOfFreedom, double denominatorDegreesOfFreedom)
        throws NotStrictlyPositiveException {
        super(numeratorDegreesOfFreedom, denominatorDegreesOfFreedom);
    }

    public FFixed(double numeratorDegreesOfFreedom, double denominatorDegreesOfFreedom,
                  double inverseCumAccuracy) throws NotStrictlyPositiveException {
        super(numeratorDegreesOfFreedom, denominatorDegreesOfFreedom, inverseCumAccuracy);
    }

    public FFixed(RandomGenerator rng, double numeratorDegreesOfFreedom,
                  double denominatorDegreesOfFreedom) throws NotStrictlyPositiveException {
        super(rng, numeratorDegreesOfFreedom, denominatorDegreesOfFreedom);
    }

    public FFixed(RandomGenerator rng, double numeratorDegreesOfFreedom,
                  double denominatorDegreesOfFreedom, double inverseCumAccuracy)
        throws NotStrictlyPositiveException {
        super(rng, numeratorDegreesOfFreedom, denominatorDegreesOfFreedom, inverseCumAccuracy);
    }

    @Override
    public double density(double x) {
        if (x < 0.0) return 0.0;
        if (x == 0.0) {
            final double d1 = getNumeratorDegreesOfFreedom();
            return d1 < 2.0 ? Double.POSITIVE_INFINITY : (d1 == 2.0 ? 1.0 : 0.0);
        }
        return super.density(x);
    }

    @Override
    public double logDensity(double x) {
        if (x < 0.0) return Double.NEGATIVE_INFINITY;
        if (x == 0.0) {
            final double d1 = getNumeratorDegreesOfFreedom();
            return d1 < 2.0 ? Double.POSITIVE_INFINITY : (d1 == 2.0 ? 0.0 : Double.NEGATIVE_INFINITY);
        }
        return super.logDensity(x);
    }
}
