package fastmath.java;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;

/**
 * Thread-safe wrapper of a Commons Math {@link RandomGenerator} which exposes the wrapped generator.
 */
public class SyncedRandomGenerator extends SynchronizedRandomGenerator {
    private final RandomGenerator delegate;

    public SyncedRandomGenerator(RandomGenerator rng) {
        super(rng);
        this.delegate = rng;
    }

    /** Returns the wrapped (not synchronised) generator. */
    public RandomGenerator getDelegate() {
        return delegate;
    }
}
