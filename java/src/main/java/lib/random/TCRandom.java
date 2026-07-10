package lib.random;

import java.util.Random;

public class TCRandom extends TCAbstractRandom {

    /** For serialization */
    private static final long serialVersionUID = 1L;

    private Random fRandom;

    public TCRandom() {
        fRandom = new Random();
    }

    public TCRandom(long seed) {
        fRandom = new Random(seed);
    }

    protected int next(int bits) {
        assert 0 < bits && bits <= 32;
        return fRandom.nextInt() >>> (32 - bits);
    }

    protected void resetSeed(long seed) {
        fRandom.setSeed(seed);
    }

    public Random getRandom() {
        return fRandom;
    }

}
