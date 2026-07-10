package sthree_lra_xnes;

import lib.matrix.TCMatrix;

/**
 * Represents a single solution in a real-valued search space.
 * This class holds the solution's representation in both the problem space (x-space)
 * and the search space (z-space), along with its evaluation value (fitness).
 */
public class TSRealSolution {

    /** The evaluation value or fitness of the solution. */
    private double fEvaluationValue;

    /** The solution vector in the z-space. */
    private TCMatrix fZVector;

    /** The solution vector in the x-space. */
    private TCMatrix fXVector;

    /**
     * Default constructor. Creates a 0-dimensional solution.
     */
    public TSRealSolution() {
        this(0);
    }

    /**
     * Constructs a solution for a given problem dimension.
     * The vectors are initialized but the evaluation value is set to NaN.
     *
     * @param dimension The dimension of the problem space.
     */
    public TSRealSolution(int dimension) {
        fEvaluationValue = Double.NaN;
        fXVector = new TCMatrix(dimension);
        fZVector = new TCMatrix(dimension);
    }

    /**
     * Copy constructor. Creates a deep copy of the source solution.
     *
     * @param src The source solution to copy.
     */
    public TSRealSolution(TSRealSolution src) {
        fEvaluationValue = src.fEvaluationValue;
        fXVector = new TCMatrix(src.fXVector);
        fZVector = new TCMatrix(src.fZVector);
    }

    /**
     * Creates and returns a deep copy of this solution.
     *
     * @return A new {@code TSRealSolution} instance with the same data.
     */
    @Override
    public TSRealSolution clone() {
        return new TSRealSolution(this);
    }

    /**
     * Copies the data from a source solution into this instance.
     *
     * @param src The source solution to copy from.
     * @return This instance after copying.
     */
    public TSRealSolution copyFrom(TSRealSolution src) {
        fEvaluationValue = src.fEvaluationValue;
        fXVector.copyFrom(src.fXVector);
        fZVector.copyFrom(src.fZVector);
        return this;
    }

    /**
     * Gets the evaluation value of the solution.
     *
     * @return The evaluation value.
     */
    public double getEvaluationValue() {
        return fEvaluationValue;
    }

    /**
     * Sets the evaluation value of the solution.
     *
     * @param value The new evaluation value.
     */
    public void setEvaluationValue(double value) {
        fEvaluationValue = value;
    }

    /**
     * Gets the solution vector in the z-space.
     *
     * @return The z-space vector as a {@code TCMatrix}.
     */
    public TCMatrix getZVector() {
        return fZVector;
    }

    /**
     * Gets the solution vector in the x-space (the primary solution vector).
     *
     * @return The x-space vector as a {@code TCMatrix}.
     */
    public TCMatrix getVector() {
        return fXVector;
    }

    /**
     * Returns a string representation of the solution.
     * Includes the evaluation value, x-vector, and z-vector.
     *
     * @return A string representation of this object.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("f(x)=").append(fEvaluationValue);
        sb.append(", x=").append(fXVector.tclone());
        sb.append(", z=").append(fZVector.tclone());
        return sb.toString();
    }
}
