package sthree_lra_xnes;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import lib.matrix.TCMatrix;
import lib.matrix.decompositions.TCEigenvalueDecomposition;
import lib.random.ICRandom;

public class S3LRAxNES {
    private static double defaultAlphaSigma(int d) {
        double power = Math.log(2.2) / Math.log(8.0);
        return 3 * Math.pow(d / 10.0, power);
    }

    private static double defaultAlphaB(int d) {
        return Math.min(2.0, 1.0 + 16.0 / (double) d);
    }

    private static double defaultRatioEtaSigmaMin(int d) {
        double power = Math.log(5.5) / Math.log(8.0);
        return 0.1 * Math.pow((double) d / 10.0, power);
    }

    private static int defaultPopulationSize(int d) {
        return 4 + (int) Math.floor(3 * Math.log(d));
    }

    private static double defaultEtaSigmaAndBOfxNES(int dimension) {
        return 0.6 * (3. + Math.log(dimension)) / (dimension * Math.sqrt(dimension));
    }

    /** problem dimension */
    private final int fDimension;

    /** population size */
    private final int fPopulationSize;

    /** mean vector */
    private TCMatrix fMean;

    /** step-size */
    private double fSigma;

    /** shape */
    private TCMatrix fB;

    /** population of candidate solutions for the current generation */
    private List<TSRealSolution> fPopulation;

    /** mean vector learning rate */
    private double fEtaMean;

    /** step-size learning rate */
    private double fEtaSigma;

    /** shape learning rate */
    private double fEtaB;

    /** weight function */
    private final double[] fWeights;

    /** random number generator */
    private ICRandom fRandom;

    /** Comparator */
    private Comparator<TSRealSolution> fComparator;

    /** natural gradient of mean */
    private TCMatrix fGMean;

    /** natural gradient of covariance matrix $\sigma^2BB^\top$ */
    private TCMatrix fGM;

    /** natural gradient of shape B */
    private TCMatrix fGB;

    /** natural gradient of step-size $\sigma$ */
    private double fGSigma;

    /** $wzz^T$ */
    private TCMatrix fWZMatrix;

    /** identity matrix */
    private TCMatrix fEye;

    /** $\sum_i w_i*I$ */
    private TCMatrix fWEye;

    /** scaling coefficients that define thresholds for directional movement */
    private final double fAlphaSigma, fAlphaMean, fAlphaB;

    /** change rates of learning rates */
    private final double fBetaSigmaForIncrease, fBetaSigmaForDecrease, fBetaMean, fBetaB;

    /** exponential moving average weights */
    private final double fExponentialMovingAverageWeightSigma, fExponentialMovingAverageWeightMean,
            fExponentialMovingAverageWeightB;

    /** minimum and maximum values of learning rates */
    private final double fEtaMeanMin, fEtaMeanMax, fEtaBMin, fEtaBMax, fEtaSigmaMin, fEtaSigmaMax;

    /** $\mu_w$ */
    private final double fMuW;

    /** expected value of features under random function */
    private double fGammaSigma, fGammaMean, fGammaB;

    /** p_B and p_m */
    private TCMatrix fEvolutionPathB, fEvolutionPathMean;

    /** temporary matrices */
    private TCMatrix fMTmp, fVecTmp, fVTmp;

    /** temporary matrices */
    private TCMatrix fBBTBeforeUpdate, fBBeforeUpdate, fInvSqrtBBT, fNormalizedNewBBT,
            fSubFromNormalizedNewBBTToEye;

    /**
     * approximated value of expected Euclidean norm of the sample of a standard normal
     * distribution, \epsilon_d.
     */
    private final double fSqrtChiN;

    /** directional movement features */
    private double fDirectionalMovementFeatureOfMean, fDirectionalMovementFeatureOfB,
            fDirectionalMovementFeatureOfStepsize;

    /** object for performing eigenvalue decomposition */
    private TCEigenvalueDecomposition fRes;

    /**
     * Constructs an S3LRAxNES instance with the default hyperparameter settings and initializes its parameters.
     *
     * @param dimension The dimension of the problem.
     * @param mean The initial mean vector.
     * @param sigma The initial step-size.
     * @param random The random number generator.
     */
    public S3LRAxNES(int dimension, TCMatrix mean, double sigma, ICRandom random) {
        this(dimension, defaultPopulationSize(dimension), mean, sigma, random, //
                2.0, // alphaMean
                defaultAlphaSigma(dimension), // alphaSigma
                defaultAlphaB(dimension), // alphaB,
                0.08, // betaMean,
                0.02, // betaSigmaForIncrease,
                0.1, // betaSigmaForDecrease,
                0.1, // betaB,
                0.01, // exponentialMovingAverageWeightMean,
                0.01, // exponentialMovingAverageWeightSigma,
                0.01, // exponentialMovingAverageWeightB,
                1.0, // initEtaMean,
                defaultEtaSigmaAndBOfxNES(dimension), // initEtaSigma,
                defaultEtaSigmaAndBOfxNES(dimension), // initEtaB,
                0.01, // lowerEtaMean,
                defaultRatioEtaSigmaMin(dimension) * defaultEtaSigmaAndBOfxNES(dimension), // lowerEtaSigma,
                0.2 * defaultEtaSigmaAndBOfxNES(dimension), // lowerEtaB,
                1.0, // upperEtaMean,
                1.0, // upperEtaSigma,
                1.5 * defaultEtaSigmaAndBOfxNES(dimension)// upperEtaB
        );
    }

    /**
     * Constructs an S3LRAxNES instance and initializes its parameters.
     *
     * @param dimension The dimension of the problem.
     * @param populationSize The number of individuals in the population.
     * @param mean The initial mean vector.
     * @param sigma The initial step-size.
     * @param random The random number generator.
     * @param alphaMean The scaling coefficient for the mean's directional movement.
     * @param alphaSigma The scaling coefficient for the step-size's directional movement.
     * @param alphaB The scaling coefficient for the shape's directional movement.
     * @param betaMean The change rate for the mean's learning rate.
     * @param betaSigmaForIncrease The change rate for increasing the step-size's learning rate.
     * @param betaSigmaForDecrease The change rate for decreasing the step-size's learning rate.
     * @param betaB The change rate for the shape's learning rate.
     * @param exponentialMovingAverageWeightMean The exponential moving average weight for the
     *        mean's evolution path p_m.
     * @param exponentialMovingAverageWeightSigma The exponential moving average weight for the
     *        step-size's directional movement feature l_sigma.
     * @param exponentialMovingAverageWeightB The exponential moving average weight for the shape's
     *        evolution path p_B.
     * @param initEtaMean The initial learning rate for the mean.
     * @param initEtaSigma The initial learning rate for the step-size.
     * @param initEtaB The initial learning rate for the shape.
     * @param lowerEtaMean The lower bound for the mean's learning rate.
     * @param lowerEtaSigma The lower bound for the step-size's learning rate.
     * @param lowerEtaB The lower bound for the shape's learning rate.
     * @param upperEtaMean The upper bound for the mean's learning rate.
     * @param upperEtaSigma The upper bound for the step-size's learning rate.
     * @param upperEtaB The upper bound for the shape's learning rate.
     */
    private S3LRAxNES(int dimension, int populationSize, TCMatrix mean, double sigma,
            ICRandom random, double alphaMean, double alphaSigma, double alphaB, double betaMean,
            double betaSigmaForIncrease, double betaSigmaForDecrease, double betaB,
            double exponentialMovingAverageWeightMean, double exponentialMovingAverageWeightSigma,
            double exponentialMovingAverageWeightB, double initEtaMean, double initEtaSigma,
            double initEtaB, double lowerEtaMean, double lowerEtaSigma, double lowerEtaB,
            double upperEtaMean, double upperEtaSigma, double upperEtaB) {
        fDimension = dimension;
        fPopulationSize = populationSize;
        fSigma = sigma;

        fAlphaMean = alphaMean;
        fAlphaSigma = alphaSigma;
        fAlphaB = alphaB;

        fBetaMean = betaMean;
        fBetaSigmaForIncrease = betaSigmaForIncrease;
        fBetaSigmaForDecrease = betaSigmaForDecrease;
        fBetaB = betaB;

        fExponentialMovingAverageWeightMean = exponentialMovingAverageWeightMean;
        fExponentialMovingAverageWeightSigma = exponentialMovingAverageWeightSigma;
        fExponentialMovingAverageWeightB = exponentialMovingAverageWeightB;

        fEtaMean = initEtaMean;
        fEtaSigma = initEtaSigma;
        fEtaB = initEtaB;

        fEtaMeanMin = lowerEtaMean;
        fEtaSigmaMin = lowerEtaSigma;
        fEtaBMin = lowerEtaB;

        fEtaMeanMax = upperEtaMean;
        fEtaSigmaMax = upperEtaSigma;
        fEtaBMax = upperEtaB;
        fRandom = random;
        fWeights = new double[fPopulationSize];
        double[] wTmp = new double[fPopulationSize];
        double wSum = 0.;
        for (int i = 0; i < fPopulationSize; ++i) {
            wTmp[i] = Math.max(0., Math.log((fPopulationSize / 2.) + 1.) - Math.log(i + 1.0));
            wSum += wTmp[i];
        }
        for (int i = 0; i < fPopulationSize; ++i) {
            fWeights[i] = (wTmp[i] / wSum) - (1. / fPopulationSize);
        }
        wSum = 0.;
        for (int i = 0; i < fPopulationSize; ++i) {
            wSum += fWeights[i];
        }

        fEye = new TCMatrix(fDimension, fDimension).eye();
        fMean = mean.clone();
        fB = new TCMatrix(fDimension, fDimension).eye();
        fComparator = new TSEvaluationValueComparator(true);
        fPopulation = new ArrayList<>();
        for (int i = 0; i < fPopulationSize; ++i) {
            fPopulation.add(new TSRealSolution(fDimension));
        }

        fGMean = new TCMatrix(fDimension);
        fGM = new TCMatrix(fDimension, fDimension);
        fGB = new TCMatrix(fDimension, fDimension);
        fWEye = new TCMatrix(fDimension, fDimension).eye().times(wSum);
        fWZMatrix = new TCMatrix(fDimension, fDimension);
        fMTmp = new TCMatrix(fDimension, fDimension);
        fVecTmp = new TCMatrix(fDimension);
        fVTmp = new TCMatrix(fDimension, fDimension);
        fGSigma = 0.0;

        double wSquareSum = 0.0;
        for (int i = 0; i < fPopulationSize; i++) {
            wSquareSum += fWeights[i] * fWeights[i];
        }
        fMuW = 1. / wSquareSum;

        fSqrtChiN = Math.sqrt(fDimension)
                * (1.0 - 1.0 / (4 * fDimension) + 1.0 / (21.0 * fDimension * fDimension));
        fEvolutionPathB = new TCMatrix(fDimension, fDimension);
        fEvolutionPathMean = new TCMatrix(fDimension);
        fInvSqrtBBT = new TCMatrix(fDimension, fDimension);
        fBBTBeforeUpdate = new TCMatrix(fDimension, fDimension);
        fNormalizedNewBBT = new TCMatrix(fDimension, fDimension);
        fBBeforeUpdate = new TCMatrix(fDimension, fDimension);
        fDirectionalMovementFeatureOfMean = 0;
        fDirectionalMovementFeatureOfB = 0;
        fDirectionalMovementFeatureOfStepsize = 0;
        fRes = null;
        fSubFromNormalizedNewBBTToEye = new TCMatrix(fDimension, fDimension);

        fGammaMean = 0.0;
        fGammaB = 0.0;
        fGammaSigma = 0.0;
    }

    /**
     * Gets the population size
     *
     * @return The population size.
     */
    public int getPopulationSize() {
        return fPopulationSize;
    }

    /**
     * Gets the random number generator.
     *
     * @return The random number generator.
     */
    public ICRandom getRandom() {
        return fRandom;
    }

    /**
     * Gets the current mean vector of the search distribution.
     *
     * @return The mean vector.
     */
    public TCMatrix getMeanVector() {
        return fMean;
    }

    /**
     * Gets the current shape B.
     *
     * @return The shape B.
     */
    public TCMatrix getB() {
        return fB;
    }

    /**
     * Gets the current step-size sigma.
     *
     * @return The step-size sigma.
     */
    public double getSigma() {
        return fSigma;
    }

    /**
     * Calculates and returns the eigenvalues of the full covariance matrix Sigma = sigma^2 * B *
     * B^T.
     *
     * @return An array of eigenvalues of the covariance matrix.
     */
    public double[] calculateEigenvaluesOfCovariance() {
        TCMatrix cov = fB.clone().times(fB.tclone()).times(fSigma * fSigma);
        return cov.eig().getRealEigenvalues();
    }

    /**
     * Calculates and returns the eigenvalues of the matrix B*B^T.
     *
     * @return An array of eigenvalues of B*B^T.
     */
    public double[] calculateEigenvaluesOfBBT() {
        TCMatrix bbt = fB.clone().times(fB.tclone());
        return bbt.eig().getRealEigenvalues();
    }

    /**
     * Samples a new candidate solution from the current distribution N(mean, sigma^2 * B * B^T).
     * The new solution vector is stored in the i-th member of the population.
     *
     * @param i The index of the candidate solution to resample.
     */
    public TSRealSolution sampleCandidateSolution(int i) {
        // Equation (1)
        TSRealSolution s = fPopulation.get(i);
        s.getZVector().randn(fRandom);
        s.getVector().times(fB, s.getZVector()).times(fSigma).add(fMean);
        return s;
    }

    /**
     * Gets the i-th solution from the population.
     *
     * @param i The index of the solution.
     * @return The i-th TSRealSolution object.
     */
    public TSRealSolution getTSRealSolution(int i) {
        return fPopulation.get(i);
    }

    /**
     * Gets the current learning rate for the step-size.
     *
     * @return The step-size learning rate.
     */
    public double getEtaSigma() {
        return fEtaSigma;
    }

    /**
     * Gets the current learning rate for the shape.
     *
     * @return The shape learning rate.
     */
    public double getEtaB() {
        return fEtaB;
    }

    /**
     * Gets the current learning rate for the mean vector.
     *
     * @return The mean vector learning rate.
     */
    public double getEtaMean() {
        return fEtaMean;
    }

    /**
     * Executes one generation of the S3LRA-xNES algorithm.
     * This involves calculating gradients, updating parameters, updating evolution paths,
     * and adapting the learning rates.
     *
     * @throws Exception if the covariance matrix has one or more negative eigenvalues,
     *         indicating it is not positive semi-definite.
     */
    public void nextGeneration() throws Exception {

        fBBTBeforeUpdate.times(fB, fMTmp.tcopyFrom(fB));
        fBBeforeUpdate.copyFrom(fB);

        calculateSqrtInvBBT();

        calculateNaturalGradients();
        updateParameters();

        calculateNormalizedNewBBT();
        calculateSubMatrix();

        updateEvolutionPathB();
        updateEvolutionPathMean();

        calculateDirectionalMovementFeatures();
        updateGammas();
        updateEtas();
    }

    /**
     * Gets the directional movement feature for the mean vector.
     *
     * @return The directional movement feature of the mean.
     */
    public double getDirectionalMovementFeatureOfMean() {
        return fDirectionalMovementFeatureOfMean;
    }

    /**
     * Gets the directional movement feature for the step-size.
     *
     * @return The directional movement feature of the step-size.
     */
    public double getDirectionalMovementFeatureOfStepsize() {
        return fDirectionalMovementFeatureOfStepsize;
    }

    /**
     * Gets the directional movement feature for the shape B.
     *
     * @return The directional movement feature of B.
     */
    public double getDirectionalMovementFeatureOfB() {
        return fDirectionalMovementFeatureOfB;
    }

    /**
     * Sets the step-size sigma.
     *
     * @param sigma The new step-size.
     */
    public void setSigma(double sigma) {
        fSigma = sigma;
    }

    /**
     * Gets the expected directional movement feature for the mean.
     *
     * @return The gamma value for the mean.
     */
    public double getGammaMean() {
        return fGammaMean;
    }

    /**
     * Gets the expected directional movement feature for the shape B.
     *
     * @return The gamma value for B.
     */
    public double getGammaB() {
        return fGammaB;
    }

    /**
     * Gets the expected directional movement feature for the step-size.
     *
     * @return The gamma value for the step-size.
     */
    public double getGammaStep() {
        return fGammaSigma;
    }

    /**
     * Gets the sum of squared weights (mu_w).
     *
     * @return The value of mu_w.
     */
    public double getMuW() {
        return fMuW;
    }

    /**
     * Gets the natural gradient for the step-size sigma.
     *
     * @return The natural gradient g_sigma.
     */
    public double getGSigma() {
        return fGSigma;
    }

    /**
     * Copies the full covariance matrix Sigma = sigma^2 * B * B^T into a destination matrix.
     *
     * @param dest The destination matrix.
     */
    public void copyCov(TCMatrix dest) {
        dest.times(fB, fMTmp.tcopyFrom(fB)).times(fSigma * fSigma);
    }

    /**
     * Copies the matrix B*B^T into a destination matrix.
     *
     * @param dest The destination matrix.
     */
    public void copyBBT(TCMatrix dest) {
        dest.times(fB, fMTmp.tcopyFrom(fB));
    }

    /**
     * Gets the natural gradient for the mean vector.
     *
     * @return The natural gradient for the mean.
     */
    public TCMatrix getGMean() {
        return fGMean;
    }

    /**
     * Gets the result of the subtraction of the identity matrix from the normalized new B*B^T.
     *
     * @return The resulting matrix.
     */
    public TCMatrix getSubFromNormalizedNewBBTToEye() {
        return fSubFromNormalizedNewBBTToEye;
    }

    /**
     * Gets the evolution path for the shape B.
     *
     * @return The evolution path p_B.
     */
    public TCMatrix getEvolutionPathB() {
        return fEvolutionPathB;
    }

    /**
     * Calculates the directional movement features for the mean, step-size, and shape.
     */
    private void calculateDirectionalMovementFeatures() {
        // Equation (14)
        fDirectionalMovementFeatureOfB = calcLengthOfMatrixEvolutionPath(fEvolutionPathB);
        // Equation (9), with the approximation in (13).
        // Equation(10)--(12) are calculated theoretically, as well as Equation (33) in the
        // supplementary document.
        fDirectionalMovementFeatureOfStepsize =
                (1 - fExponentialMovingAverageWeightSigma)
                        * (1 - fExponentialMovingAverageWeightSigma)
                        * fDirectionalMovementFeatureOfStepsize
                        + fExponentialMovingAverageWeightSigma
                                * (2.0 - fExponentialMovingAverageWeightSigma)
                                * (Math.exp(2 * fEtaSigma * fGSigma)
                                        - 2 * Math.exp(fEtaSigma * fGSigma) + 1.0)
                                * fDimension / 2.0 * fMuW / (fEtaSigma * fEtaSigma);
        // Equation (21)
        fDirectionalMovementFeatureOfMean = fEvolutionPathMean.normL2();
    }

    /**
     * Updates the expected values (gamma) of the directional movement features.
     */
    private void updateGammas() {
        // Equations (25)--(27)
        fGammaSigma = (1. - fExponentialMovingAverageWeightSigma)
                * (1. - fExponentialMovingAverageWeightSigma) * fGammaSigma
                + fExponentialMovingAverageWeightSigma
                        * (2.0 - fExponentialMovingAverageWeightSigma);
        fGammaMean = (1. - fExponentialMovingAverageWeightMean)
                * (1. - fExponentialMovingAverageWeightMean) * fGammaMean
                + fExponentialMovingAverageWeightMean * (2.0 - fExponentialMovingAverageWeightMean);
        fGammaB = (1. - fExponentialMovingAverageWeightB) * (1. - fExponentialMovingAverageWeightB)
                * fGammaB
                + fExponentialMovingAverageWeightB * (2.0 - fExponentialMovingAverageWeightB);
    }

    /**
     * Calculates the "length" of a matrix evolution path, defined as trace(mat * mat) / 2.
     *
     * @param mat The matrix evolution path.
     * @return The calculated length.
     */
    private double calcLengthOfMatrixEvolutionPath(TCMatrix mat) {
        return fMTmp.times(mat, mat).trace() / 2.0;
    }

    /**
     * Adapts the learning rates for the mean, step-size, and shape based on the
     * comparison between directional movement features divided by corresponding scaling coefficient
     * and directional movement feature's expected value.
     */
    private void updateEtas() {
        // Equations (28)
        fEtaMean *= Math.exp(fBetaMean
                * (fDirectionalMovementFeatureOfMean / fAlphaMean - Math.sqrt(fGammaMean)));
        fEtaMean = clip(fEtaMean, fEtaMeanMin, fEtaMeanMax);

        // Equations (29)
        double pow = fDirectionalMovementFeatureOfStepsize / fAlphaSigma - fGammaSigma;
        // Change-rate switching mechanism for step-size learning rate, equation (31)
        if (pow > 0) {
            fEtaSigma *= Math.exp(fBetaSigmaForIncrease * pow);
        } else {
            fEtaSigma *= Math.exp(fBetaSigmaForDecrease * pow);
        }
        fEtaSigma = clip(fEtaSigma, fEtaSigmaMin, fEtaSigmaMax);

        // Equations (30)
        fEtaB *= Math.exp(fBetaB * (fDirectionalMovementFeatureOfB / fAlphaB - fGammaB));
        fEtaB = clip(fEtaB, fEtaBMin, fEtaBMax);

    }

    /**
     * Clips a value to be within a specified minimum and maximum range.
     *
     * @param rawValue The value to clip.
     * @param min The minimum allowed value.
     * @param max The maximum allowed value.
     * @return The clipped value.
     */
    private double clip(double rawValue, double min, double max) {
        return Math.min(Math.max(rawValue, min), max);
    }

    /**
     * Updates the evolution path for the mean vector (p_m).
     */
    private void updateEvolutionPathMean() {
        // Equation (22) with the approximation in equation (24)
        // Equation (23) is calculated theoretically as the same process as the equation (54) to
        // equation (55) in the supplementary document.
        fEvolutionPathMean.times(1.0 - fExponentialMovingAverageWeightMean)
                .add(fVecTmp
                        .times(fGMean,
                                Math.sqrt(fExponentialMovingAverageWeightMean
                                        * (2.0 - fExponentialMovingAverageWeightMean) * fMuW)
                                        / fSqrtChiN));
    }

    /**
     * Calculates the difference between the normalized new B*B^T and the identity matrix.
     *
     * @requires fNormalizedNewBBT
     */
    private void calculateSubMatrix() {
        fSubFromNormalizedNewBBTToEye.sub(fNormalizedNewBBT, fEye);
    }

    /**
     * Updates the evolution path for the shape B (p_B).
     *
     * @requires fSubFromNormalizedNewBBTToEye
     */
    private void updateEvolutionPathB() {
        // Equation (18)
        double approxCovB =
                (fEtaB * fEtaB / 2.0 * ((long) fDimension * (long) fDimension + fDimension - 2.0))
                        / fMuW;
        // Equation (15)
        fEvolutionPathB.times(1.0 - fExponentialMovingAverageWeightB)
                .add(fMTmp.times(fSubFromNormalizedNewBBTToEye,
                        Math.sqrt(fExponentialMovingAverageWeightB
                                * (2.0 - fExponentialMovingAverageWeightB) / approxCovB)));
    }

    private void calculateNormalizedNewBBT() {
        fMTmp.times(fB, fVTmp.tcopyFrom(fB));
        fNormalizedNewBBT.times(fInvSqrtBBT, fMTmp.times(fInvSqrtBBT));
    }

    /**
     * Calculates $\sqrt{BB^\top}^{-1}$.
     *
     * @throws Exception if $\sqrt{BB^\top}^{-1}$ has one or more negative eigenvalues.
     */
    private void calculateSqrtInvBBT() throws Exception {
        fRes = fBBTBeforeUpdate.eig();
        fMTmp.copyFrom(fRes.getD()).sqrt();
        if (fMTmp.isNan()) {
            throw new Exception("Negative eigenvalue(s)");
        }
        fMTmp.inverse();
        fVTmp.copyFrom(fRes.getV());
        fMTmp.times(fVTmp.tclone());
        fInvSqrtBBT.times(fVTmp, fMTmp);
    }

    /** Updates the distribution parameters, mean vector $m$, step-size $\sigma$, and shape B. */
    private void updateParameters() {
        fMean.add(fVecTmp.times(fB, fGMean).times(fEtaMean * fSigma));// Equation (6)
        fSigma *= Math.exp(fEtaSigma / 2.0 * fGSigma);// Equation (7)
        fB.times(fMTmp.copyFrom(fGB).times(fEtaB / 2.0).expm(true));// Equation (8)
    }

    /**
     * Calculates the natural gradients for the mean, step-size, and shape.
     * This involves sorting the population by evaluation values.
     */
    private void calculateNaturalGradients() {
        Collections.sort(fPopulation, fComparator);
        fGMean.fill(0.0);
        fGM.fill(0.0);
        TCMatrix wz = new TCMatrix(fDimension);
        TCMatrix zT = new TCMatrix(1, fDimension);
        for (int i = 0; i < fPopulationSize; ++i) {
            TCMatrix zi = fPopulation.get(i).getZVector();
            wz.copyFrom(zi).times(fWeights[i]);
            fGMean.add(wz);// Equation (2)
            fGM.add(fWZMatrix.times(wz, zT.tcopyFrom(zi)));// Equation (3)
        }
        fGM.sub(fWEye);// Equation (3)
        fGM.add(fMTmp.tcopyFrom(fGM)).div(2.);
        fGSigma = fGM.trace() / fDimension;// Equation (4)
        fGB.copyFrom(fGM).sub(fEye.times(fGSigma));// Equation (5)
        fEye.eye();
    }

}
