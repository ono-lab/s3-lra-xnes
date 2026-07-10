package main;

import lib.matrix.TCMatrix;
import lib.random.ICRandom;
import lib.random.TCRandom;
import sthree_lra_xnes.S3LRAxNES;
import sthree_lra_xnes.TSRealSolution;

/**
 * Sample code of optimization with the proposed method, S3-LRA-xNES.
 */
public class Main {

    static double Sphere(TCMatrix x) {
        double eval = 0;
        for (int i = 0; i < x.getDimension(); ++i) {
            double xi = x.getValue(i);
            eval += xi * xi;
        }
        return eval;
    }

    static double Rastrigin(TCMatrix x) {
        int n = x.getDimension();
        double f = 10.0 * n;
        for (int i = 0; i < n; i++) {
            double xi = x.getValue(i);
            f += xi * xi - 10.0 * Math.cos(2.0 * Math.PI * xi);
        }
        return f;
    }

    static double Ellipsoid(TCMatrix x) {
        int n = x.getDimension();
        double f = 0.0;
        for (int i = 0; i < n; i++) {
            double xi = x.getValue(i);
            double fi = Math.pow(1000.0, (double) i / ((double) (n - 1.0))) * xi;
            f += fi * fi;
        }
        return f;
    }

    static double kTablet(TCMatrix x, int k) {
        int n = x.getDimension();
        double f = 0.0;
        for (int i = 0; i < k; i++) {
            double xi = x.getValue(i);
            f += xi * xi;
        }
        for (int i = k; i < n; i++) {
            double xi = x.getValue(i);
            f += (100.0 * xi) * (100.0 * xi);
        }
        return f;
    }

    static double Bohachevsky(TCMatrix x) {
        int n = x.getDimension();
        double f = 0;
        for (int i = 0; i < n - 1; i++) {
            double xi = x.getValue(i);
            double xi1 = x.getValue(i + 1);
            f += xi * xi + 2.0 * xi1 * xi1 + 0.7
                    - 0.3 * Math.cos(3.0 * Math.PI * xi)
                    - 0.4 * Math.cos(4.0 * Math.PI * xi1);
        }
        return f;
    }

    static double Schaffer(TCMatrix x) {
        int n = x.getDimension();
        double f = 0;
        for (int i = 0; i < n - 1; i++) {
            double xi = x.getValue(i);
            double xi1 = x.getValue(i + 1);
            double yi = xi * xi + xi1 * xi1;
            f += Math.pow(yi, 0.25) * (Math.pow(Math.sin(50.0 * Math.pow(yi, 0.1)), 2.0) + 1.0);
        }
        return f;
    }
    // Add here any benchmark problems that receives a solution-vector as a TCMatrix
    // and returns its evaluation value as a double.

    public static void main(String[] args) {
        /** The problem dimension. **/
        int dimension = 20;
        /** The initial mean vector. */
        double meanElem = 3;
        /** The initial step-size. */
        double sigma = 2;
        /** The random seed. */
        long seed = 0l;
        /** The evaluation budget */
        long maxEval = dimension == 80 ? (long) 2e7 : (long) 1e7;
        /** The target function value. */
        double criterion = 1e-8;

        ICRandom random = new TCRandom(seed);
        TCMatrix mean = new TCMatrix(dimension).fill(meanElem);

        double best = Double.POSITIVE_INFINITY;
        long noOfEvals = 0;
        try {
            S3LRAxNES proposed = new S3LRAxNES(dimension, mean, sigma, random);
            /** is the optimization successful? */
            boolean successFlag = false;
            int populationSize = proposed.getPopulationSize();

            // Stopping criterion: Evaluation budget
            while (noOfEvals < maxEval) {

                // Sample and evaluate candidate solutions.
                for (int i = 0; i < populationSize; ++i) {
                    TSRealSolution s = proposed.sampleCandidateSolution(i);
                    // Change the objective function if you want:
                    double value = Sphere(s.getVector());
                    // double value = Rastrigin(s.getVector()); // like this,
                    // double value = Ellipsoid(s.getVector()); // or like this.
                    s.setEvaluationValue(value);
                    best = Math.min(best, value);
                }
                noOfEvals += populationSize;

                // Executes one generation of the S3LRA-xNES algorithm.
                proposed.nextGeneration();

                // Stopping criterion: target function value
                if (best < criterion) {
                    successFlag = true;
                    break;
                }

            }
            System.out.println("Success:" + successFlag);
            System.out.println("NoOfEvals:" + noOfEvals);
            System.out.println("Best:" + best);
        } catch (Exception e) {
            System.err.println("Error: " + e);
        }
    }
}
