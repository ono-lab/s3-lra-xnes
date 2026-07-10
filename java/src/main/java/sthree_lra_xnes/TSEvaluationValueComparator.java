package sthree_lra_xnes;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Compares solutions based on the objective function of the single-objective problem.
 */
public class TSEvaluationValueComparator implements Comparator<TSRealSolution>, Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Is this a minimization problem?
     * true: yes, false:no (, the problem is a maximization problem.)
     */
    private boolean fMinimization;

    /**
     * @param minimization Is this a minimization problem? true: yes, false:no (a maximization one.)
     */
    public TSEvaluationValueComparator(boolean minimization) {
        fMinimization = minimization;
    }

    @Override
    public int compare(TSRealSolution a, TSRealSolution b) {
        int sgn = 0;
        if (a.getEvaluationValue() - b.getEvaluationValue() < 0) {
            sgn = -1;
        } else if (a.getEvaluationValue() - b.getEvaluationValue() > 0) {
            sgn = 1;
        } else {
            sgn = 0;
        }
        return fMinimization ? sgn : -sgn;

    }

}
