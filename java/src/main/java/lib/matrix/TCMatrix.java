package lib.matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import lib.matrix.decompositions.TCCholeskyDecomposition;
import lib.matrix.decompositions.TCEigenvalueDecomposition;
import lib.matrix.decompositions.TCLUDecomposition;
import lib.matrix.decompositions.TCQRDecomposition;
import lib.matrix.decompositions.TCSingularValueDecomposition;
import lib.random.ICRandom;

public class TCMatrix implements Serializable, Cloneable {

    /** Serial version ID. */
    private static final long serialVersionUID = 1L;

    /** Array storing the matrix elements in row-major order. **/
    private final double[] fElements;

    /** Number of rows (height). **/
    private final int fNoOfRow;

    /** Number of columns (width). **/
    private final int fNoOfColumn;

    /**
     * Creates an noOfRow x noOfColumn matrix.
     * Each element of the matrix is initialized to 0.0.
     *
     * @param noOfRow The height of the matrix.
     * @param noOfColumn The width of the matrix.
     */
    public TCMatrix(final int noOfRow, final int noOfColumn) {
        fElements = new double[noOfColumn * noOfRow];
        fNoOfRow = noOfRow;
        fNoOfColumn = noOfColumn;
    }

    /**
     * Creates a column vector of the specified height.
     * Each element is initialized to 0.0.
     *
     * @param column The height of the vector.
     */
    public TCMatrix(final int column) {
        this(column, 1);
    }

    /**
     * Copy constructor.
     *
     * @param mat The source matrix to copy.
     */
    public TCMatrix(final TCMatrix mat) {
        fElements = new double[mat.fElements.length];
        fNoOfRow = mat.getRowDimension();
        fNoOfColumn = mat.getColumnDimension();
        copyFrom(mat);
    }

    /**
     * Initializes from a 1D array.
     * A deep copy of the array's contents is performed.
     *
     * @param array
     */
    public TCMatrix(final double[] array) {
        this(array.length, 1);
        System.arraycopy(array, 0, fElements, 0, array.length);
    }

    /**
     * Initializes from a 2D array.
     * A deep copy of the array's contents is performed.
     *
     * @param array
     */
    public TCMatrix(final double[][] array) {
        this(array.length, array[0].length);
        int iStepIndex = 0;
        for (int i = 0; i < fNoOfRow; i++) {
            System.arraycopy(array[i], 0, fElements, iStepIndex, fNoOfColumn);
            iStepIndex += fNoOfColumn;
        }
    }

    /**
     * Performs a deep copy of the specified source matrix into this matrix.
     *
     * @param src The source matrix to copy.
     * @return This matrix.
     */
    public TCMatrix copyFrom(final TCMatrix src) {
        assert (fNoOfRow == src.fNoOfRow && fNoOfColumn == src.fNoOfColumn);

        System.arraycopy(src.fElements, 0, fElements, 0, fElements.length);
        return this;
    }

    /**
     * Deep copies the transpose of the source matrix 'src' into this matrix.
     *
     * @param src The source matrix.
     * @return This matrix.
     */
    public TCMatrix tcopyFrom(final TCMatrix src) {
        assert (fNoOfRow == src.fNoOfColumn && fNoOfColumn == src.fNoOfRow);

        for (int i = 0; i < fNoOfColumn; i++) {
            for (int j = 0; j < fNoOfRow; j++) {
                setValue(j, i, src.getValue(i, j));
            }
        }
        return this;
    }

    /**
     * Clones this matrix in a transposed state.
     *
     * @return Transposed clone of this matrix.
     */
    public TCMatrix tclone() {
        return new TCMatrix(fNoOfColumn, fNoOfRow).tcopyFrom(this);
    }

    /**
     * Copies a submatrix from 'src' to a specified location in this matrix.
     * The submatrix is defined by the source row and column ranges.
     * The indices for the start and end of the submatrix are inclusive.
     *
     * @param src The source matrix.
     * @param srcFromRow The starting row index of the submatrix in the source.
     * @param srcToRow The ending row index of the submatrix in the source.
     * @param srcFromColumn The starting column index of the submatrix in the source.
     * @param srcToColumn The ending column index of the submatrix in the source.
     * @param dstFromRow The starting row index of the destination in this matrix.
     * @param dstFromColumn The starting column index of the destination in this matrix.
     * @return this
     */
    public TCMatrix copySubmatrixFrom(final TCMatrix src, final int srcFromRow, final int srcToRow,
            final int srcFromColumn, final int srcToColumn, final int dstFromRow,
            final int dstFromColumn) {

        final int submatrixColumnDimension = srcToColumn - srcFromColumn + 1;
        final int submatrixRowDimension = srcToRow - srcFromRow + 1;

        assert (0 <= dstFromRow && dstFromRow + submatrixRowDimension <= fNoOfRow && // Check if the
                                                                                     // submatrix
                                                                                     // height
                                                                                     // exceeds the
                                                                                     // bounds of
                                                                                     // the
                                                                                     // destination
                                                                                     // matrix.
                0 <= dstFromColumn && dstFromColumn + submatrixColumnDimension <= fNoOfColumn && // Check
                                                                                                 // if
                                                                                                 // the
                                                                                                 // submatrix
                                                                                                 // width
                                                                                                 // exceeds
                                                                                                 // the
                                                                                                 // bounds
                                                                                                 // of
                                                                                                 // the
                                                                                                 // destination
                                                                                                 // matrix.
                0 <= srcFromRow && srcFromRow <= srcToRow && srcToRow < src.fNoOfRow && // Check if
                                                                                        // the
                                                                                        // submatrix
                                                                                        // height
                                                                                        // exceeds
                                                                                        // the
                                                                                        // bounds of
                                                                                        // the
                                                                                        // source
                                                                                        // matrix.
                0 <= srcFromColumn && srcFromColumn <= srcToColumn
                && srcToColumn < src.fNoOfColumn); // Check if the submatrix width exceeds the
                                                   // bounds of the source matrix.

        for (int i = 0; i < submatrixRowDimension; i++) {
            for (int j = 0; j < submatrixColumnDimension; j++) {
                setValue(i + dstFromRow, j + dstFromColumn,
                        src.getValue(i + srcFromRow, j + srcFromColumn));
            }
        }
        return this;
    }

    /**
     * Returns the number of columns of this matrix (width).
     *
     * @return
     */
    public int getColumnDimension() {
        return fNoOfColumn;
    }

    /**
     * Returns the number of rows of this matrix (height).
     *
     * @return
     */
    public int getRowDimension() {
        return fNoOfRow;
    }

    /**
     * Returns the total number of elements in the matrix (rows * columns).
     *
     * @return
     */
    public int getDimension() {
        return fElements.length;
    }

    /**
     * Returns the element at the specified row and column.
     *
     * @param r
     * @param c
     * @return
     */
    public double getValue(final int r, final int c) {
        assert (0 <= r && r < fNoOfRow && 0 <= c && c < fNoOfColumn);

        return fElements[r * fNoOfColumn + c];
    }

    /**
     * Returns the element at a single index 'r' from the underlying 1D array.
     *
     * @param r
     * @return
     */
    public double getValue(final int r) {
        assert (0 <= r && r < fElements.length);
        return fElements[r];
    }

    /**
     * For a 1x1 matrix, this is the same as getValue(0, 0).
     *
     * @return
     */
    public double getValue() {
        assert (fElements.length == 1);
        return fElements[0];
    }

    /**
     * Sets the element at the specified row and column to a value.
     *
     * @param r
     * @param c
     * @param val
     * @return
     */
    public TCMatrix setValue(final int r, final int c, final double val) {
        assert (0 <= r && r < fNoOfRow && 0 <= c && c < fNoOfColumn);

        fElements[r * fNoOfColumn + c] = val;
        return this;
    }

    /**
     * Sets the element at a single index 'r' in the underlying 1D array to a value.
     *
     * @param r
     * @param val
     * @return
     */
    public TCMatrix setValue(final int r, final double val) {
        assert (0 <= r && r < fElements.length);
        fElements[r] = val;
        return this;
    }

    /**
     * For a 1x1 matrix, this is the same as setValue(0, 0, val).
     *
     * @return
     */
    public TCMatrix setValue(final double val) {
        assert (fElements.length == 1);
        fElements[0] = val;
        return this;
    }

    /**
     * Performs this = this + m.
     * The dimensions of this matrix and m must be the same.
     *
     * @param m
     * @return this
     */
    public TCMatrix add(final TCMatrix m) {
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] += m.fElements[i];
        }
        return this;
    }

    /**
     * Performs this = m1 + m2.
     * The dimensions of this matrix, m1, and m2 must be the same.
     *
     * @param m1
     * @param m2
     * @return
     */
    public TCMatrix add(final TCMatrix m1, final TCMatrix m2) {
        assert (fNoOfRow == m1.fNoOfRow && fNoOfRow == m2.fNoOfRow && fNoOfColumn == m1.fNoOfColumn
                && fNoOfColumn == m2.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m1.fElements[i] + m2.fElements[i];
        }
        return this;
    }

    /**
     * Adds a scalar value to each element of this matrix.
     *
     * @param val
     * @return this
     */
    public TCMatrix add(final double val) {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] += val;
        }
        return this;
    }

    /**
     * Sets this matrix to the result of adding a scalar value to each element of matrix m.
     *
     * @param m
     * @param val
     * @return
     */
    public TCMatrix add(final TCMatrix m, final double val) {
        assert (val != 0.0);
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m.fElements[i] + val;
        }
        return this;
    }

    /**
     * Performs this = this - m.
     * The dimensions of this matrix and m must be the same.
     *
     * @param m
     * @return
     */
    public TCMatrix sub(final TCMatrix m) {
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);

        for (int i = 0; i < fElements.length; i++) {
            fElements[i] -= m.fElements[i];
        }
        return this;
    }

    /**
     * Performs this = m1 - m2.
     * The dimensions of this matrix, m1, and m2 must be the same.
     *
     * @param m1
     * @param m2
     * @return
     */
    public TCMatrix sub(final TCMatrix m1, final TCMatrix m2) {
        assert (fNoOfRow == m1.fNoOfRow && fNoOfRow == m2.fNoOfRow && fNoOfColumn == m1.fNoOfColumn
                && fNoOfColumn == m2.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m1.fElements[i] - m2.fElements[i];
        }
        return this;
    }

    /**
     * Subtracts a scalar value from each element of this matrix.
     *
     * @param val
     * @return this
     */
    public TCMatrix sub(final double val) {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] -= val;
        }
        return this;
    }

    /**
     * Sets this matrix to the result of subtracting a scalar value from each element of matrix m.
     *
     * @param m
     * @param val
     * @return
     */
    public TCMatrix sub(final TCMatrix m, final double val) {
        assert (val != 0.0);
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m.fElements[i] - val;
        }
        return this;
    }

    /**
     * Performs this = m * m2 and returns this.
     * Classified into three patterns for optimization.
     *
     * @param m
     * @param m2
     * @return this
     */
    public TCMatrix times(final TCMatrix m, final TCMatrix m2) {
        assert (m.fNoOfColumn == m2.fNoOfRow && fNoOfRow == m.fNoOfRow
                && fNoOfColumn == m2.fNoOfColumn);

        if (m2.fNoOfColumn > 1) {
            if (fElements.length > 600) { // If the matrix size is large, an ikj loop order is
                                          // faster due to better cache utilization.
                Arrays.fill(fElements, 0.0);
                int mIStepIndex = 0;
                for (int i = 0; i < fElements.length; i += fNoOfColumn) {
                    int m2KStepIndex = 0;
                    for (int k = 0; k < m.fNoOfColumn; k++) {
                        for (int j = 0; j < fNoOfColumn; j++) {
                            fElements[i + j] +=
                                    m.fElements[mIStepIndex + k] * m2.fElements[m2KStepIndex + j];
                        }
                        m2KStepIndex += fNoOfColumn;
                    }
                    mIStepIndex += m.fNoOfColumn;
                }
            } else { // If the matrix size is small, caching the value to be written is faster.
                int mIStepIndex = 0;
                for (int i = 0; i < fElements.length; i += fNoOfColumn) {
                    for (int j = 0; j < fNoOfColumn; j++) {
                        int m2KStepIndex = 0;
                        double val = 0.0;
                        for (int k = 0; k < m.fNoOfColumn; k++) {
                            val += m.fElements[mIStepIndex + k] * m2.fElements[m2KStepIndex + j];
                            m2KStepIndex += fNoOfColumn;
                        }
                        fElements[i + j] = val;
                    }
                    mIStepIndex += m.fNoOfColumn;
                }
            }
        } else { // In the case of matrix-vector multiplication, a double loop is faster.
            Arrays.fill(fElements, 0.0);
            int iStepIndex = 0;
            for (int i = 0; i < fNoOfRow; i++) {
                for (int j = 0; j < m.fNoOfColumn; j++) {
                    fElements[i] += m.fElements[iStepIndex + j] * m2.fElements[j];
                }
                iStepIndex += m.fNoOfColumn;
            }
        }
        return this;
    }

    /**
     * Performs this = this * m. The result overwrites this matrix.
     * Both 'this' and 'm' must be square matrices of the same size.
     * Otherwise, using times(m1, m2) is recommended.
     *
     * @param m
     * @return
     */
    public TCMatrix times(final TCMatrix m) {

        assert (fNoOfRow == fNoOfColumn && fNoOfColumn == m.fNoOfRow
                && m.fNoOfRow == m.fNoOfColumn);

        final int size = fNoOfRow;
        final double[] ithRow = new double[size];

        for (int i = 0; i < fElements.length; i += size) {
            System.arraycopy(fElements, i, ithRow, 0, size);
            Arrays.fill(fElements, i, i + size, 0.0);
            int mKIndex = 0;
            for (int k = 0; k < size; k++) {
                for (int j = 0; j < size; j++) {
                    fElements[i + j] += ithRow[k] * m.fElements[mKIndex + j];
                }
                mKIndex += size;
            }
        }
        return this;
    }

    /**
     * Multiplies every element of this matrix by a scalar value.
     *
     * @param val The scalar value to multiply by.
     * @return this
     */
    public TCMatrix times(final double val) {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] *= val;
        }
        return this;
    }

    /**
     * Sets this matrix to the result of multiplying each element of matrix m by a scalar value.
     *
     * @param m
     * @param val
     * @return
     */
    public TCMatrix times(final TCMatrix m, final double val) {
        assert (val != 0.0);
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m.fElements[i] * val;
        }
        return this;
    }

    /**
     * Performs this = m * this. The result overwrites this matrix.
     * 'this' and 'm' must be n x n square matrices of the same size.
     * This may be slightly slower due to the creation of a new double[n] array.
     * If performance is critical, using times(m1, m2) is recommended.
     *
     * @param m The matrix to multiply from the left.
     * @return this
     */
    public TCMatrix timesLeft(final TCMatrix m) {
        assert (fNoOfRow == fNoOfColumn && fNoOfColumn == m.fNoOfRow
                && m.fNoOfRow == m.fNoOfColumn);

        final int size = fNoOfRow;
        final double[] jthColumn = new double[size];

        for (int j = 0; j < size; j++) {
            int iStepIndex = 0;
            for (int i = 0; i < size; i++) {
                jthColumn[i] = fElements[iStepIndex + j];
                fElements[iStepIndex + j] = 0;
                iStepIndex += size;
            }
            for (int i = 0; i < fElements.length; i += size) {
                for (int k = 0; k < size; k++) {
                    fElements[i + j] += m.fElements[i + k] * jthColumn[k];
                }
            }
        }
        return this;
    }

    /**
     * Performs this = this .* m (element-wise multiplication).
     *
     * @param m
     * @return
     */
    public TCMatrix timesElement(final TCMatrix m) {
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);

        for (int i = 0; i < fElements.length; i++) {
            fElements[i] *= m.fElements[i];
        }

        return this;
    }

    /**
     * Performs this = m1 .* m2 (element-wise multiplication).
     *
     * @param m1
     * @param m2
     * @return
     */
    public TCMatrix timesElement(final TCMatrix m1, final TCMatrix m2) {
        assert (fNoOfRow == m1.fNoOfRow && fNoOfRow == m2.fNoOfRow && fNoOfColumn == m1.fNoOfColumn
                && fNoOfColumn == m2.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m1.fElements[i] * m2.fElements[i];
        }
        return this;
    }

    /**
     * Divides every element of this matrix by a scalar value.
     *
     * @param val
     * @return
     */
    public TCMatrix div(final double val) {
        assert (val != 0.0);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] /= val;
        }
        return this;
    }

    /**
     * Sets this matrix to the result of dividing each element of matrix m by a scalar value.
     *
     * @param m
     * @param val
     * @return
     */
    public TCMatrix div(final TCMatrix m, final double val) {
        assert (val != 0.0);
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m.fElements[i] / val;
        }
        return this;
    }

    /**
     * Performs this = this ./ m (element-wise division).
     *
     * @param m
     * @return
     */
    public TCMatrix divElement(TCMatrix m) {
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);

        for (int i = 0; i < fElements.length; i++) {
            fElements[i] /= m.fElements[i];
        }
        return this;
    }

    /**
     * Performs this = m1 ./ m2 (element-wise division).
     *
     * @param m1
     * @param m2
     * @return
     */
    public TCMatrix divElement(final TCMatrix m1, final TCMatrix m2) {
        assert (fNoOfRow == m1.fNoOfRow && fNoOfRow == m2.fNoOfRow && fNoOfColumn == m1.fNoOfColumn
                && fNoOfColumn == m2.fNoOfColumn);
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = m1.fElements[i] / m2.fElements[i];
        }
        return this;
    }

    public double trace() {
        final int length = Math.min(fNoOfColumn, fNoOfRow);
        double sum = 0.0;
        int iStepIndex = 0;
        for (int i = 0; i < length; i++) {
            sum += fElements[iStepIndex + i];
            iStepIndex += fNoOfColumn;
        }
        return sum;
    }

    /**
     * Do cholesky decomposition.
     *
     * @return TCCholeskyDecomposition Class of this matrix
     */
    public TCCholeskyDecomposition chol() {
        return new TCCholeskyDecomposition(this);
    }

    /**
     * Do LU decomposition.
     *
     * @return TCLUDecomposition Class of this matrix
     */
    public TCLUDecomposition lu() {
        return new TCLUDecomposition(this);
    }

    /**
     * Do QR decomposition.
     *
     * @return TCQRDecomposition Class of this matrix
     */
    public TCQRDecomposition qr() {
        return new TCQRDecomposition(this);
    }

    /**
     * Do Eigenvalue decomposition.
     *
     * @return TCEigenvalueDecomposition Class of this matrix
     */
    public TCEigenvalueDecomposition eig() {
        return new TCEigenvalueDecomposition(this);
    }

    /**
     * Do SingularValue decomposition.
     *
     * @return TCSingularValuDecomposition Class of this matrix.
     */
    public TCSingularValueDecomposition svd() {
        return new TCSingularValueDecomposition(this);
    }

    /**
     * Returns determinant of this matrix.
     *
     * @return the determinant of this matrix
     */
    public double det() {
        return new TCLUDecomposition(this).det();
    }

    /**
     * Returns TCMatrix X such that satisfies A*X = B, where A is this matrix.
     *
     * @param B
     * @return The solution of A*X = B.
     */
    public TCMatrix solve(TCMatrix B) {
        return (getRowDimension() == getColumnDimension() ? (new TCLUDecomposition(this)).solve(B)
                : (new TCQRDecomposition(this)).solve(B));
    }

    /**
     * Returns an inverse of this matrix.
     * Updates this object.
     *
     * @return this matrix
     */
    public TCMatrix inverse() {
        final TCMatrix identityMatrix = new TCMatrix(this).eye();
        this.copyFrom(solve(identityMatrix));
        return this;
    }

    /**
     * Returns rank of this matrix.
     *
     * @return the rank of this matrix
     */
    public int rank() {
        return new TCSingularValueDecomposition(this).rank();
    }

    /**
     * Returns condition number (singular value decomposition) of this matrix.
     *
     * @return the condition number of this matrix
     */
    public double cond() {
        return new TCSingularValueDecomposition(this).cond();
    }

    public double normF() {
        double sum = 0.0;
        for (final double val : fElements) {
            sum += val * val;
        }
        return Math.sqrt(sum);
    }

    /**
     * Calculates the L2 norm of this matrix.
     *
     * @return L2 norm
     */
    public double normL2() {
        assert (getColumnDimension() == 1);

        double sum = 0.0;
        for (final double val : fElements) {
            sum += val * val;
        }
        return Math.sqrt(sum);
    }

    /**
     * Calculates the Euclidean distance between this matrix and another vector.
     *
     * @param vector The vector to calculate the distance to.
     * @return Euclidean distance
     */
    public double distanceL2(final TCMatrix vector) {
        assert (fNoOfColumn == 1 && vector.fNoOfColumn == 1);

        double sum = 0.0;
        for (int i = 0; i < fElements.length; i++) {
            final double diff = fElements[i] - vector.fElements[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    /**
     * Normalizes the current vector to a unit vector.
     * This matrix must be a column vector.
     *
     * @return this
     */
    public TCMatrix enforceToUnitVector() {
        assert (fNoOfColumn == 1);

        final double norm = normL2();
        div(norm);
        return this;
    }

    /**
     * Calculates the inner product with a vector.
     * Both 'this' and 'vector' must be column vectors of the same dimension.
     *
     * @param vector
     * @return <this, vector>
     */
    public double innerProduct(final TCMatrix vector) {
        assert (fNoOfColumn == 1 && vector.fNoOfColumn == 1 && fNoOfRow == vector.fNoOfRow);
        double sum = 0.0;
        for (int i = 0; i < fElements.length; i++) {
            sum += fElements[i] * vector.fElements[i];
        }
        return sum;
    }

    /**
     * Overwrites all elements of the current matrix with a specified value.
     *
     * @param val
     * @return this
     */
    public TCMatrix fill(final double val) {
        Arrays.fill(fElements, val);
        return this;
    }

    /**
     * Initializes the current matrix with uniformly distributed random numbers in the range [0, 1].
     *
     * @param random The random number generator.
     * @return this
     */
    public TCMatrix rand(final ICRandom random) {
        for (int i = 0; i < getDimension(); i++) {
            fElements[i] = random.nextDouble();
        }
        return this;
    }

    /**
     * Overwrites the current matrix with random numbers from a standard normal distribution.
     *
     * @param random The random number generator.
     * @return this
     */
    public TCMatrix randn(final ICRandom random) {
        for (int i = 0; i < getDimension(); i++) {
            fElements[i] = random.nextGaussian();
        }
        return this;
    }

    /**
     * Replaces the current matrix with an identity matrix.
     *
     * @return this
     */
    public TCMatrix eye() {
        Arrays.fill(fElements, 0.0);
        final int size = Math.min(fNoOfRow, fNoOfColumn);
        for (int i = 0; i < size; i++) {
            fElements[i * fNoOfColumn + i] = 1.0;
        }
        return this;
    }

    /**
     * Replaces this matrix with the matrix exponential of this, expm(this).
     *
     * @param sym True if the current matrix is symmetric.
     * @return this
     */
    public TCMatrix expm(final boolean sym) {
        assert (fNoOfColumn == fNoOfRow);

        final TCEigenvalueDecomposition eig = eig();
        final TCMatrix expD = eig.getD();
        final TCMatrix v = eig.getV();

        for (int i = 0; i < fNoOfColumn; i++) {
            expD.setValue(i, i, Math.expm1(expD.getValue(i, i)) + 1.0);
        }

        if (sym) {
            final TCMatrix vT = v.tclone();
            this.copyFrom(v.times(expD).times(vT));
        } else {
            final TCMatrix vI = v.clone().inverse();
            this.copyFrom(v.times(expD).times(vI));
        }
        return this;
    }

    /**
     * Replaces each element x of the matrix with exp(x).
     *
     * @return this
     */
    public TCMatrix exp() {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = Math.expm1(fElements[i]) + 1.0;
        }
        return this;
    }

    /**
     * Replaces each element x of the matrix with sin(x).
     *
     * @return this
     */
    public TCMatrix sin() {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = Math.sin(fElements[i]);
        }
        return this;
    }

    /**
     * Replaces each element x of the matrix with cos(x).
     *
     * @return this
     */
    public TCMatrix cos() {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = Math.cos(fElements[i]);
        }
        return this;
    }

    /**
     * Replaces each element x of the matrix with tan(x).
     *
     * @return this
     */
    public TCMatrix tan() {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = Math.tan(fElements[i]);
        }
        return this;
    }

    /**
     * Replaces each element x of the matrix with sqrt(x).
     *
     * @return this
     */
    public TCMatrix sqrt() {
        for (int i = 0; i < fElements.length; i++) {
            fElements[i] = Math.sqrt(fElements[i]);
        }
        return this;
    }

    /**
     * Compares each element of this matrix with the corresponding element of matrix m,
     * and replaces it with the larger of the two.
     *
     * @param m matrix
     * @return this
     */
    public TCMatrix max(final TCMatrix m) {
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);

        for (int i = 0; i < fElements.length; i++) {
            if (fElements[i] < m.fElements[i]) {
                fElements[i] = m.fElements[i];
            }
        }
        return this;
    }

    /**
     * Replaces each element a_(i,j) of the matrix with max(a_(i,j), val).
     *
     * @param val The value to compare with.
     * @return this
     */
    public TCMatrix max(final double val) {
        for (int i = 0; i < fElements.length; i++) {
            if (fElements[i] < val) {
                fElements[i] = val;
            }
        }
        return this;
    }

    /**
     * Finds the maximum value in an m-by-n matrix and its corresponding index(es).
     * If there are multiple occurrences of the maximum value, all their indices are retrieved.
     *
     * @param indices A list to store the 'index' of the maximum value(s).
     *        Note that the 'index' for an element (i,j) in an m-by-n matrix is represented as i*n +
     *        j.
     * @return max The maximum value.
     */
    public double max(final ArrayList<Integer> indices) {
        indices.clear();
        double max = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < fElements.length; i++) {
            if (max < fElements[i]) {
                indices.clear();
                indices.add(i);
                max = fElements[i];
            } else if (max == fElements[i]) {
                indices.add(i);
            }
        }
        return max;
    }

    /**
     * Compares each element of this matrix with the corresponding element of matrix m,
     * and replaces it with the smaller of the two.
     *
     * @param m matrix
     * @return this
     */
    public TCMatrix min(final TCMatrix m) {
        assert (fNoOfRow == m.fNoOfRow && fNoOfColumn == m.fNoOfColumn);

        for (int i = 0; i < fElements.length; i++) {
            if (fElements[i] > m.fElements[i]) {
                fElements[i] = m.fElements[i];
            }
        }
        return this;
    }

    /**
     * Replaces each element a_(i,j) of the matrix with min(a_(i,j), val).
     *
     * @param val The value to compare with.
     * @return this
     */
    public TCMatrix min(final double val) {
        for (int i = 0; i < fElements.length; i++) {
            if (val < fElements[i]) {
                fElements[i] = val;
            }
        }
        return this;
    }

    /**
     * Finds the minimum value in an m-by-n matrix and its corresponding index(es).
     * If there are multiple occurrences of the minimum value, all their indices are retrieved.
     *
     * @param indices A list to store the 'index' of the minimum value(s).
     *        Note that the 'index' for an element (i,j) in an m-by-n matrix is represented as i*n +
     *        j.
     * @return min The minimum value.
     */
    public double min(final ArrayList<Integer> indices) {

        indices.clear();
        double min = Double.MAX_VALUE;
        for (int i = 0; i < fElements.length; i++) {
            if (fElements[i] < min) {
                indices.clear();
                indices.add(i);
                min = fElements[i];
            } else if (fElements[i] == min) {
                indices.add(i);
            }
        }
        return min;
    }

    @Override
    public TCMatrix clone() {
        return new TCMatrix(this);
    }

    /**
     * Creates and returns a clone of the submatrix corresponding to the rows from src_fromRow
     * to src_toRow and columns from src_fromColumn to src_toColumn. All endpoint indices are
     * inclusive.
     *
     * @param srcFromRow first row index
     * @param srcToRow last row index
     * @param srcFromColumn first column index
     * @param srcToColumn last column index
     * @return A clone of the submatrix.
     */
    public TCMatrix cloneSubmatrix(final int srcFromRow, final int srcToRow,
            final int srcFromColumn, final int srcToColumn) {

        final int subRowDimension = srcToRow - srcFromRow + 1;
        final int subColumnDimension = srcToColumn - srcFromColumn + 1;
        final TCMatrix result = new TCMatrix(subColumnDimension, subRowDimension);
        result.copySubmatrixFrom(this, srcFromRow, srcToRow, srcFromColumn, srcToColumn, 0, 0);
        return result;
    }

    /**
     * Returns a clone of the row vector at the specified row index.
     *
     * @param rowIndex The row index.
     * @return The row vector.
     */
    public TCMatrix cloneRowVector(final int rowIndex) {
        final TCMatrix result = new TCMatrix(1, fNoOfColumn);
        final int columnIndex = rowIndex * fNoOfColumn;
        System.arraycopy(fElements, columnIndex, result.fElements, 0, fNoOfColumn);
        return result;
    }

    /**
     * Returns a clone of the column vector at the specified column index.
     *
     * @param columnIndex The column index.
     * @return The column vector.
     */
    public TCMatrix cloneColumnVector(final int columnIndex) {
        final TCMatrix result = new TCMatrix(fNoOfRow);
        int iStepIndex = 0;
        for (int i = 0; i < fNoOfRow; ++i) {
            result.setValue(i, fElements[iStepIndex + columnIndex]);
            iStepIndex += fNoOfColumn;
        }
        return result;
    }

    /**
     * Checks if the matrix contains any NaN elements.
     *
     * @return true: contains NaN, false: does not contain NaN
     */
    public boolean isNan() {
        for (final double val : fElements) {
            if (Double.isNaN(val)) {
                return true;
            }
        }
        return false;
    }

    public String toString() {

        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < fNoOfRow; i++) {
            sb.append("[");
            sb.append(getValue(i, 0));
            for (int j = 1; j < fNoOfColumn; j++) {
                sb.append(", " + getValue(i, j));
            }
            if (i == fNoOfRow - 1) {
                sb.append("]");
            } else {
                sb.append("], \n");
            }
        }
        sb.append("]\n");

        return sb.toString();
    }

    /**
     * Copies a matrix to a specified row of this matrix.
     *
     * @param m The matrix to copy. Its number of columns must be the same as this matrix.
     * @param row The starting row to copy to. (row + m.fNoOfRow) must be within the row bounds of
     *        this matrix.
     * @return This matrix after copying.
     */
    public TCMatrix copyAtRow(final TCMatrix m, final int row) {
        assert (fNoOfColumn == m.fNoOfColumn && 0 <= row && row + m.fNoOfRow <= fNoOfRow);

        final int additionalIndex = row * fNoOfColumn;
        for (int i = 0; i < m.fElements.length; i += fNoOfColumn) {
            for (int j = 0; j < m.fNoOfColumn; ++j) {
                fElements[i + additionalIndex + j] = m.fElements[i + j];
            }
        }
        return this;
    }

    /**
     * Copies a matrix to a specified column of this matrix.
     *
     * @param m The matrix to copy. Its number of rows must be the same as this matrix.
     * @param column The starting column to copy to. (column + m.fNoOfColumn) must be within the
     *        column bounds of this matrix.
     * @return This matrix after copying.
     */
    public TCMatrix copyAtColumn(final TCMatrix m, final int column) {
        assert (fNoOfRow == m.fNoOfRow && 0 <= column && column + m.fNoOfColumn <= fNoOfColumn);

        int mIStepIndex = 0;
        for (int i = 0; i < fElements.length; i += fNoOfColumn) {
            for (int j = 0; j < m.fNoOfColumn; j++) {
                fElements[i + j + column] = m.fElements[mIStepIndex + j];
            }
            mIStepIndex += m.fNoOfColumn;
        }
        return this;
    }

    /**
     * Enforces symmetry on this matrix. The upper triangle is copied to the lower triangle.
     *
     * @return The resulting symmetric matrix (this).
     */
    public TCMatrix enforceSymmetry() {
        assert (fNoOfColumn == fNoOfRow);

        int iStepIndex = 0;
        for (int i = 0; i < fNoOfRow; i++) {
            int jStepIndex = 0;
            for (int j = 0; j < i; j++) {
                fElements[jStepIndex + i] = fElements[iStepIndex + j];
                jStepIndex += fNoOfColumn;
            }
            iStepIndex += fNoOfColumn;
        }
        return this;
    }

    @Override
    public boolean equals(final Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof TCMatrix)) {
            return false;
        }
        final TCMatrix o_ = (TCMatrix) o;
        return getColumnDimension() == o_.getColumnDimension()
                && getRowDimension() == o_.getRowDimension()
                && Arrays.equals(fElements, o_.fElements);
    }

    @Override
    public int hashCode() {
        final int multiplier = 37;
        int result = 17;
        for (int i = 0; i < getRowDimension(); i++) {
            for (int j = 0; j < getColumnDimension(); j++) {
                final long l = Double.doubleToLongBits(getValue(i, j));
                result = multiplier * result + (int) (l ^ (l >>> 32));
            }
        }
        return result;
    }
}
