/**
 *
 * @author Ben Bohannon
 */
public class Matrix {

    /**
     * Backing array for the Matrix.
     */
    private double[][] array;

    /**
     * Initializes Matrix with given rows and columns, with zeros for elements.
     * @param rows Number of rows
     * @param columns Number of columns
     */
    public Matrix(final int rows, final int columns) {
        array = new double[rows][columns];
    }

    /**
     * Creates a Matrix with the input double array.
     * @param input Data to initialize Matrix.
     */
    public Matrix(final double[][] input) {
        if (input == null || input.length == 0 || input[0].length == 0) {
            throw new java.lang.IllegalArgumentException("Invalid input to "
                    + "create Matrix.");
        }

        array = new double[input.length][input[0].length];
        int numOfColumns = input[0].length;

        for (int i = 0; i < input.length; i++) {
            for (int j = 0; j < numOfColumns; j++) {
                if (input[i].length != numOfColumns) {
                    throw new java.lang.IllegalArgumentException("Input data is"
                            + "not rectangular!");
                }

                array[i][j] = input [i][j];
            }
        }

    }

    /**
     * Copy Constructor for a Matrix.
     * @param A Matrix to copy.
     */
    public Matrix(Matrix A) {
        array = new double[A.getRows()][A.getColumns()];

        for (int i = 0; i < array.length; i++) {
            System.arraycopy(A.array[i], 0, array[i], 0, array.length);
        }
    }

    /**
     * Return the element at the input row and column.
     * @param row Row to index.
     * @param column Column to index.
     * @return Element at requested position.
     */
    public final double get(final int row, final int column) {
        if (row < 0 || column < 0) {
            throw new java.lang.IllegalArgumentException("Cannot get negative "
                    + "index!");
        } else if (row > array.length - 1 || column > array[row].length - 1) {
            throw new java.lang.IllegalArgumentException("Indexed outside of "
                    + "Matrix.");
        }

        return array[row][column];
    }

    /**
     * Sets the element at the input position to the input data.
     * @param row Row to index.
     * @param column Column to index.
     * @param data Data to place.
     * @return Data that was replaced.
     */
    public final double set(final int row, final int column,
            final double data) {

        if (row < 0 || column < 0) {
            throw new java.lang.IllegalArgumentException("Cannot get negative "
                    + "index!");
        } else if (row > array.length || column > array[row].length) {
            throw new java.lang.IllegalArgumentException("Indexed outside of "
                    + "Matrix.");
        }

        double temp = array[row][column];
        array[row][column] = data;
        return temp;
    }

    /**
     * Returns the number of Rows in this Matrix.
     * @return number of Rows.
     */
    public final int getRows() {
        return array.length;
    }

    /**
     * Returns the number of Columns in this Matrix.
     * @return number of columns.
     */
    public final int getColumns() {
        return array[0].length;
    }

    public final double[] getRow(int row) {
        if (row < 0 || row > array.length) {
            throw new java.lang.IllegalArgumentException("Invalid row "
                    + "requested.");
        }

        double[] temp = new double[array[0].length];

        System.arraycopy(array[row], 0, temp, 0, temp.length);

        return temp;
    }

    public final double[] getColumn(int col) {
        if (col < 0 || col > array[0].length) {
            throw new java.lang.IllegalArgumentException("Invalid column "
                    + "requested.");
        }

        double[] temp = new double[array.length];

        for (int i = 0; i < temp.length; i++) {
            temp[i] = array[i][col];
        }

        return temp;
    }

    /**
     * Adds two Matrices together and returns their summed Matrix.
     * Does not change either Matrix.
     * @param other Matrix to add to.
     * @return The sum of the two Matrices.
     */
    public final Matrix add(final Matrix other) {
        if (this.getRows() != other.getRows()
                || this.getColumns() != other.getColumns()) {
            throw new java.lang.ArithmeticException("Can't add Matrices of "
                    + "different sizes.");
        }

        Matrix temp = new Matrix(this.getRows(), this.getColumns());

        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                temp.array[i][j] = this.array[i][j] + other.array[i][j];
            }
        }

        return temp;
    }

    /**
     * Subtracts the argument Matrix from this Matrix.
     * @param other other Matrix to subtract by
     * @return Result of the subtraction.
     */
    public final Matrix subtract(final Matrix other) {
        if (this.getRows() != other.getRows()
                || this.getColumns() != other.getColumns()) {
            throw new java.lang.ArithmeticException("Can't add Matrices of "
                    + "different sizes.");
        }

        Matrix temp = new Matrix(this.getRows(), this.getColumns());

        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                temp.array[i][j] = this.array[i][j] - other.array[i][j];
            }
        }

        return temp;
    }

    /**
     * Multiplies the Matrix by the input scalar value.
     * @param scalar value to multiply by.
     * @return The resultant Matrix.
     */
    public final Matrix multiply(final double scalar) {
        Matrix temp = new Matrix(this.getRows(), this.getColumns());

        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                temp.array[i][j] = this.array[i][j] * scalar;
            }
        }

        return temp;
    }

    /**
     * Multiplies this Matrix by other Matrix.
     * @param other Other matrix to multiply
     * @return The resultant Matrix.
     */
    public final Matrix multiply(final Matrix other) {
        if (this.getColumns() != other.getRows()) {
            throw new java.lang.ArithmeticException("Matrix dimensions cannot"
                    + " be multiplied!");
        }

        Matrix temp = new Matrix(this.getRows(), other.getColumns());

        for (int i = 0; i < this.getRows(); i++) {
            for (int j = 0; j < other.getColumns(); j++) {
                double sum = 0;
                //Get row i of this, and dot with column j of other.
                for (int pos = 0; pos < other.getRows(); pos++) {
                    sum += this.array[i][pos] * other.array[pos][j];
                }

                temp.array[i][j] = sum;
            }
        }

        return temp;
    }

    /**
     * Creates and returns the transpose of this Matrix.
     * @return transpose of this Matrix.
     */
    public final Matrix transpose() {
        Matrix temp = new Matrix(getColumns(), getRows());

        for (int i = 0; i < this.getRows(); i++) {
            for (int j = 0; j < this.getColumns(); j++) {
                temp.array[j][i] = this.array[i][j];
            }
        }

        return temp;
    }

    @Override
    public final boolean equals(final Object obj) {
        if (!(obj instanceof Matrix)) {
            return false;
        }
        if (obj == this) {
            return true;
        }
        Matrix other = (Matrix) obj;

        if (this.getRows() != other.getRows()
                || this.getColumns() != other.getColumns()) {
            return false;
        }

        for (int i = 0; i < this.getRows(); i++) {
            for (int j = 0; j < this.getColumns(); j++) {
                if (this.get(i, j) != other.get(i, j)) {
                    return false;
                }
            }
        }

        return true;
    }

    @Override
    public final int hashCode() {
        return 0;
    }

    @Override
    public final String toString() {
        String finalString = "";

        for (int i = 0; i < this.getRows(); i++) {
            String rowString = "[ ";
            for (int j = 0; j < this.getColumns(); j++) {
                //rowString += (array[i][j] + ", ");
                rowString += (String.format("%.15g", array[i][j]) + ", ");
            }
            rowString += "]\n";
            finalString += rowString;
        }

        return finalString;
    }

    /**
     * Gets a square Identity Matrix of size n by n.
     * @param n Dimension of Matrix.
     * @return Identity Matrix.
     */
    public static Matrix getIdentity(final int n) {
        Matrix temp = new Matrix(n, n);

        for (int i = 0; i < n; i++) {
            temp.set(i, i, 1);
        }

        return temp;
    }

    /**
     * Gets a square Hilbert Matrix of size n by n.
     * @param n Dimension of Matrix.
     * @return Hilbert Matrix.
     */
    public static Matrix getHilbert(final int n) {
        Matrix temp = new Matrix(n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp.array[i][j] = 1.0 / (1.0 + i + j);
            }
        }

        return temp;
    }
}
