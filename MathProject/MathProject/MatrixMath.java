
/**
 *
 * @author Ben Bohannon
 */
public final class MatrixMath {

    /**
     * Private constructor, so this class cannot be instantiated.
     */
    private MatrixMath() {
        //MatrixMath is a static class. It should neither be able to be
        //instantiated, nor have any instance methods/variables.
    }

    /**
     * Decomposes the input Matrix into L, lower triangular, and U,
     * upper triangular, Matrices.
     * @param A Matrix to decompose.
     * @return The L and U Matrix in that order in an array.
     */
    public static Matrix[] lu_fact(Matrix A) {

        if (A.getRows() != A.getColumns()) {
            throw new java.lang.IllegalArgumentException("Only Square matrices"
                    + " can be decomposed!");
        }

        int rows = A.getRows();
        int cols = A.getColumns();

        double[][] LData = new double[rows][cols];
        double[][] UData = new double[rows][cols];

        //The main diagonal of L is ones.
        for (int i = 0; i < rows; i++) {
            LData[i][i] = 1;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (i <= j) {
                    //solve a U
                    UData[i][j] = A.get(i, j);
                    for (int k = 0; k < i; k++) {
                        UData[i][j] -= UData[k][j] * LData[i][k];
                    }
                } else {
                    //solve a L
                    LData[i][j] = A.get(i, j);
                    for (int k = 0; k < j; k++) {
                        LData[i][j] -= UData[k][j] * LData[i][k];
                    }
                    LData[i][j] /= UData[j][j];
                }
            }
        }

        Matrix[] array = new Matrix[2];
        array[0] = new Matrix(LData);
        array[1] = new Matrix(UData);

        return array;

    }

    /**
     * Factorizes the input Matrix into Q, orthogonal Matrix, and R, upper
     * triangular Matrix using Givens Rotations.
     * @param A Matrix to factorize.
     * @return Matrix Q and R in that order in an array.
     */
    public static Matrix[] qr_fact_givens(Matrix A) {

        int rows = A.getRows();
        int cols = A.getColumns();

        Matrix R = new Matrix(A);
        Matrix Q = Matrix.getIdentity(cols);

        for (int j = 0; j < cols; j++) {
            for (int i = rows - 1; i > j; i--) {
               //Make the identity matrix to hold our Rotation Matrix.
                Matrix G = Matrix.getIdentity(cols);

                //Find the values of cos and sin.
                double[] XY = givensRotation(R.get(i - 1, j),
                        R.get(i, j));

//                System.out.println("Dealing with elements "
//                      + R.getElement(i - 1, j)
//                      + ", " + R.getElement(i, j));
//                System.out.println("Cos = " + XY[0] + "\nSin = " + XY[1]);

                //Create the G rotation matrix.
                G.set(i - 1, i - 1, XY[0]);    //cos
                G.set(i - 1, i, -1 * XY[1]);   //-sin
                G.set(i, i - 1, XY[1]);        //sin
                G.set(i, i, XY[0]);            //cos

//                System.out.println("G matrix:\n" + G);

                //Apply G to find Q and R.
                R = G.multiply(R);
                Q = G.multiply(Q);

//                System.out.println("new R matrix:\n" + R);
//                System.out.println("new Q matrix:\n" + Q);
            }
        }

        Matrix[] temp = new Matrix[2];
        temp[0] = Q.transpose();
        temp[1] = R;
        return temp;
    }

    /**
     * Factorizes the input Matrix into Q, orthogonal Matrix, and R, upper
     * triangular Matrix using HouseHolder Reflections.
     * @param A Matrix to factorize.
     * @return Matrix Q and R in that order in an array.
     */
    public static Matrix[] qr_fact_househ(Matrix A) {

        Matrix Q = Matrix.getIdentity(A.getRows());
        Matrix R = new Matrix(A);

//        System.out.println("Initial Matrix Q:\n" + Q);
//        System.out.println("Initial Matrix R:\n" + R);

        int rows = R.getRows();
        int cols = R.getColumns();

        for (int i = 0; i < cols - 1; i++) {
            //Do HOUSEHOLDER!
            Matrix I = Matrix.getIdentity(R.getRows());

            //Get u
            double[] u = R.getColumn(i);
            for (int k = 0; k < i; k++) {
                u[k] = 0;
            }
            u[i] -= getMagnitude(u);

//            System.out.println("U vector is: " + Arrays.toString(u));

            //Get v
            double mag = getMagnitude(u);
            for (int k = 0; k < u.length; k++) {
                u[k] /= mag;
            }

//            System.out.println("V vector is: " + Arrays.toString(u)
//                    + "\n");

            Matrix H = I.subtract(v_Vtranspose(u).multiply(2.0));

//            System.out.println("H matrix is:\n" + H);

            R = H.multiply(R);
            Q = Q.multiply(H);

//            System.out.println("Middle Matrix Q:\n" + Q);
//            System.out.println("Middle Matrix R:\n" + R);
        }

        Matrix[] temp = new Matrix[2];
        temp[0] = Q;
        temp[1] = R;
        return temp;
    }

    /**
     * Solves the equation LUx = b for x.
     * @param L Lower triangular matrix
     * @param U Upper triangular matrix
     * @param b vector.
     * @return solution vector "x".
     */
    public static double[] solve_lu_b(final Matrix L, final Matrix U,
            final double[] b) {

        double[] y = new double[b.length];

        //First solve Ly = b
        for (int i = 0; i < y.length; i++) {
            double temp = b[i];

            for (int k = 0; k < i; k++) {
                temp -= L.get(i, k) * y[k];
            }
            temp /= L.get(i, i);

            y[i] = temp;
        }

//        System.out.println("Y vector:\n" + Arrays.toString(y));
//        double[][] temptemp = new double[y.length][1];
//        for (int i = 0; i < y.length; i++) {
//            temptemp[i][0] = y[i];
//        }
//        System.out.println("LY Matrix:\n" + L.multiply(new Matrix(temptemp)));

        double[] x = new double[b.length];

        //Now solve Ux = y
        for (int i = b.length - 1; i > -1; i--) {
            double temp = y[i];

            for (int k = i + 1; k < b.length; k++) {
                temp -= U.get(i, k) * x[k];
            }

            temp /= U.get(i, i);

            x[i] = temp;
        }

//        System.out.println("X vector:\n" + Arrays.toString(x));
//        double[][] temp = new double[x.length][1];
//        for (int i = 0; i < x.length; i++) {
//            temp[i][0] = x[i];
//        }
//        System.out.println("Ux Matrix:\n" + U.multiply(new Matrix(temp)));

        return x;
    }

    /**
     * Solves the equation QRx = b for x.
     * @param Q Orthogonal Matrix.
     * @param R Upper triangular Matrix.
     * @param b vector.
     * @return solution vector "x".
     */
    public static double[] solve_qr_b(Matrix Q, Matrix R, double[] b) {
        //QRx = b OR Rx = Q^t * b = y
        //Use backward sub for Rx = y

        //Make b into a Matrix.
        double[][] temp = new double[b.length][1];
        for (int i = 0; i < b.length; i++) {
            temp[i][0] = b[i];
        }

        //Get the right side of our equation.
        Matrix y = Q.transpose().multiply(new Matrix(temp));

        //Backwards sub to solve for x.
        double[] x = new double[b.length];

        for (int i = b.length - 1; i > -1; i--) {
            double temp2 = y.get(i, 0);

            for (int k = i + 1; k < b.length; k++) {
                temp2 -= R.get(i, k) * x[k];
            }

            temp2 /= R.get(i, i);

            x[i] = temp2;
        }

        return x;
    }

    /**
     * Finds the maximum norm of a Matrix.
     * @param A Matrix to calculate.
     * @return Maximum Norm of the Matrix.
     */
    public static double maximumNorm(final Matrix A) {
        double error = 0;
        int rows = A.getRows();
        int cols = A.getColumns();

        for (int i = 0; i < rows; i++) {
            double nextError = 0;
            for (int j = 0; j < cols; j++) {
                nextError += Math.abs(A.get(i, j));
            }
            if (nextError > error) {
                error = nextError;
            }
        }

        return error;
    }

    /**
     * Encodes the input stream.
     * @param inputStream
     * @return The encoded stream.
     */
    public static double[] encode(double[] inputStream) {

        Matrix A0 = new Matrix(inputStream.length + 3, inputStream.length + 3);
        double[] input = new double[inputStream.length + 3];

        for (int i = 0; i < inputStream.length; i++) {
            input[i] = inputStream[i];
        }


        for (int i = 0; i < A0.getRows(); i++) {
            A0.set(i, i, 1);
            if (i - 2 < 0) {
                A0.set(i, i + A0.getColumns() - 2, 1);
            } else {
                A0.set(i, i - 2, 1);
            }
            if (i - 3 < 0) {
                A0.set(i, i + A0.getColumns() - 3, 1);
            } else {
                A0.set(i, i - 3, 1);
            }
        }

        Matrix A1 = new Matrix(inputStream.length + 3, inputStream.length + 3);

        for (int i = 0; i < A1.getRows(); i++) {
            A1.set(i, i, 1);
            if (i - 1 < 0) {
                A1.set(i, i + A1.getColumns() - 1, 1);
            } else {
                A1.set(i, i - 1, 1);
            }
            if (i - 3 < 0) {
                A1.set(i, i + A1.getColumns() - 3, 1);
            } else {
                A1.set(i, i - 3, 1);
            }
        }

        double[][] temp = new double[input.length][1];
        for (int i = 0; i < input.length; i++) {
            temp[i][0] = input[i];
        }
        Matrix x = new Matrix(temp);

        Matrix y0 = A0.multiply(x);
        Matrix y1 = A1.multiply(x);

        double[] encode = new double[y0.getRows()];
        for (int i = 0; i < y0.getRows(); i++) {
            encode[i] = 10 * (y0.get(i, 0) % 2) + y1.get(i, 0) % 2;
        }

        return encode;
    }

    public static double[][] power_method(Matrix A, double tol,
            double[] initial) {
        double[][] returnData = new double[3][];
        returnData[1] = new double[1]; //EigenValue
        returnData[2] = new double[1]; //Number of iterations needed.

        //Turn input double[] into a Matrix.
        double[][] temp = new double[initial.length][1];
        for (int i = 0; i < initial.length; i++) {
            temp[i][0] = initial[i];
        }
        Matrix x = new Matrix(temp);
        Matrix xOld;
        int iterations = 0;
        do {
            iterations++;
            xOld = new Matrix(x);

            Matrix temp2 = A.multiply(x);
            x = temp2.multiply(1.0 / xOld.get(0, 0));

        } while(iterations < 5000
                && Math.abs(x.get(0, 0) - xOld.get(0, 0)) > tol);

        if (iterations == 5000) {
            System.out.println("Power method could not converge within"
                    + " the given tolerance in <5000 iterations");
            return null; //If Power method does not converge, return null.
        }

        returnData[2][0] = iterations;
        returnData[1][0] = x.get(0, 0);
        returnData[0] = x.getColumn(0);

        return returnData;
    }

    /**
     * Gets the magnitude of the input vector.
     * @param vector Vector to calculate.
     * @return Magnitude of the vector.
     */
    private static double getMagnitude(final double[] vector) {
        double mag = 0;

        for (double d : vector) {
            mag += d * d;
        }

        return java.lang.Math.sqrt(mag);
    }


    /**
     * Multiplies the input vector v by its transpose, and returns the
     * resultant matrix.
     * @param v Vector to multiply.
     * @return Resultant Matrix.
     */
    private static Matrix v_Vtranspose(final double[] v) {
        Matrix mat = new Matrix(v.length, v.length);

        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v.length; j++) {
                mat.set(i, j, v[i] * v[j]);
            }
        }

//        System.out.println("V * V^t :\n" + mat);

        return mat;
    }

    /**
     * Finds and returns the cos and sin needed for a givens rotation.
     * @param x Element to not zero out?
     * @param y Element to zero out.
     * @return cos and sin in an array in that order.
     */
    private static double[] givensRotation(double x, double y) {
        double[] values = new double[2];

        if (y == 0) {
            values[0] = 1;
            values[1] = 0;
            return values;
        } else {
            double mag = Math.sqrt(x * x + y * y);
            values[0] = x / mag;
            values[1] = -y / mag;
        }
        return values;
    }

    public static double[][] jacobi(Matrix mat, double[] b, double[] initial,
        double tol) {

        //returns tolerance & x
        //make lu and d matrices
        Matrix lu = findL(mat).add(findU(mat));
        Matrix d = findD(mat);
        Matrix dInv = new Matrix(d);


        //inverse of d
        for (int i = 0; i < mat.getRows(); i++) {
            for (int j = 0; j < mat.getColumns(); j++) {
                if (i == j) {
                    if (dInv.get(i, j) != 0) {
                        dInv.set(i, j, -1 / dInv.get(i, j));
                    }
                }
            }
        }

        //System.out.println("Inverse of D: " + dInv);

        //make b a matrix
        double[][] temp = new double[b.length][1];
        for (int i = 0; i < b.length; i++) {
            temp[i][0] = b[i];
        }

        Matrix bo = new Matrix(temp);

        temp = new double[b.length][1];
        for (int i = 0; i < b.length; i++) {
            temp[i][0] = initial[i];
        }
        Matrix xo = new Matrix(temp);

        //right side of formula
        int count = 0;
        double tolTemp;
        Matrix x;
        do {
            x = dInv.multiply((lu.multiply(xo)).add(bo));
            tolTemp = getMagnitude(x.subtract(xo).getColumn(0));
            count++;
            xo = new Matrix(x);
            System.out.println(x);
        } while (count < 500 && Math.abs(tolTemp) > tol);

        if (count == 500) {
            System.out.println("Solution did not converge within the tolerance"
                    + " in 500 iterations.");
            return null;
        }

        double[] iteration = new double[1];
        iteration[0] = count;
        double[][] returnVals = {x.getColumn(0), iteration};

        return returnVals;
    }

    public static double[][] gauss_seidel(Matrix mat, double[] b,
            double[] initial, double tol) {
        Matrix ld = findL(mat).add(findD(mat));
        Matrix u = findU(mat);

        //make b a matrix
        double[][] temp2 = new double[b.length][1];
        for (int i = 0; i < b.length; i++) {
            temp2[i][0] = b[i];
        }
        Matrix bo = new Matrix(temp2);

        temp2 = new double[b.length][1];
        for (int i = 0; i < b.length; i++) {
            temp2[i][0] = initial[i];
        }
        Matrix xo = new Matrix(temp2);

        Matrix right;

        double[] x = new double[b.length];
        int count = 0;
        Matrix x1;
        double tolTemp;
        do {
            //Solve (L+D)x = right
            right = (u.multiply(xo).add(bo).multiply(-1));
            count++;
            double temp;
            for (int i = 0; i < x.length; i++) {
                temp = right.get(i, 0);

                for (int k = 0; k < i; k++) {
                    temp -= ld.get(i, k) * x[k];
                }
                temp /= ld.get(i, i);
                x[i] = temp;
            }

            double[][] temp3 = new double[x.length][1];
            for (int i = 0; i < x.length; i++) {
                temp3[i][0] = x[i];
            }

            x1 = new Matrix(temp3);

            tolTemp = getMagnitude(x1.subtract(xo).getColumn(0));
            count++;
            xo = x1;
        } while (count < 500 && tolTemp > tol);

        if (count == 500) {
            System.out.println("Solution did not converge within the tolerance"
                    + " in 500 iterations.");
            return null;
        }
        //make matrix x into an array
        double[][] xFinal = new double[1][x1.getRows()];
        for (int i = 0; i < x1.getRows(); i++) {
            xFinal[0][i] = x1.get(i, 0);
        }
        double[] iteration = new double[1];
        iteration[0] = count;
        double[][] returnVals = {xFinal[0], {count}};

        return returnVals;

    }

    private static Matrix findL(Matrix mat) {
        int row = mat.getRows();
        int col = mat.getColumns();
        Matrix l = new Matrix(mat);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                if (i <= j) {
                    l.set(i, j, 0);
                }
            }
        }

        //System.out.println("L Matrix: " + l);

        return l;
    }

    private static Matrix findU(Matrix mat) {
        int row = mat.getRows();
        int col = mat.getColumns();
        Matrix u = new Matrix(mat);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                if (i >= j) {
                    u.set(i, j, 0);
                }
            }
        }

        //System.out.println("U Matrix: " + u);

        return u;
    }

    private static Matrix findD(Matrix mat) {
        int row = mat.getRows();
        int col = mat.getColumns();
        Matrix d = new Matrix(mat.getRows(), mat.getColumns());

        for (int i = 0; i < row; i++) {
            for(int j = 0; j <col; j++) {
                if (i == j) {
                    d.set(i, j, mat.get(i, j));
                }
            }
        }

        //System.out.println("D Matrix: " + d);

        return d;
    }

}
