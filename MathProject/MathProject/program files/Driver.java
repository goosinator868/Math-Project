
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Scanner;
import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;


/**
 *
 * @author Ben Bohannon and Sydney Young
 */
public final class Driver {

    //TODO: Write code/get data for part 3.
    //TODO: Fix/Test Gauss-Seidel and jacobi methods.

    //TODO: Do writing component for part 1, part 2, part 3.

    private static BufferedReader reader;

    /**
     * Private constructor to prevent instantiation.
     */
    private Driver() {
        //Nobody instantiate a Driver!
    }

    public static void main(String[] args) {
        reader = new BufferedReader(new InputStreamReader(System.in));


        System.out.println("Welcome to the Matrix Calculator!\n"
                + "What would you like to do?\n");

        boolean isFinished = false;
        while (!isFinished) {
            System.out.println("1) LU Decomposition\n"
                    + "2) QR Decomposition\n"
                    + "3) Solve Ax = b\n"
                    + "4) Power Method\n"
                    + "5) Hilbert Matrix Solutions and Error\n"
                    + "6) Encoder/Decoder\n"
                    + "7) Exit\n");

            int choice = getIntInput(1, 7);

            switch (choice) {
                case 1: //LU DECOMPOSITION------------------------------------
                    System.out.println("Input name of file containing the"
                            + " Matrix to factorize.");
                    System.out.println("**Note: The file must be in the"
                            + " same directory as the Driver.class file**");

                    try {
                        String path = reader.readLine();

                        Matrix A = parser(path);

                        Matrix[] LU = MatrixMath.lu_fact(A);

                        System.out.println("L Matrix:\n" + LU[0]);
                        System.out.println("U Matrix:\n" + LU[1]);

                        System.out.println("\n LU Error: "
                                + MatrixMath.maximumNorm(
                                        LU[0].multiply(LU[1]).subtract(A)));

                    } catch (Exception e) {
                        System.out.println("Failed to read file.");
                        System.out.println(e.toString());
                    }

                    break;

                case 2: //QR DECOMPOSITION-------------------------------------
                    System.out.println("1) Use Givens\n"
                            + "2) Use HouseHolder");

                    int choice2 = getIntInput(1, 2);

                    switch(choice2) {
                        case 1: //GIVENS
                            System.out.println("Input name of file containing"
                                    + "the Matrix to decompose.");
                            System.out.println("**Note: The file must be in the"
                                    + " same directory as the Driver.class file"
                                    + "**");

                            try {
                                String path = reader.readLine();

                                Matrix A = parser(path);

                                Matrix[] QR = MatrixMath.qr_fact_givens(A);

                                System.out.println("Q Matrix:\n" + QR[0]);
                                System.out.println("R Matrix:\n" + QR[1]);

                                System.out.println("\n QR Error: "
                                + MatrixMath.maximumNorm(
                                        QR[0].multiply(QR[1]).subtract(A)));

                            } catch (Exception e) {
                                System.out.println("Failed to read file.");
                            }
                            break;

                        case 2: //HOUSEHOLDER
                            System.out.println("Input name of file containing"
                                    + "the Matrix to decompose.");
                            System.out.println("**Note: The file must be in the"
                                    + " same directory as the Driver.class file"
                                    + "**");

                            try {
                                String path = reader.readLine();

                                Matrix A = parser(path);

                                Matrix[] QR = MatrixMath.qr_fact_givens(A);

                                System.out.println("Q Matrix:\n" + QR[0]);
                                System.out.println("R Matrix:\n" + QR[1]);

                                System.out.println("\n QR Error: "
                                + MatrixMath.maximumNorm(
                                        QR[0].multiply(QR[1]).subtract(A)));

                            } catch (Exception e) {
                                System.out.println("Failed to read file.");
                            }
                            break;
                        default:
                            System.out.println("Something went wrong.");
                    }
                    break;

                case 3: //Ax = b--------------------------------------------
                    System.out.println("1) Solve using LU Decomposition\n"
                            + "2) Solve using QR Decomposition\n"
                            + "3) Solve using Jacobi Method\n"
                            + "4) Solve using Gauss-Seidel Method\n"
                            + "5) Go Back");
                    System.out.println("NOTE: These methods need an augmented "
                            + "matrix in order to function.\n");

                    int choice3 = getIntInput(1, 5);

                    Matrix Ab;
                    double[] b;
                    Matrix A;
                    double[] solution;

                    switch (choice3) {
                        case 1: //LU SOLVE
                            System.out.println("Input name of file containing"
                                    + "the augmented Matrix to decompose.");
                            System.out.println("**Note: The file must be in the"
                                    + " same directory as the Driver.class file"
                                    + "**");

                            try {
                                String path = reader.readLine();
                                Ab = parser(path);

                                b = Ab.getColumn(Ab.getColumns() - 1);
                                A = new Matrix(Ab.getRows(),
                                        Ab.getColumns() - 1);
                                for(int i = 0; i < Ab.getRows(); i++) {
                                    for (int j = 0; j < Ab.getColumns() - 1; j++) {
                                        A.set(i, j, Ab.get(i, j));
                                    }
                                }

                                Matrix[] LU = MatrixMath.lu_fact(A);

                                solution = MatrixMath.solve_lu_b(LU[0], LU[1], b);

                                //Convert solution to a Matrix.
                                double[][] temp = new double[solution.length][1];
                                for (int i = 0; i < solution.length; i++) {
                                    temp[i][0] = solution[i];
                                }
                                Matrix sol = new Matrix(temp);
                                //Convert b vector to Matrix.
                                temp = new double[b.length][1];
                                for (int i = 0; i < b.length; i++) {
                                    temp[i][0] = b[i];
                                }
                                Matrix bMat = new Matrix(temp);

                                //Print solution
                                System.out.println("Solution (transpose):\n"
                                        + Arrays.toString(solution));

                                //FIND ERROR and PRINT IT.
                                System.out.println("Solution error: "
                                    + MatrixMath.getMagnitude(
                                            A.multiply(sol)
                                             .subtract(bMat)
                                             .getColumn(0)));

                            } catch (Exception e) {
                                System.out.println("Failed to read file.");
                            }

                            break;
                        case 2: //QR SOLVE
                            System.out.println("Input name of file containing"
                                    + "the augmented Matrix to decompose.");
                            System.out.println("**Note: The file must be in the"
                                    + " same directory as the Driver.class file"
                                    + "**");

                            try {
                                String path = reader.readLine();

                                Ab = parser(path);

                                b = Ab.getColumn(Ab.getColumns() - 1);
                                A = new Matrix(Ab.getRows(),
                                        Ab.getColumns() - 1);
                                for(int i = 0; i < Ab.getRows(); i++) {
                                    for (int j = 0; j < Ab.getColumns() - 1;
                                            j++) {
                                        A.set(i, j, Ab.get(i, j));
                                    }
                                }

                                Matrix[] QR = MatrixMath.qr_fact_givens(A);

                                solution = MatrixMath.solve_qr_b(QR[0],
                                        QR[1], b);

                                //Convert solution to a Matrix.
                                double[][] temp = new double[solution.length][1];
                                for (int i = 0; i < solution.length; i++) {
                                    temp[i][0] = solution[i];
                                }
                                Matrix sol = new Matrix(temp);
                                //Convert b vector to Matrix.
                                temp = new double[b.length][1];
                                for (int i = 0; i < b.length; i++) {
                                    temp[i][0] = b[i];
                                }
                                Matrix bMat = new Matrix(temp);

                                //Print solution
                                System.out.println("Solution (transpose):\n"
                                        + Arrays.toString(solution));

                                //FIND ERROR and PRINT IT.
                                System.out.println("Solution error: "
                                    + MatrixMath.getMagnitude(
                                            A.multiply(sol)
                                             .subtract(bMat)
                                             .getColumn(0)));

                            } catch (Exception e) {
                                System.out.println("Failed to Read File.");
                            }
                            break;
                        case 3: //JACOBI SOLVE
                            System.out.println("Input name of file containing"
                                    + "the augmented Matrix to decompose.");
                            System.out.println("**Note: The file must be in the"
                                    + " same directory as the Driver.class file"
                                    + "**");

                            try {
                                String path = reader.readLine();

                                Ab = parser(path);

                                b = Ab.getColumn(Ab.getColumns() - 1);
                                A = new Matrix(Ab.getRows(),
                                        Ab.getColumns() - 1);
                                for (int i = 0; i < Ab.getRows(); i++) {
                                    for (int j = 0; j < Ab.getColumns() - 1;
                                            j++) {
                                        A.set(i, j, Ab.get(i, j));
                                    }
                                }
                                System.out.println("Input initial "
                                        + "guess vector");
                                double[] initial = vectorParser(A.getRows());

                                System.out.println("Input error tolerance.");
                                double tol = getDoubleInput(false);


                                double[][] info = MatrixMath.jacobi(A, b,
                                        initial, tol, false);

                                if (info != null) {
                                    solution = info[0];
                                    //Convert solution to a Matrix.
                                    double[][] temp
                                            = new double[solution.length][1];
                                    for (int i = 0; i < solution.length; i++) {
                                        temp[i][0] = solution[i];
                                    }
                                    Matrix sol = new Matrix(temp);
                                    //Convert b vector to Matrix.
                                    temp = new double[b.length][1];
                                    for (int i = 0; i < b.length; i++) {
                                        temp[i][0] = b[i];
                                    }
                                    Matrix bMat = new Matrix(temp);

                                    //Print solution
                                    System.out.println("Solution (transpose):\n"
                                            + Arrays.toString(solution));

                                    //FIND ERROR and PRINT IT.
                                    System.out.println("Solution error: "
                                        + MatrixMath.getMagnitude(
                                                A.multiply(sol)
                                                 .subtract(bMat)
                                                 .getColumn(0)));
                                }

                            } catch (Exception e) {
                                System.out.println("Failed to Read File.");
                            }
                            break;
                        case 4: //GAUSS-SEIDEL SOLVE
                            System.out.println("Input name of file containing"
                                    + "the augmented Matrix to decompose.");
                            System.out.println("**Note: The file must be in the"
                                    + " same directory as the Driver.class file"
                                    + "**");

                            try {
                                String path = reader.readLine();

                                Ab = parser(path);

                                b = Ab.getColumn(Ab.getColumns() - 1);
                                A = new Matrix(Ab.getRows(),
                                        Ab.getColumns() - 1);
                                for (int i = 0; i < Ab.getRows(); i++) {
                                    for (int j = 0; j < Ab.getColumns() - 1;
                                            j++) {
                                        A.set(i, j, Ab.get(i, j));
                                    }
                                }
                                System.out.println("Input initial "
                                        + "guess vector");
                                double[] initial = vectorParser(A.getRows());

                                System.out.println("Input error tolerance.");
                                double tol = getDoubleInput(false);


                                double[][] info = MatrixMath.gauss_seidel(A, b,
                                        initial, tol, false);

                                if (info != null) {
                                    solution = info[0];
                                    //Convert solution to a Matrix.
                                    double[][] temp = new double[solution.length][1];
                                    for (int i = 0; i < solution.length; i++) {
                                        temp[i][0] = solution[i];
                                    }
                                    Matrix sol = new Matrix(temp);
                                    //Convert b vector to Matrix.
                                    temp = new double[b.length][1];
                                    for (int i = 0; i < b.length; i++) {
                                        temp[i][0] = b[i];
                                    }
                                    Matrix bMat = new Matrix(temp);

                                    //Print solution
                                    System.out.println("Solution (transpose):\n"
                                            + Arrays.toString(solution));

                                    //FIND ERROR and PRINT IT.
                                    System.out.println("Solution error: "
                                        + MatrixMath.getMagnitude(
                                                A.multiply(sol)
                                                 .subtract(bMat)
                                                 .getColumn(0)));
                                }
                            } catch (Exception e) {
                                System.out.println("Failed to Read File.");
                            }

                            break;
                        case 5:
                            break;
                        default:
                            System.out.println("Something went wrong.");
                    }
                    break;
                case 4: //POWER METHOD--------------------------------------
                    System.out.println("Input name of file containing the"
                            + " Matrix to calculate.");
                    System.out.println("**Note: The file must be in the"
                            + " same directory as the Driver.class file**");

                    try {
                        String path = reader.readLine();

                        A = parser(path);

                        System.out.println("Input initial "
                                        + "eigenvector guess");
                        double[] initial = vectorParser(A.getRows());

                        System.out.println("Input error tolerance.");
                        double tol = getDoubleInput(false);

                        double[][] powerStuff
                                = MatrixMath.power_method(A, tol, initial);

                        if (powerStuff != null) {
                            System.out.println("Eigenvector (transpose):"
                                    + Arrays.toString(powerStuff[0]));
                            System.out.println("Eigenvalue: "
                                    + powerStuff[1][0]);
                            System.out.println("Number of iterations "
                                    + powerStuff[2][0]);
                        }


                    } catch (Exception e) {
                        System.out.println("Failed to read file.");
                    }
                    break;
                case 5: //HILBERT MATRIX THINGS-----------------------
                    doHilbert();
                    break;
                case 6: //ENCODER DECODER-----------------------------

                    System.out.println("1) Use Jacobi\n"
                            + "2) Use Gauss Seidel\n"
                            + "3) Go Back");

                    int choice4 = getIntInput(1, 3);

                    switch (choice4) {
                        case 1:
                            jacobiEncode();
                            break;
                        case 2:
                            gaussSeidelEncode();
                            break;
                        case 3:
                            break;
                        default:
                            System.out.println("Something went wrong.");
                            break;
                    }

                    break;
                case 7: //QUIT----------------------------------------
                    isFinished = true;
                    break;
                default:
                    System.out.println("Something went wrong.");
                    break;
            }
        }

        System.out.println("Closing.");

    }

    private static boolean doHilbert() {
        System.out.println("\nWhat would you like to do with"
                + " Hilbert Matrices?\n");

        System.out.println(
                "1) Solve Hx = b for n = 2,3...20? (Problem 1.d)\n"
                + "2) Create txt file of errors.\n"
                + "3) Go Back.\n");

        int choice = getIntInput(1, 3);

        int n;
        switch (choice) {
                case 1:
                    System.out.println("Input n, max size of Hilbert Matrix");
                    n = getIntInput(2, Integer.MAX_VALUE - 1);
                    writeHilbertData(n + 1);
                    return false;
                case 2:
                    System.out.println("Input n, max size of Hilbert Matrix");
                    n = getIntInput(2, Integer.MAX_VALUE);
                    writeErrorData(n + 1);
                    return false;
                case 3:
                    return false;
                default:
                    System.out.println("Something went wrong.");
        }

        return false;
    }

    private static void writeHilbertData(int n) {
        try {
            PrintWriter writer = new PrintWriter("HilbertOutput.txt", "UTF-8");

            System.out.print("Writing..");

            for (int i = 2; i < n; i++) {
                //Print Hilbert Matrix.
                Matrix h = Matrix.getHilbert(i);
                writer.println("----- " + i + "x" + i
                        + " Hilbert Matrix -------\n");
                writer.println("Matrix:\n" + h);

                //Print b vector.
                double[] b = new double[i];
                double value = Math.pow(0.1, i / 3.0);
                for (int j = 0; j < i; j++) {
                    b[j] = value;
                }
                writer.println("b vector:\n" + Arrays.toString(b));

                //Do the calculations.
                Matrix[] lu = MatrixMath.lu_fact(h);
                double[] x = MatrixMath.solve_lu_b(lu[0], lu[1], b);

                //Print x vector
                writer.println("\nSolution x vector:\n" + Arrays.toString(x));

                //Print LU Decomp Error.
                writer.println("\nLU error: "
                        + MatrixMath.maximumNorm(
                                lu[0].multiply(lu[1]).subtract(h)));

                //Print X vector Error.
                double[][] temp = new double[x.length][1];
                double[][] temp2 = new double[b.length][1];
                for (int j = 0; j < x.length; j++) {
                    temp[j][0] = x[j];
                    temp2[j][0] = b[j];
                }
                writer.println("Solution error: "
                        + MatrixMath.getMagnitude(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2)).getColumn(0))
                        + "\n");
            }

            writer.close();
            System.out.println("Done writing! Output to HilbertOutput.txt");
        } catch (Exception ex) {
            System.out.println("Writing failed!\n" + ex.toString());
        }
    }

    private static void gaussSeidelEncode() {
        try {
            System.out.println("Input length of random binary stream.");

            int len = getIntInput(1, Integer.MAX_VALUE);

            double[] x = xGenerator(len);

            System.out.println("Binary stream:\n" + Arrays.toString(x));

            System.out.println("Press enter to encode.");

            reader.readLine();

            double[] y = MatrixMath.encode(x);

            System.out.println("Encoded stream:\n"
                    + Arrays.toString(y));

            System.out.println("Press enter to decode.");

            reader.readLine();

            double[] initial = new double[y.length];
            for (int i = 0; i < initial.length; i++) {
                initial[i] = 1;
            }

            for (int i = 0; i < y.length; i++) {
                y[i] = (int) y[i] / 10; //Get y0
            }

            Matrix A0 = new Matrix(y.length, y.length);
            for (int i = 0; i < A0.getRows(); i++) {
                A0.set(i, i, 1);
                if (i - 2 > 0) {
                    A0.set(i, i - 2, 1);
                }
                if (i - 3 > 0) {
                    A0.set(i, i - 3, 1);
                }
            }

            double[][] solutionStuff
                    = MatrixMath.gauss_seidel(A0, y,
                            initial, 1e-9, true);

            double[] sol = new double[solutionStuff[0].length - 3];
            for (int i = 0; i < sol.length; i++) {
                sol[i] = solutionStuff[0][i];
            }

            System.out.println("Solution stream:\n"
                    + Arrays.toString(sol));
            System.out.println("Number of iterations: "
                    + solutionStuff[1][0]);

        } catch (Exception e) {
            System.out.println("Error encoding/decoding.");
        }
    }

    private static void jacobiEncode() {
        try {
            System.out.println("Input length of random binary stream.");

            int len = getIntInput(1, Integer.MAX_VALUE);

            double[] x = xGenerator(len);

            System.out.println("Binary stream:\n" + Arrays.toString(x));

            System.out.println("Press enter to encode.");

            reader.readLine();

            double[] y = MatrixMath.encode(x);

            System.out.println("Encoded stream:\n"
                    + Arrays.toString(y));

            System.out.println("Press enter to decode.");

            reader.readLine();

            double[] initial = new double[y.length];
            for (int i = 0; i < initial.length; i++) {
                initial[i] = 1;
            }

            for (int i = 0; i < y.length; i++) {
                y[i] = (int) y[i] / 10; //Get y0
            }

            Matrix A0 = new Matrix(y.length, y.length);
            for (int i = 0; i < A0.getRows(); i++) {
                A0.set(i, i, 1);
                if (i - 2 > 0) {
                    A0.set(i, i - 2, 1);
                }
                if (i - 3 > 0) {
                    A0.set(i, i - 3, 1);
                }
            }

            double[][] solutionStuff
                    = MatrixMath.jacobi(A0, y,
                            initial, 1e-9, true);

            double[] sol = new double[solutionStuff[0].length - 3];
            for (int i = 0; i < sol.length; i++) {
                sol[i] = solutionStuff[0][i];
            }

            System.out.println("Solution stream:\n"
                    + Arrays.toString(sol));
            System.out.println("Number of iterations: "
                    + solutionStuff[1][0]);

        } catch (Exception e) {
            System.out.println("Error encoding/decoding.");
        }
    }

    private static void writeErrorData(int n) {
        try {
            final int iterations = n;
            PrintWriter writer = new PrintWriter("ErrorOutput.txt", "UTF-8");

            System.out.print("Writing..");

            writer.println("LU error for increasing n:");
            for (int i = 2; i < iterations; i++) {
                //Get Hilbert Matrix.
                Matrix h = Matrix.getHilbert(i);

                //Do calculations.
                Matrix[] lu = MatrixMath.lu_fact(h);

                //Print LU Decomp Error.
                writer.println(MatrixMath.maximumNorm(
                        lu[0].multiply(lu[1]).subtract(h)));
            }
            writer.println("\nLU solution error for increasing n:");
            for (int i = 2; i < iterations; i++) {
                //Get Hilbert Matrix.
                Matrix h = Matrix.getHilbert(i);

                double[] b = new double[i];
                double value = Math.pow(0.1, i / 3.0);
                for (int j = 0; j < i; j++) {
                    b[j] = value;
                }

                //Do calculations.
                Matrix[] lu = MatrixMath.lu_fact(h);
                double[] x = MatrixMath.solve_lu_b(lu[0], lu[1], b);

                double[][] temp = new double[x.length][1];
                double[][] temp2 = new double[b.length][1];
                for (int j = 0; j < x.length; j++) {
                    temp[j][0] = x[j];
                    temp2[j][0] = b[j];
                }
                writer.println(MatrixMath.getMagnitude(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2)).getColumn(0)));
            }

            writer.println("\nQR Givens error for increasing n:");
            for (int i = 2; i < iterations; i++) {
                //Get Hilbert Matrix.
                Matrix h = Matrix.getHilbert(i);

                //Do calculations.
                Matrix[] qrGivens = MatrixMath.qr_fact_givens(h);

                //Print QR-Givens Error
                writer.println(MatrixMath.maximumNorm(
                        qrGivens[0].multiply(qrGivens[1]).subtract(h)));
            }
            writer.println("\nQR Givens solution error for increasing n:");
            for (int i = 2; i < iterations; i++) {
                //Get Hilbert Matrix.
                Matrix h = Matrix.getHilbert(i);

                double[] b = new double[i];
                double value = Math.pow(0.1, i / 3.0);
                for (int j = 0; j < i; j++) {
                    b[j] = value;
                }

                //Do calculations.
                Matrix[] lu = MatrixMath.qr_fact_givens(h);
                double[] x = MatrixMath.solve_qr_b(lu[0], lu[1], b);

                double[][] temp = new double[x.length][1];
                double[][] temp2 = new double[b.length][1];
                for (int j = 0; j < x.length; j++) {
                    temp[j][0] = x[j];
                    temp2[j][0] = b[j];
                }
                writer.println(MatrixMath.getMagnitude(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2)).getColumn(0)));
            }

            writer.println("\nQR householder error for increasing n:");
            for (int i = 2; i < iterations; i++) {
                //Get Hilbert Matrix.
                Matrix h = Matrix.getHilbert(i);

                //Do calculations.
                Matrix[] qrHouse = MatrixMath.qr_fact_househ(h);

                //Print QR-Householder Error
                writer.println(MatrixMath.maximumNorm(
                        qrHouse[0].multiply(qrHouse[1]).subtract(h)));
            }
            writer.println("\nQR HouseHolder solution error for increasing n:");
            for (int i = 2; i < iterations; i++) {
                //Get Hilbert Matrix.
                Matrix h = Matrix.getHilbert(i);

                double[] b = new double[i];
                double value = Math.pow(0.1, i / 3.0);
                for (int j = 0; j < i; j++) {
                    b[j] = value;
                }

                //Do calculations.
                Matrix[] lu = MatrixMath.qr_fact_househ(h);
                double[] x = MatrixMath.solve_qr_b(lu[0], lu[1], b);

                double[][] temp = new double[x.length][1];
                double[][] temp2 = new double[b.length][1];
                for (int j = 0; j < x.length; j++) {
                    temp[j][0] = x[j];
                    temp2[j][0] = b[j];
                }
                writer.println(MatrixMath.getMagnitude(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2)).getColumn(0)));
            }

            writer.close();
            System.out.println("Done writing! Output to ErrorOutput.txt");
        } catch (Exception ex) {
            System.out.println("Writing failed!\n" + ex.toString());
        }
    }

    private static Matrix parser(String fileName) {
        Matrix m = null;
        try {
            File file = new File(fileName);
            int r = 0;
            ArrayList<ArrayList> mat3 = new ArrayList<ArrayList>();
            double[] mat2;
            Scanner kb = new Scanner(file);
            while (kb.hasNextLine()) {
                String line = kb.nextLine();
                ArrayList<Double> mat = new ArrayList<Double>();
                while (line.indexOf(' ') != -1) {
                    double num = Double.parseDouble(line.substring(0,
                            line.indexOf(' ')));
                    line = line.substring(line.indexOf(' ') + 1, line.length());
                    mat.add(num);
                }
                double num = Double.parseDouble(line.substring(0,
                        line.length()));
                mat.add(num);
                //mat2 = mat.toArray();
                mat3.add(mat);
                r++;
            }
            double[][] mat4 = new double[r][];
            for (int i = 0; i < r; i++) {
                ArrayList<Double> mat = mat3.get(i);
                mat4[i] = new double[mat.size()];
                for (int j = 0; j < mat.size(); j++) {
                    mat4[i][j] = mat.get(j);
                }
            }
            m = new Matrix(mat4);
        } catch (FileNotFoundException e) {
            System.out.println("File not found \n" + e.toString());
        }
        return m;
    }

    private static double[] vectorParser(int size) {

        double[] vector = new double[size];
        int current = 0;
        while (current < size) {
            try {
                System.out.println("Input " + (size - current) + " numbers,"
                        + " divided by spaces.");
                String input = reader.readLine();

                int i = 0;
                String num = "";
                while (i < input.length() && current < size) {
                    char nextChar = input.charAt(i);
                    if (nextChar == '.' || Character.isDigit(nextChar)) {
                        num += nextChar;
                    } else if (nextChar == ' ' || nextChar == '\n') {
                        vector[current++] = Double.parseDouble(num);
                        num = "";
                    }
                    i++;
                }
                if (!num.equals("")) {
                    vector[current++] = Double.parseDouble(num);
                    num = "";
                }
            } catch (Exception e) {
                System.out.println("Invalid data entered.");
            }
        }

        return vector;

//        //[8, 3, 9, 5, 1]
//        ArrayList<Double> mat = new ArrayList<Double>();
//        input = input.substring(input.indexOf("[") + 1, input.length() - 1);
//        while (input.indexOf(",") != -1) {
//            System.out.println(input.indexOf(","));
//            mat.add(Double.parseDouble(input.substring(0, input.indexOf(","))));
//            input = input.substring(input.indexOf(",") + 2, input.length());
//        }
//        mat.add(Double.parseDouble(input));
//        double[] v = new double[mat.size()];
//        for(int i = 0; i < mat.size(); i++) {
//            v[i] = mat.get(i);
//        }
//        return v;
    }

    public static double[] xGenerator(int n) {
        double[] x = new double[n];
        for (int i = 0; i < x.length; i++) {
            x[i] = Math.round(Math.random());
        }
        return x;
    }

    private static int getIntInput(int lower, int upper) {
        boolean validInput = false;
            String input;
            int choice = 0;
            while (!validInput) {
                try {
                    input = reader.readLine();
                    choice = Integer.parseInt(input);
                    if (choice > upper || choice < lower) {
                        throw new RuntimeException();
                    }
                    validInput = true;
                } catch(Exception e) {
                    System.out.println("Invalid input!");
                }
            }
        return choice;
    }

    private static double getDoubleInput(boolean allowNegative) {
        String input;
        double data = 0;
        boolean validInput = false;
        while (!validInput) {
            try {
                input = reader.readLine();
                data = Double.parseDouble(input);
                if (!allowNegative && data < 0) {
                    throw new Exception("No negatives");
                }
                validInput = true;
            } catch (Exception e) {
                System.out.println("Invalid input!");
            }
        }
        return data;
    }


}