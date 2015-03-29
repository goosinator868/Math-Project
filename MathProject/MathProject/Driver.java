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
 * @author Ben Bohannon
 */
public final class Driver {

    /**
     * Private constructor to prevent instantiation.
     */
    private Driver() {
        //Nobody instantiate a Driver!
    }

    public static void main(String[] args) {
        System.out.println("Welcome to the Matrix Calculator!\n"
                + "What would you like to do?");

        boolean isFinished = false;
        while (!isFinished) {
            System.out.println("\n1) Manually enter a Matrix.\n"
                    + "2) Read Matrices from a txt file.\n"
                    + "3) Use Hilbert Matrices.\n"
                    + "4) Exit.\n");

            BufferedReader reader =
                    new BufferedReader(new InputStreamReader(System.in));

            boolean validInput = false;
            String input = "";
            int choice = 0;
            while (!validInput) {
                try {
                    input = reader.readLine();
                    choice = Integer.parseInt(input);
                    if (choice > 4 || choice < 0) {
                        throw new RuntimeException();
                    }
                    validInput = true;
                } catch(Exception e) {
                    System.out.println("Invalid input!");
                }
            }

            switch (choice) {
                case 1:
                    System.out.print("[0, 1, 2, 3]");
                    System.out.println(Arrays.toString(inputParser("[0, 1, 2, 3]")));
                    break;
                case 2:
                    //System.out.println("Not available yet!");
                    System.out.println(parser(input));
                    break;
                case 3:
                    isFinished = doHilbert(reader);
                    break;
                case 4:
                    isFinished = true;
                    break;
                default:
                    System.out.println("Something went wrong.");
                    break;
            }
        }

        System.out.println("Closing.");

    }

    private static boolean doHilbert(BufferedReader reader) {
        System.out.println("\nWhat would you like to do with"
                + " Hilbert Matrices?\n");

        System.out.println(
                "1) Solve Hx = b for n = 2,3...20? (Problem 1.d)\n"
                + "2) Create txt file of errors.\n"
                + "3) Go Back.\n");

        boolean validInput = false;
        String input;
        int choice = 0;
        while (!validInput) {
            try {
                input = reader.readLine();
                choice = Integer.parseInt(input);
                if (choice > 3 || choice < 0) {
                    throw new RuntimeException();
                }
                validInput = true;
            } catch(Exception e) {
                System.out.println("Invalid input!");
            }
        }

        switch (choice) {
                case 1:
                    writeHilbertData();
                    return false;
                case 2:
                    writeErrorData();
                    return false;
                case 3:
                    return false;
                default:
                    System.out.println("Something went wrong.");
        }

        return false;
    }

    private static void writeHilbertData() {
        try {
            PrintWriter writer = new PrintWriter("HilbertOutput.txt", "UTF-8");

            System.out.print("Writing..");

            for (int i = 2; i < 21; i++) {
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
                        + MatrixMath.maximumNorm(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2))) + "\n");
            }

            writer.close();
            System.out.println("Done writing! Output to HilbertOutput.txt");
        } catch (Exception ex) {
            System.out.println("Writing failed!\n" + ex.toString());
        }
    }

    private static void writeErrorData() {
        try {
            final int iterations = 21;
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
                writer.println(MatrixMath.maximumNorm(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2))));
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
                writer.println(MatrixMath.maximumNorm(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2))));
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
                writer.println(MatrixMath.maximumNorm(
                                h.multiply(new Matrix(temp))
                                 .subtract(new Matrix(temp2))));
            }

            writer.close();
            System.out.println("Done writing! Output to ErrorOutput.txt");
        } catch (Exception ex) {
            System.out.println("Writing failed!\n" + ex.toString());
        }
    }

    private static Matrix parser(String fileName) {
        Matrix m = new Matrix(1, 1);
        try {
            File file = new File(fileName);
            int r = 1;
            ArrayList<ArrayList> mat3 = new ArrayList<ArrayList>();
            double[] mat2;
            Scanner kb = new Scanner(file);
            while (kb.hasNextLine()) {
                String line = kb.nextLine();
                ArrayList<Double> mat = new ArrayList<Double>();
                while (line.indexOf(' ') != -1) {
                    double num = Double.parseDouble(line.substring(0, line.indexOf(' ')));
                    line = line.substring(line.indexOf(' ') + 1, line.length());
                    mat.add(num);
                }
                double num = Double.parseDouble(line.substring(0, line.length()));
                mat.add(num);
                //mat2 = mat.toArray();
                mat3.add(mat);
                r++;
            }
            double[][] mat4 = new double[r][1];
            for (int i = 0; i < r; i++) {
                ArrayList<Double> mat = mat3.get(i);
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

    private static double[] inputParser(String input) {
        //[8, 3, 9, 5, 1]
        ArrayList<Double> mat = new ArrayList<Double>();
        input = input.substring(input.indexOf("[") + 1, input.length() - 1);
        while (input.indexOf(",") != -1) {
            System.out.println(input.indexOf(","));
            mat.add(Double.parseDouble(input.substring(0, input.indexOf(","))));
            input = input.substring(input.indexOf(",") + 2, input.length());
        }
        mat.add(Double.parseDouble(input));
        double[] v = new double[mat.size()];
        for(int i = 0; i < mat.size(); i++) {
            v[i] = mat.get(i);
        }
        return v;
    }

    public static double[] class xGenerator() {
        int rand = Math.round(Math.random()*7) + 3;
        double[] x = new double[rand];
        for (int i = 0; i < x.length; i++) {
            x[i] = Math.round(Math.random());
        }
    }


}