import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import org.junit.Test;
import org.junit.Before;

/**
 *
 * @author Ben Bohannon
 */
public class LeslieMatrixCalc {
    double[] initial = {210000, 190000, 180000, 210000, 200000, 170000, 120000,
        90000, 50000};
    double[][] mat = {{0, 1.2, 1.1, 0.9, 0.1, 0, 0, 0, 0},
                      {0.7, 0, 0, 0, 0, 0, 0, 0, 0},
                      {0, 0.85, 0, 0, 0, 0, 0, 0, 0},
                      {0, 0, 0.9, 0, 0, 0, 0, 0, 0},
                      {0, 0, 0, 0.9, 0, 0, 0, 0, 0},
                      {0, 0, 0, 0, 0.88, 0, 0, 0, 0},
                      {0, 0, 0, 0, 0, 0.8, 0, 0, 0},
                      {0, 0, 0, 0, 0, 0, 0.77, 0, 0},
                      {0, 0, 0, 0, 0, 0, 0, 0.4, 0}};
    Matrix A;
    double[] x;

    @Before
    public void setUp() {
        A = new Matrix(mat);
        x = initial;
    }

    @Test
    public void outputMatrix() {
        try {
            PrintWriter writer = new PrintWriter(new File("PopulationOutput.txt"));

            //Make x into a Matrix.
            double[][] temp = new double[x.length][1];
            for (int i = 0; i < x.length; i++) {
                temp[i][0] = x[i];
            }
            Matrix xMat = new Matrix(temp);
            double oldPop = 0;

            for (int i = 2001; i < 2051; i++) {
                xMat = A.multiply(xMat);

                double[] temp2 = xMat.getColumn(0);
                double total = 0;
                for (int j = 0; j < temp2.length; j++) {
                    total += temp2[j];
                }

                if (i % 10 == 0) {
                    writer.println("Year " + i + " distribution:\n"
                            + Arrays.toString(xMat.getColumn(0)));

                    writer.println("Total: " + total);
                    if (oldPop != 0) {
                        writer.println("Change: " + total / oldPop + "\n");
                    }
                }

                oldPop = total;

            }
            double[] ini = {1, 1, 1, 1, 1, 1, 1, 1, 1};
            double[][] eigenInfo = MatrixMath.power_method(A, 1e-9, ini);

            writer.println("\nEigenvalue: " + eigenInfo[1][0] + "\n");
            writer.println("-------Halved in 2020 Data----------");

            xMat = new Matrix(temp);
            for (int i = 2001; i < 2051; i++) {
                if (i == 2020) {
                    A.set(0, 1, A.get(0, 1) / 2.0);
                    writer.println("\nBirth rate halved!\n");
                }
                xMat = A.multiply(xMat);

                double[] temp2 = xMat.getColumn(0);
                double total = 0;
                for (int j = 0; j < temp2.length; j++) {
                    total += temp2[j];
                }

                if (i % 10 == 0) {
                    writer.println("Year " + i + " distribution:\n"
                            + Arrays.toString(xMat.getColumn(0)));

                    writer.println("Total: " + total);
                    if (oldPop != 0) {
                        writer.println("Change: " + total / oldPop + "\n");
                    }
                }

                oldPop = total;

            }

            eigenInfo = MatrixMath.power_method(A, 1e-9, ini);

            writer.println("\nEigenvalue: " + eigenInfo[1][0] + "\n");


            writer.close();

        } catch (Exception e) {

        }
    }


}
