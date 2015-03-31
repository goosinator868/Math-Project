import java.util.Arrays;
import org.junit.Test;
import org.junit.Before;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import static org.junit.Assert.assertTrue;

/**
 *
 * @author Ben Bohannon
 */
public class MatrixMathTests {

    Matrix A, B, C, I, H;

    @Before
    public void setUp() {
        double[][] dataIdentity = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; //3x3
        double[][] dataA = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; //3x3
        double[][] dataB = {{1, 0, 1}}; //1x3
        double[][] dataC = {{1, 2}, {3, 4}};

        double[][] dataH =
            {{1.0,  1.0/2, 1.0/3, 1.0/4},
            {1.0/2, 1.0/3, 1.0/4, 1.0/5},
            {1.0/3, 1.0/4, 1.0/5, 1.0/6},
            {1.0/4, 1.0/5, 1.0/6, 1.0/7}};

        A = new Matrix(dataA);
        B = new Matrix(dataB);
        C = new Matrix(dataC);
        I = new Matrix(dataIdentity);
        H = new Matrix(dataH);
    }


    @Test(timeout = 1000)
    public void LUDecompTest() {
        Matrix[] LU = MatrixMath.lu_fact(A);

        System.out.println("Matrix L:\n" + LU[0]);
        System.out.println("Matrix U:\n" + LU[1]);

        assertEquals(A, LU[0].multiply(LU[1]));

        LU = MatrixMath.lu_fact(H);

        System.out.println("Matrix L:\n" + LU[0]);
        System.out.println("Matrix U:\n" + LU[1]);
        
        System.out.println("LU Error: " +
                MatrixMath.maximumNorm(LU[0].multiply(LU[1]).subtract(H)));
    }

    @Test(timeout = 2000)
    public void householderFactorize() {
        System.out.println("\nStarting HouseHolder QR Factorization!\n");
        
        Matrix[] QR = MatrixMath.qr_fact_househ(H);

        System.out.println("Matrix Q:\n" + QR[0]);
        System.out.println("Matrix R:\n" + QR[1]);

        System.out.println("Matrix QR:\n" + QR[0].multiply(QR[1]));
        
        assertTrue(checkOrthogonalness(QR[0]));
        
        System.out.println("HouseHolder Error: " +
                MatrixMath.maximumNorm(QR[0].multiply(QR[1]).subtract(H)));
    }

    @Test(timeout = 1000)
    public void givensFactorize() {
        System.out.println("\nStarting Givens QR Factorization!\n");
        Matrix[] QR = MatrixMath.qr_fact_givens(H);

        System.out.println("Matrix Q:\n" + QR[0]);
        System.out.println("Matrix R:\n" + QR[1]);

        System.out.println("Matrix QR:\n" + QR[0].multiply(QR[1]));

        assertTrue(checkOrthogonalness(QR[0]));
        
        System.out.println("Givens Error: " +
                MatrixMath.maximumNorm(QR[0].multiply(QR[1]).subtract(H)));
    }

    @Test(timeout = 1000)
    public void luSolveTest() {
        System.out.println("\nStarting LU solution Test!\n");

        Matrix[] LU = MatrixMath.lu_fact(H);
        double[] test = {0.0464159, 0.0464159, 0.0464159, 0.0464159};

        double[] solution = MatrixMath.solve_lu_b(LU[0], LU[1], test);

        System.out.println("solution to Hilbert Matrix:\n"
                + Arrays.toString(solution) + "\n");
    }

    @Test(timeout = 1000)
    public void qrSolveTest() {
        System.out.println("\nStarting QR solution Test!\n");

        Matrix[] QR = MatrixMath.qr_fact_givens(H);
        double[] test = {0.0464159, 0.0464159, 0.0464159, 0.0464159};

        double[] solution = MatrixMath.solve_qr_b(QR[0], QR[1], test);

        System.out.println("solution to Hilbert Matrix:\n"
                + Arrays.toString(solution) + "\n");

        //Temporary test
        double[][] temp = new double[solution.length][1];
        for (int i = 0; i < solution.length; i++) {
            temp[i][0] = solution[i];
        }
        Matrix sol = new Matrix(temp);

        System.out.println("Plugging solution back in:\n" +
                QR[0].multiply(QR[1].multiply(sol)) + "\n");
    }

    private boolean checkOrthogonalness(Matrix Q) {
        double[][] vectors = new double[Q.getColumns()][Q.getRows()];

        for (int i = 0; i < vectors.length; i++) {
            vectors[i] = Q.getColumn(i);
        }

        double maximumDot = 0;

        for (double[] current : vectors) {
            for (double[] dotter : vectors) {
                if (current != dotter) {
                    double dot = 0;

                    for (int i = 0; i < dotter.length; i++) {
                        dot += current[i] * dotter[i];
                    }

                    double absDot = java.lang.Math.abs(dot);
                    if (absDot > 1e-10) {
                        System.out.println("Dot was: " + dot);
                        return false;
                    } else if (absDot > maximumDot) {
                        maximumDot = absDot;
                    }
                }
            }
        }

        System.out.println("Maximum dot was: " + maximumDot);

        return true;
    }

}
