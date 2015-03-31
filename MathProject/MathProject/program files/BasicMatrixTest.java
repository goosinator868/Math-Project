import org.junit.Test;
import org.junit.Before;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 *
 * @author Ben Bohannon
 */
public class BasicMatrixTest {

    Matrix A, B, I;

    @Before
    public void setUp() {
        double[][] dataIdentity = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; //3x3
        double[][] dataA = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; //3x3
        double[][] dataB = {{1, 0, 1}}; //1x3

        A = new Matrix(dataA);
        B = new Matrix(dataB);
        I = new Matrix(dataIdentity);

    }

    @Test
    public void toStringTest() {
        System.out.println("Matrix A:\n" + A);
        System.out.println("Matrix B:\n" + B);
        System.out.println("Matrix I:\n" + I);
    }

    @Test(timeout = 1000)
    public void addTest() {
        double[][] compareData = {{2, 2, 3}, {4, 6, 6}, {7, 8, 10}};
        Matrix compare = new Matrix(compareData);

        Matrix addition = A.add(I);

        assertEquals(addition, compare);

        System.out.println("A + I:\n" + addition);
    }

    @Test(timeout = 1000)
    public void subtractTest() {
        double[][] compareData = {{0, 2, 3}, {4, 4, 6}, {7, 8, 8}};
        Matrix compare = new Matrix(compareData);

        Matrix subtraction = A.subtract(I);

        assertEquals(subtraction, compare);

        System.out.println("A-I:\n" + subtraction);
    }

    @Test(timeout = 1000)
    public void multiplyScalarTest() {
        double[][] compareData = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
        Matrix compare = new Matrix(compareData);

        Matrix multiply = A.multiply(2);

        assertEquals(compare, multiply);

        System.out.println("2 * A:\n" + multiply);
    }

    @Test(timeout = 1000)
    public void multiplyMatrixTest() {
        double[][] compareData = {{8, 10, 12}};
        Matrix compare = new Matrix(compareData);

        Matrix multiply = B.multiply(A);

        assertEquals(compare, multiply);

        System.out.println("B * A:\n" + multiply);

        try {
            A.multiply(B);
            fail();
        } catch (java.lang.Exception e) {
            System.out.println("Caught illegal multiply!");
        }
    }

    @Test(timeout = 200)
    public void transposeTest() {
        double[][] compareData = {{1, 4, 7}, {2, 5, 8}, {3, 6, 9}};
        Matrix compare = new Matrix(compareData);

        Matrix transpose = A.transpose();

        assertEquals(transpose, compare);

        System.out.println("Transpose of A:\n" + transpose);
    }
    
    @Test(timeout = 200) 
    public void hilbertTest() {
        System.out.println("Hilbert Matrix:\n" + Matrix.getHilbert(5));
    }
}
