import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import java.util.Arrays;


/**
 * Created by Sebas.Hollow on 15.11.27.
 */

public class Main {
    public final static boolean DEBUG = true;
    public final static double PRECISION = 0.00001;
    public final static int ITERATIONS = 500;
    static double[][] d = {{1.342,   .202,  -.599,   .432},
                           { .202,  1.342,   .202,  -.599},
                           {-.599,   .202,  1.342,   .202},
                           { .432,  -.599,   .202,  1.342}};
    static double[][] c = {{0.36,  0,      0,      0},
                           {0,     0.36,   0,      0},
                           {0,     0,      0.36,   0},
                           {0,     0,      0,      0.36}};
    static double[] b = {1.941, -.230, -1.941, -.230};
    public static RealMatrix aMatrix = MatrixUtils.createRealMatrix(c).add(MatrixUtils.createRealMatrix(d));
    public static RealVector bVector = MatrixUtils.createRealVector(b);


    public static void main(String[] args) {
   /*     solveSOR(); System.out.println(); */
        //solveConjugate(); System.out.println();

    new InverseIteration().solveInverseIteration(6);
    }

    public static void log (String s){
        if (DEBUG) System.out.println(s);
    }

    public static void testSOR(){
        double[][] a = {{4, -1, -1},
                        {6,  8,  0},
                        {-5, 0, 12}};
        double[] b = {-2, 45, 80};
        aMatrix = MatrixUtils.createRealMatrix(a);
        bVector = MatrixUtils.createRealVector(b);
        SORMethod method = SORMethod.create(aMatrix, bVector);
        if (method == null)
            return;
        double[] xArray = method.solve(100, 1.1);
        System.out.println(Arrays.toString(xArray));
    }

    public static void testConjugateGradientMethod(){
        double[][] arrayA = {{2, 1, 0.95},
                             {1, 2, 1},
                             {0.95, 1, 2}};
        double[] arrayB = {3.95, 4, 3.95};
        aMatrix = MatrixUtils.createRealMatrix(arrayA);
        bVector = MatrixUtils.createRealVector(arrayB);
        ConjugateGradientMethod method = new ConjugateGradientMethod(aMatrix, bVector);
        method.solve(10);
    }

    public static void solveConjugate(){
        ConjugateGradientMethod method = new ConjugateGradientMethod(aMatrix, bVector);
        double[] answer = method.solve(ITERATIONS);
        if (answer == null)
            System.out.println("The method did not converge. Increasing amount of iterations might solve this problem.");
    }

    public static void solveSOR(){
        SORMethod method = SORMethod.create(aMatrix, bVector);
        double[] answer = method.solve(ITERATIONS, 1.1);
        if (answer == null)
            System.out.println("The method did not converge. Increasing amount of iterations might solve this problem.");
    }
}
