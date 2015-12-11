import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import java.text.MessageFormat;
import java.util.Arrays;

/**
 * Created by Sebas.Hollow on 28/11/2015.
 */
public class SORMethod {
    public RealMatrix aMatrix;
    public RealVector bVector;
    public double[] xs;
    public int length;
    public double omega;

    private SORMethod (RealMatrix a, RealVector b){
        aMatrix = a;
        bVector = b;
        length = a.getRowDimension();
        xs = new double[length];
    }

    public static SORMethod create(RealMatrix a, RealVector b){
        if (!a.equals(a.transpose())) {
            System.out.println("Warning: Matrix is not symmetrical. Might not converge.");
        }

        if (!isConvergent(a))
            System.out.println("Warning: Convergence criteria is not met. Might not converge.");

        return new SORMethod(a, b);
    }


    public double[] solve(int iterations, double omega) {
        this.omega = omega;
        double[] prev;
        for (int i = 1; i <= iterations; i++){
            prev = xs.clone();
            xs = iterate(xs);

            System.out.println(MessageFormat.format("SORMethod.solve {0}: {1}", i, Arrays.toString(xs)));
            if (Main.PRECISION > getMaxDifference(xs, prev))
                return xs;
        }
        return null;
    }

    private double[] iterate(double[] xs){
        for (int j = 0; j < length; j++)
            xs[j] = calculateX(xs, j);
        return xs;
    }

    private double calculateX(double[] xs, int index){
        double x = bVector.getEntry(index);
        double[] row = aMatrix.getRow(index);

        for (int i = 0; i < length; i++){
            if (index != i)
                x -= (row[i] * xs[i]);
        }
        x /= row[index];
        x = xs[index] * (1 - omega) + omega * (x);
        return x;
    }

    private static boolean isConvergent(RealMatrix m){
        for (int i = 0; i < m.getRowDimension(); i++)
            if (!rowConvergent(m.getRow(i), i))
                return false;
        return true;
    }

    private static boolean rowConvergent(double[] row, int diagonalValueIndex){
        double rowSum = 0;
        for (double aRow : row)
            rowSum += Math.abs(aRow);
        return 2*row[diagonalValueIndex] > rowSum;
    }

    private double getMaxDifference (double[] array, double[] array2){
        double maxDiff = 0;
        for (int i = 0; i < array.length; i++){
            double diff = Math.abs((array[i] - array2[i]));
            if (maxDiff < diff)
                maxDiff = diff;
        }
        return maxDiff;
    }
}
