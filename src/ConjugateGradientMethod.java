import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import java.text.MessageFormat;
import java.util.Arrays;

/**
 * Created by Sebastianas.Hollow on 28/11/2015.
 */
public class ConjugateGradientMethod {
    public static RealMatrix aMatrix;
    public static RealVector bVector;
    public static RealVector vectorZ;
    public static RealVector prevZ;
    public static RealVector vectorP;
    public static RealVector vectorR;
    public static RealVector x;
    public static double beta;
    public static double alpha;

    ConjugateGradientMethod(RealMatrix a, RealVector b){
        aMatrix = a;
        bVector = b;
        x = MatrixUtils.createRealVector(new double[aMatrix.getColumnDimension()]);
    }

    public double[] solve(int iterations) {
        //0 Iteration
        vectorZ = aMatrix.preMultiply(x).subtract(bVector);
        vectorP = vectorZ.copy();

        //1+ Iteration
        for (int i = 0; i <= iterations; i++){
            vectorR = calculateR(vectorP);
            alpha = calculateAlpha(vectorZ, vectorP, vectorR);
            x = calculateX(x, alpha, vectorZ);
            prevZ = vectorZ;
            vectorZ = calculateZ(prevZ, alpha, vectorR);
            //Error checking
            printInfo(i);
            if (getError(vectorZ) < Math.pow(Main.PRECISION, 2))
                break;
            beta = calculateBeta(prevZ, vectorZ);
            vectorP = calculateP(vectorZ, beta, vectorP);
        }
        return x.getData();
    }

    private double calculateAlpha (RealVector zVector, RealVector pVector, RealVector rVector){
        RealMatrix z = vToMatrix(zVector);
        RealMatrix r = vToMatrix(rVector);
        RealMatrix p = vToMatrix(pVector);
        double alphaTop = z.transpose().preMultiply(p).getNorm();
        double alphaLow = r.transpose().preMultiply(p).getNorm();
        double alpha = alphaTop/alphaLow;
        Main.log ("Alpha: " + alpha);
        return alpha;
    }
    private double calculateBeta (RealVector z0, RealVector z1){
        RealMatrix prevZ = vToMatrix(z0);
        RealMatrix z = vToMatrix(z1);
        double betaTop = z.transpose().preMultiply(z).getNorm();
        double betaLow = prevZ.transpose().preMultiply(prevZ).getNorm();
        double beta = betaTop/betaLow;
        Main.log ("Beta: " + beta);
        return beta;
    }

    private RealVector calculateR (RealVector p){
        RealVector r = aMatrix.preMultiply(p);
        Main.log ("R vector: " + Arrays.toString(r.getData()));
        return r;
    }

    private RealVector calculateP (RealVector z, double beta, RealVector prevP){
        RealVector p = z.add(prevP.mapMultiply(beta));
        Main.log (MessageFormat.format("P vector: {0} ", Arrays.toString(z.getData())));
        return p;
    }

    private RealVector calculateZ(RealVector prevZ, double a, RealVector r){
        RealVector z = prevZ.subtract(r.mapMultiply(a));
        Main.log (MessageFormat.format("Z vector: {0} ", Arrays.toString(z.getData())));
        return z;
    }

    private RealVector calculateX (RealVector x, double a, RealVector z){
        x = x.subtract(z.mapMultiply(a));
        Main.log (MessageFormat.format("X vector: {0} ", Arrays.toString(x.getData())));
        return x;
    }


    private RealMatrix vToMatrix (RealVector v){
        return MatrixUtils.createRealMatrix(new double[][]{v.getData()});
    }

    private double getError (RealVector z){
        RealMatrix matrix = vToMatrix(z);
        double numbers = matrix.transpose().preMultiply(matrix).getNorm();
        //Main.log("Error " + numbers + "\n");
        return numbers;
    }


    private void printInfo (int i){
        double[] data = vToMatrix(x).getRow(0);
        String format = "ConjugateGradient {0}: {1}\n";
        System.out.println(MessageFormat.format(format, i, Arrays.toString(data)));
    }
}
