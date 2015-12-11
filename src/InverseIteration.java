import org.apache.commons.math.linear.*;

import java.text.MessageFormat;
import java.util.Arrays;

/**
 * Created by Sebas.Hollow on 11/29/15.
 */
public class InverseIteration {
    double[][] matrixData = {{3.24, 2,     0,      0},
                             { 2,   3.24,  2,      0},
                             { 0,    2,    3.24,   2},
                             { 0,    0,     2,     3.24}};
    double[] initialEigenvector = {0,0,0,1};
    RealMatrix m = MatrixUtils.createRealMatrix(matrixData);
    RealMatrix realIdentityMatrix = MatrixUtils.createRealIdentityMatrix(4);
    xEigenVectorPair xVec = new xEigenVectorPair(MatrixUtils.createRealVector(initialEigenvector));

    public void solveInverseIteration(double eigenValue){
        double newValue = 0.0;
        RealVector yVec;

        do {
            //FixMe: get rid of this ugliness
            try {
                yVec = getYVec(eigenValue, xVec.get());
            } catch (SingularMatrixException e){
                return;
            }
            xVec.push(getXVec(yVec));
            System.out.println("EigenVector: " + xVec.get());

            eigenValue = newValue;
            newValue = getEigenValueApproximation(xVec.get());

        } while(Math.max(Math.abs(newValue - eigenValue), xVec.difference()) > Main.PRECISION);
    }

    private double getEigenValueApproximation(RealVector rv){
        rv = m.preMultiply(rv).ebeMultiply(rv);
        double sum = 0;
        for (double value : rv.getData())
            sum += value;
        System.out.println("EigenValue = " + sum + "\n");
        return sum;
    }

    private RealVector getXVec(RealVector vec){
        vec = vec.mapDivide(vec.getNorm());
        return vec;
    }

    private RealVector getYVec(double eigenValue, RealVector xVector){
        RealMatrix m = matrixMinusEigen(eigenValue);
        RealVector vector = inverse(m).preMultiply(xVector);
        System.out.println(MessageFormat.format("Y vector data: {0}", Arrays.toString(vector.getData())));
        return vector;
    }

    private RealMatrix inverse(RealMatrix matrixToInverse){
        return new LUDecompositionImpl(matrixToInverse).getSolver().getInverse();
    }

    private RealMatrix matrixMinusEigen(double eigenValue){
        return m.subtract(realIdentityMatrix.scalarMultiply(eigenValue));
    }

    private static class xEigenVectorPair {
        RealVector vec;
        RealVector prevVec;

        xEigenVectorPair(RealVector rv) {
            vec = rv;
            prevVec = null;
        }

        public void push (RealVector rv){
            prevVec = vec;
            vec = rv;
        }

        public RealVector get (){
            return vec;
        }

        public double difference(){
            return vec.subtract(prevVec).getNorm();
        }
    }

}
