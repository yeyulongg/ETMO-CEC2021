package etmo.metaheuristics.dmoea_lem;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class test {
    public  static void main(String [] args){
        double testmatrix[][] = new double[2][2];
        testmatrix[0][0] = 1;
        testmatrix[0][1] = 2;
        testmatrix[1][0] = 3;
        testmatrix[1][1] = 4;
        RealMatrix testmatrix1 = new Array2DRowRealMatrix(testmatrix);

        RealMatrix inversetestMatrix = inverseMatrix(testmatrix1);
        System.out.println("两个矩阵相乘后的结果为：\t"+testmatrix1.multiply(inversetestMatrix));
    }
    //求逆函数
    public static RealMatrix inverseMatrix(RealMatrix A) {
        RealMatrix result = new LUDecomposition(A).getSolver().getInverse();
        return result;
    }

}
