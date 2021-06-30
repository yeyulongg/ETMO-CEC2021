package etmo.metaheuristics.dmoea_lem;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.comparators.CrowdingComparator;
import org.apache.commons.math3.genetics.Population;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.SymmLQ;

public class Autoencoding {

    public SolutionSet ae_Prediction(SolutionSet population,SolutionSet curr_NDS, SolutionSet his_NDS, int NP, Problem problem, int generations) throws JMException {

        int m = problem.getNumberOfObjectives();    //the number of objectives
        int N_curr = curr_NDS.size();               //the number of solution in curr_NDS
        int N_his = his_NDS.size();                 //the number of solution in his_NDS
        int d = curr_NDS.get(0).numberOfVariables();//the dimension of variable
        int N_l=0;
        if(N_curr>N_his)
            N_l = N_his;
        else
            N_l = N_curr;
        Distance distance = new Distance();
        distance.crowdingDistanceAssignment(his_NDS,problem.getNumberOfObjectives());
        his_NDS.sort(new CrowdingComparator());

        distance.crowdingDistanceAssignment(curr_NDS,problem.getNumberOfObjectives());
        curr_NDS.sort(new CrowdingComparator());

        SolutionSet his2_NDS = new SolutionSet();
        SolutionSet curr2_NDS = new SolutionSet();
        double Array_his[][] = new double[N_l][d];
        double Array_curr[][] = new double[N_l][d];
        for(int k = 0;k < N_l; k++){
           his2_NDS.add(his_NDS.get(k));
           curr2_NDS.add(curr_NDS.get(k));
        }
        //获取两个时间窗的非支配个体的矩阵形式种群
        for(int i=0;i<N_l;i++){
            for(int j=0;j<d;j++){
                Array_his[i][j] = his2_NDS.get(i).getDecisionVariables()[j].getValue();
                Array_curr[i][j] = curr2_NDS.get(i).getDecisionVariables()[j].getValue();
            }
        }
        RealMatrix Maxtri_his = new Array2DRowRealMatrix(Array_his);
        RealMatrix Maxtri_curr = new Array2DRowRealMatrix(Array_curr);
        RealMatrix Q = Maxtri_his.multiply(Maxtri_his.transpose());
        RealMatrix P = Maxtri_curr.multiply(Maxtri_his.transpose());
        double lambda = 1e-5;

        //生成单位矩阵,乘lambda
        double ident[][] = new double[N_l][N_l];
        for(int i=0;i<N_l;i++){
            ident[i][i] = lambda;
        }
        ident[N_l-1][N_l-1] = 0;
        RealMatrix reg = new Array2DRowRealMatrix(ident);
        RealMatrix inver = inverseMatrix(Q.add(reg));
        RealMatrix M = P.multiply(inver);//the learned matrix M

        RealMatrix matrix_varM = M.multiply(Maxtri_his);
        double Array_varM[][] = matrix_varM.getData();
        double Array_var[][] = new double[N_l][d];
        for(int i=0;i<N_l;i++){
            for(int j=0;j<d;j++){
                Array_var[i][j] = Math.pow((Array_curr[i][j]-Array_varM[i][j]),2);
            }
        }
        double mean[] = new double[N_l];
        double sum_row=0.0;
        double sum_column=0.0;
        for (int i=0;i<N_l;i++){
            //计算每一行的均值
            for(int j=0;j<d;j++){
                sum_row+=Array_var[i][j];
            }
            mean[i] = sum_row/d;
            //获取每行均值后，对各行均值求均值
            sum_column+= mean[i];
        }
        double v = sum_column/N_l;

        double Array_v[][] = new double[N_l][d];

        //对v进行处理成矩阵形式
        for(int i=0;i<N_l;i++){
            for(int j=0;j<d;j++){
                Array_v[i][j] = v;
            }
        }
        RealMatrix matrix_v = new Array2DRowRealMatrix(Array_v);
        //预测种群的矩阵形式
        RealMatrix matrix_pre_solution = (M.multiply(Maxtri_curr)).add(matrix_v);
        double Array_pre_solution[][] = matrix_pre_solution.getData();//数组形式

        //将预测的矩阵形式的种群转换成 SolutionSet,初始化
        SolutionSet init_Pop = new SolutionSet();


        //若预测种群个体数大于 种群数一半，则只取前一半
        if(N_l>NP/2){

            for(int i=0;i<NP/2;i++){
                init_Pop.add(curr_NDS.get(i));
            }
            //获取两个时间窗的非支配个体的矩阵形式种群
            for(int i=0;i<NP/2;i++){
                for(int j=0;j<d;j++){
                    init_Pop.get(i).getDecisionVariables()[j].setValue(Array_pre_solution[i][j]);
                }
            }
        }
        else{
            //初始化的时候让其等同于curr_NDS,避免初始为空中群
            for(int i=0;i<N_l;i++){
                init_Pop.add(curr_NDS.get(i));
            }

            for(int i=0;i<N_l;i++){
                for(int j=0;j<d;j++){
                    init_Pop.get(i).getDecisionVariables()[j].setValue(Array_pre_solution[i][j]);
                }
            }
        }
        int left = NP - init_Pop.size();
        for(int i=0;i<left;i++){
            if(i<curr2_NDS.size())
            init_Pop.add(curr_NDS.get(i));
        }
        //检查预测种群的决策变量范围
        for(int i=0;i<init_Pop.size();i++){
            for(int j=0;j<init_Pop.get(0).numberOfVariables();j++){
                double lower = curr_NDS.get(0).getDecisionVariables()[j].getLowerBound();
                double upper = curr_NDS.get(0).getDecisionVariables()[j].getUpperBound();
                if(init_Pop.get(i).getDecisionVariables()[j].getValue()<lower)
                    init_Pop.get(i).getDecisionVariables()[j].setValue(lower);
                if(init_Pop.get(i).getDecisionVariables()[j].getValue()>upper)
                    init_Pop.get(i).getDecisionVariables()[j].setValue(upper);
            }
        }
        int init_Pop_size = init_Pop.size();
        //int popSize = population.size();
        int popSize = NP;
        int diff =  popSize - init_Pop_size;
        if(diff>0){
            //预测种群数未达到原种群数,则原种群中删除预测种群数个个体，空出位置放预测个体。
            for(int i=0;i<init_Pop_size;i++){
                population.remove(i);
                population.add(init_Pop.get(i));
            }
        }
        return population;
    }
    //求逆函数
    public static RealMatrix inverseMatrix(RealMatrix A) {
            A.transpose();
            double[][] B = A.getData();
            if(B.length!=B[1].length){
                System.out.println("singur");
            }
            RealMatrix result = new LUDecomposition(A).getSolver().getInverse();
            return result;
    }

}
