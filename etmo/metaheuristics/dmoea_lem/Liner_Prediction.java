package etmo.metaheuristics.dmoea_lem;

import etmo.core.SolutionSet;
import etmo.util.JMException;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Liner_Prediction extends Predition {
    private List<double[]> centers_;

    private List<double[]> vars_;

    private List<double[]> predict_Centors_;
    public SolutionSet population_;
    private static final double Gaussian_radio = 2000.0;
    public Liner_Prediction(SolutionSet population){
        super(population);
        centers_ = new ArrayList<double[]>();
        vars_ = new ArrayList<double[]>();
        predict_Centors_ = new ArrayList<double[]>();
    }
    public SolutionSet prediton_exe(SolutionSet population)throws JMException ,ClassNotFoundException{
        updateCenterSet(population);
        int length=centers_.size();
        if(length!=1)
        {
            Random random =new Random();
            //返回两代质点的差，构成向量
            double[] vector = Support.DVector(centers_.get(length-1),centers_.get(length-2));
            for(int n=0;n<population.size();n++)
            {
                for(int i=0;i<numberOfVariables;i++)
                {
                    double temp=population.get(n).getDecisionVariables()[i].getValue()+vector[i];
                    temp+=random.nextGaussian()/Gaussian_radio;//将原种群的个体加上向量值，从而得到预测种群
                    double lower = population.get(n).getDecisionVariables()[i].getLowerBound();
                    double upper = population.get(n).getDecisionVariables()[i].getUpperBound();
                    /*
                     * max,min
                     */
                    if(temp<lower)
                        temp=lower;
                    if(temp>upper)
                        temp = upper;
                    //更新个体的决策变量
                    population.get(n).getDecisionVariables()[i].setValue(temp);
                }
            }
        }
        //将当前种群的非支配个体中心点加入集合
        double[] preCenter=Support.compCenter_nondomin(population);
        predict_Centors_.add(preCenter);
        return population;
    }
    public void updateCenterSet(SolutionSet population) throws JMException {
        double[] center=Support.compCenter_nondomin(population);
        centers_.add(center);
    }
}
