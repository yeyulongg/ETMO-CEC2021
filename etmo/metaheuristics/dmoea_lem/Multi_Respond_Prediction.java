package etmo.metaheuristics.dmoea_lem;

import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.Ranking;
import etmo.util.comparators.CrowdingComparator;

import java.util.ArrayList;

public class Multi_Respond_Prediction {
    private double corr_;
    public SolutionSet multi_Respond_Prediction(ArrayList<SolutionSet> dynamicPopulationSets,SolutionSet population)throws JMException,ClassNotFoundException{

        int size = dynamicPopulationSets.size();
        if(size!=1)
        {
            SolutionSet population_t = dynamicPopulationSets.get(size-1);
            SolutionSet population_t_1 = dynamicPopulationSets.get(size-2);
            //选择N/10个非支配个体算相关系数
            SolutionSet nondominated_population_t = select_Nondominated(population_t);
            SolutionSet nondominated_population_t_1= select_Nondominated(population_t_1);
            Similar(nondominated_population_t,nondominated_population_t_1);//计算相关系数

            //计算两个环境的非支配解的质心
            double[] center_t = Support.compCenter_nondomin(population_t);
            double[] center_t_1 = Support.compCenter_nondomin(population_t_1);

            //质心向量
            double[] vector1 = norVector(center_t,center_t_1);//正向
            double[] vector2 = Support.reVector(vector1);//反向

            if(corr_>=0.7){
               //System.out.println(corr_);
                population = new ModelB().modelB(population_t,vector1);
                //population = new MultiDirection_Prediction().multiDirection_Prediction(dynamicPopulationSets,population);
            }
            else if (corr_<0.2){
               //System.out.println(corr_);
                population = new ModelA().modelA(population_t);
            }
            else {
                //System.out.println(corr_);
                population = new ModelC().modelC(population_t,vector1,vector2);
            }

        }

        return population;
    }

    public SolutionSet select_Nondominated(SolutionSet population){     //用于 选择 (检测环境是否变化的) 解   非支配选择
        SolutionSet nondominated_Population = new SolutionSet();
        Distance distance = new Distance();
        Ranking ranking = new Ranking(population);
        int remain = population.size()/10;      //选N/10个
        int index = 0;
        SolutionSet front = null;
        // Obtain the next front
        front = ranking.getSubfront(index);

        while ((remain > 0) && (remain >= front.size())) {
            // Assign crowding distance to individuals
            distance.crowdingDistanceAssignment(front, population.get(0).getNumberOfObjectives());
            // Add the individuals of this front
            for (int k = 0; k < front.size(); k++) {
                nondominated_Population.add(front.get(k));
            } // for
            // Decrement remain
            remain = remain - front.size();
            // Obtain the next front
            index++;
            if (remain > 0) {
                front = ranking.getSubfront(index);
            } // if
        } // while
        // Remain is less than front(index).size, insert only the best one
        if (remain > 0) { // front contains individuals to insert
            distance.crowdingDistanceAssignment(front,population.get(0).getNumberOfObjectives());
            front.sort(new CrowdingComparator());
            for (int k = 0; k < remain; k++) {
                nondominated_Population.add(front.get(k));
            } // for
            remain = 0;
        } // if
        return nondominated_Population;
    }
    public void Similar(SolutionSet nondominated_population_t,SolutionSet nondominated_population_t_1){
        int numberOfObjectives = nondominated_population_t.get(0).getNumberOfObjectives();
        double x[][] = new double[numberOfObjectives][nondominated_population_t.size()];
        double y[][] = new double[numberOfObjectives][nondominated_population_t_1.size()];
        for(int i=0;i<nondominated_population_t_1.size();i++)
            for(int j=0;j<numberOfObjectives;j++)
            {
                x[j][i]=nondominated_population_t.get(i).getObjective(j);
                y[j][i]=nondominated_population_t_1.get(i).getObjective(j);
            }
        double[] radio= Support.oldRadio(x);
        double[] dsOld= Support.Ds(radio, x);
        double[] dsNew= Support.Ds(radio, y);
        corr_=Math.abs(Support.corr(dsOld, dsNew));
       // System.out.println(corr_);
    }
    public double[] norVector(double[] center_t,double[] center_t_1) {
        double[] vector =Support.DVector(center_t_1,center_t);
        return vector;
    }

}
