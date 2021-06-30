package etmo.metaheuristics.dmoea_lem;

import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.Ranking;
import etmo.util.comparators.CrowdingComparator;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.ToDoubleBiFunction;

public class MultiDirection_Prediction {
    private static final double Gaussian_radio = 100.0;
    private double[] zideal_; //ideal point
    private double[] znadir_;//Nadir point
    public SolutionSet multiDirection_Prediction(List<SolutionSet> DynamicPopulationSets, SolutionSet population) throws JMException,ClassNotFoundException{
        zideal_ = new double[DynamicPopulationSets.get(0).get(0).getNumberOfObjectives()];
        znadir_ = new double[DynamicPopulationSets.get(0).get(0).getNumberOfObjectives()];
        int size = DynamicPopulationSets.size();
        List<SolutionSet> list1 = new <SolutionSet>ArrayList();
        List<SolutionSet> list2 = new <SolutionSet>ArrayList();


            SolutionSet st1 = DynamicPopulationSets.get(size-1);
            SolutionSet st2 = DynamicPopulationSets.get(size-2);
            //环境t聚类
            for (int i=0;i<st1.size();i++){
                SolutionSet sols = new SolutionSet();
                sols.add(st1.get(i));
                list1.add(sols);
            }
            //聚类前数据准备
            estimateIdealPoint(st1);//计算最小的目标函数值
            estimateNadirPoint(st1);//计算最大的目标函数值
            normalizationObjective(st1);//归一化
            computeDistanceToIdealPoint(st1);
            list1 = new HierarchicalClustering1(list1).clusteringAnalysis(5);

            //环境t-1聚类
            for (int i=0;i<st2.size();i++){
                SolutionSet sols = new SolutionSet();
                sols.add(st2.get(i));
                list2.add(sols);
            }
            estimateIdealPoint(st2);//计算最小的目标函数值
            estimateNadirPoint(st2);//计算最大的目标函数值
            normalizationObjective(st2);//归一化
            computeDistanceToIdealPoint(st2);
            list2 = new HierarchicalClustering1(list2).clusteringAnalysis(5);
            //list1为t-1环境的聚类（包含4个SolutionSet) list2为t环境的聚类（包含4个SolutionSet)
            //分别对t-1环境和t环境的每个聚类计算中心点
            ArrayList<double[]> NonDom_Clus_centerSet_t = new ArrayList<double[]>();
            ArrayList<double[]> NonDom_Clus_centerSet_t_1 = new ArrayList<double[]>();

            NonDom_Clus_centerSet_t = comp_Nondom_center(list1);
            NonDom_Clus_centerSet_t_1 = comp_Nondom_center(list2);

            //ArrayList<double[]> Clus_centerSet_t_1 = new ArrayList<double[]>();
            //ArrayList<double[]> Clus_centerSet_t = new ArrayList<double[]>();

            //Clus_centerSet_t_1 = comp_center(list2);
            //Clus_centerSet_t = comp_center(list1);

            //找出当前类与上一个环境中哪个类最近，跟近的类做向量，则预测种群为当前类加上向量
            //list1.get(i)为种群，list1.get(i).get(j)为个体j
            int popSize = population.size();
            population.clear();
            ArrayList<double[]> dirta = comp_DirtaT(NonDom_Clus_centerSet_t_1,NonDom_Clus_centerSet_t);
            Random random = new Random();
            for(int i=0;i<list1.size();i++){
                int[] flag = Support.nonDominateOprate(list1.get(i));
                for(int j=0;j<list1.get(i).size();j++){
                    if(flag[j] == 1){
                        for(int k=0;k<list1.get(i).get(j).numberOfVariables();k++){
                            double temp = list1.get(i).get(j).getDecisionVariables()[k].getValue()+dirta.get(i)[k];
                            //temp+=random.nextGaussian()/Gaussian_radio;
                            double lower = list1.get(i).get(j).getDecisionVariables()[k].getLowerBound();
                            double upper = list1.get(i).get(j).getDecisionVariables()[k].getUpperBound();
                            if(temp<lower)
                                temp=lower;
                            if(temp>upper)
                                temp = upper;
                            //将原种群的个体加上向量值，从而得到预测种群
                            list1.get(i).get(j).getDecisionVariables()[k].setValue(temp);
                        }
                    }
                    else {
                        list1.get(i).get(j).initialize();
                    }
                }
            }

            /*
            int Clus_Size = list1.size();
            for (int i=0;i<Clus_Size;i++){
                int solutionSet_Size = list1.get(i).size();
                for (int j=0;j<solutionSet_Size;j++){
                    int numberofVariables = list1.get(i).get(0).numberOfVariables();
                    for (int k=0;k<numberofVariables;k++){
                        double temp = list1.get(i).get(j).getDecisionVariables()[k].getValue()+dirta.get(i)[k];
                        temp+=random.nextGaussian()/Gaussian_radio;
                        double lower = list1.get(i).get(j).getDecisionVariables()[k].getLowerBound();
                        double upper = list1.get(i).get(j).getDecisionVariables()[k].getUpperBound();
                        if(temp<lower)
                            temp=lower;
                        if(temp>upper)
                            temp = upper;
                        //将原种群的个体加上向量值，从而得到预测种群
                        list1.get(i).get(j).getDecisionVariables()[k].setValue(temp);
                    }
                }
            }*/
            //将每个类中的个体添加到种群中

            population.clear();
            for (int i=0;i<list1.size();i++){
                for (int j=0;j<list1.get(i).size();j++)
                {
                    population.add(list1.get(i).get(j));
                }
            }
        return population;
    }
    //计算每一聚类的非支配质心
    public static ArrayList<double[]> comp_Nondom_center(List<SolutionSet> list)throws  JMException{
        ArrayList<double[]> Nondom_centerList = new ArrayList<double[]>();
        for(int i=0;i<list.size();i++){
            SolutionSet population = list.get(i);
            Ranking ranking = new Ranking(population);
            SolutionSet Nondom_population = new SolutionSet();
            for(int j=0;j<ranking.getSubfront(0).size();j++){
                Nondom_population.add(ranking.getSubfront(0).get(j));
            }
            double[] result = new double[population.get(0).numberOfVariables()];//存放中心点决策变量
            double sum;
            int count;
            for(int z=0;z<Nondom_population.get(0).numberOfVariables();z++){
                sum = 0.0;
                count = 0;
                for(int k=0;k<Nondom_population.size();k++){
                    sum+=Nondom_population.get(k).getDecisionVariables()[z].getValue();
                    count++;
                }
                result[z] = sum/count;
            }
            Nondom_centerList.add(result);
        }
        return Nondom_centerList;
    }

    //计算每一聚类的中心
    public static ArrayList<double[]> comp_center(List<SolutionSet> list) throws JMException {
        ArrayList<double[]> centerList = new ArrayList<double[]>();
        for (int i=0;i<list.size();i++){
            SolutionSet population = list.get(i);
            double[] result = new double[population.get(0).numberOfVariables()];//存放中心点决策变量
            double sum;
            int count;
            for(int j=0;j<population.get(0).numberOfVariables();j++){
                sum = 0.0;
                count = 0;
                for(int k=0;k<population.size();k++){
                    sum+=population.get(k).getDecisionVariables()[j].getValue();
                    count++;
                }
                result[j] = sum/count;
            }
            centerList.add(result);
        }
        return centerList;
    }


    //找出当前类与上一个环境中哪个类最近，跟近的类做向量dirtaT
    public ArrayList<double[]> comp_DirtaT(ArrayList<double[]> Clus_centerSet_t_1,ArrayList<double[]> Clus_centerSet_t){
        ArrayList<double[]> dirList = new ArrayList<double[]>();//存放dirt向量
        int len = Clus_centerSet_t.size();
        double [][] value = new double[len][len];
        for(int i=0;i<Clus_centerSet_t.size();i++){
            for (int j=0;j<Clus_centerSet_t_1.size();j++){
                double [] center_t = Clus_centerSet_t.get(i);
                double [] center_t_1 = Clus_centerSet_t_1.get(j);
                double result=0.0;
                for(int k=0;k<center_t.length;k++){
                    result+=Math.pow((center_t[k]-center_t_1[k]),2);
                }
                value[i][j] = Math.sqrt(result);
            }
        }
        for(int i=0;i<value.length;i++){
            double [] dirta = new double[value.length];
            for(int j=0;j<value[0].length;j++){
                double min = value[i][0];
                dirta = DVector(Clus_centerSet_t.get(i),Clus_centerSet_t_1.get(0));
                if (value[i][j]<min){
                    min=value[i][j];
                    dirta = DVector(Clus_centerSet_t.get(i),Clus_centerSet_t_1.get(j));
                }
            }
            dirList.add(dirta);
        }
        return dirList;
    }

    //考虑相减是谁减谁，计算dirta
    public static double[] DVector(double[] center_t,double[] center_t_1) {
        int length = center_t.length;
        double[] vector = new double[length];
        for (int n = 0; n < length; n++)
            vector[n] = center_t[n]-center_t_1[n];
        return vector;
    }
    /*
     * Estimate the Ideal Point
     */
    public void estimateIdealPoint(SolutionSet solutionSet){
        for(int i=0; i<solutionSet.get(0).getNumberOfObjectives();i++){
            zideal_[i] = 1.0e+30;
            for(int j=0; j<solutionSet.size();j++){
                if(solutionSet.get(j).getObjective(i) < zideal_[i]){
                    zideal_[i] = solutionSet.get(j).getObjective(i);
                }
            }

        }
    }

    /*
     * Estimate the Nadir Point
     */
    public void estimateNadirPoint(SolutionSet solutionSet){
        for(int i=0; i<solutionSet.get(0).getNumberOfObjectives();i++){
            znadir_[i] = -1.0e+30;
            for(int j=0; j<solutionSet.size();j++){
                if(solutionSet.get(j).getObjective(i) > znadir_[i]){
                    znadir_[i] = solutionSet.get(j).getObjective(i);
                }
            }

        }
    }

    /*
     * Normalization
     */
    public void normalizationObjective(SolutionSet solutionSet){
        for(int i=0; i<solutionSet.size(); i++){
            Solution sol = solutionSet.get(i);

            for(int j=0; j<solutionSet.get(0).getNumberOfObjectives(); j++){
                double val = 0.0;
                val = (sol.getObjective(j) - zideal_[j])/(znadir_[j]-zideal_[j]);
                //val = (sol.getObjective(j) - zideal_[j]);
                sol.setNormalizedObjective(j, val);
            }
        }
    }

    /*
     * Compute the Convergence Distance of each Solutions Which use the distance of
     * each solution to the Ideal Point
     */
    //HCM计算两个个体角度公式的分母中的一个
    public void computeDistanceToIdealPoint(SolutionSet solutionSet){
        for(int i=0; i<solutionSet.size(); i++){
            Solution sol = solutionSet.get(i);
            double normDistance = 0.0;
            double sumValue = 0.0;
            for(int j=0; j<solutionSet.get(0).getNumberOfObjectives(); j++){
                normDistance += sol.getNormalizedObjective(j) * sol.getNormalizedObjective(j);
                sumValue +=  sol.getNormalizedObjective(j);
            }
            normDistance = Math.sqrt(normDistance);

            sol.setDistanceToIdealPoint(normDistance);
            sol.setSumValue(sumValue);
        }
    }

}
