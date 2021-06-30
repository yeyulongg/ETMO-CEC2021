package etmo.metaheuristics.dmoea_lem;
import etmo.core.*;
import etmo.metaheuristics.dmoea_lem.service.DMOEA_LEM;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.Configuration;
import etmo.util.JMException;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

public class AE_DMOEA_Multi_main  {
    public static Logger logger_; // Logger object
    public static FileHandler fileHandler_; // FileHandler object

    public static void main(String[] args) throws JMException,
            SecurityException, IOException,ClassNotFoundException{
        logger_ = Configuration.logger_;
        fileHandler_ = new FileHandler("MOEAD_Dynamic.log");
        logger_.addHandler(fileHandler_);

        ProblemSet problemSet; // The problem to solve
        DynamicAlgorithm algorithm; // The algorithm to use
        Operator crossover; // Crossover operator
        Operator mutation; // Mutation operator
        Operator selection;

        HashMap parameters; // Operator parameters
        int fc = 20;
        int sc = 10;
        int nc = 30;
        int populationSize=0;
        for(int fun=33;fun<=40;fun++)
        {
            problemSet = null;
            switch (fun){
                case 33:
                    problemSet = ETMOF33.getProblem(fc,sc);
                    break;
                case 34:
                    problemSet = ETMOF34.getProblem(fc,sc);
                    break;
                case 35:
                    problemSet = ETMOF35.getProblem(fc,sc);
                    break;
                case 36:
                    problemSet = ETMOF36.getProblem(fc,sc);
                    break;
                case 37:
                    problemSet = ETMOF37.getProblem(fc,sc);
                    break;
                case 38:
                    problemSet = ETMOF38.getProblem(fc,sc);
                    break;
                case 39:
                    problemSet = ETMOF39.getProblem(fc,sc);
                    break;
                case 40:
                    problemSet = ETMOF40.getProblem(fc,sc);
                    break;
                default:break;
            }
            int taskNumber = problemSet.size();
            String[][] pf = new String[taskNumber][nc+1];
            for(int j=0;j<taskNumber;j++){
                for(int i=0;i<nc+1;i++){
                    pf[j][i] = "PF/DynamicPF/" + problemSet.get(j).getName() + "/POF_Tt=" + i + ".txt";
                }
            }
                System.out.println("taskNumber = "+taskNumber);

                algorithm = new AE_DMOEA_Multi(problemSet);
                if(problemSet.get(0).getNumberOfObjectives()==2){
                    populationSize =100;
                }else{
                    populationSize =105;
                }
                algorithm.setInputParameter("populationSize", populationSize);
                //algorithm.setInputParameter("maxEvaluations", populationSize*fc*(nc+1));
                algorithm.setInputParameter("maxEvaluations", taskNumber*populationSize*fc*(nc+2));
                algorithm.setInputParameter("dataDirectory","F:/weight");//MOEAD 权重目录
                algorithm.setInputParameter("T", 20) ;
                algorithm.setInputParameter("delta", 0.9) ;
                algorithm.setInputParameter("nr", 2) ;

                // Crossover operator DE
                parameters = new HashMap();
                parameters.put("CR", 1.0);
                parameters.put("F", 0.5);
                crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

                // Mutation operator
                parameters = new HashMap();
                parameters.put("probability", 1.0 / problemSet.get(0).getNumberOfObjectives());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                // Selection Operator
                parameters = null ;
                selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters) ;

                // Add the operators to the algorithm
                algorithm.addOperator("crossover", crossover);
                algorithm.addOperator("mutation", mutation);
                algorithm.addOperator("selection", selection);

                DecimalFormat form = new DecimalFormat("#.####E0");
                System.out.println("RunID\t" + "IGD for "+problemSet.get(0).getName()+" to "+problemSet.get(taskNumber-1).getName());
                QualityIndicator indicator;
                int times = 21;
                double[] aveIGD = new double[taskNumber];
                for (int i = 1; i <= times; i++) {
                    double[] MIGD = new double[taskNumber];
                    List<SolutionSet> dynamicPopulations = algorithm.execute();

                    //动态环境中，dynamicPopulations存了三十个环境变化中中每个环境的最后一代，由于多任务混合了两个任务的个体
                    SolutionSet[] resPopulation = new SolutionSet[taskNumber];
                    for(int j=0;j<taskNumber;j++)
                        resPopulation[j] = new SolutionSet();
                    for (int k = 0; k < dynamicPopulations.size(); k++) {
                        //每个环境中存储的最后一代为多个任务的个体，故需要做一个划分,对第k个环境划分成tasknumber部分，分别存在resPopulation数组中
                        for (int z = 0; z < dynamicPopulations.get(k).size(); z++) {
                            Solution sol = dynamicPopulations.get(k).get(z);
                            int pid = z / (dynamicPopulations.get(k).size() / taskNumber);
                            int start = problemSet.get(pid).getStartObjPos();
                            int end = problemSet.get(pid).getEndObjPos();
                            Solution newSolution = new Solution(end - start + 1);
                            for (int q = start; q <= end; q++)
                                newSolution.setObjective(q - start, sol.getObjective(q));
                            resPopulation[pid].add(newSolution);
                        }
                        //对第K个环境的每一个任务计算IGD值
                        for(int q=0;q<taskNumber;q++){
                            indicator = new QualityIndicator(problemSet.get(q), pf[q][k]);
                             resPopulation[q].printObjectivesToFile("DMOEA_LEM_Dynamic/DMOEA_LEM_"+problemSet.get(q).getNumberOfObjectives()+"Obj_"+
                                problemSet.get(q).getName()+ "_" + problemSet.get(q).getNumberOfVariables() + "D_run"+i+"_environment"+k+".txt");
                        MIGD[q] += indicator.getIGD(resPopulation[q]);
                        resPopulation[q].clear();
                        }
                    }
                    for(int p=0;p<taskNumber;p++){
                        MIGD[p] = MIGD[p]/dynamicPopulations.size();
                        try {
                            File file = new File("D:\\data\\"+problemSet.get(p).getName()+"task"+p+".txt");
                            FileOutputStream fos = null;
                            if(!file.exists()){
                                file.createNewFile();//文件不存在，创建
                                fos = new FileOutputStream(file);//首次获取写入
                            }else{
                                fos = new FileOutputStream(file,true);//续写
                            }
                            OutputStreamWriter osw = new OutputStreamWriter(fos,"UTF-8");
                            osw.write(MIGD[p]+"\n");
                            osw.flush(); // 把缓存区内容压入文件
                        } catch (IOException e) {
                            e.printStackTrace();
                        }

                        aveIGD[p]+=MIGD[p];
                        System.out.println(i + "\t" + form.format(MIGD[p]));
                    }

                }
                for(int i=0;i<taskNumber;i++){
                    System.out.println("Average MIGD for " + problemSet.get(i).getName() + ": " + form.format(aveIGD[i] / times));
                }
                System.out.println();
        }
    }
}
