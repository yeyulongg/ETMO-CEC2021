package etmo.metaheuristics.Multitask_MOEAD;

import etmo.core.*;
import etmo.operators.crossover.CrossoverFactory;
import etmo.operators.mutation.MutationFactory;
import etmo.operators.selection.SelectionFactory;
import etmo.problems.benchmarks_ETMO.*;
import etmo.qualityIndicator.QualityIndicator;
import etmo.util.Configuration;
import etmo.util.JMException;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

public class MultiTask_AE_MOEAD_Main {
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
            int taskNumber = 0;
            int populationSize = 0;

            for(int fun=33;fun<=40;fun++)
            {
                problemSet = null;
                switch (fun){
                    case 33:
                        problemSet = ETMOF33.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 100*taskNumber;
                        break;
                    case 34:
                        problemSet = ETMOF34.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 100*taskNumber;
                        break;
                    case 35:
                        problemSet = ETMOF35.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 100*taskNumber;
                        break;
                    case 36:
                        problemSet = ETMOF36.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 100*taskNumber;
                        break;
                    case 37:
                        problemSet = ETMOF37.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 150*taskNumber;
                        break;
                    case 38:
                        problemSet = ETMOF38.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 150*taskNumber;
                        break;
                    case 39:
                        problemSet = ETMOF39.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 100*taskNumber;
                        break;
                    case 40:
                        problemSet = ETMOF40.getProblem(fc,sc);
                        taskNumber = problemSet.size();
                        populationSize = 100*taskNumber;
                        break;
                    default:break;
                }
                System.out.println("taskNumber = "+taskNumber);

                algorithm = new MultiTask_AE_MOEAD(problemSet);

                algorithm.setInputParameter("populationSize", populationSize);
                //algorithm.setInputParameter("maxEvaluations", populationSize*fc*(nc+1));
                algorithm.setInputParameter("maxEvaluations", populationSize*fc*(nc+1));
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
                    parameters.put("probability", 1.0 / problemSet.getMaxDimension());
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
                    System.out.println("RunID\t" + "MIGD for " + problemSet.get(0).getName()+" to "+problemSet.get(taskNumber-1).getName());


                    String[] pf = new String[nc+1];
                    for(int j=0;j<taskNumber;j++) {
                        for (int i = 0; i < nc + 1; i++) {
                            pf[i] = "PF/DynamicPF/" + problemSet.get(j).getName() + "/POF_Tt=" + i + ".txt";
                        }
                    }
                    QualityIndicator indicator;
                    int times = 5;
                    double aveIGD = 0;
                    for (int t = 1; t <= times; t++) {
                        double MIGD = 0;
                        List<SolutionSet> dynamicPopulations = algorithm.execute();
                        List<SolutionSet>[] resPopulation = new List[problemSet.size()];
                        for (int j = 0; j < dynamicPopulations.size();j++)
                        {
                            SolutionSet[] solutionSets = new SolutionSet[taskNumber];
                            for(int k=0;k<dynamicPopulations.get(j).size();k++)
                            {
                                Solution sol = dynamicPopulations.get(j).get(k);
                                int pid = sol.getSkillFactor();
                                solutionSets[pid].add(sol);
                            }
                            for(int p=0;p<problemSet.size();p++) {
                                            resPopulation[p].add(solutionSets[p]);
                            }

                        }
                        for(int j=0;j<taskNumber;j++){
                            for(int p=0;p<dynamicPopulations.size();p++){
                                indicator = new QualityIndicator(problemSet.get(j), pf[p]);
                                resPopulation[j].get(p).printObjectivesToFile("MultiTask_AEMOEAD/MultiTask_AEMOEAD_"+problemSet.get(j).getNumberOfObjectives()+"Obj_"+
                                        problemSet.get(j).getName()+ "_" + problemSet.get(j).getNumberOfVariables() + "D_run"+t+"_"+p+".txt");
                                MIGD += indicator.getIGD(resPopulation[j].get(p));
                            }
                            MIGD = MIGD/dynamicPopulations.size();
                            aveIGD += MIGD;
                            System.out.println(t + "\t" + form.format(MIGD));
                        }
                    }
                    for(int i=0;i<taskNumber;i++) {
                        System.out.println("Average MIGD for " + problemSet.get(i).getName() + ": " + form.format(aveIGD / times));
                    }
                    System.out.println();
                }
            }
}
