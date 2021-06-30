package etmo.metaheuristics.dmoea_lem.service;

import etmo.core.*;
import etmo.metaheuristics.dmoea_lem.service.DMOEA_LEM;
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

public class DMOEA_LEM_main {
    public static Logger logger_; // Logger object
    public static FileHandler fileHandler_; // FileHandler object

    public static void main(String[] args) throws JMException,
            SecurityException, IOException,ClassNotFoundException{
        logger_ = Configuration.logger_;
        fileHandler_ = new FileHandler("MOEAD_Dynamic.log");
        logger_.addHandler(fileHandler_);

        ProblemSet problemSet1; // The problem to solve
        ProblemSet problemSet2;// multitask中的某一个task的Problem
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
            problemSet1 = null;
            switch (fun){
                case 33:
                    problemSet1 = ETMOF33.getProblem(fc,sc);
                    populationSize = 100;
                    break;
                case 34:
                    problemSet1 = ETMOF34.getProblem(fc,sc);
                    populationSize = 100;
                    break;
                case 35:
                    problemSet1 = ETMOF35.getProblem(fc,sc);
                    populationSize = 100;
                    break;
                case 36:
                    problemSet1 = ETMOF36.getProblem(fc,sc);
                    populationSize = 100;
                    break;
                case 37:
                    problemSet1 = ETMOF37.getProblem(fc,sc);
                    populationSize = 105;
                    break;
                case 38:
                    problemSet1 = ETMOF38.getProblem(fc,sc);
                    populationSize = 105;
                    break;
                case 39:
                    problemSet1 = ETMOF39.getProblem(fc,sc);
                    populationSize = 100;
                    break;
                case 40:
                    problemSet1 = ETMOF40.getProblem(fc,sc);
                    populationSize = 100;
                    break;
                default:break;
            }
            int taskNumber = problemSet1.size();
            System.out.println("taskNumber = "+taskNumber);
            for (int tsk=0;tsk<taskNumber;tsk++) {

                problemSet2 = problemSet1.getTask(tsk);
                algorithm = new DMOEA_LEM(problemSet2);

                algorithm.setInputParameter("populationSize", populationSize);
                //algorithm.setInputParameter("maxEvaluations", populationSize*fc*(nc+1));
                algorithm.setInputParameter("maxEvaluations", populationSize*fc*(nc+2));
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
                parameters.put("probability", 1.0 / problemSet2.getMaxDimension());
                parameters.put("distributionIndex", 20.0);
                mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

                // Selection Operator
                parameters = null ;
                selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters) ;

                // Add the operators to the algorithm
                algorithm.addOperator("crossover", crossover);
                algorithm.addOperator("mutation", mutation);
                algorithm.addOperator("selection", selection);

                System.out.println("RunID\t" + "MIGD for " + problemSet2.get(0).getName());
                DecimalFormat form = new DecimalFormat("#.####E0");

                String[] pf = new String[nc+1];
                for(int i=0;i<nc+1;i++){
                    pf[i] = "PF/DynamicPF/" + problemSet2.get(0).getName() + "/POF_Tt=" + i + ".txt";
                }

                QualityIndicator indicator;
                int times = 1;
                double aveIGD = 0;
                for (int i = 1; i <= times; i++) {
                    double MIGD = 0;
                    List<SolutionSet> dynamicPopulations = algorithm.execute();
                    for(int p=0;p<dynamicPopulations.size();p++){
                        indicator = new QualityIndicator(problemSet2.get(0), pf[p]);
                        dynamicPopulations.get(p).printObjectivesToFile("DMOEA_LEM_Dynamic/DMOEA_LEM_"+problemSet2.get(0).getNumberOfObjectives()+"Obj_"+
                                problemSet2.get(0).getName()+ "_" + problemSet2.get(0).getNumberOfVariables() + "D_run"+i+"_"+p+".txt");
                        MIGD += indicator.getIGD(dynamicPopulations.get(p));
                    }
                    MIGD = MIGD/dynamicPopulations.size();
                    aveIGD += MIGD;
                    System.out.println(i + "\t" + form.format(MIGD));
                }
                System.out.println("Average MIGD for " + problemSet2.get(0).getName() + ": " + form.format(aveIGD / times));

                System.out.println();
            }
        }
    }
}
