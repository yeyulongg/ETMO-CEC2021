package etmo.metaheuristics.Multitask_MOEAD;

import etmo.core.*;
import etmo.metaheuristics.dmoea_lem.Autoencoding;
import etmo.metaheuristics.dmoea_lem.Utils;
import etmo.util.Distance;
import etmo.util.JMException;
import etmo.util.PseudoRandom;
import etmo.util.Ranking;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

public class MultiTask_AE_MOEAD extends DynamicAlgorithm {


        /**
         * Constructor
         *
         * @param problemSet
         * Problem to solve
         */

        private static final double EPS = 1.0e-8;
        int populationSize;
        int maxEvaluations;
        int evaluations;
        int generations;
        // population repository
        SolutionSet population;
        SolutionSet tempPopulation;
        List<SolutionSet> dynamicPopulationSets;

        SolutionSet offspringPopulation;
        SolutionSet union;
        Distance distance = new Distance();
        double[][] rmpMatrix;
        int[] numVars;


        Operator mutationOperator;
        Operator crossoverOperator;
        Operator selectionOperator;

        double[] z_;//Z向量，理想点
        double[][] lambda_;//lambda权重向量
        int T_;//邻居个数
        int[][] neighborhood_;
        double delta_;

        //nr: maximal number of solutions replaced by each child solution

        int nr_;
        Solution[] indArray_;//存放m个个体，这m个个体中，每个个体的某一个目标函数值是最小的，也就是m个最小值构成理想点，这些最小值存在Z_当中
        String functionType_;//函数类型，以下使用的是TCHE1 为切比雪夫
        String dataDirectory_;

        public MultiTask_AE_MOEAD(ProblemSet problemSet) {
            super(problemSet);
            functionType_ = "_TCHE1";
            //	System.out.println("sup: " + problemSet_.get(0).getHType());
        } // DMOEA_LEM

        @Override
        public List<SolutionSet> execute() throws JMException, ClassNotFoundException,NullPointerException {


            /*
             * fc indicates the frequency of change
             * nc number of changes in each run
             */
            int fc = 20;
            int nc = 30;

            //Read the parameters
            maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
            populationSize = ((Integer) this.getInputParameter("populationSize")).intValue();
            dataDirectory_ = this.getInputParameter("dataDirectory").toString();
            T_ = ((Integer) this.getInputParameter("T")).intValue();
            nr_ = ((Integer) this.getInputParameter("nr")).intValue();
            delta_ = ((Double) this.getInputParameter("delta")).doubleValue();
            neighborhood_ = new int[populationSize][T_];//初始化邻居空间
            z_ = new double[problemSet_.get(0).getNumberOfObjectives()];//初始化理想点空间
            lambda_ = new double[populationSize][problemSet_.get(0).getNumberOfObjectives()];

            //Initialize the variables
            population = new SolutionSet(populationSize);
            dynamicPopulationSets = new ArrayList<SolutionSet>();
            problemSet_.get(0);
            indArray_ = new Solution[problemSet_.get(0).getNumberOfObjectives()];
            // Read the operators
            mutationOperator = operators_.get("mutation");
            crossoverOperator = operators_.get("crossover");
            selectionOperator = operators_.get("selection");


            numVars = new int[problemSet_.size()];
            for(int i=0; i<problemSet_.size(); i++) {
                numVars[i] = problemSet_.get(i).getNumberOfVariables();
            }
            evaluations = 0;
            generations = 1;

            // STEP 1. Initialization
            // STEP 1.1. Compute euclidean distances between weight vectors and find T
            initUniformWeight();
            initNeighborhood();
            //STEP 1.2. Initialize population
            initPopulation();
            // STEP 1.3. Initialize z_
            initIdealPoint();
            // STEP 2. Update
            // Generations
            while (evaluations < maxEvaluations) {
                //保留每个环境的最后一代种群
                if((generations) % fc == 0 ){
                    SolutionSet solutionSet = new SolutionSet();
                    for(int i=0;i<population.size();i++){
                        solutionSet.add(population.get(i));
                    }
                    dynamicPopulationSets.add(solutionSet);
                }

                generations++;

                if((generations)%fc==1&&generations>fc){
                    int size = dynamicPopulationSets.size();
                    if(size>=2){
                        //AE-Prediction
                        Ranking ranking = new Ranking(dynamicPopulationSets.get(size-1));
                        SolutionSet curr_NDS = ranking.getSubfront(0);
                        Ranking ranking1 = new Ranking(dynamicPopulationSets.get(size-2));
                        SolutionSet his_NDS = ranking1.getSubfront(0);
                        population = new Autoencoding().ae_Prediction(population,curr_NDS,his_NDS,populationSize,problemSet_.get(0),generations);
                        //population = new SVR_DMOEA().SVR_MOEAD(dynamicPopulationSets,problemSet_.get(0),population);
                        //多响应预测模型预测下一个环境的种群
                        //population = new Multi_Respond_Prediction().multi_Respond_Prediction((ArrayList<SolutionSet>) dynamicPopulationSets,population);
                        //多方向预测模型
                        //population = new MultiDirection_Prediction().multiDirection_Prediction(dynamicPopulationSets,population);

                        for(int i=0;i<population.size();i++){
                            problemSet_.get(0).dynamicEvaluate(population.get(i),generations);
                            evaluations++;
                        }
                    }
                    else{
                        population.clear();
                        initPopulation();
                    }
                }

                int[] permutation = new int[populationSize];
                Utils.randomPermutation(permutation, populationSize);
                for (int i = 0; i < populationSize; i++) {
                    int n = permutation[i];
                    int type;
                    double rnd = PseudoRandom.randDouble();

                    // STEP 2.1. Mating selection based on probability
                    if (rnd < delta_) // if (rnd < realb)
                    {
                        type = 1;   // neighborhood
                    } else {
                        type = 2;   // whole population
                    }
                    Vector<Integer> p = new Vector<Integer>();
                    matingSelection(p, n, 2, type);//选交配池，选出放入P中

                    // STEP 2.2. Reproduction
                    Solution child;
                    Solution[] parents = new Solution[3];
                    parents[0] = population.get(p.get(0));
                    parents[1] = population.get(p.get(1));
                    parents[2] = population.get(n);
                    int[] skillfactors = new int[3];
                    skillfactors[0] = parents[0].getSkillFactor();
                    skillfactors[1] = parents[1].getSkillFactor();
                    skillfactors[2] = parents[2].getSkillFactor();


                        // Apply DE crossover
                        child = (Solution) crossoverOperator.execute(new Object[]{population.get(n), parents});
                        // Apply mutation
                        mutationOperator.execute(child);
                        // Evaluation-算出新个体函数值
                    if(skillfactors[2]==skillfactors[1]&&skillfactors[2]==skillfactors[0]){
                        problemSet_.get(skillfactors[2]).dynamicEvaluate(child,generations);
                        evaluations++;
                        child.setSkillFactor(skillfactors[2]);
                    }else {

                    }

                        // STEP 2.3. Repair. Not necessary
                        // STEP 2.4. Update z_
                        updateReference(child);
                        // STEP 2.5. Update of solutions
                        updateProblem(child, n, type);//更新种群，子代代替父代

                }//for
            }
            if(dynamicPopulationSets.size() != nc){
                System.out.println("The real number of changes is " + dynamicPopulationSets.size() +
                        " which is not the same with the preset " + nc);
            }
            return dynamicPopulationSets;
        }
        public void initUniformWeight () {
            if ((problemSet_.get(0).getNumberOfObjectives() == 2) && (populationSize <= 300)) {
                for (int n = 0; n < populationSize; n++) {
                    double a = 1.0 * n / (populationSize - 1);
                    lambda_[n][0] = a;
                    lambda_[n][1] = 1 - a;
                } // for
            } // if
            else {
                String dataFileName;
                dataFileName = "W" + problemSet_.get(0).getNumberOfObjectives() + "D_" +
                        populationSize + ".dat";
                try {
                    // Open the file
                    FileInputStream fis = new FileInputStream(dataDirectory_ + "/" + dataFileName);
                    InputStreamReader isr = new InputStreamReader(fis);
                    BufferedReader br = new BufferedReader(isr);
                    int numberOfObjectives = 0;
                    int i = 0;
                    int j = 0;
                    String aux = br.readLine();
                    while (aux != null) {
                        StringTokenizer st = new StringTokenizer(aux);
                        j = 0;
                        numberOfObjectives = st.countTokens();
                        while (st.hasMoreTokens()) {
                            double value = (new Double(st.nextToken())).doubleValue();
                            lambda_[i][j] = value;
                            //System.out.println("lambda["+i+","+j+"] = " + value) ;
                            j++;
                        }
                        aux = br.readLine();
                        i++;
                    }
                    br.close();
                } catch (Exception e) {
                    System.out.println("initUniformWeight: failed when reading for file: " + dataDirectory_ + "/" + dataFileName);
                    e.printStackTrace();
                }
            } // else
            //System.exit(0) ;
        }
        public void initNeighborhood () {
            double[] x = new double[populationSize];
            int[] idx = new int[populationSize];
            for (int i = 0; i < populationSize; i++) {
                // calculate the distances based on weight vectors 计算向量间的欧式距离
                for (int j = 0; j < populationSize; j++) {
                    x[j] = Utils.distVector(lambda_[i], lambda_[j]);
                    idx[j] = j;
                }
                // find 'niche' nearest neighboring subproblems
                Utils.minFastSort(x, idx, populationSize, T_);
                System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
            }
        }
        public void initPopulation () throws JMException, ClassNotFoundException {
            population = new SolutionSet(populationSize);
            for(int i = 0;i < populationSize; i++){
                Solution newSolution = new Solution(problemSet_);
                int id = i % problemSet_.size();
                problemSet_.get(id).dynamicEvaluate(newSolution,generations);
                evaluations++;
                newSolution.setSkillFactor(id);
                population.add(newSolution);
            }
        }
        void initIdealPoint () throws JMException, ClassNotFoundException {
            for (int i = 0; i < problemSet_.get(0).getNumberOfObjectives(); i++) {
                z_[i] = 1.0e+30;
                indArray_[i] = new Solution(problemSet_);
                problemSet_.get(0).dynamicEvaluate(indArray_[i],generations);
                evaluations++;
            }
            for (int i = 0; i < populationSize; i++) {
                updateReference(population.get(i));
            }
        }
        void updateReference (Solution individual){
            for (int n = 0; n < problemSet_.get(0).getNumberOfObjectives(); n++) {
                if (individual.getObjective(n) < z_[n]) {
                    z_[n] = individual.getObjective(n);
                    indArray_[n] = individual;
                }
            }
        }
        public void matingSelection(Vector<Integer> list, int cid, int size, int type) {
            // list : the set of the indexes of selected mating parents
            // cid  : the id of current subproblem
            // size : the number of selected mating parents
            // type : 1 - neighborhood; otherwise - whole population
            int ss;
            int r;
            int p;
            ss = neighborhood_[cid].length;
            while (list.size() < size) {
                if (type == 1) {
                    r = PseudoRandom.randInt(0, ss - 1);
                    p = neighborhood_[cid][r];//从cid编号的向量中选邻居r
                    //p = population[cid].table[r];
                } else {
                    p = PseudoRandom.randInt(0, populationSize - 1);
                }
                boolean flag = true;
                for (int i = 0; i < list.size(); i++) {
                    if (list.get(i) == p) // p is in the list
                    {
                        flag = false;
                        break;
                    }
                }
                //if (flag) list.push_back(p);
                if (flag) {
                    list.addElement(p);
                }
            }
        }
        void updateProblem(Solution indiv, int id, int type) {
            // indiv: child solution
            // id:   the id of current subproblem
            // type: update solutions in - neighborhood (1) or whole population (otherwise)
            int size;
            int time;
            time = 0;
            if (type == 1) {
                size = neighborhood_[id].length;
            } else {
                size = population.size();
            }
            int[] perm = new int[size];
            Utils.randomPermutation(perm, size);
            for (int i = 0; i < size; i++) {
                int k;
                if (type == 1) {
                    k = neighborhood_[id][perm[i]];
                } else {
                    k = perm[i];      // calculate the values of objective function regarding the current subproblem
                }
                double f1, f2;
                f1 = fitnessFunction(population.get(k), lambda_[k]);
                f2 = fitnessFunction(indiv, lambda_[k]);
                if (f2 < f1) {
                    population.replace(k, new Solution(indiv));
                    time++;
                }
                // the maximal number of solutions updated is not allowed to exceed 'limit'
                if (time >= nr_) {
                    return;
                }
            }
        }
        double fitnessFunction(Solution individual, double[] lambda) {
            double fitness;
            fitness = 0.0;
            if (functionType_.equals("_TCHE1")) {
                double maxFun = -1.0e+30;
                for (int n = 0; n < problemSet_.get(0).getNumberOfObjectives(); n++) {
                    double diff = Math.abs(individual.getObjective(n) - z_[n]);
                    double feval;
                    if (lambda[n] == 0) {
                        feval = 0.0001 * diff;
                    } else {
                        feval = diff * lambda[n];
                    }
                    if (feval > maxFun) {
                        maxFun = feval;
                    }
                }
                fitness = maxFun;
            }
            else {
                System.out.println("MOEAD.fitnessFunction: unknown type " + functionType_);
                System.exit(-1);
            }
            return fitness;
        }
}
