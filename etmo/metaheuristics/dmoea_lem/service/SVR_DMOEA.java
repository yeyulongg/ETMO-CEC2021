package etmo.metaheuristics.dmoea_lem.service;

import etmo.core.Problem;
import etmo.core.SolutionSet;
import etmo.util.JMException;
import libsvm.*;

import java.text.DecimalFormat;
import java.util.List;

public class SVR_DMOEA {
    public SolutionSet SVR_MOEAD(List<SolutionSet> DynamicPopulationSet, Problem problem, SolutionSet population) throws JMException, ClassNotFoundException {

        int popSize = population.size();
        //先对种群重新初始化
        for(int i=0;i<popSize;i++){
            population.get(i).initialize();
        }
        int d = population.get(0).numberOfVariables();
        int T = DynamicPopulationSet.size();
        int od = 4;
        DecimalFormat df = new DecimalFormat("0.0000");
        for(int i=0;i<popSize;i++){
            for(int j=0;j<d;j++){
                svm_node[][] datas = new svm_node[T-od][od];//训练集向量表,初始为空
                double[] labels = new double[T-od];

                for(int jj=0;jj<T-od;jj++){
                    int train_iter = 0;
                    for(int kk=jj;kk<=jj+od;kk++){
                        if(kk!=jj+od){
                            svm_node node = new svm_node();
                            node.index = train_iter+1;
                            //将数据小数保留到四位
                            double p = DynamicPopulationSet.get(kk).get(i).getDecisionVariables()[j].getValue();
                            node.value = Double.parseDouble(df.format(p));
                            datas[jj][train_iter] = node;
                            train_iter++;
                        }
                        else {
                            double p = DynamicPopulationSet.get(kk).get(i).getDecisionVariables()[j].getValue();
                            labels[jj] = Double.parseDouble(df.format(p));
                        }
                    }
                }
                //定义svm_problem对象
                svm_problem svr_problem = new svm_problem();
                svr_problem.l = T-od;//向量个数
                svr_problem.x = datas;//训练集向量表
                svr_problem.y = labels;//对应的label
                //定义svm_problem对象
                svm_parameter parameter = new svm_parameter();
                parameter.svm_type = svm_parameter.EPSILON_SVR;//svm类型为支持向量回归
                parameter.kernel_type = svm_parameter.RBF;//核函数为rbf
                parameter.cache_size = 200;
                parameter.eps = 0.05;
                parameter.C = 1000;
                parameter.shrinking = 0;
                svm_model model = svm.svm_train(svr_problem,parameter);// svm.svm_train()训练出SVM分类模型

                svm_node[][] testdatas = new svm_node[1][od]; // 测试集的向量表，初始为空
                int test_iter=0;
                for(int kk=T-od;kk<=T-1;kk++){
                    svm_node node = new svm_node();
                    node.index = test_iter+1;
                    double p = DynamicPopulationSet.get(kk).get(i).getDecisionVariables()[j].getValue();
                    node.value = Double.parseDouble(df.format(p));
                    testdatas[0][test_iter] = node;
                    test_iter++;
                }

                double predictValue = svm.svm_predict(model,testdatas[0]);
                double lower = population.get(i).getDecisionVariables()[j].getLowerBound();
                double upper = population.get(i).getDecisionVariables()[j].getUpperBound();
                if(predictValue<lower)
                    population.get(i).getDecisionVariables()[j].setValue(lower);
                if(predictValue>upper)
                    population.get(i).getDecisionVariables()[j].setValue(upper);
                population.get(i).getDecisionVariables()[j].setValue(predictValue);
            }
        }
        return population;
    }
}
