package etmo.metaheuristics.dmoea_lem.service;

import libsvm.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

public class SVR_test {
    public static void main(String[] args)
    {
        List<Double> label = new ArrayList<Double>();
        List<svm_node[]> nodeSet = new ArrayList<svm_node[]>();
        getData(nodeSet,label,"C:\\Users\\Oliver\\Desktop\\SVM_Test\\src\\file\\train.txt");
        int dataRange = nodeSet.get(0).length;

        //训练集向量表
        svm_node[][] datas = new svm_node[nodeSet.size()][dataRange];//训练集的向量表,初始为空
        for(int i=0;i<datas.length;i++){
            for(int j=0;j<dataRange;j++){
                datas[i][j] = nodeSet.get(i)[j];//将文件读取到的训练赋给datas
            }
        }
        //训练集向量表对应的labelvalue
        double[] labels = new double[label.size()];//训练集向量表对应的label数组,初始为空
        for(int i=0;i<labels.length;i++){
            labels[i] = label.get(i);//将文件读取到的label复制给labels
        }

        //定义svm_problem对象
        svm_problem problem = new svm_problem();
        problem.l = nodeSet.size();//向量个数
        problem.x = datas;//训练集向量表
        problem.y = labels;//对应的label
        //定义svm_problem对象
        svm_parameter parameter = new svm_parameter();
        parameter.svm_type = svm_parameter.EPSILON_SVR;//svm类型为支持向量回归
        parameter.kernel_type = svm_parameter.RBF;//线性回归
        parameter.cache_size = 100;
        parameter.eps = 0.00001;
        parameter.C = 1.9;

        //训练SVM分类模型
        System.out.println(svm.svm_check_parameter(problem, parameter));// 如果参数没有问题，则返回null,否则返回error描述。
        svm_model model = svm.svm_train(problem, parameter);// svm.svm_train()训练出SVM分类模型

        // 获取测试数据
        List<Double> testlabel = new ArrayList<Double>();
        List<svm_node[]> testnodeSet = new ArrayList<svm_node[]>();
        getData(testnodeSet, testlabel, "C:\\Users\\Oliver\\Desktop\\SVM_Test\\src\\file\\test.txt");

        svm_node[][] testdatas = new svm_node[testnodeSet.size()][dataRange]; // 测试集的向量表，初始为空
        for (int i = 0; i < testdatas.length; i++) {
            for (int j = 0; j < dataRange; j++) {
                testdatas[i][j] = testnodeSet.get(i)[j];//读取文件复制给testdatas
            }
        }
        double[] testlables = new double[testlabel.size()]; // 测试集对应的label(真实值)
        for (int i = 0; i < testlables.length; i++) {
            testlables[i] = testlabel.get(i);
        }

        // 预测测试数据的lable
        double err = 0.0;
        for (int i = 0; i < testdatas.length; i++) {
            double truevalue = testlables[i];
            System.out.print(truevalue + " ");
            double predictValue = svm.svm_predict(model, testdatas[i]);//通过回归模型预测的label
            System.out.println(predictValue);
            err += Math.abs(predictValue - truevalue);
        }
        System.out.println("err=" + err / datas.length);
    }
    public static void getData(List<svm_node[]> nodeSet, List<Double> label, String filename) {
        try {
            FileReader fr = new FileReader(new File(filename));
            BufferedReader br = new BufferedReader(fr);
            String line = null;
            while ((line = br.readLine()) != null) {
                String[] datas = line.split(",");
                svm_node[] vector = new svm_node[datas.length - 1];
                for (int i = 0; i < datas.length - 1; i++) {
                    svm_node node = new svm_node();
                    node.index = i + 1;
                    node.value = Double.parseDouble(datas[i]);
                    vector[i] = node;
                }
                nodeSet.add(vector);
                double lablevalue = Double.parseDouble(datas[datas.length - 1]);
                label.add(lablevalue);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
