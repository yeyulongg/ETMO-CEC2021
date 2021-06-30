package etmo.metaheuristics.dmoea_lem;

import etmo.core.SolutionSet;
import etmo.util.JMException;

import java.util.Random;

public class ModelC {
    private static final double Gaussian_radio = 200.0;

    public SolutionSet modelC(SolutionSet population, double[] vector1, double[] vector2) throws JMException, ClassNotFoundException {

        Random random =new Random();
        int[] flag = Support.nonDominateOprate(population);
        for(int n=0;n<population.size();n++)
        {
            Random ran =new Random();
            double r = ran.nextDouble();
            double[] vector=new double[vector1.length];
            for(int i=0;i<vector.length;i++)
                vector[i]=vector2[i]+r*(vector1[i]-vector2[i]);
            for(int i=0;i<population.get(0).numberOfVariables();i++)
            {
                double temp=population.get(n).getDecisionVariables()[i].getValue()+vector[i];
                temp+=random.nextGaussian()/Gaussian_radio;
                double lower = population.get(n).getDecisionVariables()[i].getLowerBound();
                double upper = population.get(n).getDecisionVariables()[i].getUpperBound();
                if(temp<lower)
                    temp=lower;
                if(temp>upper)
                    temp = upper;
                population.get(n).getDecisionVariables()[i].setValue(temp);
            }
        }
        return population;
    }
}
