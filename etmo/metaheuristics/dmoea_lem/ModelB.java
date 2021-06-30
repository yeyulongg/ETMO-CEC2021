package etmo.metaheuristics.dmoea_lem;

import etmo.core.SolutionSet;
import etmo.util.JMException;

import java.util.Random;

public class ModelB {
    private static final double Gaussian_radio = 200.0;
    public SolutionSet modelB(SolutionSet population,double[] vector)throws JMException {

        Random random =new Random();
        for(int n=0;n<population.size();n++)
        {
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
