package etmo.metaheuristics.dmoea_lem;

import etmo.core.SolutionSet;
import etmo.encodings.variable.Real;
import etmo.util.JMException;

import java.util.Random;

public class ModelA {
    private static final double Gaussian_radio = 200.0;

    public SolutionSet modelA(SolutionSet population) throws ClassNotFoundException, JMException {    //A.size=1 or corr<0.2
        int[] flag = Support.nonDominateOprate(population);
        for (int n = 0; n < population.size(); n++) {
            if (flag[n] == 1) continue;
            else
                population.get(n).initialize();
        }
        Random random = new Random();
        for (int n = 0; n < population.size(); n++) {
            for (int i = 0; i < population.get(0).numberOfVariables(); i++) {
                double temp = population.get(n).getDecisionVariables()[i].getValue();
                temp += random.nextGaussian() / Gaussian_radio;
                double lower = population.get(n).getDecisionVariables()[i].getLowerBound();
                double upper = population.get(n).getDecisionVariables()[i].getUpperBound();

                if (temp < lower)
                    temp = lower;
                if (temp > upper)
                    temp = upper;
                population.get(n).getDecisionVariables()[i].setValue(temp);
                if ((temp < lower) || (temp > upper))
                    population.get(n).getDecisionVariables()[i] = new Real(lower, upper);
            }
        }
        return population;
    }
}
