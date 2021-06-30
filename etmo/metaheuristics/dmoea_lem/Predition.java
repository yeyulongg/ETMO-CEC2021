package etmo.metaheuristics.dmoea_lem;

import etmo.core.SolutionSet;

import javax.management.JMException;
import java.io.IOException;

public class Predition {
    protected int numberOfObjectives;

    protected int numberOfVariables;

    protected String problemName;

    protected int index;

    protected int change;

    public Predition(SolutionSet population){
        if(population.size()<1)
            System.out.println("Size error-Prediction!");
        else
        {
            this.numberOfObjectives=population.get(0).getNumberOfObjectives();
            this.numberOfVariables=population.get(0).numberOfVariables();
            this.problemName=population.get(0).getProblemSet().get(0).getName();
        }
    }
    public void printIndex()
    {
        System.out.println(this.index);
    }

    public void dynamicResponse(SolutionSet population) throws JMException, ClassNotFoundException{};
    public void dynamicResponse(SolutionSet population,SolutionSet population_1) throws JMException, ClassNotFoundException{};
    public void dynamicResponse(SolutionSet population,SolutionSet population_1,SolutionSet population_2) throws JMException, ClassNotFoundException{};    //ira产生的population_,new_,old_

    public void update(SolutionSet population_) throws JMException, etmo.util.JMException {};

    public void printCenter() throws IOException {};

    public void printVars() throws IOException{};
}
