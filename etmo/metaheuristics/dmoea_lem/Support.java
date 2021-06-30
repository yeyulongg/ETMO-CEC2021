package etmo.metaheuristics.dmoea_lem;

import java.awt.Menu;
import java.util.List;
import java.util.Random;

import etmo.core.SolutionSet;
import etmo.core.SolutionSet;
import etmo.util.JMException;

public class Support {

    public static double distVector(double[] vector1, double[] vector2) {
        int dim = vector1.length;
        double sum = 0;
        for (int n = 0; n < dim; n++)
            sum += Math.pow((vector1[n] - vector2[n]),2);
        return Math.sqrt(sum);
    }

    public static void minSort(double x[], int index[])
    {
        int[] isMin = new int[x.length];
        for(int n=0;n<x.length;n++)
            isMin[n] = 0;
        double min = 100000000.0;
        for (int i = 0; i < index.length; i++)
        {
            min = 100000000.0;
            for (int j = 0; j < x.length; j++)
            {
                if(isMin[j]==1) continue;
                else if(x[j]<min)
                {
                    min=x[j];
                    index[i]=j;
                }
            }
            isMin[index[i]]=1;
        } // for

    } // minFastSort


    public static int[] randomSort(int n) {
        int result[] = new int[n];
        for(int i=0;i<n;i++)
            result[i]=i;
        Random rand =new Random();
        for(int m=n;m>0;m--)
        {
            int k=rand.nextInt(m);
            if(k==(m-1))
                continue;
            else
            {
                int temp=result[m-1];
                result[m-1]=result[k];
                result[k]=temp;
            }
        }
        return result;
    }
    //计算种群所有个体的质点
    public static double[] compCenter(SolutionSet population_) throws JMException
    {
        double[] result= new double[population_.get(0).numberOfVariables()];
        double sum=0.0;
        for(int j =0;j<population_.get(0).numberOfVariables();j++)
        {
            sum=0.0;
            for(int i=0;i<population_.size();i++)
                sum+=population_.get(i).getDecisionVariables()[j].getValue();
            result[j] = sum/population_.get(0).numberOfVariables();
        }
        return result;
    }
    //计算种群中非支配个体的质点
    public static double[] compCenter_nondomin(SolutionSet population_) throws JMException
    {
        double[] result= new double[population_.get(0).numberOfVariables()];
        double sum=0.0;
        int count=0;
        int flag[] = nonDominateOprate(population_);
        for(int j =0;j<population_.get(0).numberOfVariables();j++)
        {
            sum=0.0;
            count=0;
            for(int i=0;i<population_.size();i++){
                if(flag[i]==0) continue;
                else{
                    sum+=population_.get(i).getDecisionVariables()[j].getValue();
                    count++;
                }
            }
            result[j] = sum/count;
        }
        return result;
    }

    public static int[] nonDominateOprate(SolutionSet population_)
    {
        int length = population_.size();
        int dim=population_.get(0).getNumberOfObjectives();
        int[] flag= new int[length];
        int judgment=0;
        for(int i=0;i<length;i++)
            flag[i]=1;
        for(int i=0;i<length;i++)
        {
            if(flag[i]==0) continue;
            for(int j=0;j<length;j++)
            {
                if((i==j)||(flag[j]==0)) continue;
                else
                {
                    judgment=0;
                    for(int k=0;k<dim;)
                    {
                        if(population_.get(i).getObjective(k)<population_.get(j).getObjective(k))
                            judgment++;
                        k++;
                        if(judgment<k)
                            break;
                    }
                    if(judgment==dim)
                        flag[j]=0;
                }
            }
        }
        return flag;
    }

    public static double verticalDistance(double[] vector,double[] weight)
    {
        int size = vector.length;
        double sum1=0.0,sum2=0.0;
        for(int n=0;n<size;n++)
        {
            sum1+=weight[n]*vector[n];
            sum2+=weight[n]*weight[n];
        }
        double value = sum1/sum2;
        double[] temp=new double[size];
        sum1=0.0;
        for(int n=0;n<size;n++)
        {
            temp[n]=vector[n]-weight[n]*value;
            sum1+=temp[n]*temp[n];
        }
        return Math.sqrt(sum1);
    }

    public static double[] DVector(double[] vector1,double[] vector2)
    {
        int length = vector1.length;
        double[] vector = new double[length];
        for (int n = 0; n < length; n++)
            vector[n] = vector2[n] - vector1[n];
        return vector;
    }

    public static double[] reVector(double[] vector_)   //反向的
    {
        int length = vector_.length;
        double[] vector = new double[length];
        for (int n = 0; n < length; n++)
            vector[n] =  -1.0*vector_[n];
        return vector;
    }

    public static double[][] comp(int p,double s)
    {
        double[][] A = new double[p][p];
        for(int i=0;i<p;i++)
            for(int j=0;j<p;j++)
            {
                if(i==j) A[j][i]=1.0;
                else if(i>j) A[j][i]=0.0;
                else
                    A[j][i]=(double)Cfunction(i, j)*Math.pow(s, j-i);
            }
        return A;
    }

    public static int Cfunction(int i,int j)
    {
        if(i==0)
            return 1;
        else if(i==1)
            return j;
        else
        {
            int a=Afunction(i, j);
            int b=Afunction(i, i);
            int c=a/b;
            return c;
        }
    }

    public static int Afunction(int i,int j)
    {
        int a=i,b=j;
        int c=b;
        while(a>1)
        {
            b--;
            c*=b;
            a--;
        }
        return c;
    }

    public static double[] oldRadio(double old[][])
    {
        double[] result = new double[old.length];
        int[] index = new int[old.length];
        for(int i=0;i<old.length;i++)
        {
            index[i]=maxIndex(old[i]);
        }
        if(old.length==2)
        {
            double x1=old[0][index[0]];
            double x2=old[0][index[1]];
            double y1=old[1][index[0]];
            double y2=old[1][index[1]];
            result[0]=1.0;
            result[1]=-(x2-x1)/(y2-y1);
        }
        else if(old.length==3)
        {
            double x1=old[0][index[0]];
            double x2=old[0][index[1]];
            double x3=old[0][index[2]];
            double y1=old[1][index[0]];
            double y2=old[1][index[1]];
            double y3=old[1][index[2]];
            double z1=old[2][index[0]];
            double z2=old[2][index[1]];
            double z3=old[2][index[2]];
            double a = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
            double b = (z2-z1)*(x3-z1)-(x2-x1)*(z3-z1);
            double c = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
            result[0] = 1.0;
            result[1] = b/a;
            result[2] = c/a;

        }

        return result;
    }

    public static double corr(double x[],double y[])
    {
        double Ex=meanV(x);
        double Ey=meanV(y);
        double dx=ssDs(x);
        double dy=ssDs(y);
        double z[]=new double[x.length];
        for(int i=0;i<x.length;i++)
        {
            x[i]=x[i]-Ex;
            y[i]=y[i]-Ey;
            z[i]=x[i]*y[i];
        }
        double Ez=meanV(z);
        return Ez/(dx*dy);
    }

    public static double[] Ds(double radio[],double x[][])
    {
        double[] ds = new double[x[0].length];
        for(int i=0;i<ds.length;i++)
        {
            double temp = 0.0;
            for(int j=0;j<x.length;j++)
                temp+=x[j][i]*radio[j];
            ds[i]=-1.0*temp;
        }

        return ds;
    }

    public static double[] staDs(double radio[],double x[][])
    {
        double[] ds = new double[x[0].length];
        for(int i=0;i<ds.length;i++)
        {
            double temp = 0.0;
            for(int j=0;j<x.length;j++)
                temp+=x[j][i]*radio[j];
            ds[i]=-1.0*temp;
        }
        //sta
        double mean = meanV(ds);
        for(int i=0;i<ds.length;i++)
            ds[i]-=mean;
        int max = maxIndex(ds);
        int min = minIndex(ds);
        double len = ds[max]-ds[min];
        for(int i=0;i<ds.length;i++)
            ds[i]=(ds[i]-ds[min]);

        return ds;
    }

    public static double fFunc(double[] x,double[] y)
    {

        double sumX1=0.0;
        double sumX2=0.0;
        double sumY1=0.0;
        double sumY2=0.0;

        int m=x.length;

        for(int i=0;i<m;i++)
        {
            sumX1+=x[i];
            sumX2+=x[i]*x[i];
            sumY1+=y[i];
            sumY2+=y[i]*y[i];
        }

        double a0= (sumX2-sumX1*sumX1/m+sumY2-sumY1*sumY1/m)/(2.0*(m-1));
        double a1=sumX2/m+sumY2/m-Math.pow((sumX1/m+sumY1/m), 2.0);
        return m*a1/a0;
    }

    public static double ssDs(double ds[])
    {
        double mean = meanV(ds);
        double sum=0.0;
        for(int i=0;i<ds.length;i++)
            sum+=Math.pow(ds[i]-mean,2.0);
        return Math.sqrt(sum/ds.length);

    }

    public static double lengthDs(double ds[])
    {
        int max = maxIndex(ds);
        int min = minIndex(ds);
        double len = ds[max]-ds[min];
        return len;
    }

    public static int maxIndex(double x[])
    {
        int index=-1;
        double max=-1.0*1000000;
        for(int i=0;i<x.length;i++)
            if(x[i]>max)
            {
                index=i;
                max=x[i];
            }
        return index;

    }

    public static int minIndex(double x[])
    {
        int index=-1;
        double min=10E8;
        for(int i=0;i<x.length;i++)
            if(x[i]<min)
            {
                index=i;
                min=x[i];
            }
        return index;

    }

    public static double meanV(double x[])
    {

        double sum=0.0;
        for(int i=0;i<x.length;i++)
            sum+=x[i];
        return sum/x.length;

    }

    public static double[] radioFunc(double[][] x,double[][] y)
    {
        double[] result=new double[x[0].length];
        for(int i=0;i<x[0].length;i++)
        {
            if((y[0][i]-x[0][i])==0)
                result[i]=1.0;
            else
                result[i]=(y[1][i]-x[1][i])/(y[0][i]-x[0][i]);
            System.out.println(y[1][i]-x[1][i]);
            System.out.println(y[0][i]-x[0][i]);
        }
        double min=minIndex(result);
        double max=maxIndex(result);
        for(int i=0;i<x[0].length;i++)
        {
            result[i]=(result[i]-min)/(max-min);
//			System.out.println(result[i]);
        }

        return result;
    }

    public static double[] comLinePara(List<double[]> centers_)
    {
        int var=centers_.get(0).length;
        int len=centers_.size();
        double[] para = new double[2];
        double[] center=new double[var];

        double sumx;
        double sumy;
        double sumxy;
        double sumxx;

        for(int i=0;i<var;i++)
        {
            sumx=0.0;
            sumy=0.0;
            sumxy=0.0;
            sumxx=0.0;

            for(int n=0;n<len;n++)
            {
                sumx+=(double)n;
                sumxx+=(double)n*n;
                sumxy+=(double)n*centers_.get(n)[i];
                sumy+=centers_.get(n)[i];
            }
            para[0]=(sumxx*sumy-sumx*sumxy)/((double)len*sumxx-Math.pow(sumx, 2.0));
            para[1]=((double)len*sumxy-sumx*sumy)/((double)len*sumxx-sumx);
            center[i]=para[0]+para[1]*(double)(len+1.0);
        }

        return center;
    }
}

