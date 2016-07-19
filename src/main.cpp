#include <iostream>
#include <random>
#include <time.h>
#include <stdlib.h>
#include "scattering.h"
#include "sample.h"


using namespace std;




int main()
{

    Scattering_2to3 A(1.3, 0.2);
    double E1, temp, s;

    temp = 0.4;
    s = 15.0;
    E1 = 15.;
// okay, let us test the result of sampling p3
    Sample_2to3 B(A,E1, temp);
    for (int i=0; i< 10000; i++)
        B.Sample_p3_Qq(A, s, temp);

/*
    Scattering_2to2 A(1.3, 0.2);
    A.Sigma_Qq_tabulate();
    A.Gamma_Qq_tabulate();

    A.Sigma_Qg_tabulate();
    A.Gamma_Qg_tabulate();

    double temp = 0.4;
    double E1 = 15.;
    Sample_2to2 B(A, E1, temp);
 
 
*/


 
    return 0;
}






// okay, this part is just for a test, if we combining importance sample and the rejection sampling, 
//see if we can increase our efficiency of sampling by using an exponential distribution of x and y, instead of the uniform distribution.
// but the result turns out there are not much different between those two. while for uniform distribution, the efficiency is 5440/100000, to sample 10000 (E2,s), takes 5.583s
// for the exponential distribution, the efficiency is 4298/100000, to sample 10000 (E2, s), takes 5.434 s

/*
double Sample_E2ANDs_max2(Scattering_2to2 A, double E1, double temp)
{
    double alpha = 2.*(E1 -
    
    
    sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));
    double dE2 = 0.01;
    double E2 = 10* temp;

    double Fmax = exp(E2/temp)/(exp(E2/temp) + 1) * beta*E2 * A.Sigma_Qq(A.Mass2_ + beta*E2, temp);
    return Fmax;
}


int Valid_xy2(double x, double y, double E1, double temp, double M2)
{
    int ctl1 = x+y>0 ? 1:0;
    int ctl2 = x+y<10 ? 1:0;
    double xmin = -0.24/(temp/M2 * sqrt(E1*E1 - M2));
    int ctl3 = x > xmin ? 1:0;
    double ymin = -xmin;
    int ctl4 = y > ymin ? 1 : 0;
    return ctl1*ctl2*ctl3*ctl4;
}

int Sample_E2ANDs2(Scattering_2to2 A, double E1, double temp)
{
    double alpha = 2.*(E1 - sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));

    double Fmax = Sample_E2ANDs_max2(A, E1, temp);

    std::default_random_engine generator;
    std::exponential_distribution<double> distribution(1.0);
   
    int i = 0;
    while(i < 10000)
//    for (int i=0; i< 100000; i++)
    {
    int status = 0;
    double x, y;
    while ( !status)
    {
        x = distribution(generator);
        y = distribution(generator);
        x = -log(x);
        y = -log(y);
        status = Valid_xy2(x, y, E1, temp, A.Mass2_);
    }


    double s = temp*(alpha*x + beta*y);
    double E2 = temp*(x+y);
    double result = exp(E2/temp) / (exp(E2/temp + 1)) * (s-A.Mass2_) * A.Sigma_Qq(s, temp);

    double ctl = ((double) rand() / RAND_MAX);
    if (ctl < result/Fmax)  i=i+1;
       //out << x << " " << y << " " << E2 << " " << s << " " << result << endl;
    
    }

    cout << i << endl;
    return 0;
}
*/


