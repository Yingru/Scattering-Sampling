#include "sample.h"

using namespace std;

Sample_2to2::Sample_2to2(Scattering_2to2& A, double E1, double temp)
{
    srand(std::random_device{}());
    E1_ = E1;
    temp_ = temp;
   
//    double SR_Qq = A.Gamma_Qq(E1_, temp_); // in the future, if I need to speed up, I can change
    					   // this into A.Get_Gamma_Qq()
//    double SR_Qg = A.Gamma_Qg(E1_, temp_);

    double SR_Qq = A.Get_gamma_Qq(E1_, temp_);
    double SR_Qg = A.Get_gamma_Qg(E1_, temp_);

    double ctl_sample = ((double) rand()/RAND_MAX);
   // cout << SR_Qq << " " << SR_Qg << " " <<  SR_Qq/(SR_Qq + SR_Qg) << endl;
    if (ctl_sample < SR_Qq/(SR_Qq+SR_Qg))
    {
       // cout << "Qq->Qq will soon start!" << endl;
	int status_s = Sample_E2ANDs_Qq(A, E1_, temp_);
	int status_t = Sample_t_Qq(A, s_, temp_);
	cout << 1 << " "<< ctl_sample << " " << E2_ << " " << s_ << " " << t_ << endl;
    }
    else
    {
       // cout << "Qg->Qg will start soon!" << endl;
	int status_s = Sample_E2ANDs_Qg(A, E1_, temp_);
	int stauts_t = Sample_t_Qg(A, s_, temp_);
	cout << 2 << " " <<ctl_sample << " " << E2_ << " " << s_ << " " << t_ << endl;
    }
}





double* Sample_2to2::Sample_vector_CoM(Scattering_2to2& A, double s, double t)
{
    double p40 = (s - A.Mass2_) / (2.*sqrt(s));
    double cos_theta4 = - t/(2*p40*p40) - 1;
    double sin_theta4 = sqrt(1 - cos_theta4*cos_theta4);
    double phi4 = 2.*M_PI*((double) rand()/RAND_MAX);
    double p4x = p40 * sin_theta4 * cos(phi4);
    double p4y = p40 * sin_theta4 * sin(phi4);
    double p4z = p40 * cos_theta4;

    double p30_ = sqrt(p40*p40 + A.Mass2_);
    double p3x_ = -p4x;
    double p3y_ = -p4y;
    double p3z_ = -p4z;

    double p3[] = {p30_, p3x_, p3y_, p3z_};
    return p3;
}

    


// right now think about that, wonder if I can leave Alpha as a out_of_class functiom
double Sample_2to2::Alpha(double Q2)
{
    double alpha_s;
    double Lambda2 = 0.04;

    if (Q2 < 0)
        alpha_s = 1.396/(log(-Q2/Lambda2));
    else if (Q2 > 0)
        alpha_s = 1.396*(0.5 - atan(log(Q2/Lambda2)/M_PI)/M_PI);
    else
        alpha_s = 1.;

    if (alpha_s > 1 || alpha_s < 0)
        alpha_s = 1.;

    return alpha_s;
}


double Sample_2to2::DebyeMass2(double alpha_s, double temp)
{
    return 15.2789*alpha_s*temp*temp;
}













//--------------- for sample Qq -> Qq ---------------------------------





// for F(E1, temp) = f_j(E2, temp) * (s-M*M) * sigma(s, temp)
// this function returns MAX(F)(E1, temp). and which in fact, is when s = smax for each E2
// and it makes our life much easier, because in that way, for each (E1, temp), we can just loop over all E2* ~(0, 10*temp), then the maximum is when s = smax = beta*E2
// and in fact, what makes our life even better, is that, instead of loop E2 from 0, we can start from E2=temp, and it is based on our knowledge of the shape of F
double Sample_2to2::Max_E2ANDs_Qq(Scattering_2to2& A, double E1, double temp)
{
    double alpha = 2.*(E1 - sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));
    double dE2 = 0.1;
    double E2 = temp;

//    double Fmax = 1.0/(exp(E2/temp) + 1) * beta*E2 * A.Sigma_Qq(A.Mass2_+beta*E2, temp);
// to speed up
    double Fmax = 1.0/(exp(E2/temp) + 1) * beta*E2 * A.Get_sigma_Qq(A.Mass2_+beta*E2, temp);

    double Fmax_next = Fmax;

    while (Fmax_next >= Fmax)
    {
        Fmax = Fmax_next;
	E2 += dE2 * temp;
//	Fmax_next = 1.0/(exp(E2/temp) + 1) * beta*E2 * A.Sigma_Qq(A.Mass2_ + beta*E2, temp);
        Fmax_next = 1.0/(exp(E2/temp) + 1) * beta*E2 * A.Get_sigma_Qq(A.Mass2_+beta*E2, temp);
    }

/*
// check whether our calculation of Fmax is correct?
    for (int iE2 = 0; iE2 < 100; iE2 ++)      
    {
        E2 = iE2 * 0.01 * temp + 1e-6;
	double result = 10/(exp(E2/temp) + 1) * beta*E2 * A.Sigma_Qq(A.Mass2_+beta*E2, temp); 
	cout << E2 << " " << result << endl;
    }
*/

    return 1.2*Fmax;
}





int Sample_2to2::Sample_E2ANDs_Qq(Scattering_2to2& A, double E1, double temp)
{
    double alpha = 2.*(E1 - sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));
    
    double Fmax = Max_E2ANDs_Qq(A, E1, temp);


// to sample E2 and s, where we have 0 < E2 < 10*temp; smin < s < smax
// we can do the transformation: E2/temp = x+y; s/M*M = temp/M*M*(alpha*x + beta*y)
// then 0 < x+y < inf (or in this case: 10)
// x > = -0.25 /(temp/(M*M) * sqrt(E1*E1 -M*M))
// y >= 0.25/(temp/(M*M) * sqrt(E1*E1-M*M))
// and this is easy to achieve by constructing: x ~ random.uniform(0,1), y~random.uniform(0,1), 
//then x = x**0.5, y = y*x, x = 1-x, x = xmin + (xmax-xmin)*x, y = ymin + (ymax-ymin)*y, then the final (x, y) is in the triangle when we desired

    
    int ctl = 0;
    int count = 0;
    double x, y, E2, s;
    double xmin = -0.25/(temp/(A.Mass2_) * sqrt(E1*E1 - A.Mass2_));
    double ymin = -xmin;
    
    while (!ctl)
    {
        x = ((double) rand()/RAND_MAX);
	y = ((double) rand()/RAND_MAX);
	x = sqrt(x);
	y *= x;
	x = 1-x;
	x = xmin + 10*x; // here it is not that general actually, because in this case, I set
			 // E2 ~ (0, 10*temp), which in fact, should be (0, inf)
	y = ymin + 10*y;

	s = temp*(alpha*x + beta*y);
	E2 = temp*(x+y);

//	double result = 1.0/(exp(E2/temp) + 1) * (s-A.Mass2_) * A.Sigma_Qq(s, temp);
        double result = 1.0/(exp(E2/temp) + 1) *(s-A.Mass2_) * A.Get_sigma_Qq(s, temp);
// for check purpose
        if (result > Fmax)  cout << "something is WRONG in Max_E2ANDs_Qq()!!" << endl;
        double ctl_sample = ((double) rand()/RAND_MAX);
        ctl = ctl_sample < result/Fmax ? 1 : 0;
	count += 1;
    }

    // cout << count << " " << ctl << " " << x << " " << y << " " << E2 << " " << s << endl;

    s_ = s;
    E2_ = E2;

    return ctl;
}



// okay, right now I am going to sample t. in fact, if you plot M2(t), you will find that it distributed as an exponential funciton.
// therefore, for sample M2(t) in Qq->Qq scattering, we are using the combination of importance and rejection sampling...
// here F(t, s) = M2(t,s), MAX(F) = F(t=tmax=0)
// and if we sample p(t) = a exp(-xa), (here a = 100 would be good enough), then F_prime(t,s) = M2(t, s) * exp(xa)  -- note: I tested in python, where the exponential contruction in numpy is different from c++, in python, it is 1/b exp(-x/b), and b = 1/a
// checked, when combining the importance sampling and rejection sampling, the efficiency is 49085/100000

int Sample_2to2::Sample_t_Qq(Scattering_2to2& A, double s, double temp)
{
//    int status_E2ANDs = Sample_E2ANDs_Qq(A, E1, temp);
// in constructor, if I already add this step, I don't think I need to do this one anymore. Better think of a way where should I leave this
// Or should I say, Sample_t_Qq(double s, double temp)? After all, it is related to s instead of E1
// should add a control condition to make sure there is no error here....
    double tmin = -pow(s-A.Mass2_, 2)/s;
    double tmax = 0;

    double scale = 10.;  // in fact, here scale>10. should be sufficient, after all, we dont want our sample to be too dramaric

    double Fmax = A.M2_Qq(tmax, s, temp) * exp(scale*tmax);
    Fmax = 1.1*Fmax;
// check if my assumption of the maximum value is correct?
// note: checked M2_Qg, and it turns out the maximum value is more complex than M2_Qq, For M2_Qq, as t increase, M2_Qq increase. 
// but for M2_Qg, it is s dependent. And the safe way is to calculate M2(tmin) and M2(tmax) and compare the result
/*
    double dt = (tmax - tmin) / 100;
    for (int i=0; i< 100; i++)
    {
        double t = tmin + i * dt;
        double result = A.M2_Qq(t, s, temp) * exp(scale*t);
        cout << t << " " << result << " "<< Fmax << endl;
    }
*/


    std::default_random_engine generator(std::random_device{}());
    std::exponential_distribution<double> distribution(scale);

    double t;
    int ctl = 0;
    int count = 0;
    while (!ctl)
    {
        t = -distribution(generator);
	while (!(t>tmin && t<tmax)) t = -distribution(generator);

	double result = A.M2_Qq(t, s, temp) * exp(scale*t);
	double ctl_sample = ((double) rand()/RAND_MAX);
	ctl = (ctl_sample < result/Fmax ? 1 : 0);
	count += 1;
    }

    t_ = t;
    //cout << count << " " << ctl << " " << s_ < " " << t_ << endl;



    return ctl;
}

























// ---------------- for sample Qg -> Qg ---------------------------
double Sample_2to2::Max_E2ANDs_Qg(Scattering_2to2& A, double E1, double temp)
{
    double alpha= 2.*(E1-sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1+sqrt(E1*E1 - A.Mass2_));
    double dE2 = 0.1;
    double E2 = 1e-6;

//    double Fmax = 1.0/(exp(E2/temp) + 1) * beta*E2 * A.Sigma_Qg(A.Mass2_+beta*E2, temp);
    double Fmax = 1.0/(exp(E2/temp)+1) * beta*E2 *A.Get_sigma_Qg(A.Mass2_+beta*E2, temp);
    double Fmax_next = Fmax;

    while (Fmax_next >= Fmax)
    {
        Fmax = Fmax_next;
	E2 += dE2 * temp;
	//Fmax_next = 1.0/(exp(E2/temp) - 1) *beta*E2 * A.Sigma_Qg(A.Mass2_+beta*E2, temp);
        Fmax_next = 1.0/(exp(E2/temp) - 1) *beta*E2 * A.Get_sigma_Qg(A.Mass2_+beta*E2, temp);
    }

// check whether the calculation of Fmax is correct??
/*
    for (int iE2 = 0; iE2< 100; iE2 ++)
    {
        E2 = iE2*dE2*temp + 1e-6;
	double result = 1.0/(exp(E2/temp) - 1) * beta*E2 * A.Sigma_Qg(A.Mass2_+beta*E2, temp);
	cout << E2 <<  " " << result << " " << Fmax << endl;
    }
*/
    return 1.2*Fmax;
}




int Sample_2to2::Sample_E2ANDs_Qg(Scattering_2to2& A, double E1, double temp)
{
    double alpha = 2.*(E1 - sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));
    double Fmax = Max_E2ANDs_Qg(A, E1, temp);

    int ctl = 0;
    int count = 0;
    double x, y, E2, s;
    double xmin = -0.25/(temp/(A.Mass2_) * sqrt(E1*E1 - A.Mass2_));
    double ymin = -xmin;

    while (!ctl)
    {
        x = ((double) rand()/RAND_MAX);
	y = ((double) rand()/RAND_MAX);
	x = sqrt(x);
	y *= x;
	x = 1-x;
	x = xmin + 10*x;
	y = ymin + 10*y;

	s = temp*(alpha*x+beta*y);
	E2 = temp*(x+y);
//	double result = 1.0/(exp(E2/temp) -1) * (s-A.Mass2_) * A.Sigma_Qg(s, temp);
        double result = 1.0/(exp(E2/temp) -1) * (s-A.Mass2_) * A.Get_sigma_Qg(s, temp);
// for check purpose
        if (result > Fmax) cout << "something is WRONG in Max_E2ANDs_Qg()!!!"<< endl; 
	double ctl_sample = ((double) rand()/RAND_MAX);

	ctl = ctl_sample < result/Fmax ? 1 : 0;
	count += 1;
    }

    s_ = s;
    E2_ = E2;

    return ctl;
}



int Sample_2to2::Sample_t_Qg(Scattering_2to2& A, double s, double temp)
{
    double tmin = -pow(s-A.Mass2_, 2)/s;
    double tmax = 0;

    double scale = 10.;
    double Fmax1, Fmax2, Fmax;
    Fmax1 = A.M2_Qg(tmax, s, temp) * exp(scale*tmax);
    Fmax2 = A.M2_Qg(tmin, s, temp) * exp(scale*tmin);
    Fmax = 1.1 * std::max(Fmax1, Fmax2);
// check check
/*
    double dt = (tmax - tmin)/100;
    for (int i=0; i<100; i++)
    {
        double t = tmin + i*dt;
	double result = A.M2_Qg(t, s, temp) * exp(scale*t);
        cout << t << " " << result << " " << Fmax << endl;
    }
*/
    std::default_random_engine generator(std::random_device{}());
    std::exponential_distribution<double> distribution(scale);

    double t;
    int ctl = 0;
    int count = 0;
    while (!ctl)
    {
        t = -distribution(generator);
	while (!(t>tmin && t<tmax)) t = -distribution(generator);

	double result = A.M2_Qg(t, s, temp) * exp(scale*t);
	double ctl_sample = ((double) rand()/RAND_MAX);
	ctl = (ctl_sample < result/Fmax ? 1 : 0);
	count += 1;
    }

    t_ = t;
    // cout << count << " " << ctl << " " << s << " " << t_ << endl;
    return ctl;
}
