#include "sample.h"

using namespace std;


double* Rotate(double, double, double*);


Sample_2to3::Sample_2to3(Scattering_2to3& A, double E1, double temp)
{
    srand(std::random_device{}());
    E1_ = E1;
    temp_ = temp;


/*
    double SR_Qq = A.Get_gamma_Qq(E1_, temp_);
    double SR_Qg = A.Get_gamma_Qg(E1_, temp_);

    double ctl_sample = ((double) rand()/RAND_MAX);
    
    if (ctl_sample < SR_Qq/(SR_Qq + SR_Qg))
    {
        int status_s = Sample_E2ANDs_Qq(A, E1_, temp_);
	int status_p = Sample_p3(A, s_, temp_);
	cout << 1 << " " << ctl_sample << " " << E2_ << " " << s_ << " " << endl;
    }
    else
    {
        int status_s = Sample_E2ANDs_Qg(A, E1_, temp_);
	int stauts_p = Sample_p3(A, s_, temp_);
	cout << 2 << " " << ctl_sample << " " << E2_ << " " << s_ << " " << endl;
    }
*/
}

/*
double* Rotate(double theta, double phi, double *p)
{
    double ct=cos(theta), st=sin(theta), cp=cos(phi), sp=sin(phi);
    double R00 = ct*cp;
    double R01 = -sp;
    double R02 = st*cp;
    double R10 = ct*sp;
    double R11 = cp;
    double R12 = st*sp;
    double R20 = -st;
    double R21 = 0.;
    double R22 = ct;

    double p_prime[] = {0, 0, 0};
    p_prime[0] = R00*p[0] + R01*p[1] + R02*p[2];
    p_prime[1] = R10*p[0] + R11*p[1] + R12*p[2];
    p_prime[2] = R20*p[0] + R21*p[1] + R22*p[2];

    return p_prime;
}
*/

double* Sample_2to3::Sample_vector_CoM(Scattering_2to3& A, double s, double *k)
{
   double p40 = k[0];
   double k0 = k[1];
   double theta4 = k[2];
   double phi4k = k[3];

   double phi4 = 2.*M_PI*((double) rand()/RAND_MAX);
   double vec_p4[] = {p40, p40*sin(theta4)*cos(phi4), p40*sin(theta4)*sin(phi4), p40*cos(theta4)};

   double cos_theta_4k = (s-A.Mass2_ - 2*sqrt(s)*(p40+k0))/ (2*p40*k0) + 1.0;
   double sin_theta_4k = sqrt(1-cos_theta_4k*cos_theta_4k);
   double vec_k_prime[] = {k0, k0*sin_theta_4k*cos(phi4k), k0*sin_theta_4k*sin(phi4k), k0*cos_theta_4k};
   double *vec_k;

   vec_k = Rotate(theta4, phi4, vec_k_prime);

   p30_ = sqrt(s) - vec_p4[0] - vec_k[0];
   p3x_ = - vec_p4[1] - vec_k[1];
   p3y_ = -vec_p4[2] - vec_k[2];
   p3z_ = -vec_p4[3] - vec_k[3];

   double vec_p3[] = {p30_, p3x_, p3y_, p3z_};
   return vec_p3;


}

double Sample_2to3::Alpha(double Q2)
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

double Sample_2to3::DebyeMass2(double alpha_s, double temp)
{
    return 15.2789*alpha_s*temp*temp;
}

















//------------------ sample Q + q -> Q + q + g ---------------------
// the good part is, to sample E2 and s, the 2->2 process should be the same as 2->3 process
double Sample_2to3::Max_E2ANDs_Qq(Scattering_2to3& A, double E1, double temp)
{
    double alpha = 2.*(E1 - sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));
    double dE2 = 0.1;
    double E2 = temp;

    double Fmax = 1.0/(exp(E2/temp) + 1) * beta*E2 *A.Get_sigma_Qq(A.Mass2_+beta*E2, temp);
    double Fmax_next = Fmax;

    while (Fmax_next >= Fmax)
    {
        Fmax = Fmax_next;
	E2 += dE2 * temp;
	Fmax_next = 1.0/(exp(E2/temp) + 1) * beta*E2 * A.Get_sigma_Qq(A.Mass2_+beta*E2, temp);
    }

    return 1.2*Fmax;
}



int Sample_2to3::Sample_E2ANDs_Qq(Scattering_2to3& A, double E1, double temp)
{
    double alpha = 2.*(E1 - sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));
    double Fmax = Max_E2ANDs_Qq(A, E1, temp);

// sample E2 and s, where we have 0< E2 < 10*temp; smin < s < smax
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

	s = temp*(alpha*x + beta*y);
	E2 = temp*(x+y);

	double result = 1.0/(exp(E2/temp) + 1.0) *(s-A.Mass2_) * A.Get_sigma_Qq(s, temp);
	if (result > Fmax) cout << "something is WRONG in Max_E2ANDs_Qq() !!" << endl;
	double ctl_sample = ((double) rand()/RAND_MAX);
	ctl = ctl_sample < result/Fmax ? 1 : 0;
	count += 1;
    }

    E2_ = E2;
    s_ = s;

    return ctl;
}




// okay, right now we have s, and now we need to sample p3
// since sigma = Integrate(p40, k0, theta4, phi_4k) * M2_2to3, and M2_2to3 = M2_2to2 * PD, PD < 1, therefore M2_2to3 < M2_2to2 <= M2_2to2(t=0)
// okay, right now the problem is the efficiency is too low

int Sample_2to3::Sample_p3_Qq(Scattering_2to3& A, double s, double temp)
{
    double tmax = 0;
    double Fmax = 48*M_PI*0.3*A.M2_Qq_2to2(tmax, s, temp);

    Fmax = 0.05*Fmax;

    double p4max = (s - A.Mass2_) / (2.*sqrt(s));
    int ctl = 0;
    int count = 0;
    double x, y, p40, k0, theta4, phi4, phi4k;
    
    double p4k, result;
    while (!ctl)
    {
        x = ((double) rand()/RAND_MAX);
	y = ((double) rand()/RAND_MAX);
	x = sqrt(x);
	y *=x;
	y = 1 -y;
	p40 = p4max * x;
	k0 = p4max * y;
        p4k = (p40-0.5*sqrt(s)) * (k0 - 0.5*sqrt(s));
	while (p4k < 0.25*A.Mass2_)
	{
	    x = ((double) rand()/RAND_MAX);
	    y = ((double) rand()/RAND_MAX);
	    x = sqrt(x);
	    y *= x;
	    y = 1 - y;
	    p40 = p4max *x;
	    k0 = p4max * y;
	    p4k = (p40 - 0.5*sqrt(s)) * (k0 - 0.5*sqrt(s));
	}
	
	theta4 = M_PI * ((double) rand()/RAND_MAX);
	phi4 = 2.*M_PI * ((double) rand()/RAND_MAX);
	phi4k = 2.*M_PI*((double) rand()/RAND_MAX);
	double Integrand[] = {p40, k0, theta4, phi4k};
        //cout << x << " " << y << " "<< p40 << " " << k0 << " " << endl;
	result = A.M2_Qq(Integrand, s, temp);
	double ctl_sample = ((double) rand()/RAND_MAX);
	if (result > Fmax) cout << "Same is Wrong in sample_p3_Qq!! " << endl;
	ctl = ctl_sample < result/Fmax ? 1 : 0;
	count += 1;
    }

    p40_ = p40;
    k0_ = k0;
    theta4_ = theta4;
    phi4_ = phi4;
    phi4k_ = phi4k;
    cout << ctl << " " << count << " " << p40_ << " "<< k0_ << " " << theta4_ << " " << phi4_ << " " << phi4k_ <<" " <<  result << " " <<Fmax << endl;
    return ctl;


}






















// ---------- sample Q + g -> Q + g + g ---------------------------
double Sample_2to3::Max_E2ANDs_Qg(Scattering_2to3& A, double E1, double temp)
{
    double alpha = 2.*(E1 - sqrt(E1*E1 - A.Mass2_));
    double beta = 2.*(E1 + sqrt(E1*E1 - A.Mass2_));
    double dE2 = 0.1;
    double E2 = 1e-6;

    double Fmax = 1.0/(exp(E2/temp) + 1) * beta*E2 * A.Get_sigma_Qg(A.Mass2_+beta*E2, temp);
    double Fmax_next = Fmax;

    while (Fmax_next >= Fmax)
    {
        Fmax = Fmax_next;
	E2 += dE2 * temp;
	Fmax_next = 1.0/(exp(E2/temp) - 1) *beta*E2 * A.Get_sigma_Qg(A.Mass2_+beta*E2, temp);
    }

    return 1.2*Fmax;
}



int Sample_2to3::Sample_E2ANDs_Qg(Scattering_2to3& A, double E1, double temp)
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
	x = 1 -x;
	x = xmin + 10*x;
	y = ymin + 10*y;

	s = temp*(alpha*x + beta*y);
	E1 = temp*(x+y);
	double result = 1.0/(exp(E2/temp) - 1) * (s-A.Mass2_) * A.Get_sigma_Qg(s, temp);

	if (result < Fmax) cout << "something is WRONG in Max_E2ANDs_Qq() !!!! " << endl;
	double ctl_sample = ((double) rand()/RAND_MAX);

	ctl = ctl_sample < result/Fmax ? 1:0;
	count += 1;

    }
     
    s_ = s;
    E2_ = E2;

    return ctl;
}





int Sample_2to3::Sample_p3_Qg(Scattering_2to3& A, double s, double temp)
{
    double tmax = 0;
    double Fmax = 48.*M_PI*0.3*A.M2_Qg_2to2(tmax, 2.*A.Mass2_ - s, s, temp);

    Fmax = 0.05*Fmax;

    double p4max = (s-A.Mass2_) / (2.*sqrt(s));
    int ctl = 0;
    int count = 0;

    double x, y, p40, k0, theta4, phi4, phi4k;

    double p4k, result;

    while(!ctl)
    {
        x = ((double) rand()/RAND_MAX);
	y = ((double) rand()/RAND_MAX);
	x = sqrt(x);
	y *= x;
	y = 1-y;
	p40 = p4max * x;
	k0 = p4max * y;
	p4k = (p40 - 0.5*sqrt(s)) * (k0 - 0.5*sqrt(s));
	while (p4k < 0.25 * A.Mass2_)
	{
	    x = ((double) rand()/RAND_MAX);
	    y = ((double) rand()/RAND_MAX);
	    x = sqrt(x);
	    y *= x;
	    y = 1 - y;
	    p40 = p4max * x;
	    k0 = p4max * y;
	    p4k = (p40 - 0.5*sqrt(s)) * (k0 - 0.5*sqrt(s));
	}

	theta4 = M_PI * ((double) rand()/RAND_MAX);
	phi4 = 2.*M_PI*((double) rand()/RAND_MAX);
	phi4k = 2.*M_PI*((double) rand()/RAND_MAX);

	double Integrand[] = {p40, k0, theta4, phi4};
	result = A.M2_Qg(Integrand, s, temp);
	double ctl_sample = ((double) rand()/RAND_MAX);
	if (result > Fmax) cout << "Something is WRONG in Sample_p3_Qg !!!" << endl;
        ctl = ctl_sample < result/Fmax ? 1 : 0;
	count += 1;
    }

    p40_ = p40;
    k0_ = k0;
    theta4_ = theta4;
    phi4_ = phi4;
    phi4k_ = phi4k;
    cout << ctl << " " << count << " " << p40_ << " " << k0_ << " " << theta4_ << " " << phi4_ << " " << phi4k_ << " " << result << " " << Fmax << " " << endl;
}
