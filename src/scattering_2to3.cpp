#include "scattering.h"

using namespace std;

Scattering_2to3::Scattering_2to3(double M, double kfactor)
:    T0_(0.1),
     E10_(1.691),
     ds_(0.2),
     dT_(0.02),
     dE1_(2.),
     Ns_MAX_(41),
     NT_MAX_(51),
     NE_MAX_(41),
     sigma_Qq_(boost::extents[Ns_MAX_][NT_MAX_]),
     sigma_Qg_(boost::extents[Ns_MAX_][NT_MAX_]),
     gamma_Qq_(boost::extents[NE_MAX_][NT_MAX_]),
     gamma_Qg_(boost::extents[NE_MAX_][NT_MAX_])
{
    Mass2_ = M*M;
    kfactor_ = kfactor;
    s0_ = Mass2_ + 1e-3;
}


double Scattering_2to3::Alpha(double Q2)
{
    double alpha_s;
    if (Q2 < 0)
        alpha_s = 1.392/(log(-Q2/Lambda2));
    else if (Q2 > 0)
        alpha_s = 1.396*(0.5 - atan(log(Q2/Lambda2)/M_PI)/M_PI);
    else
        alpha_s = 1;

    if (alpha_s > 1 || alpha_s < 0)
        alpha_s = 1;

    return alpha_s;

}


double Scattering_2to3::DebyeMass2(double alpha_s, double temp)
{
    return 15.2789*alpha_s*temp*temp;
//    return 0.73;
}


// k = {p4, k, theta4, phi_4k}
int Scattering_2to3::ValidParams(double *k, double s)
{
    int status1 = (k[0] + k[1]) < sqrt(s) ? 1 : 0;
    int status2 = (k[0] + k[1]) > (s-Mass2_)/(2*sqrt(s)) ? 1 : 0;
    int status3 = (2*sqrt(s)*(k[0]+k[1]) - 4*k[0]*k[1]) < (s-Mass2_) ? 1 : 0;
    return status1*status2*status3;
}


double* Rotate(double theta, double phi, double *p)
{
    double ct = cos(theta), st = sin(theta), cp = cos(phi), sp = sin(phi);
    double R00 = ct*cp;
    double R01 = -sp;
    double R02 = st*cp;
    double R10 = ct*sp;
    double R11 = cp;
    double R12 = st*sp;
    double R20 = -st;
    double R21 = 0.;
    double R22 = ct;

    double p_prime[] = {0,0,0};
    p_prime[0] = R00*p[0] + R01*p[1] + R02*p[2];
    p_prime[1] = R10*p[0] + R11*p[1] + R12*p[2];
    p_prime[2] = R20*p[0] + R21*p[1] + R22*p[2];

    return p_prime;
}



double Scattering_2to3::M2_Qq(double *k, double s, double temp)
{
    double p40 = k[0];
    double k0 = k[1];
    double theta4 = k[2];
    double phi4 = 0.0; // phi4 can be anything between [0, 2pi)
    double phi4k = k[3];
    double vec_p4[] = {p40*sin(theta4), 0.0, p40*cos(theta4)};
    double cos_theta_k4 = (s-Mass2_ - 2*sqrt(s)*(p40+k0))/(2*p40*k0) + 1.0;
    if (cos_theta_k4 < -1 || cos_theta_k4 > 1)
        cout << "why?? Error in M2_Qq!! " << endl;
    else
    {
        double sin_theta_k4 = sqrt(1-cos_theta_k4*cos_theta_k4);
	double vec_k_prime[] = {k0*sin_theta_k4*cos(phi4k), k0*sin_theta_k4*sin(phi4k), k0*cos_theta_k4};
	double *vec_k;
	vec_k = Rotate(theta4, 0.0, vec_k_prime);
	double kx = vec_k[0];
	double ky = vec_k[1];
	double kz = vec_k[2];
	double qx = -vec_p4[0];
	double qy = -vec_p4[1];

	double y = 0.5*log((k0+kz)/(k0-kz));

        
	double kperp2 = kx*kx + ky*ky;
	double qperp2 = qx*qx + qy*qy;

        double phi_kq = (kx*qx + ky*qy)/(sqrt(kperp2) * sqrt(qperp2));
        

	double x_bar = sqrt(kperp2) * exp(fabs(y)) / sqrt(s);
	double x = sqrt(kperp2) * exp(y) / sqrt(s);

	double df_qk[] = {qx-kx, qy-ky};
	double qtkt = qperp2 + kperp2 - 2*(kx*qx + ky*qy);
	double t_prime = -qperp2/(s-Mass2_)*(qtkt/(1-x) + kperp2/x + x/(1-x)*Mass2_) - qperp2;
	double u = (1-qperp2/(s-Mass2_))*(2*Mass2_-s+qtkt/(1-x) + kperp2/x + x/(1-x)*Mass2_) - qperp2;
	double alpha_t = Alpha(t_prime);
	double m_D2 = DebyeMass2(alpha_t, temp);
	double a = kperp2 + x*x*Mass2_ + m_D2;
	double b = pow(df_qk[0], 2) + pow(df_qk[1],2) + x*x*Mass2_ + m_D2;
	double PD = pow(kx/a + df_qk[0]/b,2) + pow(ky/a + df_qk[1]/b,2);

	double M2to2 = 64./9.*M_PI*M_PI*alpha_t*alpha_t*(pow(Mass2_-u,2) + pow(s-Mass2_, 2) + 2*Mass2_*t_prime) / pow(t_prime - kfactor_*m_D2, 2);
	double alpha_s = 0.3;
	double M2_ = 12.*4.*M_PI*alpha_s*M2to2*pow(1-x_bar, 2)*PD;

        //cout << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << " " << qperp2 << " " << kperp2 << " " << y << " " << phi_kq << " " << M2_ << endl;
//	cout << qperp2 << " " << kperp2 << " " << y << " "<< phi_kq << endl;
//	cout << t_prime << " " << u << endl;
//	cout << qtkt << " " << x << endl;
//	cout << M2to2 << " " << M2_ << endl;

	double result =M2_ * sin(theta4);
	return result;
    }
}






double Wrapper_M2_Qq(double *k, size_t dim, void *params)
{
    struct gsl_params2to3 *Wp = (struct gsl_params2to3*) params;
    double s = Wp->s_;
    double temp = Wp->temp_;
    double Mass2 = Wp->pt_class_->Mass2_;

    int status_params = Wp->pt_class_->ValidParams(k, s);
    if (status_params == 0)
        return 0;
    else
    {
        double M2 = Wp->pt_class_->M2_Qq(k, s, temp);
        return M2;	
    }
}






double Scattering_2to3::Sigma_Qq(double s, double temp)
{
    double result, error;
    double xl[] = {0,0,0,0};
    double xu[] = {(s-Mass2_)/(2*sqrt(s)), (s-Mass2_)/(2*sqrt(s)), M_PI, 2*M_PI};
// if I want to use the same set up as BAMPS, which ask p4z > 0 ==> theta4 ~ (0, M_PI/2)
//    double xu[] = {(s-Mass2_)/(2*sqrt(s)), (s-Mass2_)/(2*sqrt(s)), 0.5*M_PI, 2.*M_PI};

 
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_monte_function G;

    struct gsl_params2to3 Wp = {};
    Wp.pt_class_ = this;
    Wp.temp_ = temp;
    Wp.s_ = s;

    G.f = &Wrapper_M2_Qq;
    G.dim = 4;
    G.params = &Wp;

    size_t  calls = 500000;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    {
       gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(4);
       gsl_monte_vegas_integrate (&G, xl, xu, 4, 100000, r, s, &result, &error);
       int ctl = 0;
       do{
           gsl_monte_vegas_integrate(&G, xl, xu, 4, calls/5, r, s, &result, &error);
	   ctl += 1;
       } while(fabs(gsl_monte_vegas_chisq(s) -1.0) > 0.5 && ctl < 100);

       gsl_monte_vegas_free(s);
    }

    gsl_rng_free(r);

    return result/(256.*pow(M_PI, 4) * (s-Mass2_));
}



//for most of the case, we are calculation s in a smaller range (do you want to make a histgram 
//about the distribution of s??) and in that way, it really a waste of 
//time and lose of precision to calculte such kind of evenly sliced grid.
//RIGHT, just did a histgram,for E1 = 10, temp=0.12, most of the s are in the low range (0,20)
void Scattering_2to3::Sigma_Qq_tabulate()
{
    double s, temp;
    int is, iT;
    for (is = 0; is < Ns_MAX_; ++is)
    {
        s = exp(log(s0_) + is*ds_);
	for (iT = 0; iT<NT_MAX_; iT++)
	{
	    temp = T0_ + iT*dT_;
	    sigma_Qq_[is][iT] = Sigma_Qq(s, temp);
//            cout << is << "   " << s << "    "  << temp <<  "   " << sigma_Qq_[is][iT] << endl;
         }
    
    }
    cout << "Sigma_Qq_tabulate() success! " << endl;
}



double cubeInterp(double x, double y, double A00, double A10, double A01, double A11)
{
    return A00*(1-x)*(1-y) + A10*x*(1-y) + A01*(1-x)*y + A11*x*y;
}


double Scattering_2to3::Get_sigma_Qq(double s, double temp)
{
    int is, iT;
    double delta_s, delta_T;
    double var_s, var_T;

    var_s = (log(s) - log(s0_)) / ds_;
    var_T = (temp - T0_)/dT_;
    
    is = int(var_s);
    iT = int(var_T);
    delta_T = var_T - iT;
    delta_s = var_s - is;

    if(is < 0 || is >= Ns_MAX_-1 || iT<0 || iT >= NT_MAX_-1)
    {
//        cout << "energy/temperaure out of range!\n";
//	cout << s << " " << endl;
	return sigma_Qq_[min(is, Ns_MAX_-1)][min(iT, NT_MAX_-1)];
    }
    else
    {
        if ((1-delta_s) < 1e-8) {is += 1; delta_s = 0.;}
	if ((1-delta_T) < 1e-8) {iT += 1; delta_T = 0.;}
	//double T1 = pow(T0_ + dT_ *iT, 2);
	//double T2 = pow(T1 + dT_, 2);
	double result = cubeInterp(delta_s, delta_T, sigma_Qq_[is][iT], sigma_Qq_[is+1][iT], sigma_Qq_[is][iT+1], sigma_Qq_[is+1][iT+1]);
	return result;
    }
}



double Wrapper_get_sigma_Qq(double s, void *params)
{
    struct gsl_params2to3 *Wp = (struct gsl_params2to3*) params;
    double temp = Wp->temp_;
 //   cout << s << "   " << temp << endl;
    double result = Wp->pt_class_->Get_sigma_Qq(s, temp);
    double Mass2 = Wp->pt_class_->Mass2_;
    //cout << s << " " << result << endl;
    return result *(s - Mass2);
}




double Scattering_2to3::Omega_Qq(double E1, double E2, double temp)
{
    struct gsl_params2to3 Wp={};
    Wp.pt_class_ = this;
    Wp.temp_ = temp;

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;

    gsl_function O;
    O.function = &Wrapper_get_sigma_Qq;
    O.params = &Wp;

    double smin = Mass2_ + 2.*E1*E2 - 2*E2*sqrt(E1*E1 - Mass2_);
    double smax = Mass2_ + 2.*E1*E2 + 2.*E2*sqrt(E1*E1 - Mass2_);

    double relerr = 1e-5;
    int status = 1;
    int nloop = 0;
    while(status && nloop < 10)
    {
        status = gsl_integration_qags(&O, smin, smax, 0, relerr, 1000, w, &result, &error);
	relerr *= 1.5;
	nloop += 1;
    }

    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    return result;

}


double Wrapper_omega_Qq(double E2, void *params)
{
    struct gsl_params2to3 *Wp = (struct gsl_params2to3 *) params;
    double E1 = Wp->E1_;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->Omega_Qq(E1, E2, temp);
    return result/(exp(E2/temp) + 1);
}




double Scattering_2to3::Gamma_Qq(double E1, double temp)
{
    struct gsl_params2to3 Wp = {};
    Wp.pt_class_ = this;
    Wp.E1_ = E1;
    Wp.temp_ = temp;

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function G;
    G.function = &Wrapper_omega_Qq;
    G.params = &Wp;

    double E2min = 0.;
    double E2max = 10.*temp;
    double relerr = 1e-5;

    int status =1 ;
    int nloop = 0;
    while (status && nloop < 10)
    {
        status = gsl_integration_qags(&G, E2min, E2max, 0, relerr, 1000, w, &result, &error);
	relerr *= 1.5;
	nloop += 1;
    }

    gsl_integration_workspace_free(w);
    gsl_set_error_handler(old_handler);

    double scattering_rate = result/(16.*M_PI*M_PI*E1*sqrt(E1*E1-Mass2_));
    return 36*scattering_rate;
}



void Scattering_2to3::Gamma_Qq_tabulate()
{
    double E1, temp;
    int iE, iT;
    for (iE = 0; iE < NE_MAX_; ++iE)
    {
        E1 = E10_ + iE*dE1_;
	for (iT=0; iT<NT_MAX_; ++iT)
	{
	    temp = T0_ + iT*dT_;
	    gamma_Qq_[iE][iT] = Gamma_Qq(E1, temp);
	}
    }
    cout << "Gamma_Qq_tabulate success! " << endl;
}



double Scattering_2to2::Get_gamma_Qq(double E1, double temp)
{
    int iE, iT;
    double delta_E, delta_T;
    double var_E, var_T;

    var_E = (E1 - E10_)/dE1_;
    var_T = (temp-T0_)/dT_;

    iE = int(var_E);
    iT = int(var_T);

    delta_E = var_E - iE;
    delta_T = var_T - iT;

    if (iE < 0 || iE>=NE_MAX_-1 || iT < 0 || iT >= NT_MAX_ -1)
    {
        return gamma_Qq_[min(iE, NE_MAX_-1)][min(iT, NT_MAX_-1)];
    }
    else
    {
       if ((1-delta_E) < 1e-8) {iE += 1; delta_E = 0.;}
       if ((1-delta_T) < 1e-8) {iT += 1; delta_T = 0.;}
       double result = cubeInterp(delta_E, delta_T, gamma_Qq_[iE][iT], gamma_Qq_[iE+1][iT], gamma_Qq_[iE][iT+1], gamma_Qq_[iE+1][iT+1]);

       return result;
    }
}












//------------------------------ Qg-> Qgg -----------------------------------

double Scattering_2to3::M2_Qg(double *k, double s, double temp)
{
    double p40 = k[0];
    double k0 = k[1];
    double theta4 = k[2];
    double phi4 = 0.0; // phi4 can be anything between [0, 2pi)
    double phi4k = k[3];
    double vec_p4[] = {p40*sin(theta4), 0.0, p40*cos(theta4)};
    double cos_theta_k4 = (s-Mass2_ - 2*sqrt(s)*(p40+k0))/(2*p40*k0) + 1.0;
    if (cos_theta_k4 < -1 || cos_theta_k4 > 1)
        cout << "why?? Error in M2_Qq!! " << endl;
    else
    {
        double sin_theta_k4 = sqrt(1-cos_theta_k4*cos_theta_k4);
	double vec_k_prime[] = {k0*sin_theta_k4*cos(phi4k), k0*sin_theta_k4*sin(phi4k), k0*cos_theta_k4};
	double *vec_k;
	vec_k = Rotate(theta4, 0.0, vec_k_prime);
	double kx = vec_k[0];
	double ky = vec_k[1];
	double kz = vec_k[2];
	double qx = -vec_p4[0];
	double qy = -vec_p4[1];

	double y = 0.5*log((k0+kz)/(k0-kz));
	double kperp2 = kx*kx + ky*ky;
	double qperp2 = qx*qx + qy*qy;

	double x_bar = sqrt(kperp2) * exp(fabs(y)) / sqrt(s);
	double x = sqrt(kperp2) * exp(y) / sqrt(s);

	double df_qk[] = {qx-kx, qy-ky};
	double qtkt = qperp2 + kperp2 - 2*(kx*qx + ky*qy);
	double t_prime = -qperp2/(s-Mass2_)*(qtkt/(1-x) + kperp2/x + x/(1-x)*Mass2_) - qperp2;
	double u = (1-qperp2/(s-Mass2_))*(2*Mass2_-s+qtkt/(1-x) + kperp2/x + x/(1-x)*Mass2_) - qperp2;
	double alpha_t = Alpha(t_prime);
	double m_D2t = DebyeMass2(alpha_t, temp);
	double alpha_s = Alpha(s-Mass2_);
	double m_D2s = DebyeMass2(alpha_s, temp);
	double alpha_u = Alpha(u-Mass2_);
	double m_D2u = DebyeMass2(alpha_u, temp);


	double a = kperp2 + x*x*Mass2_ + m_D2t;
	double b = pow(df_qk[0], 2) + pow(df_qk[1],2) + x*x*Mass2_ + m_D2t;
	double PD = pow(kx/a + df_qk[0]/b,2) + pow(ky/a + df_qk[1]/b,2);

        double X1 = pow(alpha_t, 2) *2.*(s-Mass2_)*(Mass2_-u)/pow(t_prime - kfactor_*m_D2t,2);
	double X2 = pow(alpha_s, 2) *4./9. *((s-Mass2_)*(Mass2_-u) + 2.*Mass2_*(s+Mass2_))/pow(s-Mass2_+m_D2s, 2);
	double X3 = pow(alpha_u, 2) *4./9. *((s-Mass2_)*(Mass2_-u) + 2.*Mass2_*(u+Mass2_))/pow(Mass2_-u+m_D2u, 2);
	double X4 = alpha_s*alpha_u*1./9.*Mass2_*(4*Mass2_-t_prime)/((s-Mass2_+m_D2s)*(Mass2_-u+m_D2u));
	double X5 = alpha_t*alpha_s *((s-Mass2_)*(Mass2_-u) + Mass2_*(s-u))/((t_prime-kfactor_*m_D2t)*(s-Mass2_+m_D2s));
	double X6 = alpha_t*alpha_s*((s-Mass2_)*(Mass2_-u) - Mass2_*(s-u))/((t_prime-kfactor_*m_D2t)*(Mass2_-u+m_D2u));
	double M2to2 = pow(4.*M_PI,2) *(X1+X2+X3+X4+X5-X6);

	double alpha_k = 0.3;
	double M2_ = 12.*4.*M_PI*alpha_k*M2to2*pow(1-x_bar, 2)*PD;


	double result =M2_ * sin(theta4);
	return result;
    }
}

double Wrapper_M2_Qg(double *k, size_t dim, void *params)
{
    struct gsl_params2to3 *Wp = (struct gsl_params2to3 *) params;
    double s = Wp->s_;
    double temp = Wp->temp_;
    double Mass2 = Wp->pt_class_->Mass2_;

    int status_params = Wp->pt_class_->ValidParams(k, s);
    if (status_params == 0)
        return 0;
    else
    {
        double M2 = Wp->pt_class_->M2_Qg(k, s, temp);
        return M2;	
    }
}


double Scattering_2to3::Sigma_Qg(double s, double temp)
{
    double result, error;
    double xl[] = {0,0,0,0};
    double xu[] = {(s-Mass2_)/(2*sqrt(s)), (s-Mass2_)/(2*sqrt(s)), M_PI, 2*M_PI};
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_monte_function G;

    struct gsl_params2to3 Wp = {};
    Wp.pt_class_ = this;
    Wp.temp_ = temp;
    Wp.s_ = s;

    G.f = &Wrapper_M2_Qg;
    G.dim = 4;
    G.params = &Wp;

    size_t  calls = 500000;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    {
       gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(4);
       gsl_monte_vegas_integrate (&G, xl, xu, 4, 100000, r, s, &result, &error);
       int ctl = 0;
       do{
           gsl_monte_vegas_integrate(&G, xl, xu, 4, calls/5, r, s, &result, &error);
	   ctl += 1;
       } while(fabs(gsl_monte_vegas_chisq(s) -1.0) > 0.5 && ctl < 100);

       gsl_monte_vegas_free(s);
    }

    gsl_rng_free(r);

    return result/(256.*pow(M_PI, 4) * (s-Mass2_));
}







void Scattering_2to3::Sigma_Qg_tabulate()
{
    double s, temp;
    int is, iT;
    for (is = 0; is < Ns_MAX_; ++is)
    {
        s = exp(log(s0_) + is *ds_);
        for (iT =0; iT<NT_MAX_; ++iT)
	{
	    temp = T0_ + iT*dT_;
	    sigma_Qg_[is][iT] = Sigma_Qq(s, temp);
	}
    }
   cout << "Sigma_Qg_tabulate success! " << endl;
}



double Scattering_2to3::Get_sigma_Qg(double s, double temp)
{
    int is, iT;
    double delta_s, delta_T;
    double var_s, var_T;

    var_s = (log(s) - log(s0_)) / ds_;
    var_T = (temp - T0_)/dT_;
    is = int(var_s);
    iT = int(var_T);
    delta_s = var_s - is;
    delta_T = var_T - iT;

    if(is < 0 || is >= Ns_MAX_-1 || iT<0 || iT >= (NT_MAX_-1))
    {
	return sigma_Qg_[min(is, Ns_MAX_-1)][min(iT, NT_MAX_-1)];
    }
    else
    {
        if ((1-delta_s) < 1e-8) {is += 1; delta_s = 0.;}
	if ((1-delta_T) < 1e-8) {iT += 1; delta_T = 0.;}

	double result = cubeInterp(delta_s, delta_T, sigma_Qg_[is][iT], sigma_Qg_[is+1][iT], sigma_Qg_[is][iT+1], sigma_Qg_[is+1][iT+1]);
	return result;
    }
}




double Wrapper_get_sigma_Qg(double s, void *params)
{
    struct gsl_params2to3 *Wp = (struct gsl_params2to3 *) params;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->Get_sigma_Qg(s, temp);
    double Mass2 = Wp->pt_class_->Mass2_;
    return result * (s-Mass2);
}



double Scattering_2to3::Omega_Qg(double E1, double E2, double temp)
{
    struct gsl_params2to3 Wp = {};
    Wp.pt_class_ = this;
    Wp.temp_ = temp;

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;

    gsl_function O;
    O.function = &Wrapper_get_sigma_Qg;
    O.params = &Wp;

    double smin = Mass2_ + 2.*E1*E2 - 2*E2*sqrt(E1*E1 - Mass2_);
    double smax = Mass2_ + 2.*E1*E2 + 2*E2*sqrt(E1*E1 - Mass2_);

    double relerr = 1e-5;
    int status = 1;
    int nloop = 0;
    while(status && nloop < 10)
    {
        status = gsl_integration_qags(&O, smin, smax, 0, relerr, 1000, w, &result, &error);
	relerr *=2;
	nloop += 1;
    }

    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);
    return result;
}





double Wrapper_omega_Qg(double E2, void *params)
{
    struct gsl_params2to3 *Wp = (struct gsl_params2to3 *) params;
    double E1 = Wp->E1_;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->Omega_Qg(E1, E2, temp);
    return result/(exp(E2/temp) - 1);
}



double Scattering_2to3::Gamma_Qg(double E1, double temp)
{
    struct gsl_params2to3 Wp = {};
    Wp.pt_class_ = this;
    Wp.E1_ = E1;
    Wp.temp_ = temp;

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;

    gsl_function G;
    G.function = &Wrapper_omega_Qg;
    G.params = &Wp;

    double E2min = 0;
    double E2max = 10.*temp;
    double relerr= 1e-5;
    int status = 1;
    int nloop = 0;

    while(status && nloop < 10)
    {
        status = gsl_integration_qags(&G, E2min, E2max, 0, relerr, 1000, w, &result, &error);
	relerr *= 2;
	nloop += 1;
    }

    gsl_integration_workspace_free(w);
    gsl_set_error_handler(old_handler);

    double scattering_rate = result/(16.*M_PI*M_PI*E1*sqrt(E1*E1 - Mass2_));
    return 16*scattering_rate;
}






