#include "scattering.h"

using namespace std;


Scattering_2to2::Scattering_2to2(double M, double kfactor)
:   T0_(0.1),
    E10_(1.6901),
    ds_(0.1),
    dT_(0.02),
    dE1_(2.0), 
    Ns_MAX_(91),
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

double Scattering_2to2::Alpha(double Q2)
{
    double alpha_s;
    if(Q2 < 0)
        alpha_s = 1.396/(log(-Q2/Lambda2));
    else if (Q2 > 0)
        alpha_s = 1.396*(0.5 - atan(log(Q2/Lambda2)/M_PI)/M_PI);
    else
        alpha_s = 1;
        
    if (alpha_s > 1 || alpha_s < 0)
        alpha_s = 1;
     
    return alpha_s;
//     return 0.3;
}


double Scattering_2to2::DebyeMass2(double alpha_s, double temp)
{
    //return 8*alpha_s/M_PI*6*temp*temp;
    return 15.2789*alpha_s*temp*temp;
    // return 0.73;
}



// >>>>>>>>>>>>>>>>>>>>> Qq -> Qq >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
double Scattering_2to2::M2_Qq(double t, double s, double temp)
{
    double u = 2.*Mass2_ - s -t;
    double alpha_t = Alpha(t);
    double m_D2 = DebyeMass2(alpha_t, temp);
//    double result = 4.0/9.0*pow(4.0*M_PI*alpha_t, 2) * (pow(Mass2_ -u,2) + pow(s-Mass2_, 2) + 2*Mass2_ *t)/pow(t - kfactor_*m_D2, 2);
    double result = 70.184*alpha_t*alpha_t* (pow(Mass2_ -u, 2)+pow(s-Mass2_, 2) + 2*Mass2_*t)/pow(t-kfactor_*m_D2, 2);
    return result;
}


double Wrapper_M2_Qq(double t, void *params)
{
    struct gsl_params2to2 *Wp = (struct gsl_params2to2 *) params;
    double s = Wp->s_;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->M2_Qq(t, s, temp);
    return result;
}

double Scattering_2to2::Sigma_Qq(double s, double temp)
{
    struct gsl_params2to2 Wp = {};
    Wp.pt_class_ = this;
    Wp.s_ = s;
    Wp.temp_ = temp;

    double result, error;
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off(); // swtich off the old gsl error handler
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &Wrapper_M2_Qq;
    F.params = &Wp;

    double tmax = 0;
    double tmin = - pow(s-Mass2_,2)/s;
    double relerr = 1e-5;
    int status = 1;
    int nloop = 0;
    while (status && nloop < 10)
    {
        status = gsl_integration_qags(&F, tmin, tmax, 0, relerr, 1000, w, &result, &error);
        relerr *= 1.5;
       // if (status) std::cout << "Increase tolerence = " << relerr << std::endl;
    }

    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

   // double cross_section = result / (16.0*M_PI*pow(s-Mass2_,2));
    double cross_section = result/(50.2655*pow(s-Mass2_,2));
    return cross_section;
}



// Okay, here comes the question, how to tabulate sigma_Qq? the range of s is about (1.69 - 650), the range of temp is (0.1, 1.0)
// tabulate in log(s) and temp!
void Scattering_2to2::Sigma_Qq_tabulate()
{
    double s = s0_;
    double temp;
    int is, iT;
    for (is = 0; is < Ns_MAX_; ++is)
    {
        s = exp(log(s0_) + is*ds_);
        for (iT = 0; iT < NT_MAX_; ++iT)
        {
	    temp = T0_ + iT*dT_;
            sigma_Qq_[is][iT] = Sigma_Qq(s, temp);
//	    cout << s << " " << temp << " " << sigma_Qq_[is][iT] << endl;
        }
    }
    cout << "Sigma_Qq_tabulate success! " << endl;
   // if (is == Ns_MAX_) std::cout << "tabulate sigma finished :).\n";
   // else std::cout << "ERROR: Scattering_rate::Sigma_tabulate failed!\n";
}



double cubeInterp(double x, double y, double A00, double A10, double A01, double A11)
{
    return A00*(1-x)*(1-y) + A10*x*(1-y) + A01*(1-x)*y + A11*x*y;
}


// here use 2D interpolation to interpolate the sigma_Qq grid (I should use GSL_2D interpolation or any other advanced interpolation)
// need to think more carefully, what is the range it should be??
double Scattering_2to2::Get_sigma_Qq(double s, double temp)
{
    int is, iT;
    double delta_s, delta_T;
    double var_s, var_T;
   
    var_s = (log(s) - log(s0_)) / ds_;
    var_T = (temp-T0_)/dT_;

    is = int(var_s);
    iT = int(var_T);
    //cout << is << " " << iT << endl;
    delta_s = var_s - is;
    delta_T = var_T - iT;

    if (is <0 || is >= (Ns_MAX_-1) || iT < 0 || iT >= (NT_MAX_-1))
    {
//        cout << "energy/temperature out of range :(" << endl;
//	cout << s << " " << temp << endl;
	return sigma_Qq_[min(is, Ns_MAX_-1)][min(iT, NT_MAX_-1)];   // need to be more careful aboutthis part.....
    }
    else
    {
        if ((1-delta_s) < 1e-8)	{ is += 1; delta_s = 0.;}
	if ((1-delta_T) < 1e-8) { iT += 1; delta_T = 0.;}
	double T1 = pow(T0_ + dT_*iT,2);
	double T2 = pow(T1 + dT_, 2);
	double result = cubeInterp(delta_s, delta_T, sigma_Qq_[is][iT]*T1, sigma_Qq_[is+1][iT]*T1, sigma_Qq_[is][iT+1]*T2, sigma_Qq_[is+1][iT+1]*T2);
	return result/(temp*temp);
    }
}

double Wrapper_get_sigma_Qq(double s, void *params)
{
    struct gsl_params2to2 *Wp = (struct gsl_params2to2 *) params;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->Get_sigma_Qq(s, temp);
    double Mass2 = Wp->pt_class_->Mass2_;
    return result *(s-Mass2);
}


double Scattering_2to2::Omega_Qq(double E1, double E2, double temp)
{
    struct gsl_params2to2 Wp = {};
    Wp.pt_class_ = this;
    Wp.temp_ = temp;

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function O;
    O.function = &Wrapper_get_sigma_Qq;
    O.params = &Wp;

    double smin = Mass2_ + 2.*E1*E2 - 2.*E2*sqrt(E1*E1-Mass2_);
    double smax = Mass2_ + 2.*E1*E2 + 2.*E2*sqrt(E1*E1-Mass2_);
//    cout << "smin, smax: " << smin << "  " << smax << " " << E1 << " " << E2 << endl;

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
    struct gsl_params2to2 *Wp = (struct gsl_params2to2 *) params;
    double E1 = Wp->E1_;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->Omega_Qq(E1, E2, temp);
    return result /(exp(E2/temp)+1);
}



double Scattering_2to2::Gamma_Qq(double E1, double temp)
{
    struct gsl_params2to2 Wp= {};
    Wp.pt_class_ = this;
    Wp.temp_ = temp;
    Wp.E1_ = E1;

    gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function G;
    G.function = &Wrapper_omega_Qq;
    G.params = &Wp;
    
    double relerr = 1e-5;
    int status = 1;
    int nloop = 0;
    double E2_min = 0;
    double E2_max = 10.*temp;	
    while (status && nloop < 10)
    {
        status = gsl_integration_qags(&G, E2_min, E2_max, 0, relerr, 1000, w, &result, &error);
        relerr *= 1.5;
//	if (status) std::cout << "Increase tolerence = " << relerr << std::endl;
        nloop += 1;
    }

    gsl_integration_workspace_free(w);
    gsl_set_error_handler(old_handler);   
    // you don't need to add normalization here! (but if you want to calculate the average quantaty, you need to normalize it!!!!)
    //double scattering_rate = result / (16*M_PI*M_PI*E1*sqrt(E1*E1-Mass2_));
    double scattering_rate = result/(157.914*E1*sqrt(E1*E1-Mass2_));
    return scattering_rate*36;
}


void Scattering_2to2::Gamma_Qq_tabulate()
{
    double E1, temp;
    int iE, iT;
    for (iE = 0; iE < NE_MAX_; ++iE)
    {
        E1 = E10_ + iE*dE1_;
	for (iT= 0; iT < NT_MAX_; ++iT)
	{
	    temp = T0_ + iT*dT_;
	    gamma_Qq_[iE][iT] = Gamma_Qq(E1, temp);
	}
    }
    cout << "Gamma_Qq_tabulate success! " << endl; 
}




// here still need some more thoughts, you should not brutely force all the out of range result = maximum value. If you plot(gamma, E1), it is not like that
double Scattering_2to2::Get_gamma_Qq(double E1, double temp)
{
    int iE, iT;
    double delta_E, delta_T;
    double var_E, var_T;

    var_E = (E1 - E10_)/dE1_;
    var_T = (temp - T0_)/dT_;
    
    iE = int(var_E);
    iT = int(var_T);

    delta_E = var_E - iE;
    delta_T = var_T - iT;

    if (iE <0 || iE >=(NE_MAX_-1) || iT<0 || iT>=(NT_MAX_-1))
    {
        return gamma_Qq_[min(iE, NE_MAX_-1)][min(iT, NT_MAX_-1)];
    }
    else
    {
        if ((1-delta_E) < 1e-8) {iE += 1; delta_E =0.;}
	if ((1-delta_T) < 1e-8) {iT += 1; delta_T =0.;}
	double result = cubeInterp(delta_E, delta_T, gamma_Qq_[iE][iT], gamma_Qq_[iE+1][iT], gamma_Qq_[iE][iT+1], gamma_Qq_[iE+1][iT+1]);
	return result;
    }
}










// >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Qg -> Qg >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
double Scattering_2to2::M2_Qg(double t, double s, double temp)
{
    double u = 2.*Mass2_ - s - t;
    double alpha_t = Alpha(t);
    double alpha_s = Alpha(s-Mass2_);
    double alpha_u = Alpha(u-Mass2_);
    double m_D2t = DebyeMass2(alpha_t, temp);
    double m_D2s = DebyeMass2(alpha_s, temp);
    double m_D2u = DebyeMass2(alpha_u, temp);
    double X1 = pow(alpha_t, 2)* 2.*(s-Mass2_)*(Mass2_-u)/pow(t-kfactor_*m_D2t,2);
    double X2 = pow(alpha_s, 2) * 4./9.*((s-Mass2_)*(Mass2_-u) + 2.*Mass2_*(s+Mass2_))/pow(s-Mass2_+m_D2s,2);
    double X3 = pow(alpha_u, 2)* 4./9.*((s-Mass2_)*(Mass2_-u) + 2*Mass2_*(u+Mass2_))/pow(Mass2_-u + m_D2u,2);
    double X4 = alpha_s*alpha_u* 1./9.*Mass2_*(4*Mass2_-t)/((s-Mass2_ + m_D2s)*(Mass2_-u + m_D2u));
    double X5 = alpha_t*alpha_s* ((s-Mass2_)*(Mass2_-u) + Mass2_*(s-u))/((t-kfactor_*m_D2t)*(s-Mass2_ + m_D2s));
    double X6 = alpha_t*alpha_u* ((s-Mass2_)*(Mass2_-u) - Mass2_*(s-u))/((t-kfactor_ * m_D2t)*(Mass2_-u + m_D2u));
    double result = pow(4.*M_PI,2)*(X1+X2+X3+X4+X5-X6);
    //cout << s << " " << t << " " <<  X1 << endl;
    return result;
}



double Wrapper_M2_Qg(double t, void *params)
{
    struct gsl_params2to2 *Wp = (struct gsl_params2to2 *) params;
    double s = Wp->s_;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->M2_Qg(t, s, temp);
    return result;
}



double Scattering_2to2::Sigma_Qg(double s, double temp)
{
    struct gsl_params2to2 Wp = {};
    Wp.pt_class_ = this;
    Wp.s_= s;
    Wp.temp_ = temp;

    double result, error;
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &Wrapper_M2_Qg;
    F.params = &Wp;

    double tmax = 0;
    double tmin = -pow(s-Mass2_, 2)/s;
    double relerr = 1e-5;
//    cout << tmax << " " << tmin << " " << Mass2_ << endl;
    int status = 1;
    int nloop = 0;
    while(status && nloop < 10)
    {
        status = gsl_integration_qags(&F, tmin, tmax, 0, relerr, 1000, w, &result, &error);
        relerr *= 1.5;
        nloop += 1;
    }
    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

    double cross_section = result / (16.0*M_PI*pow(s-Mass2_,2));
    return cross_section;
}


void Scattering_2to2::Sigma_Qg_tabulate()
{
    double s = s0_;
    double temp;
    int is, iT;
    for (is = 0; is < Ns_MAX_; ++is)
    {
        s = exp(log(s0_) + is*ds_);
        for (iT = 0; iT < NT_MAX_; ++iT)
        {
	    temp = T0_ + iT*dT_;
            sigma_Qg_[is][iT] = Sigma_Qg(s, temp);
//	    cout << s << " " << temp << " " << sigma_Qg_[is][iT] << endl;
        }
    }
    cout << "Sigma_Qg_tabulate success! " << endl;
   // if (is == Ns_MAX_) std::cout << "tabulate sigma finished :).\n";
   // else std::cout << "ERROR: Scattering_rate::Sigma_tabulate failed!\n";
}


double Scattering_2to2::Get_sigma_Qg(double s, double temp)
{
    int is, iT;
    double delta_s, delta_T;
    double var_s, var_T;
   
    var_s = (log(s) - log(s0_)) / ds_;
    var_T = (temp-T0_)/dT_;

    is = int(var_s);
    iT = int(var_T);

    delta_s = var_s - is;
    delta_T = var_T - iT;

    //cout << is << " " << iT << endl;
    if (is <0 || is >= (Ns_MAX_-1) || iT < 0 || iT >= (NT_MAX_-1))
    {
     // cout << "energy/temperature out of range :(" << endl;
//	cout << s << " " << temp << endl;
	return sigma_Qg_[min(is, Ns_MAX_ -1)][min(iT, NT_MAX_-1)];   // need to be more careful aboutthis part.....
    }
    else
    {
        if ((1-delta_s) < 1e-8)	{ is += 1; delta_s = 0.;}
	if ((1-delta_T) < 1e-8) { iT += 1; delta_T = 0.;}
	double T1 = pow(T0_ + dT_*iT,2);
	double T2 = pow(T1 + dT_, 2);
	double result = cubeInterp(delta_s, delta_T, sigma_Qg_[is][iT]*T1, sigma_Qg_[is+1][iT]*T1, sigma_Qg_[is][iT+1]*T2, sigma_Qg_[is+1][iT+1]*T2);
	return result/(temp*temp);
    }
}



double Wrapper_get_sigma_Qg(double s, void *params)
{
    struct gsl_params2to2 *Wp = (struct gsl_params2to2 *) params;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->Get_sigma_Qg(s, temp);
    double Mass2 = Wp->pt_class_->Mass2_;
    return result *(s-Mass2);
}




double Scattering_2to2::Omega_Qg(double E1, double E2, double temp)
{
    struct gsl_params2to2 Wp={};
    Wp.pt_class_ = this;
    Wp.temp_= temp;

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;
    gsl_function O;
    O.function = &Wrapper_get_sigma_Qg;
    O.params = &Wp;

    double smin = Mass2_ + 2.*E1*E2 - 2.*E2*sqrt(E1*E1-Mass2_);
    double smax = Mass2_ + 2.*E1*E2 + 2.*E2*sqrt(E1*E1-Mass2_);
    double relerr = 1e-5;
//    cout << smin << " " << smax << " " << Mass2_ << " " << E1 << " " << E2 << endl;
    int status = 1;
    int nloop = 0;

    while (status && nloop < 10)
    {
        status = gsl_integration_qags(&O, smin, smax, 0, relerr, 1000, w, &result, &error);
        relerr *= 1.5;
        nloop += 1;
    }

    gsl_set_error_handler(old_handler);
    gsl_integration_workspace_free(w);

    return result;
}



double Wrapper_omega_Qg(double E2, void *params)
{
    struct gsl_params2to2 *Wp = (struct gsl_params2to2 *) params;
    double E1 = Wp->E1_;
    double temp = Wp->temp_;
    double result = Wp->pt_class_->Omega_Qg(E1, E2, temp);
    return result /(exp(E2/temp)-1);
}


double Scattering_2to2::Gamma_Qg(double E1, double temp)
{
    struct gsl_params2to2 Wp = {};
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
    double relerr = 1e-5;

    int status = 1;
    int nloop = 0;
    while(status && nloop < 10)
    {
        status = gsl_integration_qags(&G, E2min, E2max, 0, relerr, 1000, w, &result, &error);
        relerr *= 1.5;
        nloop += 1;
    }

    gsl_integration_workspace_free(w);
    gsl_set_error_handler(old_handler);

    double scattering_rate = result / (16.*M_PI*M_PI*E1*sqrt(E1*E1-Mass2_));
    return scattering_rate*16;
}







void Scattering_2to2::Gamma_Qg_tabulate()
{
    double E1, temp;
    for (int iE=0; iE<NE_MAX_; iE++)
    {
        E1 = E10_ + iE*dE1_;
        for (int iT=0; iT<NT_MAX_; iT++)
	{
	    temp = T0_ + iT*dT_;
	    gamma_Qg_[iE][iT] = Gamma_Qg(E1, temp);
	}
    }

    cout << "Gamma_Qq_tabulate success! " << endl;
}




double Scattering_2to2::Get_gamma_Qg(double E1, double temp)
{
    int iE, iT;
    double delta_E, delta_T;
    double var_E, var_T;

    var_E = (E1 - E10_)/dE1_;
    var_T = (temp - T0_)/dT_;
    
    iE = int(var_E);
    iT = int(var_T);

    delta_E = var_E - iE;
    delta_T = var_T - iT;

    if (iE <0 || iE >=(NE_MAX_-1) || iT<0 || iT>=(NT_MAX_-1))
    {
        return gamma_Qg_[min(iE, NE_MAX_-1)][min(iT, NT_MAX_-1)];
    }
    else
    {
        if ((1-delta_E) < 1e-8) {iE += 1; delta_E =0.;}
	if ((1-delta_T) < 1e-8) {iT += 1; delta_T =0.;}
	double result = cubeInterp(delta_E, delta_T, gamma_Qg_[iE][iT], gamma_Qg_[iE+1][iT], gamma_Qg_[iE][iT+1], gamma_Qg_[iE+1][iT+1]);
	return result;
    }
}






















