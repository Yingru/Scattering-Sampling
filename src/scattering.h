#ifndef scattering_H
#define scattering_H

#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_integration.h>
#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte_vegas.h>


typedef boost::multi_array<double, 2> array_2D;

class Scattering_2to2
{
private:
    double kfactor_;
    double Lambda2 = 0.04;
    double s0_, T0_, dT_, ds_, E10_, dE1_;
    int Ns_MAX_, NT_MAX_, NE_MAX_;
    double Alpha(double Q2);
    double DebyeMass2(double alpha_s, double temp);

public:
    Scattering_2to2(double M, double kfactor);
    ~Scattering_2to2() {};
    double Mass2_;
    array_2D sigma_Qq_, sigma_Qg_;
    array_2D gamma_Qq_, gamma_Qg_;

    double M2_Qq(double t, double s, double temp);
    double Sigma_Qq(double s, double temp);
    void Sigma_Qq_tabulate();
    double Get_sigma_Qq(double s, double temp);
    double Omega_Qq(double E1, double E2, double temp);
    double Gamma_Qq(double E1, double temp);
    void Gamma_Qq_tabulate();
    double Get_gamma_Qq(double E1, double temp);


    double M2_Qg(double t, double s, double temp);
    double Sigma_Qg(double s, double temp);
    void Sigma_Qg_tabulate();
    double Get_sigma_Qg(double s, double temp);
    double Omega_Qg(double E1, double E2, double temp);
    double Gamma_Qg(double E1, double temp);
    void Gamma_Qg_tabulate();
    double Get_gamma_Qg(double E1, double temp);
};


struct gsl_params2to2
{
    double temp_;
    double s_;
    double E1_;
    Scattering_2to2 * pt_class_;
};



class Scattering_2to3
{
private:
    double kfactor_;
    double Lambda2 = 0.04;
    double s0_, T0_, E10_,  dT_, ds_, dE1_;
    int  NT_MAX_, Ns_MAX_, NE_MAX_;
    double Alpha(double Q2);
    double DebyeMass2(double alpha_s, double temp);

public:
    Scattering_2to3(double M, double kfactor);
    ~Scattering_2to3() {};
    double Mass2_;
    array_2D sigma_Qq_, sigma_Qg_;
    array_2D gamma_Qq_, gamma_Qg_;

    int ValidParams(double *k, double s);
    double Jacobi_Qq(double *k, double s);

    double M2_Qq(double *k, double s, double temp);
    double Sigma_Qq(double s, double temp);
    void Sigma_Qq_tabulate();
    double Get_sigma_Qq(double s, double temp);
    double Omega_Qq(double E1, double E2, double temp);
    double Gamma_Qq(double E1, double temp);
    void Gamma_Qq_tabulate();
    double Get_gamma_Qq(double E1, double temp);

    double M2_Qg(double *k, double s, double temp);
    double Sigma_Qg(double s, double temp);
    void Sigma_Qg_tabulate();
    double Get_sigma_Qg(double s, double temp);
    double Omega_Qg(double E1, double E2, double temp);
    double Gamma_Qg(double E1, double temp);
    void Gamma_Qg_tabulate();
    double Get_gamma_Qg(double E1, double temp);

};



struct gsl_params2to3
{
    double s_;
    double temp_;
    double E1_;
    Scattering_2to3 *pt_class_;
};

#endif
