#ifndef sample_H
#define sample_H

#include <iostream>
#include <stdlib.h>
#include <random>
#include <time.h>
#include "scattering.h"

class Sample_2to2
{
private:
    double E1_, temp_;

    double Alpha(double Q2);
    double DebyeMass2(double alpha_s, double temp);


 //  double Max_E2ANDs_Qq(const Scattering_2to2& A, double E1, double temp);
  //  int Sample_E2ANDs_Qq(const Scattering_2to2& A , double E1, double temp);
  //  int Sample_t_Qq(const Scattering_2to2& A, double s, double temp);
public:
    Sample_2to2(Scattering_2to2& A, double E1, double temp);
//    Sample_2to2(const Scattering_2to2& A, double E1, double temp);
    ~Sample_2to2(){} ;    
    double p30_, p3x_, p3y_, p3z_;
    double* Sample_vector_CoM(Scattering_2to2& A, double t, double s);

    double E2_, s_, t_;
 
    double Max_E2ANDs_Qq(Scattering_2to2& A, double E1, double temp);
    int Sample_E2ANDs_Qq(Scattering_2to2& A, double E1, double temp);
    int Sample_t_Qq(Scattering_2to2 &A, double s, double temp);
 
    double Max_E2ANDs_Qg(Scattering_2to2& A, double E1, double temp);
    int Sample_E2ANDs_Qg(Scattering_2to2& A, double E1, double temp);
    int Sample_t_Qg(Scattering_2to2& A, double s, double temp);
};


class Sample_2to3
{
private:
    double E1_, temp_;
    double Alpha(double Q2);
    double DebyeMass2(double alpha_s, double temp);

public:
    Sample_2to3(Scattering_2to3& A, double E1, double temp);
    ~Sample_2to3(){};

    double p30_, p3x_, p3y_, p3z_;
    double* Sample_vector_CoM(Scattering_2to3& A, double s, double *k);
    double E2_, s_, p40_, k0_, theta4_, phi4_, phi4k_;
    double Max_E2ANDs_Qq(Scattering_2to3& A, double E1, double temp);
    int Sample_E2ANDs_Qq(Scattering_2to3& A, double E1, double temp);
    int Sample_p3_Qq(Scattering_2to3& A, double s, double temp);

    double Max_E2ANDs_Qg(Scattering_2to3& A, double E1, double temp);
    int Sample_E2ANDs_Qg(Scattering_2to3& A, double E1, double temp);
    int Sample_p3_Qg(Scattering_2to3& A, double s, double temp);
};


#endif
