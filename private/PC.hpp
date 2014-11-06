#ifndef _PAIRCOPULANEGLL_HPP_
#define _PAIRCOPULANEGLL_HPP_

#include "VineCPP_header.hpp"

void LoadBounds(double *bounds);

void PairCopulaAIC(double *AIC, double *theta,int family, double *U,double *V,unsigned int n);
void PairCopulaAIC(double *AIC, double *theta,int family, int rotation, double *U,double *V,unsigned int n);
void PairCopulaAIC_Rotated_Obs(double *AIC, double *theta,int family, int rotation, double *U,double *V,unsigned int n);

void PairCopulaCDF(int family, const double *theta, double *U, double *V, double *p, unsigned int n);
void PairCopulaCDF(int family, int rotation, const double *theta, double *U, double *V, double *p, unsigned int n);

void PairCopulaFit(double *theta,int family,double *U,double *V,unsigned int n);
void PairCopulaFit(double *theta,int family, int rotation, double *U,double *V,unsigned int n);
void PairCopulaFit_Rotated_Obs(double *theta,int family, int rotation, double *U,double *V,unsigned int n);

void PairCopulaHfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaVfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaHfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaVfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaHfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaVfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);

void PairCopulaIndepTest(double *pValue, double *U,double *V,unsigned int n);
void PairCopulaIndepTest(double *pValue, double *U,double *V,unsigned int n,double *TestStat);
bool PairCopulaIndepTest(double *U,double *V,unsigned int n, double *tau);

void PairCopulaInvHfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvVfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvHfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvVfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvHfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvVfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);

void PairCopulaInvHfun_VecPar(int family, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvVfun_VecPar(int family, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvHfun_VecPar(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);
void PairCopulaInvVfun_VecPar(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n);

double PairCopulaNegLL(int family, const double *theta, double *U, double *V, unsigned int n);
double PairCopulaNegLL(int family, int rotation, const double *theta, double *U, double *V, unsigned int n);
double PairCopulaNegLL_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, unsigned int n);

void PairCopulaPDF(int family, const double *theta, double *U, double *V, double *p, unsigned int n);
void PairCopulaPDF(int family, int rotation, const double *theta, double *U, double *V, double *p, unsigned int n);
void PairCopulaPDF_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *p, unsigned int n);

void PairCopulaRand(int family, const double *theta, double *U, double *V, unsigned int n);
void PairCopulaRand(int family, int rotation, const double *theta, double *U, double *V, unsigned int n);

void PairCopulaSelect(int *family, double *theta, int *rotation, double *U, double *V, unsigned int n, double *familyset, int m);

double SD_Kendall_Tau(double *U, double *V, unsigned int n);
void SD_Kendall_Tau_Matrix(double *tau, double *U, unsigned int d, unsigned int n);

void VineCopulaFit(VineCopula* Vine, double *CLL, double *x0, double *U, unsigned int n);
void VineCopulaFitSeq(VineCopula* Vine, double *CLL, double *U, unsigned int n);

void VineCopulaGetPseudoObs(VineCopula* Vine, double *U, double *V, unsigned int n);
void VineCopulaGetPseudoObs(VineCopula* Vine, double *U, double *H, double *V, unsigned int n);

double VineCopulaNegLL(VineCopula* Vine, double *U, int CutOffTree, unsigned int n);
double VineCopulaNegLL(const double *Thetas, VineCopula* Vine, double *U, int CutOffTree, unsigned int n);

void VineCopulaRand(VineCopula* Vine, double *U, unsigned int n);

void VineCopulaStructureSelect(int type, double *Structure, double *Families, double *Rotations, std::vector<double>& Thetas, double *U, unsigned int d, unsigned int n, int StructuringRule, double *familyset, int m);

#endif
