/*
 *
 * Calculate the Chi-Squared P-Value
 * The functions: approx_gamma, igf and chisqr were from:
 * http://www.codeproject.com/Articles/432194/How-to-Calculate-the-Chi-Squared-P-Value
 * under the permission:
 * http://creativecommons.org/licenses/publicdomain
 *
 *
 */

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <unistd.h>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>

#include "param.h"

extern Param paramd;

/*
double gamma(double N) {
    const long double SQRT2PI = 2.5066282746310005024157652848110452530069867406099383;
    long double Z = (long double)N;
    long double Sc = powl((Z + A), (Z + 0.5));
    Sc *= expl(-1.0 * (Z + A));
    Sc /= Z;
 
    long double F = 1.0;
    long double Ck;
    long double Sum = SQRT2PI;

    for(int K = 1; K < A; K++) {
        Z++;
	Ck = powl(A - K, K - 0.5);
	Ck *= expl(A - K);
	Ck /= F;
	Sum += (Ck / Z);
	F *= (-1.0 * K);
    }
    return (double)(Sum * Sc);
}
*/

double approx_gamma(double Z) {
    const double RECIP_E = 0.36787944117144232159552377016147;  // RECIP_E = (E^-1) = (1.0 / E)
    const double TWOPI = 6.283185307179586476925286766559;  // TWOPI = 2.0 * PI

    double D = 1.0 / (10.0 * Z);
    D = 1.0 / ((12 * Z) - D);
    D = (D + Z) * RECIP_E;
    D = pow(D, Z);
    D *= sqrt(TWOPI / Z);
 
    return D;
} 

static double igf(double S, double Z) {
    if (Z < 0.0) { return 0.0; }
    
    double Sc = (1.0 / S);
    Sc *= pow(Z, S);
    Sc *= exp(-Z);

    double Sum, Nom, Denom;
    Sum = Nom = Denom = 1.0;

    for(int I = 0; I < 200; I++) {
	Nom *= Z;
	S++;
	Denom *= S;
	Sum += (Nom / Denom);
    }

    return Sum * Sc;
}

double chisqr(int Dof, double Cv) {
    if (Cv < 0 || Dof < 1) { return 0.0; }
    double K = ((double)Dof) * 0.5;
    double X = Cv * 0.5;
    if (Dof == 2) { return exp(-1.0 * X); }
    double PValue = igf(K, X);
    if (std::isnan(PValue) || std::isinf(PValue) || PValue <= 1e-8) { return 1e-14; } 
    PValue /= approx_gamma(K);
	
    return (1.0 - PValue);
}

double get_chisqr_p(unsigned short * first, unsigned short * second) {
    double CriticalValue = 0.0;
    double XSqr; 
    unsigned int NumberOfEffectiveDatapoints = 0;
    for (unsigned index = 0; index < 100; index++) {
        if (first[index] + second[index] > 0.0) NumberOfEffectiveDatapoints++;
        XSqr = first[index] - second[index];
        if (first[index] > 0.0) CriticalValue += (XSqr * XSqr) / first[index];
    }
    if (NumberOfEffectiveDatapoints < 2) NumberOfEffectiveDatapoints = 2;
    double returnP = chisqr(NumberOfEffectiveDatapoints, CriticalValue);
    std::cerr << CriticalValue << "\t" << NumberOfEffectiveDatapoints << "\t" << returnP << std::endl;

    return returnP;
}

double X2BetweenTwo(unsigned short * FirstOriginal, unsigned short * SecondOriginal, unsigned int dispots) {
    // default dispots = 100
    // delare
    //
    // coverage normalization ( in case tumor >> normal )
    //
    if (paramd.Normalization == 1) {
        double sum_first = 0;
        double sum_second = 0;
        for (int i = 0; i < dispots; i++) {
            sum_first += FirstOriginal[i];
            sum_second += SecondOriginal[i];
        }
        if ((sum_first > sum_second) && (sum_second !=0)) {
            for (int i = 0; i < dispots; i++) {
                SecondOriginal[i] = SecondOriginal[i] * (sum_first/sum_second);
            }
        }
        if ((sum_first <= sum_second) && (sum_first != 0)) {
            for (int i = 0; i < dispots; i++) {
                FirstOriginal[i] = FirstOriginal[i] * (sum_second/sum_first);
            }
        }
    }
    // normalization

    double *ExpFirst, *ExpSecond, *SumBoth;

    ExpFirst  = new double [dispots + 1];
    ExpSecond = new double [dispots + 1];
    SumBoth   = new double [dispots + 1];

    double SumFirst, SumSecond, SumTotal;
    SumFirst  = SumSecond = SumTotal  = 0.0;
    // calculation 
    for (unsigned i = 0; i < dispots; i++) {
        SumBoth[i] = FirstOriginal[i] + SecondOriginal[i];
        SumFirst += FirstOriginal[i];
        SumSecond += SecondOriginal[i];
    }
    SumTotal = SumFirst + SumSecond; 
    for (unsigned i = 0; i < dispots; i++) {
        ExpFirst[i] = SumBoth[i] * SumFirst / SumTotal;
        ExpSecond[i] = SumBoth[i] * SumSecond / SumTotal;
    }

    double result = 0.0; 
    unsigned Degree = 0;
    for (unsigned i = 0; i < dispots; i++) {
	if (FirstOriginal[i] + SecondOriginal[i] > 0.0) Degree++;
        if (ExpFirst[i]) result +=  (FirstOriginal[i] - ExpFirst[i]) * (FirstOriginal[i] - ExpFirst[i]) / ExpFirst[i];
        if (ExpSecond[i]) result +=  (SecondOriginal[i] - ExpSecond[i]) * (SecondOriginal[i] - ExpSecond[i]) / ExpSecond[i];
    }

    // release mem
    double PValue;
    
    if (Degree == 1 || result == 0)  PValue = 1.0;
    else {
	Degree--;
	PValue = chisqr(Degree, result);
    }

    if ( PValue < 0 ) PValue = PValue * (-1);

    delete [] ExpFirst;
    delete [] ExpSecond;
    delete [] SumBoth;

    return PValue;
}

