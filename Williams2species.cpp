
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include<string>
#include<iomanip>
#include <math.h>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <typeinfo>

#define pi 3.141592653589793

using namespace std;
using namespace Eigen;



//omega is the rotation
// 2q is the degree of legendre polynomial to use
int q = 5;
int depth = 6;
int gridsize = pow(2, depth)+1;
double e = pow(2, 0.5);
double mp = 1;
double me = mp / 10;
double G = 1;
double h = 10;
double c = 2 * pi * 137.035999084 * e * e / h;
double omega = 0.00;
double pc = 100;
double ec = 56.1;
// I have put an extra pi in here to match the nonrotating results : there is some difference because in michael's writing he 
// adjusts the value of the reduced planck's constant, while here I have adjusted planck's constant
double kappa = h * c / 5 * pow(3.0 / (8 * pi),(2.0 / 3));
double betap = pow(pow(me + mp, 2) * (pow(G * (me * me + mp * mp) - e * e,2) + 2 * G * G * (me * me * mp * mp - me * me * me * mp - mp * mp * mp * me) + 6 * G * mp * me * e * e),0.5);
double gammap = G * (me * me * me - mp * mp * mp) + e * e * (mp - me);
double gammam = G * (me * me * me + mp * mp * mp) - e * e * (mp + me);
double alphap = pow((gammam + betap) / (2 * kappa * me), 0.5);
double alpham = pow(abs((gammam - betap) / (2 * kappa * me)), 0.5);

string outputstr = "iteration0outputcpp.txt";
string surfacestr = "iteration0surfacecpp.txt";


double rhop(double c1List[], double c2List[], double x,  double theta) {
    double sum = 2 * pow(omega,2) / (G * (me + mp));
    for (int m = 0;m < q+1; m++) {
        if (x == 0) {
            sum += legendre(2 * m, cos(theta)) * (c1List[m] * sph_bessel(2 * m, alphap * x) + c2List[m] * pow(pi / (2 * abs(alpham) * x), 0.5) * cyl_bessel_i(2 * m + 0.5, abs(alpham) * x) * pow(-1, m));
        }
        else {
            sum+= legendre(2 * m, cos(theta)) * (c1List[m] * sph_bessel(2 * m, alphap * x) + c2List[m] * pow(pi / (2 * abs(alpham) * x), 0.5) * cyl_bessel_i(2 * m + 0.5, abs(alpham) * x) * pow(-1, m));
        }
    }
    return sum;
}

double rhoe(double c1List[], double c2List[], double x, double theta) {
    double sum = 2 * pow(omega,2) / (G * (me + mp));
    for (int m = 0;m < q+1; m++) {
        sum += legendre(2 * m, cos(theta)) * ((betap + gammap) * c1List[m] * sph_bessel(2 * m, alphap * x) + (-betap + gammap) * c2List[m] * pow(pi / (2 * abs(alpham) * x), 0.5) * cyl_bessel_i(2 * m + 0.5, abs(alpham) * x) * pow(-1, m));   
    }
    return sum / (2 * mp * (e * e + G * me * mp));
}


double pdrf(double c1List[], double c2List[], double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q; m++) {
        sum += legendre(2 * m, cos(theta)) * (c1List[m] * alphap * boost::math::sph_bessel_prime(2 * m, alphap * x) + c2List[m] * pow(-1, m) * pow(pi / (2 * abs(alpham)), 0.5) * (boost::math::cyl_bessel_i_prime(2 * m + 0.5, abs(alpham) * x) * abs(alpham) / pow(x, 0.5) - 1.0 / 2.0 / pow(x, 1.5) * cyl_bessel_i(2 * m + 0.5, abs(alpham) * x)));
    }
    return sum;
}

double pdtf(double c1List[], double c2List[],double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q; m++) {

        sum -= boost::math::legendre_p_prime(2 * m, cos(theta)) * sin(theta) / x * (c1List[m] * sph_bessel(2 * m, alphap * x) + c2List[m] * pow(-1, m) * pow(pi / (2 * abs(alpham) * x), 0.5) * cyl_bessel_i(2 * m + 0.5, abs(alpham) * x));
    } 
    return sum;
}

double  edrf(double c1List[], double c2List[],double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q; m++) {
        sum += legendre(2 * m, cos(theta)) / (2 * mp * (e * e + G * me * mp)) * ((betap + gammap) * c1List[m] * alphap * boost::math::sph_bessel_prime(2 * m, alphap * x) + (-betap + gammap) * c2List[m] * pow(-1,m) * pow(pi / (2 * abs(alpham)),0.5) * (boost::math::cyl_bessel_i_prime(2 * m + 0.5, abs(alpham) * x) * abs(alpham) / pow(x, 0.5) - 1.0 / 2.0 / pow(x, 1.5) * cyl_bessel_i(2 * m + 0.5, abs(alpham) * x)));
    }
    return sum;
}

double edtf(double c1List[], double c2List[], double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q; m++) {
        sum -= boost::math::legendre_p_prime(2 * m, cos(theta)) * sin(theta) / (x * 2 * mp * (e * e + G * me * mp)) * ((betap + gammap) * c1List[m] * sph_bessel(2 * m, alphap * x) + (-betap + gammap) * c2List[m] * pow(-1, m) * pow(pi / (2 * abs(alpham) * x), 0.5) * cyl_bessel_i(2 * m + 0.5, abs(alpham) * x));
    }
    return sum;
}

double FindpSurface(double c1List[], double c2List[], double theta) {                   
    //there may be multiple zeroes which is a problem
     // it is clear that at x = 0, sigma is positive, so we will just increase until it becomes negative
    double x1 = 0;
    double x2 = 0.3;
    double epsilon = 0.0000001;
    double s = rhop(c1List, c2List, x2, theta);
    int i = 0;
    do {
        x1 += 0.3;
        x2 += 0.3;
        s = rhop(c1List, c2List, x2, theta);
        i += 1;
        //cout << s << endl;
    } while (s > 0 && i < 100);
//now we have a zero in[x1, x2]
//implement bisection
    do {
        //cout << rhop(c1List, c2List, (x1 + x2) / 2, theta) << " "<< x1<< " "<< x2<< endl;
        if (rhop(c1List, c2List, (x1 + x2) / 2, theta) > 0) {
            x1 = (x1 + x2) / 2;
        }
        else {
            x2 = (x1 + x2) / 2;
            i += 1;
        }
    } while (abs(rhop(c1List, c2List, (x1 + x2) / 2, theta)) > epsilon && i < 1000);
return x1;
}

double FindeSurface(double c1List[], double c2List[], double theta) {
    //there may be multiple zeroes which is a problem
     // it is clear that at x = 0, sigma is positive, so we will just increase until it becomes negative
    double x1 = 0;
    double x2 = 0.3;
    double epsilon = 0.0000001;
    double s = rhoe(c1List, c2List, x2, theta);
    int i = 0;
    do {
        x1 += 0.3;
        x2 += 0.3;
        s = rhoe(c1List, c2List, x2, theta);
        i += 1;
    } while (s > 0 && i < 100);
    //now we have a zero in[x1, x2]
    //implement bisection
    do {
        if (rhoe(c1List, c2List, (x1 + x2) / 2, theta) > 0) {
            x1 = (x1 + x2) / 2;
        }
        else {
            x2 = (x1 + x2) / 2;
            i += 1;
        }
    } while (abs(rhoe(c1List, c2List, (x1 + x2) / 2, theta)) > epsilon && i < 1000);
    return x1;
}

// the next two definition are for the potentials in the bulk region where we can find both explicitly
double PsieI(double c1List[], double c2List[], double x, double theta) {
    return (4 * pi * (e * e + G * me * mp) / (pow(e * e + G * me * mp, 2) - (e * e - G * me * me) * (e * e - G * mp * mp))) * (kappa * (me / mp * rhop(c1List, c2List, x, theta) + (e * e - G * mp * mp) / (e * e + G * me * mp) * rhoe(c1List, c2List, x, theta)) - omega * omega / 3 * x * x * (1 - legendre(2, cos(theta))) * (mp + me * (e * e - G * mp * mp) / (e * e + G * me * mp)));
}

double PsipI(double c1List[],double c2List[],double x,double theta) {
    return (4 * pi * (e * e - G * me * me) / ((e * e - G * me * me) * (e * e - G * mp * mp) - pow(e * e + G * me * mp, 2))) * (-kappa * (me / mp * rhop(c1List, c2List, x, theta) + (e * e + G * mp * me) / (e * e - G * me * me) * rhoe(c1List, c2List, x, theta)) + omega * omega / 3 * x * x * (1 - legendre(2, cos(theta))) * (mp + me * (e * e + G * me * mp) / (e * e - G * me * me)));
}
//the next two definitions are for the potentials outside of that particles support
double  PsipO(double bList[], double v, double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q + 1; m++) {
        if (x > 0) {
            sum += legendre(2 * m, cos(theta)) * bList[m] / pow(x, (2 * m + 1));
        }
    }
    return (v + sum);
}

double PsieO(double aList[], double phi, double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q + 1; m++) {
        sum += legendre(2 * m, cos(theta)) * aList[m] / pow(x, 2 * m + 1);
    }
    return (phi + sum);
}

            // The next two definitions are for the potentials in the atmospheric range, we therefore use an A in their names.Note that in any iteration only one of these functions should be called.
double PsieA(double c1List[], double c2List[], double bList[], double v, double x, double theta) {
    return (-4 * pi * kappa * rhoe(c1List, c2List, x, theta) + (e * e + G * me * mp) * PsipO(bList, v, x, theta) + 4 * pi * omega * omega / 3 * x * x * (1 - legendre(2, cos(theta)))) / (e * e - G * me * me);
}

double PsipA(double c1List[], double c2List[], double aList[], double phi, double x, double theta) {
    return (-4 * pi * kappa * me / mp * rhop(c1List, c2List, x, theta) + (e * e + G * me * mp) * PsieO(aList, phi, x, theta) + 4 * pi * omega * omega / 3 * x * x * (1 - legendre(2, cos(theta)))) / (e * e - G * mp * mp);
}

double PsipOdr(double bList[], double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q + 1; m++) {
        sum -= legendre(2 * m, cos(theta)) * bList[m] * (2 * m + 1) / (pow(x,2 * m + 2));
    }
    return sum;
}

double PsipOdt(double bList[], double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q + 1; m++) {
        sum -= boost::math::legendre_p_prime(2 * m, cos(theta)) * sin(theta) * bList[m] / pow(x, (2 * m + 2));
    }
    return sum;
}

double PsieOdr(double aList[], double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q + 1; m++) {
        sum -= legendre(2 * m, cos(theta)) * aList[m] * (2 * m + 1) / pow(x ,(2 * m + 2));
    }
    return sum;
}

double PsieOdt(double aList[], double x, double theta) {
    double sum = 0;
    for (int m = 0;m < q + 1; m++) {
        sum -= boost::math::legendre_p_prime(2 * m,cos(theta)) * sin(theta) * aList[m] / pow(x ,(2 * m + 2));
    }
    return sum;
}

double GradpO(double c1List[], double c2List[], double bList[], double x, double theta) {
    double pdr = pdrf(c1List, c2List, x, theta);
    double pdt = pdtf(c1List, c2List, x, theta);
    double dpsidr = PsipOdr(bList, x, theta);
    double dpsidt = PsipOdt(bList, x, theta);
    return(pdr * dpsidr * (sin(pdt) * sin(dpsidt) + cos(pdt) * cos(dpsidt)));
}


double GradeO(double c1List[], double c2List[], double aList[], double x, double theta) {
    double edr = edrf(c1List, c2List, x, theta);
    double edt = edtf(c1List, c2List, x, theta);
    double dpsidr = PsieOdr(aList, x, theta);
    double dpsidt = PsieOdt(aList, x, theta);
    return(edr * dpsidr * (sin(edt) * sin(dpsidt) + cos(edt) * cos(dpsidt)));
}


double GradpI(double c1List[], double c2List[], double x, double theta) {
    double pdr = pdrf(c1List, c2List, x, theta);
    double edr = edrf(c1List, c2List, x, theta);
    double pdt = pdtf(c1List, c2List, x, theta);
    double edt = edtf(c1List, c2List, x, theta);
    double dpsidr = (4 * pi * (e * e - G * me * me) / ((e * e - G * me * me) * (e * e - G * mp * mp) - pow(e * e + G * me * mp, 2))) * (-kappa * (me / mp * pdr + (e * e + G * mp * me) / (e * e - G * me * me) * edr) + 2 * omega * omega / 3 * x * (1 - legendre(2, cos(theta))) * (mp + me * (e * e + G * me * mp) / (e * e - G * me * me)));
    double dpsidt = (4 * pi * (e * e - G * me * me) / ((e * e - G * me * me) * (e * e - G * mp * mp) - pow(e * e + G * me * mp, 2))) * (-kappa * (me / (x * mp) * pdt + (e * e + G * mp * me) / (x * (e * e - G * me * me)) * edt) + omega * omega / 3 * x * (1 - boost::math::legendre_p_prime(2, cos(theta)) * sin(theta)) * (mp + me * (e * e + G * me * mp) / (e * e - G * me * me)));
    return(pdr * dpsidr * (sin(pdt) * sin(dpsidt) + cos(pdt) * cos(dpsidt)));
}

double GradpA(double c1List[], double c2List[], double aList[], double x, double theta) {
    double pdr = pdrf(c1List, c2List, x, theta);
    double edr = PsieOdr(aList, x, theta);
    double pdt = pdtf(c1List, c2List, x, theta);
    double edt = PsieOdt(aList, x, theta);
    double dpsidr = 1 / (e * e - G * mp * mp) * (-4 * pi * kappa * me / mp * pdr + (e * e + G * me * mp) * edr + 8 * pi * omega * omega / 3 * x * (1 - legendre(2, cos(theta))));
    double dpsidt = 1 / (x * (e * e - G * mp * mp)) * (-4 * pi * kappa * me / mp * pdt + (e * e + G * me * mp) * edt + 4 * pi * omega * omega / 3 * x * x * boost::math::legendre_p_prime(2, cos(theta)) * sin(theta));
    return(pdr * dpsidr * (sin(pdt) * sin(dpsidt) + cos(pdt) * cos(dpsidt)));
}

double GradeI(double c1List[], double c2List[], double x, double theta) {
    double pdr = pdrf(c1List, c2List, x, theta);
    double edr = edrf(c1List, c2List, x, theta);
    double pdt = pdtf(c1List, c2List, x, theta);
    double edt = edtf(c1List, c2List, x, theta);
    double dpsidr = (4 * pi * (e * e + G * me * mp) / (pow(e * e + G * me * mp, 2) - (e * e - G * me * me) * (e * e - G * mp * mp))) * (kappa * (me / mp * pdr + (e * e - G * mp * me) / (e * e + G * me * mp) * edr) - 2 * omega * omega / 3 * x * (1 - legendre(2,cos(theta))) * (mp + me * (e * e - G * me * mp) / (e * e + G * mp * me)));
    double dpsidt = (4 * pi * (e * e + G * me * me) / (pow(e * e + G * me * mp, 2) - (e * e - G * me * me) * (e * e - G * mp * mp))) * (-kappa * (me / (x * mp) * pdt + (e * e - G * mp * mp) / (x * (e * e + G * mp * me)) * edt) - omega * omega / 3 * x * (1 - boost::math::legendre_p_prime(2, cos(theta)) * sin(theta)) * (mp + me * (e * e - G * mp * mp) / (e * e + G * mp * me)));
    return(edr * dpsidr * (sin(edt) * sin(dpsidt) + cos(edt) * cos(dpsidt)));
}

double GradeA(double c1List[], double c2List[], double bList[], double x, double theta) {
    double edr = edrf(c1List, c2List, x, theta);
    double pdr = PsipOdr(bList, x, theta);
    double edt = edtf(c1List, c2List, x, theta);
    double pdt = PsipOdt(bList, x, theta);
    double dpsidr = 1 / (e * e - G * me * me) * (-4 * pi * kappa * edr + (e * e + G * me * mp) * pdr + 8 * pi * omega * omega / 3 * x * (1 - legendre(2, cos(theta))));
    double dpsidt = 1 / (x * (e * e - G * me * me)) * (-4 * pi * kappa * edt + (e * e + G * me * mp) * pdt + 4 * pi * omega * omega / 3 * x * x * boost::math::legendre_p_prime(2, cos(theta)) * sin(theta));
    return(edr * dpsidr * (sin(edt) * sin(dpsidt) + cos(edt) * cos(dpsidt)));
}

//this is the integration for (f, g, h, k).It uses trapezoidal ruleand Romberg integration, as in Williams.
double Rombergf(int i, double** surfaceList, double c1List[], double c2List[], double bList[], double phi, double aList[], double v, char flag) {
    //first we create the necessary matrix
    double Rmat[depth][depth];
    double sum;
    int step;

    if (flag == 'E') {
        //compute the differences using the bulk definition
        Rmat[0][0] = pi / 2 * ((PsipI(c1List, c2List, surfaceList[0][1], surfaceList[0][0]) - PsipO(bList, v, surfaceList[0][1], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (PsipI(c1List, c2List, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0]) - PsipO(bList, v, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    else {
        //compute the differences using the atmospheric definition
        Rmat[0][0] = pi / 2 * ((PsipA(c1List, c2List, aList, phi, surfaceList[0][1], surfaceList[0][0]) - PsipO(bList, v, surfaceList[0][1], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (PsipA(c1List, c2List, aList, phi, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0]) - PsipO(bList, v, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    for (int j = 0;j < depth; j++) {

        sum = 0;
        step = int((gridsize - 1) / pow (2,(j + 1)));
        for (int k = 0;k < pow(2, j); k++) {

            if (flag == 'E') {
                sum += (PsipI(c1List, c2List, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0]) - PsipO(bList, v, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) * sin(surfaceList[(2 * k + 1) * step][0]);
            }
            else {
                sum += (PsipA(c1List, c2List, aList, phi, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0]) - PsipO(bList, v, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) * sin(surfaceList[(2 * k + 1) * step][0]);
            }    
        }
        Rmat[j + 1][0] = 0.5 * (Rmat[j][0] + pi / (2 * pow(2, j)) * sum);
        for (int k = 0;k < j + 1; k++) {
            Rmat[j + 1][k+1]=(Rmat[j + 1][k] + 1 / (pow(4, k + 1) - 1) * (Rmat[j + 1][k] - Rmat[j][k]));
        }
    }
    //for (int j = 0; j < depth; j++) {
    //    for (int k = 0; k < depth; k++) {
    //        cout<< Rmat[j][k]<<" ";
    //    }
    //    cout << endl;
    //}
    return (Rmat[depth][depth]);
}

double Rombergg(int i, double** surfaceList, double c1List[], double c2List[], double bList[], double phi, double aList[], double v, char flag) {
    double Rmat[depth][depth];
    double sum;
    int step;
    if (flag == 'E') {
        //compute the differences using the bulk definition
        Rmat[0][0] = pi / 2 * ((PsieA(c1List, c2List, bList, v, surfaceList[0][2], surfaceList[0][0]) - PsieO(aList, phi, surfaceList[0][2], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (PsieA(c1List, c2List, bList, v, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0]) - PsieO(aList, phi, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    else {
        //compute the differences using the atmospheric definition
        Rmat[0][0]=pi / 2 * ((PsieI(c1List, c2List, surfaceList[0][2], surfaceList[0][0]) - PsieO(aList, phi, surfaceList[0][2], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (PsieI(c1List, c2List, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0]) - PsieO(aList, phi, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    for (int j = 0;j < depth; j++) {

        sum = 0;
        step = int((gridsize - 1) / pow(2,(j + 1)));
        for (int k = 0;k < pow(2, j); k++) {
            if (flag == 'E') {
                sum += (PsieA(c1List, c2List, bList, v, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0]) - PsieO(aList, phi, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) * sin(surfaceList[(2 * k + 1) * step][0]);
            }
            else {
                sum += (PsieI(c1List, c2List, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0]) - PsieO(aList, phi, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) *sin(surfaceList[(2 * k + 1) * step][0]);
            }
        }
        Rmat[j + 1][0] = (0.5 * (Rmat[j][0] + pi / (2 * pow(2, j)) * sum));
        for (int k = 0;k < j + 1; k++) {
            Rmat[j + 1][k+1]=(Rmat[j + 1][k] + 1 / (pow(4, (k + 1)) - 1) * (Rmat[j + 1][k] - Rmat[j][k]));
        }
    }
    return (Rmat[depth][depth]);
}

double Rombergh(int i, double** surfaceList, double c1List[], double c2List[], double bList[], double phi, double aList[], double v, char flag) {
    double Rmat[depth][depth];
    double sum;
    int step;
    if (flag == 'E') {
        Rmat[0][0] = pi / 2 * ((GradpI(c1List, c2List, surfaceList[0][1], surfaceList[0][0]) - GradpO(c1List, c2List, bList, surfaceList[0][1], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (GradpI(c1List, c2List, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0]) - GradpO(c1List, c2List, bList, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    else {
        Rmat[0][0]=pi / 2 * ((GradpA(c1List, c2List, aList, surfaceList[0][1], surfaceList[0][0]) - GradpO(c1List, c2List, bList, surfaceList[0][1], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (GradpA(c1List, c2List, bList, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0]) - GradpO(c1List, c2List, bList, surfaceList[gridsize - 1][1], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    for (int j = 0;j < depth; j++) {
        sum = 0;
        step = int((gridsize - 1) / pow(2,(j + 1)));
        for (int k = 0;k < pow(2, j); k++) {
            if (flag == 'E') {
                sum += (GradpI(c1List, c2List, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0]) - GradpO(c1List, c2List, bList, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) * sin(surfaceList[(2 * k + 1) * step][0]);
            }
            else {
                sum += (GradpA(c1List, c2List, aList, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0]) - GradpO(c1List, c2List, bList, surfaceList[(2 * k + 1) * step][1], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) * sin(surfaceList[(2 * k + 1) * step][0]);
            }
        }
        Rmat[j + 1][0] = (0.5 * (Rmat[j][0] + pi / (2 * pow(2, j)) * sum));
        for (int k = 0;k < j + 1; k++) {
            Rmat[j + 1][k+1]=(Rmat[j + 1][k] + 1 / (pow(4,(k + 1)) - 1) * (Rmat[j + 1][k] - Rmat[j][k]));
        }
    }
    return (Rmat[depth][depth]);
}

double Rombergk(int i, double** surfaceList, double c1List[], double c2List[], double bList[], double phi, double aList[], double v, char flag) {
    double Rmat[depth][depth];
    double sum;
    int step;
    if (flag == 'E') {
        Rmat[0][0] = pi / 2 * ((GradeA(c1List, c2List, bList, surfaceList[0][2], surfaceList[0][0]) - GradeO(c1List, c2List, aList, surfaceList[0][2], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (GradeA(c1List, c2List, bList, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0]) - GradeO(c1List, c2List, aList, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    else {
        Rmat[0][0] = pi / 2 * ((GradeI(c1List, c2List, surfaceList[0][2], surfaceList[0][0]) - GradeO(c1List, c2List, aList, surfaceList[0][2], surfaceList[0][0])) * legendre(2 * i - 2, cos(surfaceList[0][0])) * (sin(surfaceList[0][0])) + (GradeI(c1List, c2List, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0]) - GradeO(c1List, c2List, aList, surfaceList[gridsize - 1][2], surfaceList[gridsize - 1][0])) * legendre(2 * i - 2, cos(surfaceList[gridsize - 1][0])) * (sin(surfaceList[gridsize - 1][0])));
    }
    for (int j = 0;j < depth; j++) {
        sum = 0;
        step = int((gridsize - 1) / pow(2 ,(j + 1)));
        for (int k = 0;k < pow(2, j); k++) {
            if (flag == 'E') {
                sum += (GradeA(c1List, c2List, bList, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0]) - GradeO(c1List, c2List, aList, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) * sin(surfaceList[(2 * k + 1) * step][0]);
            }
            else {
                sum += (GradeI(c1List, c2List, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0]) - GradeO(c1List, c2List, aList, surfaceList[(2 * k + 1) * step][2], surfaceList[(2 * k + 1) * step][0])) * legendre(2 * i - 2, cos(surfaceList[(2 * k + 1) * step][0])) * sin(surfaceList[(2 * k + 1) * step][0]);
            }
        }
        Rmat[j + 1][0] = (0.5 * (Rmat[j][0] + pi / (2 * pow(2,j)) * sum));
        for (int k = 0;k < j + 1; k++) {
            Rmat[j + 1][k+1]=(Rmat[j + 1][k] + 1 / (pow(4,(k + 1)) - 1) * (Rmat[j + 1][k] - Rmat[j][k]));
        }
    }
    return (Rmat[depth][depth]);
}

int main() {
    cout.precision(16);
    // now we can begin the iteration
    //our starting point should be a nonrotating configuration

     double thetaGrid [gridsize];
     for (int i = 0;i < gridsize; i++) {
         thetaGrid[i] = (pi / (2 * (gridsize - 1)) * i);
     }


     Vector2d iC(pc, ec);
     Matrix2d startMat;
     startMat(0, 0) = 1;
     startMat(0, 1) = 1;
     startMat(1, 0) = (betap + gammap) / (2 * mp * (e * e + G * mp * me));
     startMat(1, 1) = (-betap + gammap) / (2 * mp * (e * e + G * mp * me));

     Vector2d soln = startMat.colPivHouseholderQr().solve(iC);
     double c1List[q+1] = { soln[0], 0, 0, 0, 0, 0 };
     double c2List[q+1] = { soln[1], 0, 0, 0, 0, 0 };

     double** surfaceList; 
     surfaceList = new double* [gridsize];
     for (int i = 0; i < gridsize; i++)
         surfaceList[i] = new double[3];
     
     for (int i = 0; i < gridsize; i++) {
         surfaceList[i][0] = thetaGrid[i];
         surfaceList[i][1] = FindpSurface(c1List, c2List, thetaGrid[i]);
         surfaceList[i][2] = FindeSurface(c1List, c2List, thetaGrid[i]);
     }

     char atmosphere;
    //check to see if the surfaces intersect
     for (int i = 0; i < gridsize; i++) {
         if (surfaceList[i][0] == 0) {
             if (surfaceList[i][1] > surfaceList[i][2]) {
                 atmosphere = 'P';
             }
             else {
                 atmosphere = 'E';
             }
         }
         if (surfaceList[i][1] < surfaceList[i][2] && atmosphere == 'P') {
             cout << "Intersecting Surfaces" << endl;
         }
         if (surfaceList[i][1] > surfaceList[i][2] && atmosphere == 'E') {
             cout << "Intersecting Surfaces" << endl;
         }
     }

    //here we should solve for what our initial conditions should be for a nonrotating configurations.This is sort of intricate and depends on
    //what sort of atmosphere there is. If it is an electronic atmosphere, we have to start with the protons, bList and v. b0 is solved for with the
    // gradient condition : b0 = -dpsidr * surfaceList[0][1] * *2 where dpsidr is computed as in gradpI.Then v can be solved for with v = PsipI(c1List, c2List, surfaceList[0][1], surfaceList[0][0]) - b0 / surfaceList[0][1]
    // Then we can move on to computing the electrons, aListand phi.We solve for a0 with the gradient condition : a0 = -dpsidr * surfaceList[0][2] * *2 where dpsidr is computed as in gradeA.
    //the we can compute phi using phi = PsieA(c1List, c2List, bList, v, surfaceList[0][2], surfaceList[0][0]) - a0 / surfaceList[0][2]

    // If the atmosphere is protonic, we have to start with the electrons.a0 is solved using a0 = -dpsidr * surfaceList[0][2] * *2 where dpsidr is computed as in gradeI.
    // Then phi can be solved for using phi = PsieI(c1List, c2List, surfaceList[0][2], surfaceList[0][0]) - a0 / surfaceList[0][2]
    // b0 can be solved for using - dpsidr * surfaceList[0][1] * *2 where dpsidr is computed as in gradpA.
    // v is solved for with v = PsipA(c1List, c2List, aList, phi, surfaceList[0][1], surfaceList[0][0]) - b0 / surfaceList[0][1]

     double bList[q + 1] = { 1.24594e+06, 0, 0, 0, 0, 0 };
     double v = -57631.1;
     double  aList[q+1] = { 741059, 0, 0, 0, 0, 0 };
     double   phi = -34483;

     int iter = 0;
     ofstream logfile;

     logfile.open(surfacestr);
     for (int j = 0;j < gridsize; j++) {
         logfile << surfaceList[j][0] << " "<< surfaceList[j][1]<< " "<< surfaceList[j][2]<< "\n";
     }
     logfile.close();
     logfile.clear();

     logfile.open(outputstr);
     for (int j = 0;j < gridsize; j++) {
         logfile << sin(surfaceList[j][0])*surfaceList[j][2] << " " << cos(surfaceList[j][0]) * surfaceList[j][2]<<" " << sin(surfaceList[j][0]) * surfaceList[j][1] << " " << cos(surfaceList[j][0]) * surfaceList[j][1] << "\n";
     }
     logfile.close();
     logfile.clear();

     

     VectorXd fList(4*q+4);
     
     //add the fs
     for (int i = 1;i < q + 2; i++) {
         fList(i-1)=-Rombergf(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
     }

    //add the gs
     for (int i = 1;i < q + 2; i++) {
         fList(q+1+i-1)=-Rombergg(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
    }

    //add the hs
     for (int i = 1;i < q + 2; i++) {
         fList(2*q+2+i-1)=-Rombergh(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
     }

    //add the ks
     for (int i = 1;i < q + 2; i++) {
         fList(3 * q + 3+i-1) = -Rombergk(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
     }


     for (int i = 0; i < 4 * q + 4; i++) {
         cout << fList[i] << " ";
     }
     cout << "\n";
     cout << "\n";

     double maxf = 0;
     for (int i = 0; i < 4 * q + 4; i++) {
         if (abs(fList[i]) > maxf){
             maxf = abs(fList[i]);
         }
     }
     if (maxf < 0.01) {
         cout << "success!"<<endl;
     }

     iter = 1;
     char charcount;
     int flag2;

     VectorXd diff(4 * q + 4);

     double s = 0.01;
     MatrixXd minverse;
     MatrixXd Jacobian(4 * q + 4, 4 * q + 4);
     double hplusList[q];
     double hminusList[q];
     double hplusList2[q + 1];
     double hminusList2[q + 1];
     double h;

     double hplusplusList[q];
     double hminusminusList[q];

     /*
     for (int j = 1; j < q + 1; j++) {
         for (int k = 0; k < q + 1; k++) {
             if (k == j) {
                 h = c2List[k] * s;
                 if (abs(h)< 0.00000001) {
                     h = 0.0000001;
                 }
                 hplusList[k] = c2List[k] + h;
                 hminusList[k] = c2List[k] - h;
                 hplusplusList[k] = c2List[k] + 2*h;
                 hminusminusList[k] = c2List[k] - 2*h;

             }
             else {
                 hplusList[k] = c2List[k];
                 hminusList[k] = c2List[k];
                 hplusplusList[k] = c2List[k];
                 hminusminusList[k] = c2List[k];
             }
         }
         cout << (Rombergf(3, surfaceList, c1List, hplusList, bList, phi, aList, v, atmosphere) - Rombergf(3, surfaceList, c1List, hminusList, bList, phi, aList, v, atmosphere)) / (2 * h) << endl;
         cout << (8*Rombergf(3, surfaceList, c1List, hplusList, bList, phi, aList, v, atmosphere) - 8*Rombergf(3, surfaceList, c1List, hminusList, bList, phi, aList, v, atmosphere)- Rombergf(3, surfaceList, c1List, hplusplusList, bList, phi, aList, v, atmosphere)+ Rombergf(3, surfaceList, c1List, hminusminusList, bList, phi, aList, v, atmosphere)) / (12 * h) << endl;
     }*/
     
     
     while (maxf > 0.01 && iter < 1) {

        //now we compute G
        //creating the jacobian is a complicated grid.the first coordinate should be q + 1 entries with Rombergf, then q + 1 entries with Rombergg
        //then q + 1 entries with rombergh, then finally q + 1 entries with Rombergk.So the numbering should be(0, q), (q + 1, 2q + 1), (2q + 2, 3q + 3), (3q + 4, 4q + 3)
        // in the across direction, we need first q entries for c1List, then q entries for c2List, then q + 1 entries for bList, then phi, then q + 1 entries
        // for aList, then v.
         for (int i = 0; i < q + 1; i++) {
             //note we follow williams here on the derivative estimates, but it is unclear if these are quite right.It could be used to avoid subtraction of small numbers.
             //recall that c1[0] and c2[0] are fixed by the initial conditions
             for (int j = 1; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = c1List[k] * s;
                         if (abs(h) < 0.0000001) {
                             h = 0.01;
                         }
                         hplusList[k] = c1List[k] + h;
                         hminusList[k] = c1List[k] - h;
                     }
                     else {
                         hplusList[k] = c1List[k];
                         hminusList[k] = c1List[k];
                     }
                 }
                 Jacobian(i,j-1) = (Rombergf(i + 1, surfaceList, hplusList, c2List, bList, phi, aList, v, atmosphere) - Rombergf(i + 1, surfaceList, hminusList, c2List, bList, phi, aList, v, atmosphere)) / (2 * h);
             }

             for (int j = 1; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = c2List[k] * s;
                         if (abs(h) < 0.00000001) {
                             h = 0.0000001;
                         }
                         hplusList[k]=c2List[k] + h;
                         hminusList[k]=c2List[k] - h;
                             
                     }
                     else {
                         hplusList[k]=c2List[k];
                         hminusList[k]=c2List[k];
                     }
                 }
                 Jacobian(i,q+j-1) = (Rombergf(i + 1, surfaceList, c1List, hplusList, bList, phi, aList, v, atmosphere) - Rombergf(i + 1, surfaceList, c1List, hminusList, bList, phi, aList, v, atmosphere)) / (2 * h);
             }

             for (int j = 0; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = bList[k] * s;
                         if (abs(h) < 1) {
                             h = pow(10,k);
                         }
                         hplusList2[k]=bList[k] + h;
                         hminusList2[k]=bList[k] - h; 
                     }
                     else {
                         hplusList2[k]=bList[k];
                         hminusList2[k]=bList[k];
                     }
                 }
                 Jacobian(i,2*q+j) = (Rombergf(i + 1, surfaceList, c1List, c2List, hplusList2, phi, aList, v, atmosphere) - Rombergf(i + 1, surfaceList, c1List, c2List, hminusList2, phi, aList, v, atmosphere)) / (2 * h);
             }
 
             Jacobian(i,3*q+1) = (Rombergf(i + 1, surfaceList, c1List, c2List, bList, phi * (1 + s), aList, v, atmosphere) - Rombergf(i + 1, surfaceList, c1List, c2List, bList, phi * (1 - s), aList, v, atmosphere)) / (2 * phi * s);

             for (int j = 0; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = aList[k] * s;
                        if (abs(h) < 1) {
                            h = pow(10,k);
                        }
                         hplusList2[k]=aList[k] + h;
                         hminusList2[k]=aList[k] - h;
                     }
                     else {
                         hplusList2[k]=aList[k];
                         hminusList2[k]=aList[k];
                     }
                 }
                 Jacobian(i,3*q+2+j) = (Rombergf(i + 1, surfaceList, c1List, c2List, bList, phi, hplusList2, v, atmosphere) - Rombergf(i + 1, surfaceList, c1List, c2List, bList, phi, hminusList2, v, atmosphere)) / (2 * h);
             }
             Jacobian(i, 4 * q + 3) = (Rombergf(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 + s), atmosphere) - Rombergf(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 - s), atmosphere)) / (2 * v * s);
         }
         
        //rombergg block
         for (int i = 0; i < q + 1; i++) {
             //note we follow williams here on the derivative estimates, but it is unclear if these are quite right.It could be used to avoid subtraction of small numbers.
             //recall that c1[0] and c2[0] are fixed by the initial conditions
             for (int j = 1; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = c1List[k] * s;
                         if (abs(h) < 0.00000001) {
                             h = 0.01;
                         }
                         hplusList[k] = c1List[k] + h;
                         hminusList[k] = c1List[k] - h;
                     }
                     else {
                         hplusList[k] = c1List[k];
                         hminusList[k] = c1List[k];
                     }
                 }
                 Jacobian(i+q+1, j - 1) = (Rombergg(i + 1, surfaceList, hplusList, c2List, bList, phi, aList, v, atmosphere) - Rombergg(i + 1, surfaceList, hminusList, c2List, bList, phi, aList, v, atmosphere)) / (2 * h);
             }
             for (int j = 1; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = c2List[k] * s;
                         if (abs(h) < 0.00000001) {
                             h = 0.0000001;
                         }
                         hplusList[k] = c2List[k] + h;
                         hminusList[k] = c2List[k] - h;

                     }
                     else {
                         hplusList[k] = c2List[k];
                         hminusList[k] = c2List[k];
                     }
                 }
                 Jacobian(i+q+1, q + j - 1) = (Rombergg(i + 1, surfaceList, c1List, hplusList, bList, phi, aList, v, atmosphere) - Rombergg(i + 1, surfaceList, c1List, hminusList, bList, phi, aList, v, atmosphere)) / (2 * h);
             }
             for (int j = 0; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = bList[k] * s;
                         if (abs(h) < 1) {
                             h = pow(10, k);
                         }
                         hplusList2[k] = bList[k] + h;
                         hminusList2[k] = bList[k] - h;
                     }
                     else {
                         hplusList2[k] = bList[k];
                         hminusList2[k] = bList[k];
                     }
                 }
                 Jacobian(i+q+1, 2 * q + j) = (Rombergg(i + 1, surfaceList, c1List, c2List, hplusList2, phi, aList, v, atmosphere) - Rombergg(i + 1, surfaceList, c1List, c2List, hminusList2, phi, aList, v, atmosphere)) / (2 * h);
             }
             Jacobian(i+q+1, 3 * q + 1) = (Rombergg(i + 1, surfaceList, c1List, c2List, bList, phi * (1 + s), aList, v, atmosphere) - Rombergg(i + 1, surfaceList, c1List, c2List, bList, phi * (1 - s), aList, v, atmosphere)) / (2 * phi * s);

             for (int j = 0; j < q + 1; j++) {
                 for (int k = 0; k < q + 1; k++) {
                     if (k == j) {
                         h = aList[k] * s;
                         if (abs(h) < 1) {
                             h = pow(10, k);

                         }
                         hplusList2[k] = aList[k] + h;
                         hminusList2[k] = aList[k] - h;
                     }
                     else {
                         hplusList2[k] = aList[k];
                         hminusList2[k] = aList[k];
                     }
                 }
                 Jacobian(i+q+1, 3 * q + 2 + j) = (Rombergg(i + 1, surfaceList, c1List, c2List, bList, phi, hplusList2, v, atmosphere) - Rombergg(i + 1, surfaceList, c1List, c2List, bList, phi, hminusList2, v, atmosphere)) / (2 * h);
             }

             Jacobian(i+q+1, 4 * q + 3) = (Rombergg(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 + s), atmosphere) - Rombergg(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 - s), atmosphere)) / (2 * v * s);
         }
             // rombergh block
        for (int i = 0; i < q + 1; i++) {
            //note we follow williams here on the derivative estimates, but it is unclear if these are quite right.It could be used to avoid subtraction of small numbers.
            //recall that c1[0] and c2[0] are fixed by the initial conditions
            for (int j = 1; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = c1List[k] * s;
                        if (abs(h) < 0.001) {
                            h = 0.01;
                        }
                        hplusList[k] = c1List[k] + h;
                        hminusList[k] = c1List[k] - h;
                    }
                    else {
                        hplusList[k] = c1List[k];
                        hminusList[k] = c1List[k];
                    }
                }
                Jacobian(i + 2*q + 2, j - 1) = (Rombergh(i + 1, surfaceList, hplusList, c2List, bList, phi, aList, v, atmosphere) - Rombergh(i + 1, surfaceList, hminusList, c2List, bList, phi, aList, v, atmosphere)) / (2 * h);
            }
            for (int j = 1; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = c2List[k] * s;
                        if (abs(h) < 0.00000001) {
                            h = 0.0000001;
                        }
                        hplusList[k] = c2List[k] + h;
                        hminusList[k] = c2List[k] - h;

                    }
                    else {
                        hplusList[k] = c2List[k];
                        hminusList[k] = c2List[k];
                    }
                }
                Jacobian(i + 2*q + 2, q + j - 1) = (Rombergh(i + 1, surfaceList, c1List, hplusList, bList, phi, aList, v, atmosphere) - Rombergh(i + 1, surfaceList, c1List, hminusList, bList, phi, aList, v, atmosphere)) / (2 * h);
            }
            for (int j = 0; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = bList[k] * s;
                        if (abs(h) < 1) {
                            h = pow(10, k);
                        }
                        hplusList2[k] = bList[k] + h;
                        hminusList2[k] = bList[k] - h;
                    }
                    else {
                        hplusList2[k] = bList[k];
                        hminusList2[k] = bList[k];
                    }
                }
                Jacobian(i + 2*q + 2, 2 * q + j) = (Rombergh(i + 1, surfaceList, c1List, c2List, hplusList2, phi, aList, v, atmosphere) - Rombergh(i + 1, surfaceList, c1List, c2List, hminusList2, phi, aList, v, atmosphere)) / (2 * h);
            }
            Jacobian(i + 2*q + 2, 3 * q + 1)= (Rombergh(i + 1, surfaceList, c1List, c2List, bList, phi * (1 + s), aList, v, atmosphere) - Rombergh(i + 1, surfaceList, c1List, c2List, bList, phi * (1 - s), aList, v, atmosphere)) / (2 * phi * s);

            for (int j = 0; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = aList[k] * s;
                        if (abs(h) < 1) {
                            h = pow(10,k);
                        }
                        hplusList2[k] = aList[k] + h;
                        hminusList2[k] = aList[k] - h;
                    }
                    else {
                        hplusList2[k] = aList[k];
                        hminusList2[k] = aList[k];
                    }
                }
                Jacobian(i + 2*q + 2, 3 * q + 2 + j) = (Rombergh(i + 1, surfaceList, c1List, c2List, bList, phi, hplusList2, v, atmosphere) - Rombergh(i + 1, surfaceList, c1List, c2List, bList, phi, hminusList2, v, atmosphere)) / (2 * h);
            }

            Jacobian(i + 2*q + 2, 4 * q + 3) = (Rombergh(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 + s), atmosphere) - Rombergh(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 - s), atmosphere)) / (2 * v * s);
        }

             //rombergk block
        for (int i = 0; i < q + 1; i++) {
            //note we follow williams here on the derivative estimates, but it is unclear if these are quite right.It could be used to avoid subtraction of small numbers.
            //recall that c1[0] and c2[0] are fixed by the initial conditions
            for (int j = 1; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = c1List[k] * s;
                        if (abs(h) < 0.0001) {
                            h = 0.01;
                        }
                        hplusList[k] = c1List[k] + h;
                        hminusList[k] = c1List[k] - h;
                    }
                    else {
                        hplusList[k] = c1List[k];
                        hminusList[k] = c1List[k];
                    }
                }
                Jacobian(i + 3*q + 3,j - 1) = (Rombergk(i + 1, surfaceList, hplusList, c2List, bList, phi, aList, v, atmosphere) - Rombergk(i + 1, surfaceList, hminusList, c2List, bList, phi, aList, v, atmosphere)) / (2 * h);
            }
            for (int j = 1; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = c2List[k] * s;
                        if (abs(h) < 0.00000001) {
                            h = 0.0000001;
                        }
                        hplusList[k] = c2List[k] + h;
                        hminusList[k] = c2List[k] - h;

                    }
                    else {
                        hplusList[k] = c2List[k];
                        hminusList[k] = c2List[k];
                    }
                }
                Jacobian(i + 3*q + 3,q + j - 1) = (Rombergk(i + 1, surfaceList, c1List, hplusList, bList, phi, aList, v, atmosphere) - Rombergk(i + 1, surfaceList, c1List, hminusList, bList, phi, aList, v, atmosphere)) / (2 * h);
            }
            for (int j = 0; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = bList[k] * s;
                        if (abs(h) < 1) {
                            h = pow(10, k);
                        }
                        hplusList2[k] = bList[k] + h;
                        hminusList2[k] = bList[k] - h;
                    }
                    else {
                        hplusList2[k] = bList[k];
                        hminusList2[k] = bList[k];
                    }
                }
                Jacobian(i + 3*q + 3 ,2 * q + j) = (Rombergk(i + 1, surfaceList, c1List, c2List, hplusList2, phi, aList, v, atmosphere) - Rombergk(i + 1, surfaceList, c1List, c2List, hminusList2, phi, aList, v, atmosphere)) / (2 * h);
            }
            Jacobian(i + 3*q + 3,3 * q + 1) = (Rombergk(i + 1, surfaceList, c1List, c2List, bList, phi * (1 + s), aList, v, atmosphere) - Rombergk(i + 1, surfaceList, c1List, c2List, bList, phi * (1 - s), aList, v, atmosphere)) / (2 * phi * s);

            for (int j = 0; j < q + 1; j++) {
                for (int k = 0; k < q + 1; k++) {
                    if (k == j) {
                        h = aList[k] * s;
                        if (abs(h) < 1) {
                            h = pow(10,k);
                        }
                        hplusList2[k] = aList[k] + h;
                        hminusList2[k] = aList[k] - h;
                    }
                    else {
                        hplusList2[k] = aList[k];
                        hminusList2[k] = aList[k];
                    }
                }
                Jacobian(i + 3*q + 3, 3 * q + 2 + j) = (Rombergk(i + 1, surfaceList, c1List, c2List, bList, phi, hplusList2, v, atmosphere) - Rombergk(i + 1, surfaceList, c1List, c2List, bList, phi, hminusList2, v, atmosphere)) / (2 * h);
            }

            Jacobian(i + 3*q + 3,4 * q + 3) = (Rombergk(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 + s), atmosphere) - Rombergk(i + 1, surfaceList, c1List, c2List, bList, phi, aList, v * (1 - s), atmosphere)) / (2 * v * s);
        }
        //now solve the linear system Jdx = -f
        /*
        for (int i = 0; i < 4 * q + 4; i++) {
            for (int j = 0; j < 4 * q + 4; j++) {
                cout << Jacobian(i, j) << " ";
            }
            cout << endl;
       }*/

        diff = Jacobian.lu().solve(fList);
        for (int i = 0; i < 4 * q + 4;i++) {
            cout << diff[i] << " ";
        }
        cout << "\n\n";
 

        //update guess
        for (int i = 1; i < q + 1; i++) {
            c1List[i] += diff(i - 1);
            c2List[i] += diff(q + i-1 );
        }
        for (int i = 0; i < q + 1; i++) {
            bList[i] += diff(2*q+i);
            aList[i] += diff(3*q+2+i);
        }
        phi += diff(3 * q + 1);
        v += diff(4 * q + 3);

        for (int i = 0; i < q + 1; i++) {
            cout << c1List[i] << " ";
        }
        cout << endl;

        for (int i = 0; i < q + 1; i++) {
            cout << c2List[i] << " ";
        }
        cout << endl;

        for (int i = 0; i < q + 1; i++) {
            cout << bList[i] << " ";
        }
        cout << endl;
        cout << v << endl;

        for (int i = 0; i < q + 1; i++) {
            cout << aList[i] << " ";
        }
        cout << endl;
        cout << phi << endl;

        //compute new surfaces
        for (int i = 0; i < gridsize; i++) {
            surfaceList[i][0] = thetaGrid[i];
            surfaceList[i][1] = FindpSurface(c1List, c2List, thetaGrid[i]);
            surfaceList[i][2] = FindeSurface(c1List, c2List, thetaGrid[i]);
        }

        flag2 = 0;
        // check to see if surfaces intersect
        for (int i = 0; i < gridsize; i++) {
            if (surfaceList[i][0] == 0) {
                if (surfaceList[i][1] > surfaceList[i][2]) {
                    atmosphere = 'P';
                }
                else {
                    atmosphere = 'E';
                }
            }
            if (surfaceList[i][1] < surfaceList[i][2] && atmosphere == 'P') {
                cout << "Intersecting Surfaces" << endl;
                flag2 = 1;
                break;
            }
            if (surfaceList[i][1] > surfaceList[i][2] && atmosphere == 'E') {
                cout << "Intersecting Surfaces" << endl;
                flag2 = 1;
                break;
            }
        }

        if (flag2 == 1) {
            break;
        }

        logfile.open(surfacestr);
        for (int j = 0;j < gridsize; j++) {
            logfile << surfaceList[0][j] << " " << surfaceList[1][j] << " " << surfaceList[2][j] << "\n";
        }
        logfile.close();
        logfile.clear();

        logfile.open(outputstr);
        for (int j = 0;j < gridsize; j++) {
            logfile << sin(surfaceList[j][0]) * surfaceList[j][2] << " " << cos(surfaceList[j][0]) * surfaceList[j][2] << " " << sin(surfaceList[j][0]) * surfaceList[j][1] << " " << cos(surfaceList[j][0]) * surfaceList[j][1] << "\n";
        }
        logfile.close();
        logfile.clear();

        //the(v, bList, cList) = x in Williams' paper. In the iteration, the change in x is given by solving G(dx)=-f
        // Here I call(f, g, h, k) fList.
        // f is always proton, g is always electron, h is always proton derivatives, k is always electron derivatives
        // we need to use different sets of equations depending on what type of atmosphere we have

         //add the fs
        for (int i = 1;i < q + 2; i++) {
            fList(i-1) = -Rombergf(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
        }
        //add the gs
        for (int i = 1;i < q + 2; i++) {
            fList(q + 1 + i-1) = -Rombergg(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
        }
        //add the hs
        for (int i = 1;i < q + 2; i++) {
            fList(2 * q + 2 + i-1) = -Rombergh(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
        }
        //add the ks
        for (int i = 1;i < q + 2; i++) {
            fList(3 * q + 3+i-1) = -Rombergk(i, surfaceList, c1List, c2List, bList, phi, aList, v, atmosphere);
        }

        for (int i = 0; i < 4 * q + 4; i++) {
            cout << fList[i] << " ";
        }
        cout << "\n";

        double maxf = 0;
        for (int i = 0; i < 4 * q + 4; i++) {
            if (abs(fList[i]) > maxf) {
                maxf = abs(fList[i]);
            }
        }
        if (maxf < 0.01) {
            cout << "success!" << endl;
            flag2 = 1;
            logfile.open("coefficients.txt");
            logfile << "c1List={soln[0], " << c1List[1] << ", " << c1List[2] << ", " << c1List[3] << ", " << c1List[4] << ", " << c1List[5]  << "}\n";
            logfile << "c2List={soln[1], " << c2List[1] << ", " << c2List[2] << ", " << c2List[3] << ", " << c2List[4] << ", " << c2List[5]  << "}\n";
            logfile << "b1List={"<<bList[0]<< ", " << bList[1] << ", " << bList[2] << ", " << bList[3] << ", " << bList[4] << ", " << bList[5]  << "}\n";
            logfile << "v=" << v << "\n";
            logfile << "aList={" << aList[0] << ", " << aList[1] << ", " << aList[2] << ", " << aList[3] << ", " << aList[4] <<", "<<aList[5]  << "}\n";
            logfile << "phi=" << phi << "\n";
        }
        if (flag2 == 1) {
            break;
        }

        iter += 1;
        charcount = '0' + iter;
        surfacestr[9] = charcount;
        outputstr[9] = charcount;
     }

}