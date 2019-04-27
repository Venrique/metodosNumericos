#include <iostream>
#include <cmath>
#include <iomanip>
#include <complex>
using namespace std;

inline double funcion(double x){
    return 1+pow(sin(x),2);
}

inline double df(double x){
    return (x+1)*exp(x);
}

inline complex<double> funcionC(complex<double > x){
    return pow(x,3)-(complex<double>(2,0)*x)+complex<double>(2,0);
}

inline complex<double> dfC(complex<double> x){
    return (complex<double>(3,0)*pow(x,2))-complex<double>(2,0);
}

inline double p_hat(double p0, double p1, double p2){
    return p0 - pow( p1-p0 , 2.0)/(p2 - 2.0*p1 + p0);
}

void bisection(double a, double b, double nmax){
    //nmax log2(b-a)/tol
    double p;
    cout << "n | a | b | p | f(a) | f(b) | f(p) \n"<<endl;
    for (int i = 0; i < nmax; ++i) {
        p = (a+b)/2;

        cout << i << "\t" << a << "\t" << b << "\t" << p << "\t" << funcion(a) << "\t" << funcion(b)<< "\t" << funcion(p)<<"\n"<< endl;
        if(funcion(a)*funcion(p)<0){
            b = p;
        }else{
            a = p;
        }
    }
    cout << setprecision(10) << "El resultado es " << p << endl;
}

void pfijo(double p0, int Nmax, double T){
    //nmax log(t/max{b-p0,p0-a})/logk
    double p;

    for(int i=0; i< Nmax; i++){
        p = funcion(p0);
        cout << setprecision(10) <<  i << "\t" << p0 << "\t" << p << "\t" <<abs(p-p0) <<endl;

        if(abs(p-p0) < T) {
            cout << "Objetivo logrado con " << i+1 << " iteraciones"<< endl;
            break;
        }
        p0 = p;
    }
}

void newton(double p0,double TOL,double Nmax){

    double p;
    for(int i=0; i < Nmax; i++){
        p = p0 - funcion(p0)/df(p0);
        cout << i << setprecision(20) << "\t" << p0 << 	"\t" << p << "\t" << abs(p0 - p) << endl;
        if(abs(p - p0) < TOL){
            cout << "Objetivo logrado con " << i+1 << " iteraciones"<< endl;
            break;
        }
        p0 = p;
    }
}

void newtonComplex(complex<double> p0,double TOL,double Nmax){

    complex<double> p;
    for(int i=0; i < Nmax; i++){
        p = p0 - funcionC(p0)/dfC(p0);
        cout << i << setprecision(20) << "\t" << p0 << 	"\t" << p << "\t" << abs(p0 - p) << endl;
        if(abs(p - p0) < TOL){
            cout << "Objetivo logrado con " << i+1 << " iteraciones"<< endl;
            break;
        }
        p0 = p;
    }
}

void secant(long double p0,long double p1,long double TOL,int Nmax){

    long double p;
    for(int i=0; i < Nmax; i++){
        p = p1 - funcion(p1)*(p1 - p0 )/( funcion(p1) - funcion(p0));
        cout << i << setprecision(20) << "\t" << p0 << "\t" << p1 << "\t" << p << "\t"<< abs(p1 - p)<< endl;
        if(abs(p - p1) < TOL)	break;

        p0 = p1;
        p1 = p;

    }

}

void aitken(double p0, int Nmax, double T){

    double p0_hat; //ultimo valor en la sucesion de Aitken
    double p1_hat; //valor actual en la sucesion de Aitken

    double p1, p2; //valores en la iteracion de punto fijo

    for(int i = 0; i <= Nmax; i++){
        p1 = funcion(p0);
        p2 = funcion(p1);

        p1_hat = p_hat(p0, p1, p2);

        if(i == 0) cout << i <<setprecision(10) << "\t" << p0 << "\t" << p1 << "\t" << p2 << "\t" << p1_hat << endl;
        else cout <<i << setprecision(10) << "\t" << p0 << "\t" << p1 << "\t" << p2 << "\t" << p1_hat << "\t" << abs(p1_hat - p0_hat) <<endl;


        if(abs(p0_hat - p1_hat) < T){
            cout << "Objetivo logrado con " << i+1 << " iteraciones"<< endl;
            break;
        }

        p0 = p1;
        p0_hat = p1_hat;

    }
}

void steffensen(long double p0, int Nmax, long double T){

    long double p, p1, p2;

    for(int i = 0; i <= Nmax; i++){

        p1 = funcion(p0);
        p2 = funcion(p1);
        p = p_hat(p0, p1, p2);

        if(i == 0) cout <<i << setprecision(10) << "\t" << p0 << "\t" << p1 << "\t" << p2 << "\t"  << p << endl;
        else cout <<  i<<setprecision(10) << "\t" << p0 << "\t" << p1 << "\t" << p2 << "\t"  << p << "\t"  << abs(p0 - p) << endl;


        if(abs(p0-p) < T){
            cout << "Objetivo logrado con " << i+1 << " iteraciones"<< endl;
            break;
        }

        p0 = p;
    }
}

int main() {
    //bisection(0,1,20);
    //pfijo(1.5,5500,pow(10,-9));
    //newtonComplex(complex<double>(1,1),pow(10,-10),500);
    //secant(0.65625,0.75625,pow(10,-12),500);
    //aitken(1.5,500,pow(10,-9));
    steffensen(1.5,500,pow(10,-9));
    return 0;
}