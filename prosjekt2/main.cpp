#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;


// declare functions
double Mmax (mat &B, int &k,int &l, int n);
void rotate(mat &B, int k, int l, int n);
vec jacobi_eigenvalues(mat &B, double epsilon, int maxit, int &k, int &l, int n, int &iterations);
int rotate_test();
int Mmax_test();
int jacobi_test();


vec jacobi_eigenvalues(mat &B, double epsilon, int maxit, int &k, int &l, int n, int &iterations) {
    // Rotate the Bmatrix to find eigenvalues
    iterations = 0;
    double Bmax = 10;
    while (abs(Bmax) > epsilon && iterations < maxit) {
        // find k, l max value
        Bmax = Mmax(B,k,l,n);
        rotate(B, k, l, n);
        iterations++;

    }
    vec eigenB = sort(B.diag());
    return eigenB;
}

void rotate(mat &B, int k, int l, int n) {
    double a_kk, a_ll, a_ik, a_il;
    double c, s, t, tau;

    // First, find the value tau = cot(2x), and for t
    if (B(k,l) != 0.0) {
        tau = (B(l,l)-B(k,k))/(2*B(k,l));
        if (tau > 0.0) {
            t = 1.0/(tau+sqrt(1+tau*tau));
        } else {
            t = 1.0/(tau-sqrt(1+tau*tau));
        }
        c = 1/sqrt(1+t*t);
        s = t*c;
    } else {
        c = 1.0;
        s = 0.0;
    }
    a_kk = B(k,k);
    a_ll = B(l,l);
    // changing the matrix elements with indices k and l
    B(k,k) = c*c*a_kk - 2.0*c*s*B(k,l) + s*s*a_ll;
    B(l,l) = s*s*a_kk + 2.0*c*s*B(k,l) + c*c*a_ll;
    B(k,l) = 0.0; // hard-coding of the zeros
    B(l,k) = 0.0;

    // and then we change the remaining elements
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = B(i,k);
            a_il = B(i,l);
            B(i,k) = c*a_ik - s*a_il;
            B(k,i) = B(i,k);
            B(i,l) = c*a_il + s*a_ik;
            B(l,i) = B(i,l);
        }
    }
}

// Find the largest offdiagonal value in an symmetric matrix and return k and l coordinates

double Mmax (mat &B, int &k,int &l, int n){
    int i, j;
    double Bmax=0.;
    for(i = 0 ; i < n ; i++){
        for(j = i+1 ; j < n ; j++) {
            if (abs(Bmax) < abs(B(i,j))) {
                Bmax = B(i,j);
                k = i;
                l = j;
            }
        }
    }
    return Bmax;
}

void UnitTests () {
    cout << "Test results: " << endl;
    int a = Mmax_test();
    int b = rotate_test();
    int c = jacobi_test();
    int tot = a+b+c;
    int number_of_tests = 3;
    cout << "Successful tests: " << tot << " out of " << number_of_tests << " tests performed!" << endl;
}

int rotate_test(){
    int n = 2;
    mat C = ones(n,n);
    vec eigen_test;
    rotate(C,0,1,n);
    eigen_test = sort(C.diag());
    if (C(0,1) == 0.0 && C(1,0) == 0.0 && eigen_test(0) <= 1e-8 && eigen_test(1) - 2.0 <= 1e-8) {
        cout << "Rotate is working as advertised!" << endl;
        return 1;
    } else {
        cout << "Rotate is a failure!" << endl;
        return 0;
    }
}

int jacobi_test (){
    int n = 2;
    int l, k;
    int maxit = 100;
    int iterations;
    double epsilon = 1e-8;
    mat C = ones(n,n);
    vec eigen_test;
    eigen_test = sort(jacobi_eigenvalues(C, epsilon, maxit, k, l, n, iterations));
    if (eigen_test(0) <= epsilon && eigen_test(1) - 2.0 <= epsilon && iterations-1 == 1) {
        cout << "Jacobi_eigenvalues is working as advertised!" << endl;
        return 1;
    } else {
        cout << "Jacobi_eigenvalues is a failure!" << endl;
        return 0;
    }
}

int Mmax_test () {
    int n = 3;
    mat C = zeros<mat>(n,n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C(i,j) = n*i+j+1;
        }
    }
    int k, l;
    double Cmax = Mmax(C, k, l, n);
    if (k == 1 && l == 2 && Cmax == 6){
        cout << "Mmax is working as advertised!" << endl;
        return 1;
    } else {
        cout << "Mmax is a failure!" << endl;
        return 0;
    }
}


int main() {

    //Will solve the equation d^2/(dr^2) u = (r^2-L)u

    //define variables


    //declarations
    int k, l;
    k = 0;
    l = 0;
    int n = 150;
    int maxit = 100000;
    int iterations;
    double rmin = 0.0;
    double rmax = 4.37;
    double h = (rmax-rmin)/(n+1);
    double dh = h*h;
    double epsilon = 1e-8;
    vec a = {0.01, 0.5, 1.0, 5.0};
    vec r = linspace(rmin, rmax, n+2);
    // R matrix holding all solutions
    mat U = zeros<mat>(n,5);

    //constructing the tridiagonal matrix: d is diagonal element, e is off diagonal

    vec d = zeros(n);
    for(int i = 0 ; i<n ; i++) {
        d(i) = r(i+1)*r(i+1) + 2/dh;
    }
    vec e = -1*ones(n)/dh;

    //Create a matrix A, holding the info of the system, and the matrix B for finding solution
    mat A = zeros<mat>(n,n);
    A.diag(+1)+=e(0);
    A.diag(-1)+=e(0);
    for(int i = 0 ; i<n ; i++) {
            A(i,i) = d(i);
    }
    int maxn = 50;
    int minn = 2;

    mat to_file_iterations = zeros(maxn-minn+1,2);
    for (int n = 2; n < maxn+1; n++){

        double h = (rmax-rmin)/(n+1);
        double dh = h*h;
        vec d = zeros(n);
        vec r = linspace(rmin, rmax, n+2);
        for(int i = 0 ; i<n ; i++) {
            d(i) = r(i+1)*r(i+1) + 2/dh;
        }

        vec e = -1*ones(n)/dh;
        mat A = zeros<mat>(n,n);
        A.diag(+1)+=e(0);
        A.diag(-1)+=e(0);
        for(int i = 0 ; i<n ; i++) {
                A(i,i) = d(i);
        }

        vec eigenvalues = jacobi_eigenvalues(A, epsilon, maxit, k, l, n, iterations);
        to_file_iterations(n-minn,0) = n;
        to_file_iterations(n-minn, 1) = iterations-1;
        cout << n << endl;

    }
    to_file_iterations.save("iterations_n.dat", raw_ascii);





    vec eigenA;
    mat eigvecA;
    eig_sym(eigenA, eigvecA, A);
    eigenA = sort(eigenA);
    for(int i = 0; i < n ; i++) {
        U(i,0) = eigenA(i);
    }

    mat B = A;

    vec eigenvalues = jacobi_eigenvalues(B, epsilon, maxit, k, l, n, iterations);

    // For part c, create a new matrix A, containinga new potential V = ar^2+1/r


    mat R1 = zeros(n,n);
    mat Rlist[] = {eigvecA, R1, R1, R1, R1};

    for(int k = 0 ; k < 4 ; k++) {
        vec d2 = zeros(n);
        for(int i = 0 ; i<n ; i++) {
            d2(i) = a(k)*a(k)*r(i+1)*r(i+1) + 1/r(i+1) + 2/dh;
        }

        for(int i = 0 ; i<n ; i++) {
            A(i,i) = d2(i);
        }
        eig_sym(eigenA,Rlist[k+1],A);
        eigenA = sort(eigenA);
        for(int j = 0; j < n ; j++) {
            U(j,k+1) = eigenA(j);
        }
    }

    mat to_file1 = zeros(n+2,5);
    for(int j=0 ; j < n ; j++) {
        to_file1(j,0) = U(j,0);
        to_file1(j,1) = U(j,1);
        to_file1(j,2) = U(j,2);
        to_file1(j,3) = U(j,3);
        to_file1(j,4) = U(j,4);
    }
    to_file1.save("solutions_u.dat", raw_ascii);

    mat to_file = zeros(n+2,5);
    for(int j=0 ; j < n ; j++) {
        to_file(j,0) = Rlist[0](j,0);
        to_file(j,1) = Rlist[1](j,0);
        to_file(j,2) = Rlist[2](j,0);
        to_file(j,3) = Rlist[3](j,0);
        to_file(j,4) = Rlist[4](j,0);
    }
    to_file.save("eigvecs_u.dat", raw_ascii);


    mat C = B;
    vec eigenvalues2 = jacobi_eigenvalues(C, epsilon, maxit, k, l, n, iterations);
    cout << "C: " << eigenvalues2(0) << ",   " << eigenvalues2(1) << ",   " << eigenvalues2(2) << endl;


    cout << "Victory!" << endl;
    UnitTests();
    return 0;

}


