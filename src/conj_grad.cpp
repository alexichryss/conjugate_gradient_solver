#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>

using namespace std;
typedef long double D_TYPE;

int Q = 0; // "source term is constant and zero"
int N = 100; // number of grid points in each dimension
int NT = 200; // number of timesteps
int L = 4; // length of the domain
int T = 1; // total time
int a = 1; // thermal diffusivity of the medium (assumed constant)

int MAX_ITER = N*N; // max iterations

// helper function to print a matrix
void print_matrix(vector<vector<D_TYPE>> M) {
    size_t n = M.size();
    string filename = "CG data " + to_string(a) + ".txt";
    
    ofstream output(filename, ios::app);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            output << M[i][j] << " ";
        }
        output << "\n";
    }
    output << "\n";
    output.close();
}

// returns the input matrix multiplied by the coefficient matrix A of Crank-Nicolson xnew = Ax
vector<vector<D_TYPE>> matvec_A(vector<vector<D_TYPE>> x, D_TYPE C) {
    size_t n = x.size();
    vector<vector<D_TYPE>> xold = x;
    vector<vector<D_TYPE>> xnew(n, vector<D_TYPE>(n,0));
    D_TYPE i1, i_1, j1, j_1; // holder for values of i+1, i-1, j+1, j-1 and their boundary conditions

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            // boundary conditions
            i1  = i == (n-1) ? 0 : x[i+1][j];
            i_1 = i == 0     ? 0 : x[i-1][j];
            j1  = j == (n-1) ? 0 : x[i][j+1];
            j_1 = j == 0     ? 0 : x[i][j-1];

            xnew[i][j] = (1+4*C)*x[i][j] - C*(i1 + i_1 + j1 + j_1) + Q; 
        }
    }
    return xnew;
}

// returns the input matrix multiplied by the coefficient matrix A of Crank-Nicolson xnew = Ax
vector<vector<D_TYPE>> matvec_B(vector<vector<D_TYPE>> x, D_TYPE C) {
    size_t n = x.size();
    vector<vector<D_TYPE>> xold = x;
    vector<vector<D_TYPE>> xnew(n, vector<D_TYPE>(n,0));
    D_TYPE i1, i_1, j1, j_1; // holder for values of i+1, i-1, j+1, j-1 and their boundary conditions

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            // boundary conditions
            i1  = i == (n-1) ? 0 : x[i+1][j];
            i_1 = i == 0     ? 0 : x[i-1][j];
            j1  = j == (n-1) ? 0 : x[i][j+1];
            j_1 = j == 0     ? 0 : x[i][j-1];

            xnew[i][j] = (1-4*C)*x[i][j] + C*(i1 + i_1 + j1 + j_1) + Q;
        }
    }
    return xnew;
}

// dot product function that returns a scalar a = b · c
D_TYPE dotp(vector<vector<D_TYPE>> b, vector<vector<D_TYPE>> c) {
    D_TYPE a = 0.0;

    for (int i=0; i<b.size(); i++) {
        for (int j=0; j<b.size(); j++) {
            a += b[i][j]*c[i][j];
        }
    }
    return a;
}

// vector addition function that scales and adds two vectors to return a vector w = αw + βv
vector<vector<D_TYPE>> axpy(D_TYPE a, vector<vector<D_TYPE>> w, D_TYPE b, vector<vector<D_TYPE>> v) {
    vector<vector<D_TYPE>> temp = w;
    
    for (int i=0; i<w.size(); i++) {
        for (int j=0; j<w.size(); j++) {
            temp[i][j] = a*w[i][j] + b*v[i][j];
        }
    }

    return temp;
}


void solver() {
    // initialize calculated parameters
    D_TYPE sig_x = L/4.0; // gaussian spread in x
    D_TYPE sig_y = L/4.0; // gaussian spread in y
    D_TYPE Sx = 2.0*sig_x*sig_x; // denominator for x
    D_TYPE Sy = 2.0*sig_y*sig_y; // denominator for y
    D_TYPE dN = double(L)/(N-1); // size of change in x,y
    D_TYPE dT = double(T)/NT; // size of change in time
    D_TYPE C = a*dT/(dN*dN); // constant C
    printf("a: %d, dN: %Lf, dT: %Lf\n", a, dN, dT);
    printf("Sx: %Lf, Sy: %Lf, C: %Lf\n", Sx, Sy, C);

    // initialize all needed arrays
    vector<D_TYPE> arr(N);
    vector<vector<D_TYPE>> b(N, vector<D_TYPE>(N,0)), r, p, y;
    vector<vector<D_TYPE>> x(b);

    // linspace for our x and y coordinates
    D_TYPE front = L/-2.0;
    for (int i = 0; i < N; i++) {
        arr[i] = front + dN*i;
    }

    // initialize domain space with heat equation
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b[i][j] = a*exp(-arr[i]*arr[i]/Sx)*exp(-arr[j]*arr[j]/Sy) + Q;
        }
    }

    // initialize variables for CG
    D_TYPE rsold, rsnew, alpha;
    // initial b from b = Bx(o)
    b = matvec_B(b,C);

    clock_t start = clock();
    for (int n=0; n<NT; n++) {
        // conjugate gradient algorithm
        r = axpy(1,b,-1,matvec_A(x,C));
        rsold = dotp(r,r);
        p = r;
        for (int i=0; i<MAX_ITER; i++) {
            y = matvec_A(p,C);
            alpha = rsold / dotp(p,y);
            x = axpy(1, x, alpha, p);
            r = axpy(1, r, -1*alpha, y);
            rsnew = dotp(r,r);
            if (sqrt(rsnew) < 1E-10) {
                break;
            }
            p = axpy(1, r, rsnew/rsold, p);
            rsold = rsnew;
        }
        // get b from b = Bx
        b = matvec_B(x,C);

        if (n==1 || n==NT/3 || n==2*NT/3 || n==NT-1) {
            printf("time taken at %d: %f\n", n, (clock() - start)/double(CLOCKS_PER_SEC));
            // print_matrix(b);
        }
    }
    printf("total time taken: %f\n", (clock() - start)/double(CLOCKS_PER_SEC));
}

int main(int argc, char** argv) {
    switch(argc) {
        // if no parameters, then calls solver with defaults
        case 1 : {
            printf("Using values: N %d, NT %d, L %d, T %d, a %d\n", N, NT, L, T, a);
            solver();
            break; }
        // if demo is given, then iterate through different grid sizes
        case 2 : {
            if (strncmp(argv[1], "demo", string(argv[1]).size()) == 0) {
                vector<int> grid = {10,48,86,124,162,200};
                for (int i=0; i<grid.size(); i++) {
                    N = grid[i];
                    printf("Using values: N %d, NT %d, L %d, T %d, a %d\n", N, NT, L, T, a);
                    solver();
                }
            }
            break; }
        // if 2 parameters, then calls solver with inputs
        case 6 : {
            printf("Processing inputs...\n");
            N = stoi(argv[1]);
            NT = stoi(argv[2]);
            L = stoi(argv[3]);
            T = stoi(argv[4]);
            a = stoi(argv[5]);
            solver();
            break; }
        // otherwise displays usage message
        default : cout << "Usage: %s <learning rate> <convergence criteria> \n";
    }

    return 0;
}
