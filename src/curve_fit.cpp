#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
using namespace std;
long double A = 0.0004;
long double E = 5.0E-7;
long double H = 3.0;
long double M = 3.0;
long double B = 3.0;

// takes 2d array of numbers, 
void solver(vector<vector<long double>> points, long double a = A, long double e = E)
{   
    printf("Learning rate: %Lf, Convergence criterion: %Le\n", a, e);
    // initialize coefficients
    long double h = H, m = M, b = B;
    long double Dh = 0, Dm = 0, Db = 0;
    // vectors for data plotting
    vector<int> iterations;
    vector<long double> gradients;
    // variables for loop
    int iter = 0;
    long double norm_grad;
    double n = points.size();

    // do-while loop polynomial curve fitting for gradient descent
    do {

        iter++; // count step
        Dh = 0; Dm = 0; Db = 0; // return gradients to 0

        // loop through each point
        for (int i = 0; i < points.size(); i++)
        {
            long double x = points[i][0], y = points[i][1];
            // calculate new gradient based on new step
            Dh = Dh - 2*y*pow(x,2) + 2*h*pow(x,4) + 2*m*pow(x,3) + 2*b*pow(x,2);
            Dm = Dm - 2*y*x + 2*h*pow(x,3) + 2*m*pow(x,2) + 2*b*x;
            Db = Db - 2*y + 2*h*pow(x,2) + 2*m*x + 2*b;
        }
        // end loop
        Dh = Dh / n; Dm = Dm /n; Db = Db / n; // normalize
        h = h - a*Dh; m = m - a*Dm; b = b - a*Db; // update calculated coefficients
        norm_grad = sqrt(pow(Dh,2) + pow(Dm,2) + pow(Db,2)); // norm of the three gradients
        // record once every 1000 iterations
        if (iter % 1000 == 0) {
            iterations.push_back(iter);
            gradients.push_back(norm_grad);
        }
    } while (norm_grad > e); // end do-while when convergence criterion is met

    // print results to console
    printf("Points: %d, Iterations: %d\nh: %Lf, m: %Lf, b: %Lf\n", int(n), iter, h, m, b);
    // output results to file
    string filename = "data_" + to_string(a) + ".txt";
    ofstream output(filename);
    for (int i = 0; i<iterations.size(); i++) {
        output << iterations[i] << " " << gradients[i] << "\n";
    }
    output.close();
}

int main(int argc, char** argv)
{
    // initialize variables to hold our data
    string lines;
    vector<vector<long double>> points;

    ifstream infile("points.txt"); //read in text file

    // while loop for each line in file
    while (getline(infile, lines))
    {
        points.push_back(vector<long double>()); // add a new vector of long doubles
        istringstream ss(lines); // raw string of chars
        long double value;
        // while loop for each value on line
        while (ss >> value)
        {
            points.back().push_back(value); //add data to the end of the last vector
        }
    }

    points.erase(points.begin()); // remove header row
    
    argv++;
    switch(argc) {
        // if no parameters, then calls solver with defaults
        case 1 : {
            printf("Using values: a = %Lf, e = %Le\n", A, E);
            solver(points);
            break; }
        // if 2 parameters, then calls solver with inputs
        case 3 : {
            printf("Processing inputs...\n");
            long double a = stold(*argv);
            argv++;
            long double e = stold(*argv);
            solver(points, a, e);
            break; }
        // otherwise displays usage message
        default : cout << "Usage: %s <learning rate> <convergence criteria> \n";
    }

}
