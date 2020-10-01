#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cstring>
#include <chrono>
#include <random>
#include <iomanip>
#include <ctime>
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/SparseCore"
 
using namespace Eigen;
using namespace std;

using Mat = Matrix<double,Dynamic,Dynamic>;
using Vec = Matrix<double,Dynamic,1>;

struct Result {
    Vec end;
    Vec start;
    string status;
};

int reducedDimension(int q)
{
    return 1+(q-1)/2+( q%2==0 ? 1 : 0) ;
}

double xlog(double x)
{
    return x==0 ? 0 : x*log(abs(x));
}

Mat fullReducedInteraction(int p, int q)
{
    int q2 = reducedDimension(q);
    Mat J(q2,q2);
    J(0,0)=1;
    for(int i = 1; i<q2; ++i)
    {
        J(0,i) = 4/2.;
        J(i,0) = 4/2.;

        for(int j =i; j<q2; ++j)
        {
            J(i,j) = 4. * cos( 2. * M_PI * p / q * i * j );
            J(j,i) = 4. * cos( 2. * M_PI * p / q * i * j );
        }       
    }

    if(q%2==0)
    {
        J(q2-1,q2-1) = cos(M_PI * p * q / 2.);
        J(0,q2-1) = 2/2.;
        J(q2-1,0) = 2/2.;

        for(int j=1; j<q2-1; ++j)
        {
            J(j,q2-1) = 4. * cos( M_PI * p * j ) / 2.;
            J(q2-1,j) = 4. * cos( M_PI * p * j ) / 2.;
        }
    }

    return J /= q;
}

double gradSingleEntropy(double m)
{
    double regularizer = 1e-3;
    return 0.5 * log( 1-m+regularizer ) - 0.5 * log( 1+m+regularizer );
}

Vec gradientFreeEnergy( Vec vars, int p, int q, double t)
{
    int q2           = reducedDimension(q);
    Mat J = fullReducedInteraction(p, q);

    Matrix<double, Dynamic, 1> gradF(q2);
    gradF = 2*(J*vars);
    gradF(0) += - t * 2 * gradSingleEntropy(vars(0)); 
    for(int i=1; i<q2-1; ++i)
    {
        gradF(i) += - t * 4 * gradSingleEntropy(vars(i)); 
    }
    gradF(q2-1) += - t * (q%2==0 ? 2 : 4) * gradSingleEntropy(vars(q2-1)); 

    return gradF;
}

Result localMinimum_gradientDescent(int p, int q, double t)
{    
    int q2           = reducedDimension(q);

    Vec startingCondition(q2);
    Vec old_position(q2);
    Vec actual_position(q2);

    //see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;
    std::mt19937 gen(rd());  //here you could set the seed, but std::random_device already does that
    std::uniform_real_distribution<float> dis(-1.0, 1.0);
    startingCondition = Vec::NullaryExpr(q2,1,[&](){return dis(gen);});

    old_position = startingCondition + 0.1 * Vec::Ones(q2);
    actual_position = startingCondition;
    
    double threshold = 1e-3;
    int maxIter      = 1000000;
    double step      = 1e-3;
    string status    = "";

    Vec grad = gradientFreeEnergy(actual_position,p,q,t);
    int iter = 0;

    while(
            ( grad.norm() > threshold ) && 
            ( iter < maxIter )          && 
            ( actual_position != old_position )
         )
    {
        iter += 1;
        old_position = actual_position;
        actual_position -= step * grad;
        for(int i=0; i<q2; ++i)
        {
            if(actual_position(i) > 1)
            {
                actual_position(i) = 1;
            } 
            else if (actual_position(i) < -1) 
            {
                actual_position(i) = - 1;
            }

        }
        grad = gradientFreeEnergy(actual_position,p,q,t);
        for(int i=0; i<q2; ++i)
        {
            if(
                    (actual_position(i) == 1 && grad(i) < 0) ||
                    (actual_position(i) ==-1 && grad(i) > 0) 
              )
            {
                grad(i) = 0;
            } 
        }
    }

    if(iter == maxIter)
    {
        status = "max iter reached";
    }
    else if(actual_position == old_position)
    {
        status = "stuck in position";
    }
    else
    {
        status = "good";
    }

    Result r = {actual_position, startingCondition, status};
    return r;
}

string to_string(Vec v) {
    string s = "[";
    for(int i=0; i < v.rows(); ++i )
    {
        s.append(to_string(v(i)));
        if(i != v.rows()-1)
        {
            s.append(", ");
        }
    }
    s.append("]");
    return s;
}

string time() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream stream;
    stream << put_time(&tm, "%Y%m%d%H%M%S");
    return stream.str();
}

int main(int argc, char** argv)
{
    int p;
    int q;
    double t;
    int iter;

    if(argc == 5)
    {
        p    = atoi(argv[1]);
        q    = atoi(argv[2]);
        t    = atof(argv[3]);
        iter = atoi(argv[4]);
    } else {
        cerr << "Usage: ./exe p q t iter" << endl;
        return 1;
    }

    vector<string> v(iter);
    
    for(int i=0; i<iter; ++i)
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        Result r = localMinimum_gradientDescent(p,q,t);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        
        v[i] = "";
        v[i].append(to_string(p));
        v[i].append("|");
        v[i].append(to_string(q));
        v[i].append("|");
        v[i].append(to_string(t));
        v[i].append("|gradientCPP|");
        v[i].append(to_string(r.end));
        v[i].append("|");
        v[i].append(to_string(r.start));
        v[i].append("|");
        v[i].append(time());
        v[i].append("|");
        v[i].append(to_string(duration/1e6));
        v[i].append("|");
        v[i].append(r.status);
        cout << v[i] << endl;
    }

    return 0;
}




