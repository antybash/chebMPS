//   path to Eigen == /usr/local/include/eigen-eigen-1306d75b4a21
//   path to Eigen == ./eigen334_copy


#ifndef EXPV_H
#define EXPV_H

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#define _EXPV_PI 3.14159265358979323

std::tuple< Eigen::VectorXd, double, double >
    expv
    (
        double t,
        Eigen::MatrixXd A,
        Eigen::VectorXd v,
        double tol,
        int m
    )
{
    int n = A.rows();

    bool NONzeroMatrix = false;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            NONzeroMatrix = ( NONzeroMatrix || ( fabs( A(i,j) ) > 1e-7 ) );
    if (!NONzeroMatrix){
        return std::make_tuple( v, 0, 0 );
    }

    tol   = std::max(tol,1.0e-4);
    m     = std::min(n,30);

    double anorm = 0.0;
    for(int i = 0; i < n; i++) for(int j = 0; j < n; ++j)
        anorm = std::max(  anorm, fabs(A(i,j))  );
    
    double    btol = 1.0e-4;
    double   gamma = 0.9;
    double   delta = 1.2;
    double   t_out = fabs(t);
    double   t_new = 0.0;
    double   t_now = 0.0;
    double s_error = 0.0;
    double     eps = 2.22e-16;
    double  rndoff = anorm * eps;
    int      nstep = 0; 
    int         mb = m;
    int      mxrej = 10;

    double  avnorm;
    int         mx;
    double err_loc;
    double     err;
    double    phi1;
    double    phi2;
    Eigen::MatrixXd F;

    int k1 = 2;
    double xm = 1.0/double(m);
    double normv = v.norm();
    double beta  = normv;
    double fact = 
        ( std::pow( (double(m)+1.0)/std::exp(1.0), double(m)+1.0))
        * std::sqrt(2*_EXPV_PI*(double(m)+1.0));
    t_new = (1.0/anorm)*std::pow((fact*tol)/(4.0*beta*anorm),xm);
    double s = std::pow(10.0,double(std::floor(std::log10(t_new))-1)); 
    t_new = std::ceil(t_new/s)*s; 
    double sgn = copysign(1.0,t);
    nstep = 0;

    Eigen::VectorXd w = v;
    double hump = normv;
    while( t_now < t_out ){
        nstep = nstep + 1;
        double t_step = std::min( t_out-t_now, t_new );
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n,m+1);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m+2,m+2);

        V.col(0) = (1.0/beta)*w;
        for(int j = 0; j < m; ++j){
            Eigen::VectorXd p = A * V.col(j);
            for(int i = 0; i <= j; ++i){
                H(i,j) = p.dot( V.col(i) );
                p = p-H(i,j)*V.col(i);
            }
            double s = p.norm();
            if (s < btol){
                k1 = 0;
                mb = j+1;
                t_step = t_out - t_now;
                break;
            }
            H(j+1,j) = s;
            V.col(j+1) = (1.0/s)*p;
        }

        if (k1 != 0){
            H(m+1,m) = 1;                   //changed index (m+2,m+1) -> (m+1,m)
            Eigen::VectorXd x = A*V.col(m); //changed index     (m+1) -> m
            avnorm = x.norm();
        }

        int ireject = 0;
        while(ireject <= mxrej){
            mx = mb + k1;
            Eigen::MatrixXd preF = sgn*t_step*H.block(0,0,mx,mx);
            F = preF.exp();
            if(k1 == 0){
                err_loc = btol;
                break;
            } else {
                phi1 = fabs( beta*F(m,0) );             //changed index (m+1,1) -> (m,0)
                phi2 = fabs( beta*F(m+1,0) * avnorm );  //changed index (m+2,1) -> (m+1,1)
                if( phi1 > 10*phi2 ){
                    err_loc = phi2;
                    xm = 1.0/double(m);
                } else if (phi1 > phi2) {
                    err_loc = (phi1*phi2)/(phi1-phi2);
                    xm = 1.0/double(m);
                } else {
                    err_loc = phi1;
                    xm = 1.0/double(m-1);
                }
            }
            if( err_loc <= delta * t_step * tol )
                break;
            else {
                t_step = gamma * t_step * std::pow(t_step*tol/err_loc,xm);
                s = pow(10.0, std::floor(std::log10(t_step))-1);
                t_step = std::ceil(t_step/s) * s;
                assert(ireject != mxrej); // Requested tolerance is too high 
                ireject+=1;
            }
        }
        mx = mb + std::max( 0, k1-1 );

        w  = V.block(0,0,V.rows(),mx) * ( beta * F.col(0).head(mx) );
        beta = w.norm();
        hump = std::max(hump,beta);

        t_now = t_now + t_step;
        t_new = gamma * t_step * std::pow(t_step*tol/err_loc,xm);
        s = pow(10.0, std::floor(std::log10(t_new))-1);
        t_new = ceil(t_new/s) * s;

        err_loc  = std::max(err_loc, rndoff);
        s_error += err_loc;
    }
    err   = s_error;
    hump /= normv;
    return std::make_tuple( w, err, hump );
}

std::tuple< Eigen::VectorXd, double, double > expv ( double t, Eigen::MatrixXd A, Eigen::VectorXd v )
{
    return expv( t, A, v, 1.0e-4, 30 );
}

#endif
