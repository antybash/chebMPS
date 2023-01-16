
#ifndef COMPRESS_SINGLE_STATE_H
#define COMPRESS_SINGLE_STATE_H

#include <uni10.hpp>
#include <numeric>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <functional>
#include <stack>
#include "helper_compress2.h"
#include "tsvd.h"

using namespace uni10;
using namespace useful;

/*!
 * \class compress_SingleState
 * \brief Base class from which all compressions are induced.
 *
 */
class compress_SingleState
{
    /*******************************************************
     *
     * compress_SingleState:
     *      - compresses an MPS
     *
     ******************************************************/
    public: 
        compress_SingleState( std::vector<UniTensor>&, int, int, double, double, bool, int );
        std::vector<UniTensor> compute( bool ); 
        double accuracy;
        double relacc;
    private:
        bool svd_ran;
        int mpsL;
        int maxBondDim;
        int minBondDim;
        int num_sweeps;
        double epsilon_trunc;
        double epsilon_diff;
        double norm_refMPS;
        bool verbose;
        std::vector< UniTensor > refMPS;
        std::stack< UniTensor > L1; 
        std::stack< UniTensor > R1;
        //
        std::stack< UniTensor > L2; 
        std::stack< UniTensor > R2;
        //
        //  void make_initial_svd_guess();
        void make_svd_compress(int);
        void make_random_guess(int);
        void make_helper1(UniTensor&, std::stack<UniTensor>&, std::stack<UniTensor>&);
        //
        void localcompress(bool);
        void localcompress_h1(bool);
        void localcompress_h2();
        void localcompress_h3(bool);
        void localcompress_h4(bool);
        void localcompress_h5();
        UniTensor localcompress_t1;
        UniTensor localcompress_normalizedtensor;
        double    localcompress_trunc_error;
        double    localcompress_norm;
        //
        void globalcompress_half(bool);
        bool globalcompress_full(int);
        std::vector<UniTensor> getMPS();
        //
        void output_current_mps();
        void output_current_mps_just_diagrams();
        void output_current_mps_just_raw_elems();
        //old functions:
        void sweep_once_left_to_right();
        void sweep_once_right_to_left();
        void sweep_all_left_to_right();
        void sweep_all_right_to_left();
        void push_onto_left ( UniTensor, int );
        void push_onto_right( UniTensor, int );
        double norm_diff_target_approx();
};

/************************************************************************************************
 *                       Implementation of   compress_SingleState   class                       *
 ***********************************************************************************************/

/*! 
 * \brief Initializes reference MPS, target bond dimension, and truncation-error-epsilon.
 */
compress_SingleState::compress_SingleState 
    ( 
     std::vector<UniTensor> &mps_ref, 
     int bond_dim, 
     int min_bond_dim,
     double eps_t, 
     double eps_d, 
     bool sr,
     int sweeps
    )
{
    svd_ran = sr; // true == svd; false == random;
    num_sweeps = sweeps;
    mpsL = mps_ref.size();
    maxBondDim = bond_dim;
    minBondDim = min_bond_dim;
    epsilon_trunc = eps_t;
    epsilon_diff  = eps_d;
    refMPS = mps_ref;
    //
    int max_bond_in_refMPS = 0;
    for(int i = 0; i < refMPS.size(); ++i)
        max_bond_in_refMPS = std::max( max_bond_in_refMPS, std::max( refMPS[i].bond()[0].dim(), refMPS[i].bond()[2].dim() ) );
    //
    norm_refMPS = normMPS( refMPS );
    std::cout << " compress_SingleState::compress_SingleState: INITIALIZATION COMPLETE " << std::endl;
    std::cout << " norm of reference MPS is: "     << std::fixed    << std::setprecision(10) << norm_refMPS        << std::endl;
    std::cout << " max bond of reference MPS is: " << std::fixed    << std::setprecision(10) << max_bond_in_refMPS << std::endl;
    std::cout << " target bond dimension is: "     << maxBondDim    << std::endl;
    std::cout << " allowed truncation error: "     << epsilon_trunc << std::endl;
    std::cout << " allowed difference error: "     << epsilon_diff  << std::endl;
}

std::vector< UniTensor >
    compress_SingleState::compute
    ( bool print_out )
{
    verbose = print_out;
    int m;
    // for (int i0 = 0; i0 < 2; i0++)
    // for(int m0 = 2; m0 < maxBondDim; m0 *= 2){

    for (int i0 = 0; i0 < 1; i0++)                  //  SPEEDY COMPRESSION
    for(int m0 = minBondDim; m0 < maxBondDim; m0 *= 2){     //  SPEEDY COMPRESSION

        m = i0*maxBondDim+m0;
        if (svd_ran)
            make_svd_compress(m);
        else
            make_random_guess(m);
        /** debug **/ std::cout << "current bond size: " << m << std::endl;
        double prev_local_norm;
        globalcompress_full(1);
        prev_local_norm = localcompress_norm;
        for(int s = 0; s < num_sweeps; ++s){
            if(globalcompress_full(10)){  /** debug **/         
                auto g = getMPS(); 
                g[0] *= sqrt(localcompress_norm);
                std::cout << " compress_SingleState::compute: successful (BondDim,Sweeps)=(" << m <<  ", " << s+1 << ").";
                std::cout << " Current error (acc,relacc) = (" << accuracy << ", " << relacc << ")." << std::endl;
                std::cout << " Norm of new MPS: " << normMPS( g ) << std::endl;
                return g;
            }
            if( fabs( prev_local_norm - localcompress_norm ) < epsilon_diff ){
                /** debug **/ std::cout << " compress_SingleState::compute: converged at current bond-size, moving on from this bond dimension " << std::endl;
                break;
            }
            prev_local_norm = localcompress_norm;
        }
    }

    auto g = getMPS();
    g[0] *= sqrt(localcompress_norm);
    std::cout << " compress_SingleState::compute: \t Could not compress within (maxBondDim,maxSweeps)=(" << maxBondDim <<  ", " << num_sweeps << ")." << " Current error (acc,relacc) = (" << accuracy << ", " << relacc << ")." << " Norm of new MPS: " << normMPS( g ) << std::endl;
    std::cout << " compress_SingleState::compute: \t norm( refMPS )     = " << norm_refMPS << std::endl;
    std::cout << " compress_SingleState::compute: \t localcompress_norm = " << localcompress_norm << std::endl;
    return g;
}

std::vector<UniTensor>
    compress_SingleState::getMPS
    ()
{
    std::vector< UniTensor > curMPS;
    std::stack < UniTensor > tmp = L1;
    curMPS.resize(L1.size());
    for(int i = 0; i < L1.size(); ++i){
        curMPS[L1.size()-1-i] = tmp.top();
        tmp.pop();
    }
    return curMPS;
}

void 
    compress_SingleState::make_svd_compress
    (int m)
{
    L1 = std::stack<UniTensor>();
    L2 = std::stack<UniTensor>();
    R1 = std::stack<UniTensor>();
    R2 = std::stack<UniTensor>();
    UniTensor triv({ Bi(1), Bo(1) }, "");
    triv.setRawElem( std::vector<double>({1.0}) );
    L2.push ( triv );
    R2.push ( triv );
    std::vector<UniTensor> x = refMPS;
    left_normalize_and_truncate_mps( x, m );
    for(int i = 0; i < mpsL; ++i){
        make_helper1( x[i], L1, L2 );
    }
}

/*! 
 * \brief Initializes the guess MPS to be random.
 */
void
    compress_SingleState::make_random_guess
    (int m)
{
    L1 = std::stack<UniTensor>();
    L2 = std::stack<UniTensor>();
    R1 = std::stack<UniTensor>();
    R2 = std::stack<UniTensor>();
    UniTensor triv({ Bi(1), Bo(1) }, "");
    triv.setRawElem( std::vector<double>({1.0}) );
    L2.push ( triv );
    R2.push ( triv );
    std::vector<UniTensor> x = random_left_normalized_mps( m, 2, mpsL );
    for(int i = 0; i < mpsL; ++i){
        make_helper1( x[i], L1, L2 );
    }
}

/*! 
 * \brief Helper function for `make_random_guess`
 * Updates stack with the new tensor \param t0.
 */
void
    compress_SingleState::make_helper1
    (UniTensor &t0, std::stack<UniTensor> &left1, std::stack<UniTensor> &left2)
{
    t0.setLabel                   ({-1,-2,0});
    refMPS[left1.size()].setLabel ({-3,-2,1});
    left2.top().setLabel          ({-1,-3});
    left2.push( t0 * (left2.top() * refMPS[left1.size()]) );
    left1.push( t0 );
}

/*!
 * \brief Performs the local update for a single-site compression
 *
 * Performs the local update for a single-site compression
 */
void
    compress_SingleState::localcompress
    (bool direction)
{
                                  // if (verbose) std::cout << " compress_SingleState::localcompress BEGINS: direction = " << (direction? "left to right" : "right to left") << ", site " << L1.size() << std::endl;
    localcompress_h1(direction);  // if (verbose) std::cout << " done localcompress_h1 " << std::endl;
    localcompress_h2();           // if (verbose) std::cout << " done localcompress_h2 " << std::endl;
    localcompress_h3(direction);  // if (verbose) std::cout << " done localcompress_h3 " << std::endl;
    localcompress_h4(direction);  // if (verbose) std::cout << " done localcompress_h4 " << std::endl;
    if (verbose)
        localcompress_h5();
}

/*!
 * \brief localcompress helper function 1 
 * First prepares state:
 * * Pops out the two tensors which are going to be updated
 * Checks validity of current state:
 * * Sizes of stacks are aligned: `|L1|+|R1| = mpsLength`, `|i2| = |i1|+1`, for `i = {L,R}`
 * * Bond dimensions of to-be-contracted tensors match
 */
void
    compress_SingleState::localcompress_h1
    (bool direction)
{
    assert(L1.size() + R1.size() == mpsL);
    assert(L2.size() - L1.size() == 1);
    assert(R2.size() - R1.size() == 1);
    //
    if (direction) // true == left grows, right shrinks
    {
        R1.pop();
        R2.pop();
    } else {
        L1.pop();
        L2.pop();
    }
    //
    assert(  L2.top().bond()[1].dim() == refMPS[ int(L1.size()) ].bond()[0].dim()  );
    assert(  R2.top().bond()[1].dim() == refMPS[ int(L1.size()) ].bond()[2].dim()  );
}

/*!
 * \brief localcompress helper function 2 
 * Applies compression-operation:
 * * Contracts `L2.top` with `R2.top` and `refMPS[ L1.size() ]`
 */
void
    compress_SingleState::localcompress_h2
    ()
{
    int k = L1.size();
    refMPS[k].setLabel({-1,1,-2});
    L2.top().setLabel ({0,-1});
    R2.top().setLabel ({2,-2});
    localcompress_t1 = L2.top() * ( refMPS[k] * R2.top() );
}

/*!
 * \brief localcompress helper function 3
 * Prepares for storage by performing SVD.
 */
void
    compress_SingleState::localcompress_h3
    (bool direction)
{
    localcompress_t1.setLabel({0,1,2});
    localcompress_t1.permute({0,1,2}, (direction?2:1) ); // if true == left grows,right shrinks; then 2; else 1
    svd10 x( localcompress_t1 );
    x.compute_triple( 1e-8, maxBondDim );
    localcompress_trunc_error = x.truncation_error();
    localcompress_norm = x.norm();
    localcompress_normalizedtensor = ( direction ? x.UU() : x.VV() );
}

/*!
 * \brief localcompress helper function 4
 * Stores normalized tensor into Stack1 and Stack2
 */
void
    compress_SingleState::localcompress_h4
    (bool direction)
{
    if (direction) // true == left grows, right shrinks
    {
        L2.top().setLabel                       ({-1,-2});
        localcompress_normalizedtensor.setLabel ({-1,-3,0});
        refMPS[ L1.size() ].setLabel            ({-2,-3,1});
        L2.push( localcompress_normalizedtensor*( L2.top() *refMPS[ L1.size() ] ) );
        L1.push( localcompress_normalizedtensor );
    } else {
        R2.top().setLabel                       ({-1,-2});
        localcompress_normalizedtensor.setLabel ({0,-3,-1});
        refMPS[ L1.size() ].setLabel            ({1,-3,-2});
        //
        R1.push( localcompress_normalizedtensor );
        R2.push( localcompress_normalizedtensor*( refMPS[ L1.size() ]*R2.top() ) );
    }
}

/*!
 * \brief localcompress helper function 5
 * Prints out information post-update
 */
void
    compress_SingleState::localcompress_h5
    ()
{
    std::cout << " status of local compression at position "                  << L1.size()                        << std::endl;
    std::cout << " norm:             " << std::fixed << std::setprecision(10) << localcompress_norm               << std::endl;
    std::cout << " truncation error: " << std::fixed << std::setprecision(10) << localcompress_trunc_error        << std::endl;
}

/*!
 * \brief globalcompress_half
 * `direction`
 * * `true`  if left grows, right shrinks
 * * `false` if left shrinks, right grows 
 */
void
    compress_SingleState::globalcompress_half
    (bool direction)
{
    for(int i = 0; i < mpsL; ++i)
        localcompress(direction);
}

/*!
 * \brief globalcompress_full
 * Sweeps the system  \param sweeps  times and returns a boolean
 * \param true  `==` accuracy attained, successful compression
 * \param false `==` accuracy not attained, keep sweeping
 */
bool
    compress_SingleState::globalcompress_full
    (int local_sweeps)
{
    for(int i = 0; i < local_sweeps; ++i){
        globalcompress_half(false); // false == right grows
        globalcompress_half(true);  // true  == left grows
        accuracy = fabs( norm_refMPS - localcompress_norm );
        relacc   = accuracy / norm_refMPS;
        if ( relacc < epsilon_diff ){
      //    std::cout << " compress_SingleState::globalcompress_full: "
      //              << "Reached within desired precision in " << i 
      //              << " sweeps with acc: " << accuracy
      //              << " and relacc: " << relacc << std::endl;
      //    std::cout << " compress_SingleState::globalcompress_full "
      //              << "(norm_refMPS, localcompress_norm)"
      //              << "(" << norm_refMPS << ", " << localcompress_norm << ")" << std::endl;
           return true;
        }
    }
    std::cout << " compress_SingleState::globalcompress_full: localcompress_norm = " << localcompress_norm << std::endl;
    return false;
}



/**********************************************
 *           Output Helper functions          *
 **********************************************/
void
    compress_SingleState::output_current_mps
    ()
{
    std::stack < UniTensor > tmp1 = L1;
    int tL = tmp1.size();
    std::vector < UniTensor > t1; 
    t1.resize(tL); 
    //
    for(int i = 0; i < tL; i++){
        t1[tL-1-i] = tmp1.top();
        tmp1.pop();
    }
    //
    std::cout << "      Printing Current MPS        " << std::endl;
    for(int i = 0; i < t1.size(); ++i){
        std::cout << i << std::endl;
        std::cout << t1[i] << std::endl;
    }
    tmp1 = R1;
    while(tmp1.size() > 0){
        std::cout << tmp1.top() << std::endl;
        tmp1.pop();
    }
    std::cout << "      Ending Print of Current MPS " << std::endl;
}

void
    compress_SingleState::output_current_mps_just_diagrams
    ()
{
    std::stack < UniTensor > tmp1 = L1;
    int tL = tmp1.size();
    std::vector < UniTensor > t1; 
    t1.resize(tL); 
    //
    for(int i = 0; i < tL; i++){
        t1[tL-1-i] = tmp1.top();
        tmp1.pop();
    }
    //
    std::cout << "      Printing Current MPS        " << std::endl;
    for(int i = 0; i < t1.size(); ++i){
        std::cout << i << std::endl;
        //std::cout << t1[i] << std::endl;
        t1[i].printDiagram();
    }
    tmp1 = R1;
    while(tmp1.size() > 0){
        std::cout << tmp1.top() << std::endl;
        tmp1.pop();
    }
    std::cout << "      Ending Print of Current MPS " << std::endl;
}

void
    compress_SingleState::output_current_mps_just_raw_elems
    ()
{
    std::stack < UniTensor > tmp1 = L1;
    int tL = tmp1.size();
    std::vector < UniTensor > t1; 
    t1.resize(tL); 
    //
    for(int i = 0; i < tL; i++){
        t1[tL-1-i] = tmp1.top();
        tmp1.pop();
    }
    //
    std::cout << "      Printing Current MPS        " << std::endl;
    for(int i = 0; i < t1.size(); ++i){
        //std::cout << t1[i] << std::endl;
        std::cout << i << std::endl;
        t1[i].printDiagram();
        t1[i].printRawElem(true,10);
    }
    tmp1 = R1;
    while(tmp1.size() > 0){
        std::cout << tmp1.top() << std::endl;
        tmp1.pop();
    }
    std::cout << "      Ending Print of Current MPS " << std::endl;
}


/***********************************
 * Begin Compress Helper Functions *
 **********************************/

std::vector< UniTensor >
    compress_mps
    ( std::vector< UniTensor > mps, int maxD, int minD, double eps, bool sr, int sweeps )
{
    compress_SingleState x ( mps, maxD, minD, eps, eps, sr, sweeps );
    return x.compute( false ); // false = printout
}

std::vector< UniTensor >
    compress_sum_mps
    ( std::vector< std::vector< UniTensor > > lstMPS, int maxD, int minD, double eps, bool sr, int sweeps )
{
    std::vector< UniTensor > mps(lstMPS[0].size()); // direct sum
    std::vector< UniTensor > tmp(lstMPS.size());    // collection of mps' on one site
    for(int i = 0; i < lstMPS[0].size(); ++i){
        for(int j = 0; j < lstMPS.size(); ++j)
            tmp[j] = lstMPS[j][i];
        mps[i] = oplus(tmp, i==0, i==(lstMPS[0].size()-1));
    }
    compress_SingleState x ( mps, maxD, minD, eps, eps, sr, sweeps );
    return x.compute( false ); // 5 sweeps, false = printout
}

std::vector< UniTensor >
    compress_mpo_mps
    ( std::vector< UniTensor > initialMPO, std::vector< UniTensor > initialMPS, int maxD, int minD, double eps, bool sr, int sweeps )
{
    std::vector<UniTensor> tmp = apply_mpo_to_mps (initialMPO, initialMPS, "compress_mpo_mps");
    //               std::cout << "compress_mpo_mps start save Lprime(Psi5) uncompressed " << std::endl;
    //               for (int i = 0; i < tmp.size(); ++i) tmp[i].save("LprimePsi5uncompress/"+intstr(i));
    compress_SingleState x ( tmp, maxD, minD, eps, eps, sr, sweeps );
    auto tmp2 = x.compute( false );
    //               std::cout << "compress_mpo_mps start save Lprime(Psi5) compressed " << std::endl;
    //               for (int i = 0; i < tmp2.size(); ++i) tmp2[i].save("LprimePsi5compress/"+intstr(i));
    return tmp2;
}

/**********************************
 *  End Compress Helper Functions *
 *********************************/

#endif
