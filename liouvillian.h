
#ifndef LIOUVILLIAN_H
#define LIOUVILLIAN_H

#include <vector>
#include <uni10.hpp>

#include "helper_compress2.h"
#include "hamiltonians.h"

using namespace useful;

std::vector< UniTensor >
    HamIden
    (int fullDim, double lam)
{
    assert ( fullDim > 0 );
    assert ( fullDim % 4 == 0 );
    //
    std::vector< UniTensor > mpo(fullDim);
    mpo[0] = ham::xxz::mpo_leftmost(lam);   // ONLY TEMP FIX
    for(int i = 1; i < fullDim-1; ++i){
        if (i%2 == 0) { mpo[i] = ham::xxz::mpo_mid(lam);     }
        else          { mpo[i] = ham::xxz::mpo_aux_mid(lam); }
    }
    mpo[fullDim-1] = ham::xxz::mpo_aux_rightmost(lam);
    return mpo;
}

std::vector< UniTensor >
    negIdenHam
    (int fullDim, double lam)
{
    assert ( fullDim > 0 );
    assert ( fullDim % 4 == 0 );
    //
    std::vector< UniTensor > mpo(fullDim);
    mpo[0] = ham::xxz::mpo_aux_leftmost_negative(lam);
    for(int i = 1; i < fullDim-1; ++i){
        if (i%2 == 0) { mpo[i] = ham::xxz::mpo_aux_mid(lam); }
        else          { mpo[i] = ham::xxz::mpo_mid(lam);     }
    }
    mpo[fullDim-1] = ham::xxz::mpo_rightmost(lam);
    return mpo;
}

std::vector< UniTensor >
    Identity
    (int fullDim, double lam)
{
    return std::vector< UniTensor > ( fullDim, mpo_actual_identity(lam) );
}

std::vector< UniTensor >
    Liouvillian
    ( int fullDim, double lam )
{
    return add_many_mpo(
            std::vector< std::vector<UniTensor> >(
                {
                    HamIden(fullDim, lam),
                 negIdenHam(fullDim, lam) 
                } )
            );
}

std::vector< UniTensor >
    decoratedLiouvillian
    ( int fullDim, double lam, double alpha )
{
    // alpha * H_P \otimes I_Q - alpha * I_P \otimes H_Q
    //
    return add_many_mpo(
            std::vector< std::vector< UniTensor > > (
                {
             mpo_times( alpha,    HamIden(fullDim, lam)),
             mpo_times( alpha, negIdenHam(fullDim, lam))
                } )
       );
             // until further notice this will be the Liouvillian
             // mpo_times(  beta,   Identity(fullDim, lam))
}

std::vector< UniTensor >
    physical_parameters_Liouvillian
    ( int fullDim, double lam, double W, double Wprime )
{
    return decoratedLiouvillian( fullDim, lam, (2.0*Wprime)/W );
}



#endif
