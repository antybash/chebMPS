
#ifndef TSVD_H
#define TSVD_H

#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <uni10.hpp>

using namespace uni10;
using namespace useful;

/*!
 * \class svd10
 * \brief Performs SVD on a UniTensor
 * This class implements various methods that are useful for MPS
 * computations.
 * 
 * Explanation:
 *  * Load Tensor (possibly specifying the orientation of the proposed SVD).
 *  * Compute SVD.
 *  * Organize information in a convenient form.
 */
class svd10 {
public:
    svd10( UniTensor& );
    svd10( UniTensor&, std::vector<int>, std::vector<int> );
    void compute_triple(double, int);
    void compute_truncation_from_bond_dimension(int);
    //
    UniTensor UU();
    UniTensor US();
    UniTensor SV();
    UniTensor VV();
    //
    double norm();
    double truncation_error();
private:
    UniTensor T;  // tensor whose svd will be computed
    UniTensor ut; // U
    UniTensor st; // S
    UniTensor vt; // V^tr
    std::vector<double> full_singular_values;
    std::vector<double> truncated_singular_values;
};

/*****************
 * Implmentation *
 ****************/


svd10::svd10
    ( UniTensor &t )
{
    T = t;
}


svd10::svd10
    ( UniTensor &t, std::vector<int> row_inds, std::vector<int> col_inds )
{
    T = t;
    // union of   row_inds   and   col_inds :
    int R = row_inds.size();
    int C = col_inds.size();
    std::vector<int> permuted_indices(R+C);
    for(int i = 0; i < R; ++i)
        permuted_indices[i] = row_inds[i];
    for(int i = 0; i < C; ++i)
        permuted_indices[R+i] = col_inds[i];
    //
    T.permute(permuted_indices, (int)row_inds.size());
}

UniTensor svd10::UU () { return ut; }
UniTensor svd10::VV () { return vt; }

UniTensor svd10::US ()
{
    std::vector<int> a = ut.label();
    std::vector<int> b = st.label();
    int k = a.size()-1;
    std::iota(a.begin(),   a.end(), 0);
    std::iota(b.begin()+1, b.end(), k);
    a[k] = -1;
    b[0] = -1;
    ut.setLabel(a);
    st.setLabel(b);
    return ut*st;
}

UniTensor svd10::SV ()
{
    std::vector<int> a = st.label();
    std::vector<int> b = vt.label();
    int k = a.size()-1;
    std::iota(a.begin(),   a.end(), 0);
    std::iota(b.begin()+1, b.end(), k);
    a[k] = -1;
    b[0] = -1;
    st.setLabel(a);
    vt.setLabel(b);
    return st*vt;
}


/*!
 * \brief Computes `U, S, V^\dag`
 */
void
    svd10::compute_triple
    ( double epsilon, int maxBondDim )
{
    std::vector< Matrix > tsvd = T.getBlock().svd();
    //truncation
    int len = tsvd[1].row(); 
    int chi = std::min(len, maxBondDim);
    if(epsilon > 0)
        for (int i = 1; i < chi; ++i){
            if ( tsvd[1].at(i,i) < epsilon ){
                chi = i;
                break;
            }
        }
    //save singular values
    full_singular_values.resize(len);
    truncated_singular_values.resize(chi);
    for(int i = 0; i < len; ++i)
        full_singular_values[i]      = tsvd[1].at(i,i);
    for(int i = 0; i < chi; ++i)
        truncated_singular_values[i] = tsvd[1].at(i,i);
    // modify the outer tensor bond dimensions.
    std::vector<Bond> allBonds = T.bond();
    std::vector<Bond>::const_iterator first = allBonds.begin();
    std::vector<Bond>::const_iterator mid   = allBonds.begin() + (int)T.inBondNum();
    std::vector<Bond>::const_iterator last  = allBonds.end();
    std::vector<Bond> uBonds (first, mid );   // = U
    uBonds.push_back( Bo(chi) );              //   U-
    std::vector<Bond> vBonds(mid, last);      //   V=
    vBonds.insert( vBonds.begin(), Bi(chi) ); //  -V
    //resize
    tsvd[0].resize( (int) tsvd[0].row(), chi );
    tsvd[1].resize( chi, chi);
    tsvd[2].resize( chi, (int) tsvd[2].col() );
    //
    ut = UniTensor( uBonds,             "u of "+T.getName() );
    st = UniTensor( {Bi(chi), Bo(chi)}, "s of "+T.getName() );
    vt = UniTensor( vBonds,             "v of "+T.getName() );
    //
    ut.putBlock( tsvd[0] );
    st.putBlock( tsvd[1] );
    vt.putBlock( tsvd[2] );
}

double
    svd10::norm
    ()
{
    double sum=0.0;
    for(int i = 0; i < truncated_singular_values.size(); ++i)
        sum += truncated_singular_values[i]*truncated_singular_values[i];
    return sum;
}

double
    svd10::truncation_error
    ()
{
    double sum_trunc=0.0;
    double sum_full=0.0;
    //
    for(int i = 0; i < truncated_singular_values.size(); ++i)
        sum_trunc += truncated_singular_values[i]*truncated_singular_values[i];
    for(int i = 0; i < full_singular_values.size(); ++i)
        sum_full += full_singular_values[i]*full_singular_values[i];
    //
    return sum_full-sum_trunc;
}

/******************************************************************
 *                       Old functions                            *
 *****************************************************************/

std::vector< UniTensor >
    tensor_U_SV
    ( UniTensor T, double epsilon, int maxBondDim)
{
//       // not yet   /******************************************************
//       // not yet    * Remark: The singular values are always rescaled to 
//       // not yet    * have norm 1; ie. sum | s_k |^2 = 1.
//       // not yet    *****************************************************/
//   std::vector< UniTensor > x = tensor_svd(T, epsilon, maxBondDim);
//   // set S label to contract with V
//   x[1].setLabel({0,-1}); 
//       // not yet   double s = 0;
//       // not yet   for(int i = 0; i < x[1].getBlock().row(); ++i)
//       // not yet       s += ( x[1].getBlock().at(i,i) ) * ( x[1].getBlock().at(i,i) );
//       // not yet   x[1] *= 1.0/sqrt( s ); //normalization
//   // set V label to contract with S
//   std::vector<int> x2b = x[2].label();
//   x2b[0] = -1;  
//   x[2].setLabel(x2b);
    svd10 a(T);
    a.compute_triple( epsilon, maxBondDim );
    return { a.UU(), a.SV() };
}

std::vector< UniTensor >
    tensor_US_V
    ( UniTensor T, double epsilon, int maxBondDim)
{
//       // not yet   /******************************************************
//       // not yet    * Remark: The singular values are always rescaled to 
//       // not yet    * have norm 1; ie. sum | s_k |^2 = 1.
//       // not yet    *****************************************************/
//   std::vector< UniTensor > x = tensor_svd(T, epsilon, maxBondDim);
//   // set U label to contract with S
//   std::vector<int> x0bonds = x[0].label();
//   x0bonds.back() = -1;   
//   x[0].setLabel(x0bonds);
//   // set S label to contract with U
//   x[1].setLabel({-1,1}); 
//       // not yet   double s = 0;
//       // not yet   for(int i = 0; i < x[1].getBlock().row(); ++i)
//       // not yet       s += ( x[1].getBlock().at(i,i) ) * ( x[1].getBlock().at(i,i) );
//       // not yet   x[1] *= 1.0/sqrt( s ); //normalization
//       // not yet   std::cout << " tensor US_V: " << x[1].norm() << std::endl;
    svd10 a(T);
    a.compute_triple( epsilon, maxBondDim );
    return { a.US(), a.VV() };
}

void
    left_normalize_and_truncate_mps
    (std::vector<UniTensor> &x, int m)
{
    int L = x.size();
    for(int i = 0; i < L-1; ++i){
        x[i].setLabel({0,1,2});
        x[i].permute({0,1,2},2);
        svd10 t(x[i]);
        t.compute_triple(-1.0, m);
        x[i]   = t.UU();
        UniTensor tmp = t.SV();
        tmp.setLabel   ({0,-1});
        x[i+1].setLabel({-1,1,2});
        x[i+1] = tmp * x[i+1];
    }
}



#endif
