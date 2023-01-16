
#ifndef TIMEEVOLUTION_H
#define TIMEEVOLUTION_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <iterator>
#include <vector>
#include <complex>
#include <stack>
#include <set>
#include <functional>
#include <ctime>
#include <uni10.hpp>

#include "compress.h"
#include "tsvd.h"
#include "expv.h"
#include "hamiltonians.h"
#include "bigVec.h"
#include "bigStack.h"

class timeEvolveTwoSite{
    public:
        //timeEvolveTwoSite(int, double, double, std::vector<UniTensor>); // assume heisenberg XXZ 
        timeEvolveTwoSite(int, double, std::vector<UniTensor>, std::vector<UniTensor>); // custom hamiltonian
        void sweep ();
        void compress_mps_state();
        void compress_mps_state(int,int);
        std::vector< UniTensor > final_mps_state();
        void save (std::string);
    private:
        int mpsL;
        double epsilon;
        double dt;
        std::vector<UniTensor>  hamiltonian;
        BigVec   cur_state;
        BigStack L3;        // workhorse
        BigStack R3;        // workhorse
        UniTensor AL,AR,AC,CC;
        UniTensor KH,HH;
        /******/
        UniTensor computeHH();
        UniTensor tensor_expv(double, const UniTensor&, const UniTensor& );
        UniTensor apply_exp_to_ham1(double, UniTensor&, const UniTensor &);
        UniTensor apply_exp_to_ham2(double, UniTensor&, const UniTensor &);
        UniTensor apply_exp_to_kam (double, UniTensor&, const UniTensor &);
        void update_AC_from_CC(bool);
        void decomp_pushL3R3_defineK_AL_AR(bool);
        void decomp_pushL3R3_defineK_AL_AC_one_site(bool);
        void decomp_pushL3R3_defineK_AL_AC_two_site(bool);
        //void initializeHamiltonian(double);
        void initialization_right_normalize_cur_state();
        void prepare_state_for_time_evolution();
        /***** 2 site functions *****/
        void make_two_site_AC(bool);
        void decomp_pushL3R3_updateAC_one_site(bool);
        /****    debug params    ****/
        bool time_debug;
        int maxNN = 4000;
};

std::vector<UniTensor> QR_RQ_using_SVD(UniTensor t, bool QR_RQ)
{
    auto x = svd10(t);
    x.compute_triple(1e-6,64); // RANDOM NUMBERS
    if(QR_RQ) { return {x.UU(), x.SV()}; }
    else      { return {x.US(), x.VV()}; }
}

std::vector< UniTensor > timeEvolveTwoSite::final_mps_state()
{ 
    // NOT GOOD
    assert(0 == "could do better");
    std::vector<UniTensor> vec(mpsL);
    for (int i = 0; i < mpsL; ++i)
        vec[i] = cur_state.getElem(i);
    return vec;
}


UniTensor timeEvolveTwoSite::computeHH()
{
        //std::cout << "computeHH" << std::endl;
    assert(    (L3.size() + R3.size() == 1 + mpsL) || (L3.size() + R3.size() == mpsL)    );
    int k = L3.size()-1;
    UniTensor left  = L3.top();
    UniTensor right = R3.top();
    if (L3.size() + R3.size() == 1 + mpsL){
            //std::cout << "one-site Hamiltonian" << std::endl;
        UniTensor ham   = hamiltonian[k];
        left.setLabel ({0,-1,3}   );
        ham.setLabel  ({-1,4,1,-2});
        right.setLabel({2,-2,5}   );
        UniTensor ans = left*(ham*right);
        ans.permute({0,1,2,3,4,5},3);
        return ans;
    } else {
            //std::cout << "two-site Hamiltonian" << std::endl;
        assert( L3.size() + R3.size() == mpsL );
        UniTensor ham1 = hamiltonian[k];
        UniTensor ham2 = hamiltonian[k+1];
        left.setLabel  ({0,-1,4}   );
        ham1.setLabel  ({-1,5,1,-2});
        ham2.setLabel  ({-2,6,2,-3});
        right.setLabel ({3,-3,7}   );
        UniTensor ans = left*(ham1*(ham2*right));
        ans.permute({0,1,2,3,4,5,6,7},4);
        return ans;
    }
}

UniTensor timeEvolveTwoSite::tensor_expv(double t, const UniTensor &H, const UniTensor &current_state_vector)
{
    Matrix mH = H.getBlock();
    int N = mH.row();      assert( N == mH.col() && N == current_state_vector.elemNum() );
    Eigen::MatrixXd A(N,N);
    Eigen::VectorXd v(N);

    std::time_t t1 = std::time(nullptr), t2;

    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j)
            A(i,j) = mH.at(i,j);   //[i*N+j];
    for(int i = 0; i < N; ++i)
        v(i) = current_state_vector.getBlock()[i];

    t2 = std::time(nullptr);
    if( N > maxNN && time_debug ){ 
        std::cout << " tensor_expv --- Matrix size = " << N;
        std::cout << "; loading time = " << t2-t1 << " seconds " << std::endl;
    }
    t1 = t2;

    std::tuple<Eigen::VectorXd,double,double> output_EXPV = expv(t, A, v);
    t2 = std::time(nullptr);
    if( N > maxNN && time_debug ){ 
        std::cout << " tensor_expv --- Matrix size = " << N;
        std::cout << "; compute time = " << t2-t1 << " seconds " << std::endl;
    }
    t1 = t2;

    Eigen::VectorXd w = std::get<0>(output_EXPV);

    UniTensor ten_w = current_state_vector;
    Matrix mat = ten_w.getBlock();
    int r = mat.row();
    int c = mat.col();
    assert(N == r*c);
    for(int i = 0; i < r; ++i)
        for(int j = 0; j < c; ++j){
            mat.at(i,j) = w(i*c+j);
        }
    t2 = std::time(nullptr);
    if( N > maxNN && time_debug ){ 
        std::cout << " tensor_expv --- Matrix size = " << N;
        std::cout << "; storing time = " << t2-t1 << " seconds " << std::endl;
    }
    t1 = t2;

    ten_w.putBlock( mat );
    return ten_w;
}

UniTensor timeEvolveTwoSite::apply_exp_to_ham1(double FULL_DELTA_T_WITH_SIGN, UniTensor &H, const UniTensor &current_state_vector)
{
    H.permute({0,1,2,3,4,5},3);
    //w1.setLabel({0,1,2});                BE CAREFUL WITH THIS ONE~!
    //w1.permute ({0,1,2},3);              I DON'T REMEMBER WHY THESE TWO LINES ARE HERE
    if(H.elemNum() > 5000*5000)
        std::cout << " ham1: " << H.bond(0).dim() << " "
                               << H.bond(1).dim() << " "
                               << H.bond(2).dim() << " "
                               << H.bond(3).dim() << " "
                               << H.bond(4).dim() << " "
                               << H.bond(5).dim() << std::endl;
    return tensor_expv(FULL_DELTA_T_WITH_SIGN, H, current_state_vector);
}                                          

UniTensor timeEvolveTwoSite::apply_exp_to_ham2(double FULL_DELTA_T_WITH_SIGN, UniTensor &H, const UniTensor &vector)
{
    H.permute({0,1,2,3,4,5,6,7},4);
    //w1.setLabel({0,1,2,3});               BE CAREFUL WITH THIS ONE~!
    //w1.permute ({0,1,2,3},4);             I DON'T REMEMBER WHY THESE TWO LINES ARE HERE
    if(H.elemNum() > 5000*5000)
        std::cout << " ham2: " << H.bond(0).dim() << " "
                               << H.bond(1).dim() << " "
                               << H.bond(2).dim() << " "
                               << H.bond(3).dim() << " "
                               << H.bond(4).dim() << " "
                               << H.bond(5).dim() << " "
                               << H.bond(6).dim() << " "
                               << H.bond(7).dim() << std::endl;
    return tensor_expv(FULL_DELTA_T_WITH_SIGN, H, vector);
}


UniTensor timeEvolveTwoSite::apply_exp_to_kam(double FULL_DELTA_T_WITH_SIGN, UniTensor &K, const UniTensor &current_cc)
{
    K.permute({0,1,2,3},2);
    return tensor_expv(FULL_DELTA_T_WITH_SIGN, K, current_cc);
}

void timeEvolveTwoSite::make_two_site_AC(bool direction)
{
    // takes AC(1site) + AR -> AC(2site)
    // takes AL + AC(1site) -> AC(2site)
    if (direction) {
        AC.setLabel({0,1,-1});
        AR.setLabel({-1,2,3});
        //std::cout << " AC " << std::endl; std::cout << AC << std::endl;
        //std::cout << " AR " << std::endl; std::cout << AR << std::endl;
        AC = AC*AR;
    } else {
        AL.setLabel({0,1,-1});
        AC.setLabel({-1,2,3});
        AC = AL*AC;
    }
}

/*! \brief QR or RQ decomp, update AC */
void timeEvolveTwoSite::decomp_pushL3R3_updateAC_one_site(bool direction)
{
    UniTensor R,L,W,t,b;
    int k = (int)L3.size()-1;
    std::vector<UniTensor> x = QR_RQ_using_SVD( AC, direction );
    if (direction) { cur_state.setElem(k,  x[0]); AC = x[1]; }
    else           { cur_state.setElem(k+1,x[1]); AC = x[0]; }
    // Append to L3 or R3; 
    if (direction)
    { 
        L = L3.top();              L.setLabel({-1,-2,-3});
        W = hamiltonian[k];        W.setLabel({-2,-5,-4,1});
        t = cur_state.getElem(k);  t.setLabel({-1,-4,0});
        b = t;                     b.setLabel({-3,-5,2}); 
        //
        t = t*(L*(W*b));
        t.permute({0,1,2},2);
        L3.push( t );
        if (k+2 < cur_state.size())
            AR = cur_state.getElem(k+2);
    } else {
        R = R3.top();                R.setLabel({-1,-2,-3});
        W = hamiltonian[k+1];        W.setLabel({1,-5,-4,-2});
        t = cur_state.getElem(k+1);  t.setLabel({0,-4,-1});
        b = t;                       b.setLabel({2,-5,-3}); 
        //
        t = t*(W*(b*R));
        t.permute({0,1,2},2);
        R3.push( t );
        if (0 <= k-1)
            AL = cur_state.getElem(k-1);
    }
    assert(L3.size()+R3.size() == mpsL+1);
    /* at this point:
     *   AC == one-site AC
     *   L3 and R3 are separated by one-site
     */
}

void timeEvolveTwoSite::sweep()
{
    int N = mpsL;
    assert(L3.size() == 1);
    assert(R3.size() == mpsL);
    AC = cur_state.getElem(0);
    AR = cur_state.getElem(1);
    make_two_site_AC(true);
    R3.pop();
    for(int n = 0; n < N-2; ++n){
        HH = computeHH();
        AC = apply_exp_to_ham2(-1.0*dt/2.0, HH, AC);
        AC.setLabel({0,1,2,3});
        AC.permute( {0,1,2,3},2);
        decomp_pushL3R3_updateAC_one_site(true); // push L3
        HH = computeHH();
        AC = apply_exp_to_ham1(1.0*dt/2.0, HH, AC);
        make_two_site_AC(true);
        R3.pop();
    }
    HH = computeHH();
    AC = apply_exp_to_ham2(-1.0*dt, HH, AC); // this defines AC(N, t+dt)
    assert(L3.size() == mpsL-1);
    assert(R3.size() == 1);
    for(int n = N-3; n >= 0; --n){         //  two site:  0  1  ...  n-2   n-1   |  n  n+1  |   n+2   n+3  ...
        AC.setLabel({0,1,2,3});
        AC.permute( {0,1,2,3},2);
        decomp_pushL3R3_updateAC_one_site(false); // push R3
        HH = computeHH();
        AC = apply_exp_to_ham1(1.0*dt/2.0, HH, AC);
        make_two_site_AC(false);
        L3.pop();
        assert(L3.size() == n+1);
        HH = computeHH();
        AC = apply_exp_to_ham2(-1.0*dt/2.0, HH, AC);
    }
    AC.setLabel({0,1,2,3});
    AC.permute ({0,1,2,3},2);
    decomp_pushL3R3_updateAC_one_site(false); // push R3
    cur_state.setElem(0,AC);
}

void timeEvolveTwoSite::compress_mps_state()
{
    //vector<UniTensor> compress_mps ( vector< UniTensor > mps, int maxD, double eps, bool sr, int sweeps )
    cur_state.fill(   compress_mps( cur_state.mps_form(), 40, 32, 1e-6, true, 10)   );
    prepare_state_for_time_evolution();
}

void timeEvolveTwoSite::compress_mps_state(int maxD, int minD)
{
    //vector<UniTensor> compress_mps ( vector< UniTensor > mps, int maxD, int minD, double eps, bool sr, int sweeps )
    cur_state.fill(  compress_mps(cur_state.mps_form(), maxD, minD, 1e-6, true, 10) );
    prepare_state_for_time_evolution();
}


/*
void timeEvolveTwoSite::initializeHamiltonian(double lam)
{
    hamiltonian.resize(mpsL);
    hamiltonian[0] = ham::xxz::mpo_leftmost(lam);
    for (int i = 1; i < mpsL-1; ++i)
        hamiltonian[i] = ham::xxz::mpo_mid(lam);
    hamiltonian[mpsL-1]= ham::xxz::mpo_rightmost(lam);
}

timeEvolveTwoSite::timeEvolveTwoSite
    (int length, double lam, double delta_time_init, std::vector< UniTensor > mps_init )
{
    mpsL       = length;
    dt         = delta_time_init;
    cur_state  = mps_init;
    initializeHamiltonian(lam);
    time_debug = true;
    //
    UniTensor triv3({ Bi(1), Bi(1), Bo(1) }, "triv3");
    triv3.setRawElem( std::vector<double>({1.0}) );
    L3.push( triv3 );
    R3.push( triv3 );
    for(int i = mpsL-1; i >= 0; --i){
        UniTensor R = R3.top();         R.setLabel({-1,-2,-3});
        UniTensor W = hamiltonian[i];   W.setLabel({1,-5,-4,-2});
        UniTensor t =   cur_state[i];   t.setLabel({0,-4,-1});
        UniTensor b =   cur_state[i];   b.setLabel({2,-5,-3}); 
        UniTensor result = t*(W*(b*R));
        result.permute({0,1,2},2);
        R3.push( result );
    }
}
*/

void timeEvolveTwoSite::initialization_right_normalize_cur_state()
{
    UniTensor tmp;
    std::vector<UniTensor> x;
    for(int i = mpsL-1; i > 0; --i){
        cur_state.getElem(i,tmp);
        tmp.setLabel({0,1,2});
        tmp.permute({0,1,2},1);
        x = QR_RQ_using_SVD( tmp, false );
        cur_state.setElem(i,x[1]);
        x[0].setLabel({-1,2});
        //
        cur_state.getElem(i-1,tmp);
        tmp.setLabel({0,1,-1});
        tmp = (tmp*x[0]);
        cur_state.setElem(i-1,tmp);
    }
}

void timeEvolveTwoSite::prepare_state_for_time_evolution()
{
    UniTensor triv3({ Bi(1), Bi(1), Bo(1) }, "triv3");
    triv3.setRawElem( std::vector<double>({1.0}) );
    L3 = BigStack( mpsL, vector_string_description(mpsL,"L3_timeEvolve_") );
    R3 = BigStack( mpsL, vector_string_description(mpsL,"R3_timeEvolve_") );
    L3.push( triv3 );
    R3.push( triv3 );
    initialization_right_normalize_cur_state(); 
    UniTensor R,W,t,b;
    for(int i = mpsL-1; i > 0; --i){
        R =               R3.top();   R.setLabel({-1,-2,-3});
        W =         hamiltonian[i];   W.setLabel({1,-5,-4,-2});
        t =   cur_state.getElem(i);   t.setLabel({0,-4,-1});
        b =   cur_state.getElem(i);   b.setLabel({2,-5,-3}); 
        std::cout << " timeEvolveTwoSite::prepare_state_for_time_evolution " << i << std::endl;
        //
        std::cout << " timeEvolveTwoSite::prepare_state_for_time_evolution " << "R: "; 
        for(int tmpi = 0; tmpi < R.bondNum(); ++tmpi) std::cout << R.bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " timeEvolveTwoSite::prepare_state_for_time_evolution " << "W: "; 
        for(int tmpi = 0; tmpi < W.bondNum(); ++tmpi) std::cout << W.bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " timeEvolveTwoSite::prepare_state_for_time_evolution " << "t: "; 
        for(int tmpi = 0; tmpi < t.bondNum(); ++tmpi) std::cout << t.bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " timeEvolveTwoSite::prepare_state_for_time_evolution " << "b: "; 
        for(int tmpi = 0; tmpi < b.bondNum(); ++tmpi) std::cout << b.bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        t = t*(W*(b*R)); // the LHS is just the result of the multiplication
        t.permute({0,1,2},2);
        R3.push( t );
    }
}

timeEvolveTwoSite::timeEvolveTwoSite
    (int length, double delta_time_init, std::vector< UniTensor > mps_init, std::vector<UniTensor> mpo_init )
{
    mpsL        = length;
    dt          = delta_time_init;
    //cur_state   = mps_init;
    cur_state = BigVec( mpsL, vector_string_description(mpsL,"cur_state_timeEvolve") );
    convertVec2BigVec(mps_init,cur_state);
    //
    hamiltonian = mpo_init;
    time_debug = true;
    std::cout << " begin init time Evolution " << std::endl;
    prepare_state_for_time_evolution();
    std::cout << " done  init time Evolution " << std::endl;
}

/*******************************/
//
//   Misc. Commands 
//
/*******************************/

void timeEvolveTwoSite::save
    (std::string filename)
{
    std::cout << " ...... saving current state ...... " << std::endl;
    UniTensor tmp10;

    for(int i = 0; i < cur_state.size(); ++i){
        tmp10 = cur_state[i];
        tmp10.setName (filename+properfilename(0,0,i));
        tmp10.save    (filename+properfilename(0,0,i));
    }
}


#endif
