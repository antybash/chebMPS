
#ifndef HELPER_COMPRESS_TWO
#define HELPER_COMPRESS_TWO

#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <uni10.hpp>
#include <random>
#include <functional>
#include <algorithm>

using namespace uni10;

/*!
 * \brief Helper commands for program
 */
namespace useful {

    Bond Bi(int n) { return uni10::Bond(uni10::BD_IN,  n); }
    Bond Bo(int n) { return uni10::Bond(uni10::BD_OUT, n); }

    std::string intstr(int n)
    {
        std::stringstream in;
        in << n;
        return in.str();
    }

    std::string doublestr(double d)
    {
        std::stringstream in;
        in << std::fixed << std::setprecision(2) << d;
        return in.str();
    }

    inline std::string properfilename(double dt, double lam, int len)
    {
        std::stringstream in;
        in << std::fixed << std::setprecision(2) << dt;
        in << std::fixed << std::setprecision(2) << lam;
        in << len;
        return in.str();
    }

    std::string random_string( size_t length )
    {
        auto randchar = []() -> char
        {
            const char charset[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            return charset[ rand() % max_index ];
        };
        std::string str(length,0);
        std::generate_n( str.begin(), length, randchar );
        return str;
    }

    std::vector<std::string> random_vector_string_description(int len, std::string prologue)
    {
        std::vector<std::string> tmps(len);
        for(int i = 0; i < len; ++i)
            tmps[i] = prologue+random_string(7)+intstr(i);
        return tmps;
    }

    std::vector<std::string> vector_string_description(int len, std::string prologue)
    {
        std::vector<std::string> tmps(len);
        for(int i = 0; i < len; ++i)
            tmps[i] = prologue+intstr(i);
        return tmps;
    }


    Matrix makeMatrix(std::vector<double> elems, int r, int c)
    {
        Matrix m(r,c); 
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                m.at(i,j) = elems[i*c+j];
        return m;
    }

    UniTensor makeMPS(std::vector<double> mat_combined, std::vector<Bond> idx, int alpha, int beta)
    {
        //  MPS rank (a,s,b) consists of s matrices of size (a,b).
        //  To store MPSmatrices in a row:  (x11,x12,x21,x22, y11,y12,y21,y22)
        UniTensor t(idx,"");
        Matrix m = makeMatrix( mat_combined, 2, alpha*beta );
        t.putBlock( m );
        t.setLabel({1,0,2});
        t.permute({0,1,2},2);
        return t;
    }

    UniTensor SpinOneMakeMPS(std::vector<double> mat_combined, std::vector<Bond> idx, int alpha, int beta)
    {
        //  MPS rank (a,s,b) consists of s matrices of size (a,b).
        //  To store MPSmatrices in a row:  (x11,x12,x21,x22, y11,y12,y21,y22)
        UniTensor t(idx,"");
        Matrix m = makeMatrix( mat_combined, 3, alpha*beta );
        t.putBlock( m );
        t.setLabel({1,0,2});
        t.permute({0,1,2},2);
        return t;
    }

    UniTensor random_mps(int a, int b, int c)
    {
        UniTensor x({Bi(a), Bi(b), Bo(c)});
        x.randomize();
        return x;
    }

    UniTensor random_mpo(int a, int b, int c, int d)
    {
        UniTensor x({Bi(a), Bi(b), Bo(c), Bo(d)});
        x.randomize();
        return x;
    }

    UniTensor reflect_tensor(UniTensor mps)
    {   // pass by value, create new mps, permutes and returns fresh copy.
        mps.setLabel({0,1,2});
        mps.permute({2,1,0},1);
        return mps;
    }

    UniTensor t0() { return makeMPS({1,0,0,-1},                        {Bi(2),Bo(1),Bo(2)},1,2); }
    UniTensor t1() { return makeMPS({0,1.0/sqrt(2.0),1.0/sqrt(2.0),0}, {Bi(2),Bo(2),Bo(1)},2,1); }

    int maximum_bond_size(const std::vector<UniTensor> &vec)
    {
        int m = 0;
        for(auto x : vec)
            for(auto b : x.bond())
                m = (b.dim() > m ? b.dim() : m);
        return m;
    }

    double expect(const std::vector<UniTensor> &mps, const std::vector<UniTensor> &mpo)
    {
        UniTensor triv3({ Bi(1), Bi(1), Bo(1) }, "triv3");
        triv3.setRawElem( std::vector<double>({1.0}) );
        UniTensor tmp = triv3;
        int L = mps.size();
        assert(L == mpo.size());
        for(int i = L-1; i >= 0; --i){
            UniTensor R = tmp;    R.setLabel({-1,-2,-3});
            UniTensor W = mpo[i]; W.setLabel({1,-5,-4,-2});
            UniTensor t = mps[i]; t.setLabel({0,-4,-1});
            UniTensor b = mps[i]; b.setLabel({2,-5,-3}); 
            tmp = t*(W*(b*R));
            tmp.permute({0,1,2},2);
        }
        triv3.setLabel({0,1,2});
        tmp = triv3*tmp;
        return tmp.getBlock()[0];
    }

    double expect2(const std::vector<UniTensor> &mps, const std::vector<UniTensor> &mpo1, const std::vector<UniTensor> &mpo2)
    {
        UniTensor triv4({ Bi(1), Bi(1), Bo(1), Bo(1) }, "triv4");
        triv4.setRawElem( std::vector<double>({1.0}) );
        UniTensor tmp = triv4;
        int L = mps.size();
        assert(L == mpo1.size());
        assert(L == mpo2.size());
        for(int i = L-1; i >= 0; --i){
            UniTensor R  = tmp;     R.setLabel({-1,-2,-3,-4});
            UniTensor W1 = mpo1[i]; W1.setLabel({1,-6,-5,-2});
            UniTensor W2 = mpo2[i]; W2.setLabel({2,-7,-6,-3});
            /******************************************/
            UniTensor t  = mps[i];  t.setLabel({0,-5,-1});
            UniTensor b  = mps[i];  b.setLabel({3,-7,-4}); 
            tmp = t*(W1*(W2*(b*R)));
            tmp.permute({0,1,2,3},2);
        }
        triv4.setLabel({0,1,2,3});
        tmp = triv4*tmp;
        return tmp.getBlock()[0];
    }

    std::vector<UniTensor> QR_RQ_tensor(UniTensor t, bool QR_RQ)
    {
        std::vector<Matrix> x;
        if (QR_RQ) x = t.getBlock().qr();
        else       x = t.getBlock().rq();
        int middleBond_size   = x[0].col();
        // modify the outer tensor bond dimensions.
        std::vector<Bond> allBonds              = t.bond();
        std::vector<Bond>::const_iterator first = allBonds.begin();
        std::vector<Bond>::const_iterator mid   = allBonds.begin() + (int)t.inBondNum();
        std::vector<Bond>::const_iterator last  = allBonds.end();
        std::vector<Bond> inc(first, mid );              // = A_L       =  C      inc == "incoming bonds"
        inc.push_back( Bo(middleBond_size) );            //   A_L -        C  -   out == "outgoing bonds"
        std::vector<Bond> out (mid, last);               //    C  =  or   A_R = 
        out.insert( out.begin(), Bi(middleBond_size) );  //  - C        - A_R   
        //
        std::string nameA, nameB;
        if (QR_RQ){ nameA = "Q in QR of "; nameB = "R in QR of "; }
        else      { nameA = "R in RQ of "; nameB = "Q in RQ of "; }
        UniTensor A = UniTensor( inc, nameA +t.getName() );
        UniTensor B = UniTensor( out, nameB +t.getName() );
        //
        A.putBlock( x[0] );
        B.putBlock( x[1] );
        return std::vector<UniTensor>({A,B});
    }


    /*!
     * \brief tensor of rank(m,d,m) that is left normalized
     */
    UniTensor random_left_normalized_3tensor( int m1, int d, int m2 )
    {
        assert(d*m1 >= m2);
        UniTensor t({Bi(m1),Bi(d),Bo(d*m1)},"");
        t.randomize();
        Matrix M = t.getBlock();
        std::vector<Matrix> x = M.qr();
        x[0].resize(m1*d,m2);
        UniTensor tt({Bi(m1),Bi(d),Bo(m2)},"left-normalized");
        tt.putBlock(x[0]);
        return tt;
    }

    std::vector<UniTensor> random_left_normalized_mps(int m, int d, int L)
    {
        assert(d >= 2);
        int i = 2;
        while(i < m)
            i*=2;
        assert(i == m);
        std::vector<UniTensor> x;
        x.reserve(L);
        assert( x.size() == 0 );
        for(int i = 2; i < m+1; i*=2){
            if (x.size()-1 == L){
                x.push_back( random_left_normalized_3tensor(i/2,d,1) );
                return x;
            }
            x.push_back( random_left_normalized_3tensor(i/2, d, i) );
        }
        while(x.size() < L-1){
            x.push_back( random_left_normalized_3tensor(m,d,m) );
        }
        x.push_back(random_left_normalized_3tensor(m,d,1));
        return x;
    }

    /************************
     * MPO Helper Functions *
     ***********************/

    UniTensor matBasis(int k, int l, int r, int c)
    {
        Matrix m(r,c);
        m.at(k,l) = 1;
        UniTensor t({Bi(r),Bo(c)},"matBasis");
        t.putBlock(m);
        return t;
    }

    UniTensor sI()
    {
        Matrix m(2,2); UniTensor t({Bi(2),Bo(2)}, "sI");
        m.identity();  t.putBlock(m);
        return t;
    }
    UniTensor sIn()
    {
        Matrix m(2,2); UniTensor t({Bi(2),Bo(2)}, "sIn");
        m.identity();  t.putBlock(-1.0*m);
        return t;
    }
    UniTensor sX()
    {
        Matrix m(2,2); UniTensor t({Bi(2),Bo(2)}, "sX");
        m.at(1,0) = 0.5; m.at(0,1) = 0.5; t.putBlock(m);
        return t;
    }
    UniTensor sYp()
    {
        Matrix m(2,2); UniTensor t({Bi(2),Bo(2)}, "sYp");
        m.at(1,0) = -0.5; m.at(0,1) = 0.5; t.putBlock(m);
        return t;
    }
    UniTensor sYn()
    {
        Matrix m(2,2); UniTensor t({Bi(2),Bo(2)}, "sYn");
        m.at(1,0) = 0.5; m.at(0,1) = -0.5; t.putBlock(m);
        return t;
    }
    UniTensor sZ()
    {
        Matrix m(2,2); UniTensor t({Bi(2),Bo(2)}, "sZ");
        m.at(0,0) = 0.5; m.at(1,1) = -0.5; t.putBlock(m);
        return t;
    }
    namespace spinOne {
        double INVSQRTTWO=sqrt(1.0/2.0);
        UniTensor sI()  { Matrix m(3,3); UniTensor t({Bi(3),Bo(3)}, "sI");  m.identity();  t.putBlock(m); return t; }
        UniTensor sIn() { Matrix m(3,3); UniTensor t({Bi(3),Bo(3)}, "sIn"); m.identity();  t.putBlock(-1.0*m); return t; }
        UniTensor sX()  { Matrix m(3,3); UniTensor t({Bi(3),Bo(3)}, "sX");  m.at(1,0) = 1.0;  m.at(0,1) = 1.0; m.at(1,2) = 1.0; m.at(2,1) = 1.0; t.putBlock(INVSQRTTWO*m); return t; }
        UniTensor sYp() { Matrix m(3,3); UniTensor t({Bi(3),Bo(3)}, "sYp"); m.at(1,0) =-1.0;  m.at(0,1) = 1.0; m.at(1,2) = 1.0; m.at(2,1) =-1.0; t.putBlock(INVSQRTTWO*m); return t; }
        UniTensor sYn() { Matrix m(3,3); UniTensor t({Bi(3),Bo(3)}, "sYn"); m.at(1,0) = 1.0;  m.at(0,1) =-1.0; m.at(1,2) =-1.0; m.at(2,1) = 1.0; t.putBlock(INVSQRTTWO*m); return t; }
        UniTensor sZ()  { Matrix m(3,3); UniTensor t({Bi(3),Bo(3)}, "sZ");  m.at(0,0) = 1.0;                                    m.at(2,2) =-1.0; t.putBlock(m); return t; }
        UniTensor t0()  { return SpinOneMakeMPS({1,0,   0,0,    0,-1},                        {Bi(3),Bo(1),Bo(2)},1,2); } 
        UniTensor t1()  { return SpinOneMakeMPS({0,1.0/sqrt(2.0),   0,0,    1.0/sqrt(2.0),0}, {Bi(3),Bo(2),Bo(1)},2,1); }
    }

    UniTensor MPO1make(UniTensor virt, UniTensor phys)
    {
        UniTensor t = otimes(virt,phys);
        t.permute({0,1,3,2},2);
        t.setLabel({0,1,2,3});
        return t;
    }

    Matrix t2m(UniTensor t) { return t.getBlock(); }
    Matrix mpm(Matrix a, Matrix b) { return a+b; }
    //
    UniTensor MPO2add(std::vector< UniTensor > v)
    {
        std::vector< Matrix > m;
        m.resize(v.size());
        std::transform(v.begin(), v.end(), m.begin(), t2m);
        Matrix zero(m[0].row(), m[0].col());
        UniTensor t(v[0].bond());
        t.putBlock( std::accumulate(m.begin(), m.end(), zero, mpm) );
        return t;
    }

    UniTensor mpo_actual_identity(double lam)
    {
        return MPO1make(matBasis(0,0,1,1),sI());
    }

    /************************************************
     * Spin-Structure Factor -- Fourier Transformed *
     ***********************************************/

    UniTensor mpo_spin_leftmost(int idx, int L, double k) {
        return MPO2add({ MPO1make(matBasis(0,0,1,2),sqrt(2.0/(double (L+1)))*sin(double(idx*k))*sZ()),
                         MPO1make(matBasis(0,1,1,2),sqrt(2.0/(double (L+1)))*sI())                    });
    }
    UniTensor mpo_spin_mid(int idx, double k) {
        return MPO2add({ MPO1make(matBasis(0,0,2,2),sI()                     ),
                         MPO1make(matBasis(1,0,2,2),sin(double(idx*k))* sZ() ),
                         MPO1make(matBasis(1,1,2,2),sI()                     ) });
    }
    UniTensor mpo_spin_rightmost(int idx, double k) {
        return MPO2add({ MPO1make(matBasis(0,0,2,1),sI()                     ),
                         MPO1make(matBasis(1,0,2,1),sin(double(idx*k))* sZ()     ) });
    }

    // identity
    UniTensor mpo_spin_aux_leftmost() {
        return MPO2add({     MPO1make(matBasis(0,0,1,2),sI())    });
    }
    UniTensor mpo_spin_aux_mid() {
        return MPO2add({ MPO1make(matBasis(0,0,2,2),sI() ),
                         MPO1make(matBasis(1,1,2,2),sI() ) });
    }
    UniTensor mpo_spin_aux_rightmost() {
        return MPO2add({ MPO1make(matBasis(0,0,2,1),sI()) });
    }

    std::vector< UniTensor > 
        mpo_spin_str
        (int L, double k)
    {
        std::vector< UniTensor > lst(2*L);
        lst[0] = mpo_spin_leftmost(1,L,k);
        for(int i = 1; i < L; ++i){
            lst[2*i-1] = mpo_spin_aux_mid();
            lst[2*i] = mpo_spin_mid(i+1,k);
        }
        lst[2*L-1] = mpo_spin_aux_rightmost();
        return lst;
    }

    /*******************************************************
     * END OF Spin-Structure Factor -- Fourier Transformed *
     ******************************************************/







    /*********************************************************
     * SPIN ONE Spin-Structure Factor -- Fourier Transformed *
     ********************************************************/

    UniTensor SPINONE_mpo_spin_leftmost(int idx, int L, double k) {
        return MPO2add({ MPO1make(matBasis(0,0,1,2),sqrt(2.0/(double (L+1)))*sin(double(idx*k))*spinOne::sZ()),
                         MPO1make(matBasis(0,1,1,2),sqrt(2.0/(double (L+1)))*spinOne::sI())                    });
    }
    UniTensor SPINONE_mpo_spin_mid(int idx, double k) {
        return MPO2add({ MPO1make(matBasis(0,0,2,2),spinOne::sI()                     ),
                         MPO1make(matBasis(1,0,2,2),sin(double(idx*k))* spinOne::sZ() ),
                         MPO1make(matBasis(1,1,2,2),spinOne::sI()                     ) });
    }
    UniTensor SPINONE_mpo_spin_rightmost(int idx, double k) {
        return MPO2add({ MPO1make(matBasis(0,0,2,1),spinOne::sI()                         ),
                         MPO1make(matBasis(1,0,2,1),sin(double(idx*k))* spinOne::sZ()     ) });
    }

    // identity
    UniTensor SPINONE_mpo_spin_aux_leftmost() {
        return MPO2add({     MPO1make(matBasis(0,0,1,2),spinOne::sI())    });
    }
    UniTensor SPINONE_mpo_spin_aux_mid() {
        return MPO2add({ MPO1make(matBasis(0,0,2,2),spinOne::sI() ),
                         MPO1make(matBasis(1,1,2,2),spinOne::sI() ) });
    }
    UniTensor SPINONE_mpo_spin_aux_rightmost() {
        return MPO2add({ MPO1make(matBasis(0,0,2,1),spinOne::sI()) });
    }

    std::vector< UniTensor > 
        SPINONE_mpo_spin_str
        (int L, double k)
    {
        std::vector< UniTensor > lst(2*L);
        lst[0] = SPINONE_mpo_spin_leftmost(1,L,k);
        for(int i = 1; i < L; ++i){
            lst[2*i-1] = SPINONE_mpo_spin_aux_mid();
            lst[2*i] =   SPINONE_mpo_spin_mid(i+1,k);
        }
        lst[2*L-1] = SPINONE_mpo_spin_aux_rightmost();
        return lst;
    }

    /*******************************************************
     * END OF Spin-Structure Factor -- Fourier Transformed *
     ******************************************************/


    /*!
     *  \brief  `|| mps ||^2`
     */
    double 
        normMPS
        ( std::vector< UniTensor > &m )
    {
        UniTensor lsf({ Bi(1), Bo(1) }, "lsf");
        lsf.setRawElem( std::vector<double>({1.0}) );
        //
        UniTensor t1;
        int L = m.size();
        for(int i = 0; i < L-1; ++i){
            t1 = m[i];
            lsf.setLabel  ({-1,-2});
            t1.setLabel   ({-1,-3,0});
            m[i].setLabel ({-2,-3,1});
            lsf = (lsf*t1)*m[i];
        }
        t1 = m[L-1];
        lsf.setLabel    ({-1,-2});
        t1.setLabel     ({-1,-3,-4});
        m[L-1].setLabel ({-2,-3,-4});
        lsf = (lsf*t1)*m[L-1];
        return lsf.getBlock()[0];
    }

    std::vector<UniTensor>
        normalizeMPS
        (std::vector<UniTensor> vec)
    {
        double d = normMPS(vec);
        vec[0] *= 1.0/sqrt(d);
        return vec;
    }

    UniTensor
        oplus
        ( std::vector< UniTensor > &ts, bool leftmost, bool rightmost )
    {
        /********************************************************************
         * Input: Vector< MPS-states >                                      *
         * Assume: Each MPS state has the same physical-dimension           *
         *         (middle leg)                                             *
         *                                                                  *
         * Output: Performs   ts[0] + ts[1] + .. + ts[n]   without any      *
         * compression. This is equivalent to applying a direct sum to the  *
         * matrices associated with each spin component.                    *
         *******************************************************************/
        int n = ts.size();
        assert(n >= 1);
        if (n == 1) return ts[0];
        //
        int d = ts[0].bond()[1].dim(); // physical bond dimension
        int sx, sy;
        int alphaI, betaI;
        //
        int alphaSUM = 0;
        int betaSUM  = 0;
        for(int i = 0; i < n; ++i){
            alphaSUM += ts[i].bond()[0].dim();
            betaSUM  += ts[i].bond()[2].dim();
        }
        if(leftmost)  alphaSUM = 1;
        if(rightmost) betaSUM = 1; 
        UniTensor Total( {Bi(alphaSUM), Bi(d), Bo(betaSUM)}, "Total" );
        Total.set_zero();
        //
        for(int a = 0; a < d; ++a){
            sx = 0; 
            sy = 0;
            for(int i = 0; i < n; ++i){
                alphaI = ts[i].bond()[0].dim();
                betaI  = ts[i].bond()[2].dim();
                for(int x = 0; x < alphaI; ++x){
                    for(int y = 0; y < betaI; ++y){
                        Total.putAt( (std::vector<int>) {x+sx, a, y+sy}, 
                            ts[i].at( (std::vector<int>) {x, a, y} ) 
                                   );
                    }
                }
                if(!leftmost) sx += alphaI;
                if(!rightmost) sy += betaI;
            }
        }
        return Total;
    }

    UniTensor
        oplus_mpo
        ( std::vector< UniTensor > &ts )
    {
        int n = ts.size();
        assert(n >= 1);
        if (n == 1) return ts[0];
        //
        int d1 = ts[0].bond()[1].dim(); // physical bond dimension
        int d2 = ts[0].bond()[2].dim(); // physical bond dimension
        int sx, sy;
        int alphaI, betaI;
        //
        int alphaSUM = 0;
        int betaSUM  = 0;
        for(int i = 0; i < n; ++i){
            alphaSUM += ts[i].bond()[0].dim();
            betaSUM  += ts[i].bond()[3].dim();
        }
        UniTensor Total( {Bi(alphaSUM), Bi(d1), Bi(d2), Bo(betaSUM)}, "Total" );
        Total.set_zero();
        //
        for(int a = 0; a < d1; ++a)
        for(int b = 0; b < d2; ++b){
            sx = 0; 
            sy = 0;
            for(int i = 0; i < n; ++i){
                alphaI = ts[i].bond()[0].dim();
                betaI  = ts[i].bond()[3].dim();
                for(int x = sx; x < sx + alphaI; ++x){
                    for(int y = sy; y < sy + betaI; ++y){
                        Total.putAt( (std::vector<int>) {x, a,b, y}, 
                            ts[i].at( (std::vector<int>) {x-sx, a,b, y-sy} ) 
                                   );
                    }
                }
                sx += alphaI;
                sy += betaI;
            }
        }
        return Total;
    }


    std::vector< UniTensor >
        add_many_mpo
        ( std::vector< std::vector<UniTensor> > lst_mpo )
    {
        std::vector< UniTensor > mpo(lst_mpo[0].size()); // direct sum
        std::vector< UniTensor > tmp(lst_mpo.size());    // collection of mps' on one site
        for(int i = 0; i < lst_mpo[0].size(); ++i){
            for(int j = 0; j < lst_mpo.size(); ++j)
                tmp[j] = lst_mpo[j][i];
            mpo[i] = oplus_mpo(tmp);
        }
        return mpo;
    }


    std::vector< UniTensor >
        apply_mpo_to_mps
        ( std::vector<UniTensor> mpo1, std::vector<UniTensor> &mps1, std::string explanation )
    {
        int n = mpo1.size();
        assert (n == mps1.size());
        std::vector< UniTensor > res(n);
        UniTensor t;
        for(int i = 0; i < n; ++i){
            mpo1[i].setLabel( {0,-1,1,2} );
            mps1[i].setLabel( {3,-1,4}   );
            t = mpo1[i]*mps1[i];
            t.combineBond( (std::vector<int>) {0,3} );
            t.combineBond( (std::vector<int>) {2,4} );
            t.permute((std::vector<int>){0,1,2},2);
            res[i] = t;
            //std::cout << explanation << " apply_mpo_to_mps: done " << i << " out of " << n << std::endl;
        }
        return res;
    }

    std::vector< UniTensor >
        mpo_times
        ( double a, std::vector< UniTensor > mpo )
    {
        mpo[0] = mpo[0]*a;
        return mpo;
    }

    std::vector< UniTensor >
        mps_times
        ( double a, std::vector< UniTensor > mps )
    {
        // This is exactly the same as mpo_times
        mps[0] = mps[0]*a;
        return mps;
    }

    /*!
     * \brief Computes inner product of two states which are in matrix-product form (ie. `vector<UniTensor>`)
     */
    double
        inner_product
        ( std::vector< UniTensor > &x, std::vector< UniTensor > &y )
    {
        assert(x.size() == y.size());
        //
        UniTensor lsf({ Bi(1), Bo(1) }, "lsf");
        lsf.setRawElem( std::vector<double>({1.0}) );
        UniTensor tmp1, tmp2;
        for(int i = 0; i < x.size(); ++i){
            tmp1 = x[i];
            tmp2 = y[i];
            lsf.setLabel ({-1,-2});
            tmp1.setLabel({-1,-3,0});
            tmp2.setLabel({-2,-3,1});
            lsf = (lsf*tmp1)*tmp2;
        }
        return lsf.getBlock()[0];
    }

    std::string 
        params2string
        ( int fullDim, std::string Dstr, std::string bstr, int coef )
    {
        std::stringstream in;
        in << "L=" << fullDim/2 << "/"
           << "L=" << fullDim/2 << "_"
           << "Delta=" << Dstr << "_"
           << "beta=" << bstr << "_"
           << "coef=" << coef;
        std::cout << "Reading: " << in.str() << std::endl;
        return in.str();
    }

    double
        cheb
        ( int order, double wprime )
    {
        return std::cos( double(order)*std::acos( wprime ));
        /*
        if (order == 0)
            return 1;
        else if (order == 1)
            return wprime;
        else{
            double tn2 = 1;
            double tn1 = wprime;
            double tn;
            for(int i = 2; i <= order; ++i){
                tn = 2.0*wprime*tn1 - tn2;
                tn2 = tn1;
                tn1 = tn;
            }
            return tn;
        }
        */
    }


    std::vector< UniTensor > spin_half_Z_profile_mpo (int i, int mpoL)
    {
        UniTensor iden ({Bi(1),Bi(2),Bo(2),Bo(1)}); iden.putBlock (sI().getBlock());
        UniTensor spinZ({Bi(1),Bi(2),Bo(2),Bo(1)}); spinZ.putBlock(sZ().getBlock());


        std::vector<UniTensor> MPO(mpoL,iden);
        MPO[i] = spinZ;
        //for(int i = 0; i < mpoL; ++i) std::cout << "(" << MPO[i].bond()[0].dim() << "," << MPO[i].bond()[1].dim() << "," << MPO[i].bond()[2].dim() << "," << MPO[i].bond()[3].dim()<< ")";
        //std::cout << std::endl;
        return MPO;
    }

    std::vector< double > spin_half_Z_profile( std::vector<UniTensor> mps)
    {
        int LL = mps.size();
        //for(int i = 0; i < LL; ++i) std::cout << "(" << mps[i].bond()[0].dim() << "," << mps[i].bond()[1].dim() << "," << mps[i].bond()[2].dim() << ")";
        //std::cout << std::endl;
        std::vector<double> prof(LL,0.0);
        for(int i = 0; i < LL; ++i){
            //std::cout << "spin1/2--- Z --- profile: " << i << std::endl;
            //std::cout << "spin1/2--- Z --- profile: " << mps.size() << std::endl;
            prof[i] = expect( mps, spin_half_Z_profile_mpo(i,LL) );
        }
        return prof;
    }


}

#endif 
