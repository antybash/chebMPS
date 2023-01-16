
#ifndef GROW_H
#define GROW_H

#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <uni10.hpp>
#include <memory>
#include "csvfile.h"
#include "compress.h"
//#include "liouvillian.h"

#define PI 3.14159265358979

using namespace useful;

/// \class stem
/// \brief Computes the Chebyshev approximation for a given thermal state.
/*!
 * Description:
 *    * Input: Thermal State, k, omega.
 *    * Output: Information to be able to compute S_zz(k,omega)
 */
class stem{
public:
    stem( std::vector<UniTensor>, std::unique_ptr<HamiltonianClass>&, int, int, double, int, 
        double, double, double, int, double, double, int, std::string, bool, bool);

    double compute(double);
    void save(bool);

    int maxBondDim;
    int minBondDim;
    int num_sweeps;
    double epsilon;
    int physLen;
    int thermalLen;
    int cur_height;
    int max_height;
    double Delta;
    double beta;
    double k_vector;
    double W;
    double Wprime;
    bool first_time;
    bool savingMode;
    //int save_regime;
    std::string stoR;
    std::unique_ptr<HamiltonianClass> hamiltonian_class;
    std::vector<double> jackson_coefs;
    std::vector<double> matrix_elements;
    std::vector<double> norm_chebvecs;
    std::vector< UniTensor > mps_car;  // " car  " == first  element in list
    std::vector< UniTensor > mps_cadr; // " cadr " == second element in list
    std::vector< UniTensor > t0;
    double normT0;
    std::vector< UniTensor > liouvillian1;
    std::vector< UniTensor > liouvillian2;
    std::vector< UniTensor > observable_mpo;
    //
    void grow_once();
    void grow_full();
    //
    std::vector<double> soms;
    std::vector<double>  oms;
    void S_of_omega();
};

/*!
 *  Computes S(om) for om varying from 0 to (bandwidth = 2.0)
 *  Break 0 to W into arbitrary (50) discrete points, scale to 0->W', compute using `compute'
 */
void
    stem::S_of_omega
    ()
{
    std::vector<double> ans_soms(50);
    std::vector<double> ans_oms (50);
    for(int i = 0; i < 50; i++){
        ans_soms[i] = compute(double(2.0)*double(i)/double(50));
        ans_oms[i]  =        (double(2.0)*double(i)/double(50));
    }
    soms = ans_soms;
    oms  = ans_oms;
}

/*!
 * \brief Initializes stem class.
 *
 * The initialization consists of defining every quantity which will
 * be used the main recursion relation of the Chebyshev expansion.
 * The only objects that need to be computed are the |t_n>. 
 */
stem::stem
    ( 
      std::vector<UniTensor> thMPS,
      std::unique_ptr<HamiltonianClass> &hc,
      int mD, 
      int minD,
      double eps, 
      int int_physLen, 
      double Delt, 
      double b, 
      double k, 
      int N,      // maximum number Chebyshev basis vectors
      double ww, 
      double wwprime, 
      int sweeps,
      //,int save_reg
      std::string temp_mps_storage,
      bool firstTime,
      bool save_at_all
    )
{
    mps_cadr          = thMPS;
    hamiltonian_class = std::move(hc);
    maxBondDim        = mD;
    minBondDim        = minD;
    num_sweeps        = sweeps;
    epsilon           = eps;
    //
      physLen         = int_physLen;
   thermalLen         = physLen*2;
    // the above change was made Feb 2, 2018
    cur_height        = 0;
    max_height        = N;
    Delta             = Delt;
    beta              = b;
    k_vector          = k;
    W                 = ww;
    Wprime            = wwprime;
    stoR              = temp_mps_storage;
    first_time        = firstTime;
    savingMode        = save_at_all;
    //save_regime = save_reg;
    // fill jackson coefficients

    jackson_coefs.reserve(N);
    for (int n = 0; n < N; ++n)
        jackson_coefs.push_back( (N-n+1)*cos(PI * n / (N+1))/(N+1) + sin(PI*n/(N+1))/tan(PI/(N+1))/(N+1) );

    if (first_time){
        liouvillian1   = hamiltonian_class -> liouvillian1(2.0*(2.0*Wprime/W));   //mpo_times( 2.0* (2.0*Wprime/W),  hamiltonian_class ->    HamIden());
        liouvillian2   = hamiltonian_class -> liouvillian2(2.0*(2.0*Wprime/W));   //mpo_times( 2.0* (2.0*Wprime/W),  hamiltonian_class -> negIdenHam());

        std::cout << " physLen " << physLen << std::endl;
        observable_mpo = hamiltonian_class -> ham_mpo_spin_str(physLen,k_vector);
        std::cout << observable_mpo.size() << std::endl;

        std::cout << " stem::stem INITIALIZATION: done making and storing\n\
            Liouvillian, done making and storing spin structure factor\n\
            (observable_mpo), starting to make t0" << std::endl;

        std::cout << " ObsMPO  " << observable_mpo.size() << std::endl;
        std::cout << " Liouv1  " << liouvillian1.size()   << std::endl;
        std::cout << " Liouv2  " << liouvillian2.size()   << std::endl;
        std::cout << " MPScadr " << mps_cadr.size()       << std::endl;

        //
        std::cout << " stem::..... " << "ObsMPO: "; 
        for(int tmpi = 0; tmpi < observable_mpo[8].bondNum(); ++tmpi) std::cout << observable_mpo[8].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "Liouv1: "; 
        for(int tmpi = 0; tmpi < liouvillian1[8].bondNum(); ++tmpi) std::cout << liouvillian1[8].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "Liouv2: "; 
        for(int tmpi = 0; tmpi < liouvillian2[8].bondNum(); ++tmpi) std::cout << liouvillian2[8].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "MPScadr: "; 
        for(int tmpi = 0; tmpi < mps_cadr[8].bondNum(); ++tmpi) std::cout << mps_cadr[8].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //////////

        //
        std::cout << " stem::..... " << "ObsMPO: "; 
        for(int tmpi = 0; tmpi < observable_mpo[9].bondNum(); ++tmpi) std::cout << observable_mpo[9].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "Liouv1: "; 
        for(int tmpi = 0; tmpi < liouvillian1[9].bondNum(); ++tmpi) std::cout << liouvillian1[9].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "Liouv2: "; 
        for(int tmpi = 0; tmpi < liouvillian2[9].bondNum(); ++tmpi) std::cout << liouvillian2[9].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "MPScadr: "; 
        for(int tmpi = 0; tmpi < mps_cadr[9].bondNum(); ++tmpi) std::cout << mps_cadr[9].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //////////

        //
        std::cout << " stem::..... " << "ObsMPO: "; 
        for(int tmpi = 0; tmpi < observable_mpo[0].bondNum(); ++tmpi) std::cout << observable_mpo[0].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "Liouv1: "; 
        for(int tmpi = 0; tmpi < liouvillian1[0].bondNum(); ++tmpi) std::cout << liouvillian1[0].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "Liouv2: "; 
        for(int tmpi = 0; tmpi < liouvillian2[0].bondNum(); ++tmpi) std::cout << liouvillian2[0].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //
        std::cout << " stem::..... " << "MPScadr: "; 
        for(int tmpi = 0; tmpi < mps_cadr[0].bondNum(); ++tmpi) std::cout << mps_cadr[0].bond()[tmpi].dim() << ", "; std::cout << std::endl;
        //////////


        mps_cadr = compress_mpo_mps( observable_mpo, mps_cadr, maxBondDim, minBondDim, 1e-8, true, num_sweeps ); // 50 == maxBondDim, false == random, true == svd
        t0       = mps_cadr; // correct for beta=0.0 according to Mathematica (Jul 10, cf EverythingExact.nb)

        std::cout << " stem::stem initialization: done constructing t0, now making t1" << std::endl;
        mps_car  = compress_sum_mps( 
                std::vector< std::vector<UniTensor> > ({
                    apply_mpo_to_mps( mpo_times(0.5,liouvillian1), t0, "initializing |t1>:  (2W'/W)HxI|t0>" ),
                    apply_mpo_to_mps( mpo_times(0.5,liouvillian2), t0, "initializing |t1>: -(2W'/W)IxH|t0>" )
                                                      }),
                maxBondDim,
                minBondDim,
                epsilon
                ,true // svd
                , num_sweeps
                );

        auto t1 = mps_car;
        assert(matrix_elements.size() == 0);
        assert(  norm_chebvecs.size() == 0);
        mps_cadr = t0; matrix_elements.push_back( inner_product(t0,t0) ); norm_chebvecs.push_back( inner_product(t0,t0) );
        mps_car  = t1; matrix_elements.push_back( inner_product(t0,t1) ); norm_chebvecs.push_back( inner_product(t1,t1) );

        int tmpL = liouvillian1.size();
        assert(tmpL == liouvillian2.size());
        assert(tmpL == mps_car.size()); 
        assert(tmpL == mps_cadr.size()); 
        assert(tmpL == t0.size() );
        assert(tmpL == observable_mpo.size() );

        assert(2 == matrix_elements.size());
        assert(2 ==   norm_chebvecs.size());

        std::cout << " finished initialization of grow.h " << std::endl;


        int tmpN = jackson_coefs.size();
//      std::string saveString = (stoR+"_matrix_elements_"+".csv");
//      std::ofstream fout;
//
//       fout.open(saveString); for (int i = 0; i < tmpN; ++i) fout << matrix_elements[i] << std::endl;
//       fout.close();
//
//       saveString = (stoR+"_norm_chebvecs_"+".csv");
//       fout.open(saveString); for (int i = 0; i < tmpN; ++i) fout << norm_chebvecs[i] << std::endl;
//       fout.close();
    } else {

        std::string saveString;
        int tmpL = mps_cadr.size();

        liouvillian1.resize   (tmpL);
        liouvillian2.resize   (tmpL);
        mps_car.resize        (tmpL);
        mps_cadr.resize       (tmpL);
        observable_mpo.resize (tmpL);
        t0.resize             (tmpL);

        assert(tmpL == liouvillian2.size());
        assert(tmpL == mps_car.size()); 
        assert(tmpL == mps_cadr.size()); 
        assert(tmpL == t0.size() );
        assert(tmpL == observable_mpo.size() );

        std::ifstream fin;
        for(int i = 0; i < tmpL; ++i){
            saveString = stoR + "_mpscar_" + properfilename(0,0,i);         mps_car[i]        = UniTensor(saveString);
            saveString = stoR + "_mpscadr_" + properfilename(0,0,i);        mps_cadr[i]       = UniTensor(saveString);
            saveString = stoR + "_t0_" + properfilename(0,0,i);             t0[i]             = UniTensor(saveString);
            saveString = stoR + "_liouvillian1_" + properfilename(0,0,i);   liouvillian1[i]   = UniTensor(saveString);
            saveString = stoR + "_liouvillian2_" + properfilename(0,0,i);   liouvillian2[i]   = UniTensor(saveString);
            saveString = stoR + "_observable_mpo_" + properfilename(0,0,i); observable_mpo[i] = UniTensor(saveString); 
        }

        int tmpN = jackson_coefs.size();

        matrix_elements.resize(tmpN);
        norm_chebvecs.resize  (tmpN);

        assert(tmpN == matrix_elements.size());
        assert(tmpN ==   norm_chebvecs.size());

      // saveString = (stoR+"_matrix_elements_"+".csv");
      // fin.open(saveString); for (int i = 0; i < tmpN; ++i) fin >> matrix_elements[i];
      // fin.close();

      // saveString = (stoR+"_norm_chebvecs_"+".csv");
      // fin.open(saveString); for (int i = 0; i < tmpN; ++i) fin >> norm_chebvecs[i];
      // fin.close();

    }

    std::cout << " stem::stem grow_full... " << std::endl;
    grow_full();
}

/*!
 * After having computed |t_n> we must compute S(k,om) for many values
 * of om. This is the goal of this function.
 */
double
    stem::compute
    ( double omega_original )
{
    double omegaPrime = 2.0*Wprime/W * omega_original;
    double sum = jackson_coefs[0] * matrix_elements[0];
    for(int i = 1; i < max_height; ++i)
        sum += 2.0 * jackson_coefs[i] * matrix_elements[i] * cheb(i,omegaPrime);
    return (2.0*Wprime/W)/(PI * sqrt( 1-omegaPrime*omegaPrime ) ) * sum;
}
void
    stem::save
    (bool firstTime)
{
    /*
     *  Save mps_car, mps_cadr, current matrix_elements (can use processChebVecs script to recover Szz( om ) )
     *  Think of anything else that may be needed to recover the current state and keep going on another system.
     */
    // jackson_coefs; matrix_elements; norm_chebvecs; mps_car; mps_cadr; t0; liouvillian1; liouvillian2; observable_mpo;
 //  std::string saveString;
 //  std::string extra_file = "storage/";
 //  UniTensor tmp10;

 //  for(int i = 0; i < mps_car.size(); ++i){
 //      saveString = extra_file + stoR + "_mpscar_" + properfilename(0,0,i);
 //      tmp10 = mps_car[i];
 //      tmp10.setName (saveString);
 //      tmp10.save    (saveString);
 //  }

 //  for(int i = 0; i < mps_cadr.size(); ++i){
 //      saveString = extra_file + stoR + "_mpscadr_" + properfilename(0,0,i);
 //      tmp10 = mps_cadr[i];
 //      tmp10.setName (saveString);
 //      tmp10.save    (saveString);
 //  }

 //  if (firstTime){

 //      for(int i = 0; i < t0.size(); ++i){
 //          saveString = extra_file + stoR + "_t0_" + properfilename(0,0,i);
 //          tmp10 = t0[i];
 //          tmp10.setName (saveString);
 //          tmp10.save    (saveString);
 //      }


 //      for(int i = 0; i < liouvillian1.size(); ++i){
 //          saveString = extra_file + stoR + "_liouvillian1_" + properfilename(0,0,i);
 //          tmp10 = liouvillian1[i];
 //          tmp10.setName (saveString);
 //          tmp10.save    (saveString);
 //      }

 //      for(int i = 0; i < liouvillian1.size(); ++i){
 //          saveString = extra_file + stoR + "_liouvillian2_" + properfilename(0,0,i);
 //          tmp10 = liouvillian2[i];
 //          tmp10.setName (saveString);
 //          tmp10.save    (saveString);
 //      }

 //      for(int i = 0; i < observable_mpo.size(); ++i){
 //          saveString = extra_file + stoR + "_observable_mpo_" + properfilename(0,0,i);
 //          tmp10 = observable_mpo[i];
 //          tmp10.setName (saveString);
 //          tmp10.save    (saveString);
 //      }

 //  }

    csvfile g2(stoR+"_matrix_elements_"+  doublestr(beta)+"_"+doublestr(k_vector)+".csv"," "); g2.precision(20);
    for (int i = 0; i < matrix_elements.size(); ++i) g2 << matrix_elements[i] << endrow;

    csvfile g3(stoR+"_norm_chebvecs_"+doublestr(beta)+"_"+doublestr(k_vector)+".csv"," "); g3.precision(20);
    for (int i = 0; i < norm_chebvecs.size(); ++i)   g3 << norm_chebvecs[i] << endrow;

}

/*!
 * Compute one recursion step in the Chebyshev approximation.
 */
void 
    stem::grow_once
    ()
{
    std::vector<UniTensor> tmp = 
        compress_sum_mps( 
                 std::vector< std::vector<UniTensor> >( 
                        {
                    apply_mpo_to_mps( liouvillian1 , mps_car, "2*L|t_n-1>" ),
                    apply_mpo_to_mps( liouvillian2 , mps_car, "2*L|t_n-1>" ),
                    mps_times( -1.0, mps_cadr )
                        } ), 
                 maxBondDim, 
                 minBondDim,
                 epsilon
                 , true // svd
                 , num_sweeps
                 );
    //                      std::cout << "stem::grow_once: " << "finished compress_sum_mps" << std::endl;

    // re-adjust car and cadr
    mps_cadr = mps_car;                                   
    mps_car  = tmp;     //normalizeMPS(tmp);

    matrix_elements.push_back( inner_product( t0,      mps_car ) );
    norm_chebvecs.push_back  ( inner_product( mps_car, mps_car ) );
    //                      std::cout << "stem::grow_once: " << "finished matrix_elements" << std::endl;
}

/*!
 * Computes all recursive steps.
 */
void 
    stem::grow_full
    ()
{
    int start = matrix_elements.size(); 
    if(first_time){
        assert(start == 2);  // first defined in stem::stem
        save(first_time);
    }

    std::cout << "grow_full: (int)start          = " << start      << std::endl;
    std::cout << "grow_full: (int)maxCheb Height = " << max_height << std::endl;
    for(int i = start; i < max_height; ++i){
        grow_once();
        std::cout << "finished growth " << i << std::endl;

        if (i % 10 == 0){
            std::cout << "Printing Matrix Elements" << std::endl;
            csvfile g2(stoR+"_matrix_elements_"+  doublestr(beta)+"_"+doublestr(k_vector)+".csv"," "); g2.precision(20);
            g2 << "num" << "matelem" << endrow;
            for(int k = 0; k < matrix_elements.size(); ++k){
                std::cout << k << ", " << matrix_elements[k] << std::endl;
                g2 << k << matrix_elements[k] << endrow;
            }
            if (i%100)
                save(false); //first_time==false
        }

    }
    //S_of_omega();
}

#endif

