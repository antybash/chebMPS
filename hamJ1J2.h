
#ifndef HAM_J1J2_H
#define HAM_J1J2_H

/*************************
 * J1--J2  MPO + AUX MPO *
 ************************/

#include "helper_compress2.h"
#include "hamBASE.h"

class hamJ1J2 : public HamiltonianClass {
    private: 
        int fullDim;  // twice physical length == thermal length
        int intL;
        double lam;
        double aniso; 

        UniTensor mpo_leftmost(double lam, double pre)
        {
            return MPO2add({ 
                             MPO1make(matBasis(0,1,1,8),     pre * sX()),
                             MPO1make(matBasis(0,2,1,8),     lam * sX()),
                             MPO1make(matBasis(0,3,1,8),     pre * sYn()),
                             MPO1make(matBasis(0,4,1,8),     lam * sYn()),
                             MPO1make(matBasis(0,5,1,8),     pre * sZ()),
                             MPO1make(matBasis(0,6,1,8),     lam * sZ()),
                             MPO1make(matBasis(0,7,1,8),           sI()) });
        }

        UniTensor mpo_mid(double lam, double   pre) 
        {
            return MPO2add({ MPO1make(matBasis(0,0,8,8),         sI()    ),
                             MPO1make(matBasis(1,0,8,8),         sX()    ),
                             MPO1make(matBasis(2,1,8,8),         sI()    ),
                             MPO1make(matBasis(3,0,8,8),         sYp()   ),
                             MPO1make(matBasis(4,3,8,8),         sI()    ),
                             MPO1make(matBasis(5,0,8,8),         sZ()    ),
                             MPO1make(matBasis(6,5,8,8),         sI()    ),
                             /********************************************/
                             MPO1make(matBasis(7,1,8,8),   pre * sX()    ),
                             MPO1make(matBasis(7,2,8,8),   lam * sX()    ),
                             MPO1make(matBasis(7,3,8,8),   pre * sYn()   ),
                             MPO1make(matBasis(7,4,8,8),   lam * sYn()   ),
                             MPO1make(matBasis(7,5,8,8),   pre * sZ()    ),
                             MPO1make(matBasis(7,6,8,8),   lam * sZ()    ),
                             MPO1make(matBasis(7,7,8,8),         sI()    ) });
        }

        UniTensor mpo_rightmost(double lam, double   pre)
        {
            return MPO2add({ 
                             MPO1make(matBasis(0,0,8,1),      sI()   ),
                             MPO1make(matBasis(1,0,8,1),      sX()   ),
                             MPO1make(matBasis(3,0,8,1),      sYp()  ),
                             MPO1make(matBasis(5,0,8,1),      sZ()   )});
        }

        UniTensor mpo_aux_mid() 
        {
            return MPO2add({ MPO1make(matBasis(0,0,8,8),sI()       ),
                             MPO1make(matBasis(1,1,8,8),sI()       ),
                             MPO1make(matBasis(2,2,8,8),sI()       ),
                             MPO1make(matBasis(3,3,8,8),sI()       ),
                             MPO1make(matBasis(4,4,8,8),sI()       ),
                             MPO1make(matBasis(5,5,8,8),sI()       ),
                             MPO1make(matBasis(6,6,8,8),sI()       ),
                             MPO1make(matBasis(7,7,8,8),sI()       ) });
        }

        UniTensor mpo_aux_rightmost() 
        {
           return  MPO1make( matBasis(0,0,8,1),sI() ) ;
        }

        UniTensor mpo_aux_leftmost_negative()  // check this !!!
        {
            return MPO2add({ MPO1make(matBasis(0,7,1,8),sIn())     });
        }

        /*****************
        * Liouv J1J2 MPO *
        *****************/

        std::vector< UniTensor > HamIden () 
        {
            std::vector< UniTensor > mpo(fullDim);
            double prefactor = (0%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo[0] = mpo_leftmost(lam,prefactor);
            mpo[1] = mpo_aux_mid();
            for(int i = 1; i < intL-1; ++i){
                prefactor = (i%2==1 ? 1.0-aniso : 1.0+aniso);
                mpo[2*i]   = mpo_mid(lam,prefactor); 
                mpo[2*i+1] = mpo_aux_mid();
            }
            prefactor = ( (intL-1)%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo[2*intL-2] = mpo_mid(lam,prefactor);
            mpo[2*intL-1] = mpo_aux_rightmost();
            return mpo;
        }

        std::vector< UniTensor > negIdenHam () 
        {
            std::vector< UniTensor > mpo(fullDim);
            double prefactor = (0%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo[0] = mpo_aux_leftmost_negative();
            mpo[1] = mpo_mid(lam,prefactor);
            for(int i = 1; i < intL-1; ++i){
                prefactor = (i%2==1 ? 1.0-aniso : 1.0+aniso);
                mpo[2*i]   = mpo_aux_mid();
                mpo[2*i+1] = mpo_mid(lam,prefactor);
            }
            prefactor = ( (intL-1)%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo[2*intL-2] = mpo_aux_mid();
            mpo[2*intL-1] = mpo_rightmost(lam,prefactor);
            return mpo;
        }

    public:

        hamJ1J2(double len, double lambda, double anisotropy)
        {
            intL    = len;
            fullDim = 2*intL;
            lam     = lambda;
            aniso   = anisotropy;
        }

        std::vector<UniTensor> thermalMPS   () override
        {
            std::vector<UniTensor> mps;
            for(int i = 0; i < intL; ++i){
                mps.push_back( t0() ); 
                mps.push_back( t1() ); 
            }
            return mps;
        }

        std::vector<UniTensor> thermalMPO   () override
        {
            std::vector<UniTensor> mpo;
            double prefactor = (0%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo.push_back( mpo_leftmost(lam,prefactor) );
            mpo.push_back( mpo_aux_mid ()              ); 
            for(int i = 1; i <= intL-2; ++i){
                prefactor = (i%2==1 ? 1.0-aniso : 1.0+aniso);
                mpo.push_back( mpo_mid     ( lam,prefactor ) );
                mpo.push_back( mpo_aux_mid ()                );
            }
            prefactor = (intL%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo.push_back( mpo_mid           (lam,prefactor) );
            mpo.push_back( mpo_aux_rightmost () );
            return mpo;
        }

        std::vector<UniTensor> zerotmpMPO   () override
        {
            std::vector<UniTensor> mpo;
            double prefactor = (0%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo.push_back( mpo_leftmost(lam,prefactor) );
            for(int i = 1; i <= intL-2; ++i){
                prefactor = (i%2==1 ? 1.0-aniso : 1.0+aniso);
                mpo.push_back( mpo_mid   ( lam,prefactor ) );
            }
            prefactor = (intL%2==1 ? 1.0-aniso : 1.0+aniso);
            mpo.push_back( mpo_rightmost ( lam,prefactor ) );
            return mpo;
        }

        std::vector<UniTensor> liouvillian1 (double bandwidth) override
        {
            // prefactor = 2*  (2*Wprime/W)
            // return mpo_times( 2.0* (2.0*Wprime/W),  hamiltonian_class ->    HamIden());
            return mpo_times( bandwidth, HamIden() );
        }

        std::vector<UniTensor> liouvillian2 (double bandwidth) override
        {
            return mpo_times( bandwidth, negIdenHam() );
        }

        std::vector<UniTensor> ham_mpo_spin_str (int physLen, double kvector)
        {
            return  mpo_spin_str(physLen,kvector); // spin-1/2 
        }


}; 

#endif

