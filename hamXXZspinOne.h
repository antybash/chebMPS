
#ifndef HAM_XXZ_SPINONE_H
#define HAM_XXZ_SPINONE_H

#include "helper_compress2.h"
#include "hamBASE.h"

/***********
 * XXZ MPO *
 **********/

class hamXXZspinOne : public HamiltonianClass {

    private:
        int fullDim; // twice physical length == thermal length
        int intL;
        double lam;

        UniTensor mpo_leftmost() 
        {
            return MPO2add({ MPO1make(matBasis(0,1,1,5),spinOne::sX()),
                             MPO1make(matBasis(0,2,1,5),spinOne::sYn()),
                             MPO1make(matBasis(0,3,1,5),lam * spinOne::sZ()),
                             MPO1make(matBasis(0,4,1,5),spinOne::sI())       });
        }

        UniTensor mpo_mid() 
        {
            return MPO2add({ MPO1make(matBasis(0,0,5,5),spinOne::sI()       ),
                             MPO1make(matBasis(1,0,5,5),spinOne::sX()       ),
                             MPO1make(matBasis(2,0,5,5),spinOne::sYp()      ),
                             MPO1make(matBasis(3,0,5,5),spinOne::sZ()       ),
                             MPO1make(matBasis(4,1,5,5),spinOne::sX()       ),
                             MPO1make(matBasis(4,2,5,5),spinOne::sYn()      ),
                             MPO1make(matBasis(4,3,5,5),lam * spinOne::sZ() ),
                             MPO1make(matBasis(4,4,5,5),spinOne::sI()       ) });
        }

        UniTensor mpo_rightmost() 
        {
            return MPO2add({ MPO1make(matBasis(0,0,5,1),spinOne::sI()       ),
                             MPO1make(matBasis(1,0,5,1),spinOne::sX()       ),
                             MPO1make(matBasis(2,0,5,1),spinOne::sYp()      ),
                             MPO1make(matBasis(3,0,5,1),spinOne::sZ()       ) });
        }

        /**************
        * Aux XXZ MPO *
        **************/

        UniTensor mpo_aux_leftmost_negative() 
        {
            return MPO2add({ MPO1make(matBasis(0,4,1,5),spinOne::sIn())     });
        }

        UniTensor mpo_aux_leftmost() 
        {
            return MPO2add({ MPO1make(matBasis(0,4,1,5),spinOne::sI())      });
        }


        UniTensor mpo_aux_mid() 
        {
            return MPO2add({ MPO1make(matBasis(0,0,5,5),spinOne::sI()       ),
                             MPO1make(matBasis(1,1,5,5),spinOne::sI()       ),
                             MPO1make(matBasis(2,2,5,5),spinOne::sI()       ),
                             MPO1make(matBasis(3,3,5,5),spinOne::sI()       ),
                             MPO1make(matBasis(4,4,5,5),spinOne::sI()       ) });
        }

        UniTensor mpo_aux_rightmost() 
        {
            return  MPO1make( matBasis(0,0,5,1),spinOne::sI() ) ;
        }

        /****************
        * Liouv XXZ MPO *
        ****************/

        std::vector< UniTensor > HamIden ()
        {
            std::vector< UniTensor > mpo(fullDim);
            mpo[0] = mpo_leftmost();
            mpo[1] = mpo_aux_mid();
            for(int i = 1; i < intL-1; ++i){
                mpo[2*i] = mpo_mid();
                mpo[2*i+1] = mpo_aux_mid();
            }
            mpo[2*intL-2] = mpo_mid();
            mpo[2*intL-1] = mpo_aux_rightmost();
            return mpo;
        }

        std::vector< UniTensor > negIdenHam () 
        {
            std::vector< UniTensor > mpo(fullDim);
            mpo[0] = mpo_aux_leftmost_negative();
            mpo[1] = mpo_mid();
            for(int i = 1; i < intL-1; ++i){
                mpo[2*i] = mpo_aux_mid();
                mpo[2*i+1] = mpo_mid();
            }
            mpo[2*intL-2] = mpo_aux_mid();
            mpo[2*intL-1] = mpo_rightmost();
            return mpo;
        }

    public:

        hamXXZspinOne(double len, double lambda)
        {
            intL    = len;
            fullDim = 2*intL;
            lam     = lambda;
        }

        std::vector<UniTensor> thermalMPS   () override
        {
            std::vector<UniTensor> mps;
            for(int i = 0; i < intL; ++i){
                mps.push_back( spinOne::t0() ); 
                mps.push_back( spinOne::t1() ); 
            }
            return mps;
        }

        std::vector<UniTensor> thermalMPO   () override
        {
            std::vector<UniTensor> mpo;
            mpo.push_back( mpo_leftmost() );
            mpo.push_back( mpo_aux_mid () );
            for(int i = 1; i < intL-1; ++i){
                mpo.push_back( mpo_mid     () );
                mpo.push_back( mpo_aux_mid () );
            }
            mpo.push_back( mpo_mid           () );
            mpo.push_back( mpo_aux_rightmost () );
            return mpo;
        }

        std::vector<UniTensor> zerotmpMPO () override
        {
            std::vector<UniTensor> mpo;
            mpo.push_back( mpo_leftmost() );
            for(int i = 1; i < intL-1; ++i){
                mpo.push_back( mpo_mid   () );
            }
            mpo.push_back( mpo_rightmost () );
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
            return  SPINONE_mpo_spin_str(physLen,kvector); // spin-1
        }


};

#endif

