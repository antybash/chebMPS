
#ifndef HAM_XXZEXT_H
#define HAM_XXZEXT_H

/***************************************************
 * External Sz Field with Transverse Component MPO *
 * Model: XXZ with Delta = 1, Alternating Sz field with magnitude hz=1, and transverse component Sx
 **************************************************/
#include "helper_compress2.h"
#include "hamBASE.h"

class hamXXZext : public HamiltonianClass {
    private:
        int fullDim;  // twice physical length == thermal length
        int intL;
        double lam;

        UniTensor mpo_leftmost(double pre)
        {
            return MPO2add({ MPO1make(matBasis(0,0,1,5),sZ()+pre*sX()), 
                             MPO1make(matBasis(0,1,1,5),sX()         ),
                             MPO1make(matBasis(0,2,1,5),sYn()        ),
                             MPO1make(matBasis(0,3,1,5),sZ()         ),
                             MPO1make(matBasis(0,4,1,5),sI())       });
        }

        UniTensor mpo_mid(double pre)
        {
            return MPO2add({ MPO1make(matBasis(0,0,5,5),sI()          ),
                             MPO1make(matBasis(1,0,5,5),sX()          ),
                             MPO1make(matBasis(2,0,5,5),sYp()         ),
                             MPO1make(matBasis(3,0,5,5),sZ()          ),
                             MPO1make(matBasis(4,0,5,5),sZ()+pre*sX() ),
                             MPO1make(matBasis(4,1,5,5),sX()          ),
                             MPO1make(matBasis(4,2,5,5),sYn()         ),
                             MPO1make(matBasis(4,3,5,5),sZ()          ),
                             MPO1make(matBasis(4,4,5,5),sI())        });
        }

        UniTensor mpo_rightmost(double pre)
        {
            return MPO2add({ MPO1make(matBasis(0,0,5,1),sI()             ),
                             MPO1make(matBasis(1,0,5,1),sX()             ),
                             MPO1make(matBasis(2,0,5,1),sYp()            ),
                             MPO1make(matBasis(3,0,5,1),sZ()             ),
                             MPO1make(matBasis(4,0,5,1),sZ()+pre*sX())  });
        }

        UniTensor mpo_aux_leftmost_negative()
        {
            return MPO2add({ MPO1make(matBasis(0,4,1,5),sIn())     });
        }

        UniTensor mpo_aux_leftmost()
        {
            return MPO2add({ MPO1make(matBasis(0,4,1,5),sI())      });
        }


        UniTensor mpo_aux_mid()
        {
            return MPO2add({ MPO1make(matBasis(0,0,5,5),sI()       ),
                             MPO1make(matBasis(1,1,5,5),sI()       ),
                             MPO1make(matBasis(2,2,5,5),sI()       ),
                             MPO1make(matBasis(3,3,5,5),sI()       ),
                             MPO1make(matBasis(4,4,5,5),sI()       ) });
        }

        UniTensor mpo_aux_rightmost()
        {
           return  MPO1make( matBasis(0,0,5,1),sI() ) ;
        }

        /****************
        * Liouv XXZ MPO *
        ****************/

        std::vector< UniTensor > HamIden ()
        {
            std::vector< UniTensor > mpo(fullDim);
            mpo[0] = mpo_leftmost(lam);
            mpo[1] = mpo_aux_mid();
            for(int i = 1; i < intL-1; ++i){
                double prefactor = (i%2==1 ? -1.0 : 1.0);
                mpo[2*i]   = mpo_mid(prefactor*lam);
                mpo[2*i+1] = mpo_aux_mid();
            }
            double prefactor = ( (intL-1)%2==1 ? -1.0 : 1.0 );
            mpo[2*(intL-1)]  = mpo_mid(prefactor*lam);
            mpo[2*intL-1]    = mpo_aux_rightmost();
            return mpo;
        }

        std::vector< UniTensor > negIdenHam () 
        {
            std::vector< UniTensor > mpo(fullDim);
            mpo[0] = mpo_aux_leftmost_negative();
            mpo[1] = mpo_mid((-1.0)*lam);                           // necessary for TIME REVERSAL: cf. Tiegel thesis
            for(int i = 1; i < intL-1; ++i){
                double prefactor = (i%2==1 ? -1.0 : 1.0);
                mpo[2*i]   = mpo_aux_mid();
                mpo[2*i+1] = mpo_mid((-1.0)*prefactor*lam);         // necessary for TIME REVERSAL: cf. Tiegel thesis
            }
            double prefactor = ( (intL-1)%2==1 ? -1.0 : 1.0 );
            mpo[2*(intL-1)] = mpo_aux_mid();
            mpo[2*intL-1]   = mpo_rightmost((-1.0)*prefactor*lam);  // necessary for TIME REVERSAL: cf. Tiegel thesis
            return mpo;
        }


    public:
        hamXXZext(int len, double lambda)
        {
            intL    = len;
            fullDim = intL*2;
            lam     = lambda;
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

        std::vector<UniTensor> thermalMPO () override
        {
            std::vector<UniTensor> mpo;
            mpo.push_back( mpo_leftmost(1.0* lam) );
            mpo.push_back( mpo_aux_mid ()         );
            for(int i = 1; i <= intL-2; ++i){
                double prefactor = (i%2==1 ? -1.0 : 1.0);
                mpo.push_back( mpo_mid     (prefactor * lam) );
                mpo.push_back( mpo_aux_mid ()                );
            }
            double prefactor = ((intL-1)%2==1 ? -1.0 : 1.0);
            mpo.push_back( mpo_mid           (prefactor * lam) );
            mpo.push_back( mpo_aux_rightmost ()                );
            return mpo;
        }

        std::vector<UniTensor> zerotmpMPO () override
        {
            std::vector<UniTensor> mpo;
            mpo.push_back( mpo_leftmost(1.0* lam) );
            for(int i = 1; i <= intL-2; ++i){
                double prefactor = (i%2==1 ? -1.0 : 1.0);
                mpo.push_back( mpo_mid   (prefactor * lam) );
            }
            double prefactor = ((intL-1)%2==1 ? -1.0 : 1.0);
            mpo.push_back( mpo_rightmost (prefactor * lam) );
            return mpo;
        }


        std::vector<UniTensor> liouvillian1 (double bandwidth) override
        {
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
