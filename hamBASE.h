
#ifndef HAM_BASE_H
#define HAM_BASE_H

#include <uni10.hpp>

using namespace useful;
using namespace uni10;

class HamiltonianClass {
    public:

        virtual std::vector<UniTensor> thermalMPS   ()               = 0;
        virtual std::vector<UniTensor> thermalMPO   ()               = 0;
        virtual std::vector<UniTensor> zerotmpMPO   ()               = 0;
        virtual std::vector<UniTensor> liouvillian1 (double)         = 0;
        virtual std::vector<UniTensor> liouvillian2 (double)         = 0;
        virtual std::vector<UniTensor> ham_mpo_spin_str (int,double) = 0;

};


#endif
