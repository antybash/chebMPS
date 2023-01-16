
#ifndef BIGVEC_H
#define BIGVEC_H

#include <uni10.hpp>
#include <cstdio>
#include <vector>
#include <string>

using namespace uni10;
using namespace useful;

class BigVec {
    public:
        BigVec(int,std::vector<std::string>);
        ~BigVec();
        UniTensor getElem(int);
        UniTensor getElem(int,bool);
        void getElem(int,UniTensor&);
        void getElem(int,UniTensor&,bool);
        void setElem(int,UniTensor&);
        std::string getFilename(int);
        int size();
        void fill(std::vector<UniTensor>&);
        std::vector<UniTensor>  mps_form();
    private:
        int Len;
        std::vector<std::string> filenames;
        std::vector<bool> existence_of_tensor;
};


BigVec::BigVec(int L, std::vector<std::string> fns)
{
    Len = L;
    filenames = fns;
    existence_of_tensor = std::vector<bool>(Len,false);
}

BigVec::~BigVec()
{
    for(int k = 0; k < Len; ++k){
        if (existence_of_tensor[k]){
            remove(filenames[k].c_str());
            existence_of_tensor[k] = false;
        }
    }
}

UniTensor BigVec::getElem(int k)
{
    assert(existence_of_tensor[k]);
    return UniTensor(filenames[k]);
}
UniTensor BigVec::getElem(int k, bool OVERride)
{
    if(!OVERride) assert(existence_of_tensor[k]);
    return UniTensor(filenames[k]);
}

void BigVec::getElem(int k, UniTensor& t)
{
    assert(existence_of_tensor[k]);
    t = UniTensor(filenames[k]);
}
void BigVec::getElem(int k, UniTensor& t, bool OVERride)
{
    if(!OVERride) assert(existence_of_tensor[k]);
    t = UniTensor(filenames[k]);
}

void BigVec::setElem(int k, UniTensor& t)
{
    t.save(filenames[k]);
    existence_of_tensor[k] = true;
}

std::string BigVec::getFilename(int k)
{
    return filenames[k];
}

int BigVec::size()
{
    return Len;
}

        void fill(std::vector<UniTensor>&);
        std::vector<UniTensor>  mps_form();

void BigVec::fill(std::vector<UniTensor> &mps)
{
    assert(mps.size() == Len);
    for (int i = 0; i < Len; ++i)
        setElem(i,mps[i]);
}

std::vector<UniTensor> BigVec::mps_form()
{
    std::vector<UniTensor> mps(Len);
    for(int i = 0; i < Len; ++i){
        assert(existence_of_tensor[i]);
        getElem(i,mps[i]);
    }
}

/////////////////
//

void convertVec2BigVec(std::vector<UniTensor> &mps, BigVec &bv)
{
    for(int i = 0; i < mps.size(); ++i)
        bv.setElem(i,mps[i]);
}

#endif
