
#ifndef BIGSTACK_H
#define BIGSTACK_H

#include <uni10.hpp>
#include <vector>
#include <stack>
#include <string>
#include <cstdio>
#include "helper_compress2.h"

using namespace uni10;
using namespace useful;

class BigStack {
    public:
        BigStack(int,std::vector<std::string>);
        ~BigStack();
        UniTensor top();
        void push(UniTensor&);
        void pop();
        std::string getFilename(int);
        int size();
    private:
        int Len;
        int current;
        std::vector<std::string> filenames;
        std::vector<bool> existence_of_tensor;
};


BigStack::BigStack(int L, std::vector<std::string> fns)
{
    Len = L;
    filenames = fns;
    current = 0;
    existence_of_tensor = std::vector<bool>(Len,false);
}

BigStack::~BigStack()
{
    for(int k = 0; k < Len; ++k)
        if(existence_of_tensor[k]){
            remove(filenames[k].c_str());
            existence_of_tensor[k] = false;
        }
}

UniTensor BigStack::top()
{
    assert(current < Len);
    assert(existence_of_tensor[current]);
    return UniTensor(filenames[current]);
}

void BigStack::push(UniTensor &t)
{
    assert(current < Len);
    t.save(filenames[current]);
    existence_of_tensor[current] = true;
    current++;
}

void BigStack::pop()
{
    assert(current > 0);
    current--;
}


std::string BigStack::getFilename(int k)
{
    return filenames[k];
}

int BigStack::size()
{
    return current;
}

#endif

