/* Random number class
*/

#ifndef RNDM_H
#define RNDM_H

#include "iostream"
#include "random"


template<typename Type = double >
class RNDM
{
    // easier to use param_type
    using param_type = typename std::uniform_real_distribution<Type>::param_type;

    // store an instance of the generator/distribution in the class object
    std::mt19937 gen;
    std::uniform_real_distribution<Type> dis_dunif;

public:
    // seed generator
    RNDM(): gen(std::random_device()()) {}

    Type dunif(Type from, Type to)
    {
        return dis_dunif(gen, param_type{from, to});
    }


};


#endif
