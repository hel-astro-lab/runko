/* Random number class
*/

#ifndef RNDM_HPP
#define RNDM_HPP

#include "iostream"
#include "random"
// #include "chrono"


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


    // double Maxwellian distribution with vb
    // Type double_maxwellian(Type vb)
    // {
    //     double fmax = 0.5*(1.0 + exp(-2.0*vb*vb));
    //     double vmin = -5.0*vb;
    //     double vmax =  5.0*vb;
    //     double vf = dis_unif(gen, param_type{vmin, vmax});
    //     double f = 0.5*(exp(-(vf-vb)*(vf-vb)/2) +
    //                     exp(-(vf+vb)*(vf+vb)/2));
    //     double x = fmax*dis_unif(gen, param_type{0.0, 1.0});
    //     
    //     if(x > f){
    //         return double_maxwellian(vb);
    //     } 
    //     return vf;
    // }

};


#endif
