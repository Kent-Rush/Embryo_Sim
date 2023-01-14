#include <iostream>
#include <Eigen>
#include "embryo_sim.hpp"


int main()
{
    //Set rng seed
    srand(1245125);

    Embryo em = Embryo<100,100>();

    em.embryo = em.embryo.Random();
    em.round();

    Embryo mutant = Embryo<100,100>();

    for (int ii = 40; ii < 60; ii++)
    {
        for (int jj = 40; jj < 60; jj++)
        {
            mutant.embryo(ii,jj) = 1;
        }
    }

    //mutant.randomize(0.1);
    mutant.identity_weights.setZero();
    mutant.up_weights.setOnes();

    for (int ii = 0; ii < 100; ii++)
    {
        em.generateImage();
        em.img.save_image("frames/frame_"+std::to_string(ii)+".bmp");
        em.step();

        mutant.generateImage();
        mutant.img.save_image("frames2/frame_"+std::to_string(ii)+".bmp");
        mutant.step();
    }

    return 0;
}
