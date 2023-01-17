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

    const int rows = 30;
    const int cols = 100;
    Embryo mutant = Embryo<rows,cols>();

    for (int ii = 10; ii < 20; ii++)
    {
        for (int jj = 40; jj < 60; jj++)
        {
            mutant.embryo(ii,jj) = 1;
        }
    }


    bitmap_image img(100,100);
    bitmap_image img2(rows,cols);

    //mutant.randomize(0.1);
    mutant.identity_weights.setZero();
    mutant.up_weights.setOnes();

    for (int ii = 0; ii < 100; ii++)
    {
        em.generateImage(img);
        img.save_image("frames/frame_"+std::to_string(ii)+".bmp");
        em.step();

        mutant.generateImage(img2);
        img2.save_image("frames2/frame_"+std::to_string(ii)+".bmp");
        mutant.step();
    }

    return 0;
}
