#include <iostream>
#include <Eigen>
#include "embryo_sim.hpp"


int main()
{
    Embryo em = Embryo<100,100>();

    em.embryo = em.embryo.Random();
    
    for (int ii = 0; ii < 100; ii++)
    {
        for (int jj = 0; jj < 100; jj++)
        {
            if (em.embryo(ii,jj) > 1)
            {
                em.embryo(ii,jj) = 1;
            }

            if (em.embryo(ii,jj) < 0)
            {
                em.embryo(ii,jj) = 0;
            }

            std::round(em.embryo(ii,jj));
        }
    }

    em.generateImage();

    em.img.save_image("test.bmp");

    return 0;
}
