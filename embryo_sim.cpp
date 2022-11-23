#include <iostream>
#include <Eigen>
#include "embryo_sim.hpp"





int main()
{
    Embryo em = Embryo(100,100);

    em.embryo_vector.Random();
    
    for (int ii = 0; ii < 100*100; ii++)
    {
        if (em.embryo_vector[ii] > 1)
        {
            em.embryo_vector[ii] = 1;
        }

        if (em.embryo_vector[ii] < 0)
        {
            em.embryo_vector[ii] = 0;
        }

        std::round(em.embryo_vector[ii]);
    }

    em.generateImage();

    em.img.save_image("test.bmp");

    return 0;
}
