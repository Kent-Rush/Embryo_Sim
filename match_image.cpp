#include <iostream>
#include <Eigen>
#include <random>
#include "embryo_sim.hpp"
#include "bitmap_image.hpp"
#include <iostream>
#include <iterator>

int main()
{
    srand(124);
    const int rows = 100;
    const int cols = 100;
    //typedef Embryo<rows,cols> embryo_t;

    bitmap_image target_bmp = bitmap_image("shape.bmp");
    bitmap_image img(rows,cols);


    unsigned char r, g, b;

    Eigen::Matrix<float,rows,cols> target;
    target.setZero();
    for (int ii = 0; ii < rows; ii++)
    {
        for (int jj = 0; jj < cols; jj++)
        {
            target_bmp.get_pixel(ii,jj,r,g,b);
            if (r<126)
            {
                target(ii,jj) = 1;
            }
        }
    }

    const size_t num_generations = 1;
    const size_t num_embryos = 10;
    const size_t num_iterations = 50;

    std::vector<Embryo<rows,cols>*> embryos;

    for (size_t ii = 0; ii < num_embryos; ii++)
    {
        embryos.push_back(new Embryo<rows,cols>());
        embryos[ii]->init_dot(6,50,50);
    }

    for (size_t ii = 0; ii < num_generations; ii++)
    {
        float scores[num_embryos];

        for (size_t jj = 0; jj < num_embryos; jj++)
        {
            embryos[jj]->mutate(1.0);

            for (size_t kk = 0; kk < num_iterations; kk++)
            {
                if (jj == 0)
                {
                    embryos[jj]->generateImage(img);
                    img.save_image("frames/frame_"+std::to_string(kk)+".bmp");
                }
                embryos[jj]->round();
                embryos[jj]->step();
                
            }
            

            scores[jj] = image_score<rows,cols>(embryos[jj]->embryo, target);
            embryos[jj]->generateImage(img);

            img.save_image("embryo_"+std::to_string(jj)+".bmp");

        }

        std::vector<float> score_vec(scores,scores+num_embryos);
        std::sort(score_vec.begin(), score_vec.end());
    }

    for (size_t ii = 0; ii < num_embryos; ii++)
    {
        delete embryos[ii];
    }

}

