#include <iostream>
#include <Eigen>
#include <random>
#include "embryo_sim.hpp"
#include "bitmap_image.hpp"
#include <iostream>

int main()
{
    srand(124);
    bitmap_image target_bmp = bitmap_image("shape.bmp");

    unsigned char r, g, b;

    Eigen::Matrix<float,105,100> target;
    target.setZero();
    for (int ii = 0; ii < 100; ii++)
    {
        for (int jj = 0; jj < 100; jj++)
        {
            target_bmp.get_pixel(ii,jj,r,g,b);
            if (r<126)
            {
                target(ii,jj) = 1;
            }
        }
    }

    Embryo<3,3> A;
    Embryo<3,3> B;
    Embryo<3,3> C;

    A.left_weights.setOnes();
    A.right_weights.setOnes();
    A.up_weights.setOnes();
    A.down_weights.setOnes();
    A.identity_weights.setOnes();

    B.identity_weights.setZero();

    splice(C,A,B);

    std::cout << "down_weights:\n" << C.down_weights << std::endl;
    std::cout << "up_weights:\n" << C.up_weights << std::endl;
    std::cout << "left_weights:\n" << C.left_weights << std::endl;
    std::cout << "right_weights:\n" << C.right_weights << std::endl;
    std::cout << "identity_weights:\n" << C.identity_weights << std::endl;


}