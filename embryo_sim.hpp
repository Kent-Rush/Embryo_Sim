#include <Eigen>
#include "bitmap_image.hpp"

using Eigen::Matrix;
using Eigen::Vector;

template<int ROWS, int COLS>
class Embryo
{
public:
    Eigen::Matrix<float,ROWS,COLS> embryo;
    Eigen::Matrix<float,ROWS,COLS> identity_weights;
    Eigen::Matrix<float,COLS-1,ROWS> left_weights;
    Eigen::Matrix<float,COLS-1,ROWS> right_weights;
    Eigen::Matrix<float,COLS,ROWS-1> up_weights;
    Eigen::Matrix<float,COLS,ROWS-1> down_weights;

    //Each submatrix is COLSxCOLS
    //AND there are ROWSxROWS submatrices
    /****
    I R 0 0 | U 0 0 0
    L I R 0 | 0 U 0 0
    0 L I R | 0 0 U 0
    0 0 L I | 0 0 0 U
    -----------------
    D 0 0 0 | I R 0 0
    0 D 0 0 | L I R 0
    0 0 D 0 | 0 L I R
    0 0 0 D | 0 0 L I
    ****/

    bitmap_image img;

    Embryo():
        img(ROWS,COLS)
    {
        identity_weights.setOnes();
        left_weights.setZero();
        right_weights.setZero();
        up_weights.setZero();
        down_weights.setZero();
    }

    void step()
    {
        Eigen::Matrix<float,ROWS,COLS> new_embryo;
        new_embryo.setZero();

        for (int ii = 0; ii < ROWS; ii++)
        {
            for (int jj = 0; jj < COLS; jj++)
            {
                new_embryo(ii,jj) += identity_weights(ii,jj)*embryo(ii,jj);
                if (jj != 0)
                {
                    new_embryo(ii,jj-1) += left_weights(ii-1,jj)*embryo(ii,jj-1);
                }
                if (jj != COLS-1)
                {
                    new_embryo(ii,jj) += right_weights(ii-1,jj)*embryo(ii,jj+1);
                }
                if (ii != 0)
                {
                    new_embryo(ii,jj) += up_weights(ii-1,jj)*embryo(ii-1,jj);
                }
                if (ii != ROWS-1)
                {
                    new_embryo(ii,jj) += up_weights(ii,jj)*embryo(ii+1,jj);
                }
            }
        }

        embryo = new_embryo;
    }

    void generateImage()
    {
        for (int ii = 0; ii < ROWS; ii++)
        {
            for (int jj = 0; jj < COLS; jj++)
            {
                if (std::abs(embryo(ii,jj) < 1e-6 ))
                {
                    img.set_pixel(ii, jj,255,255,255);
                }
                else
                {
                    img.set_pixel(ii, jj,0,0,0);
                }
            }
        }
    }
};