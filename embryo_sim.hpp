#include <Eigen>
#include "bitmap_image.hpp"

using Eigen::Matrix;
using Eigen::Vector;

class Embryo
{
public:

    const int rows;
    const int columns;
    Eigen::VectorXd embryo_vector;
    Eigen::MatrixXd markov_matrix;
    int num_generations;
    int num_growth_cycles;
    int num_contestants;

    bitmap_image img;

    Embryo(const int rows_, const int columns_):
        rows(rows_),
        columns(columns_),
        img(rows_,columns_)
    {
        embryo_vector(rows_*columns_);
        markov_matrix(rows_*columns_, rows_*columns_);
        embryo_vector.setZero();
        markov_matrix.Identity();
    }

    void generateImage()
    {
        for (int ii = 0; ii < rows; ii++)
        {
            for (int jj = 0; jj < columns; jj++)
            {
                if (std::abs(embryo_vector[ii*columns+jj] - 1) < 1e-6 )
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