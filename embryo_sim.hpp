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
    Eigen::Matrix<float,ROWS,COLS-1> left_weights;
    Eigen::Matrix<float,ROWS,COLS-1> right_weights;
    Eigen::Matrix<float,ROWS-1,COLS> up_weights;
    Eigen::Matrix<float,ROWS-1,COLS> down_weights;
    float mutation_rate;

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

    

    Embryo():
        mutation_rate(1e-2f)
    {
        identity_weights.setOnes();
        left_weights.setZero();
        right_weights.setZero();
        up_weights.setZero();
        down_weights.setZero();
    }

    void randomize(float scale)
    {
        identity_weights = identity_weights.Random()*scale;
        left_weights = left_weights.Random()*scale;
        right_weights = right_weights.Random()*scale;
        up_weights = up_weights.Random()*scale;
        down_weights = down_weights.Random()*scale;
    }

    static void mutate(Embryo& embryo_out, Embryo& embryo_in)
    {
        float k = 1e-3;
        embryo_out.identity_weights = embryo_in.identity_weights + embryo_in.identity_weights.Random()*k;
        embryo_out.left_weights = embryo_in.left_weights + embryo_in.left_weights.Random()*k;
        embryo_out.right_weights = embryo_in.right_weights + embryo_in.right_weights.Random()*k;
        embryo_out.up_weights = embryo_in.up_weights + embryo_in.up_weights.Random()*k;
        embryo_out.down_weights = embryo_in.down_weights + embryo_in.down_weights.Random()*k;
    }

    void round()
    {
        for (int ii = 0; ii < ROWS; ii++)
        {
            for (int jj = 0; jj < COLS; jj++)
            {
                if (embryo(ii,jj) > 1)
                {
                    embryo(ii,jj) = 1;
                }

                if (embryo(ii,jj) < 0)
                {
                    embryo(ii,jj) = 0;
                }

                embryo(ii,jj) = std::round(embryo(ii,jj));
            }
        }
    }

    void step()
    {
        Eigen::Matrix<float,ROWS,COLS> new_embryo;
        new_embryo.setZero();

        for (int ii = 0; ii < ROWS; ii++)
        {
            for (int jj = 0; jj < COLS; jj++)
            {
                // bool none_adjacent = embryo(ii,jj) > 0.5;
                // if (ii != 0)
                // {
                //     none_adjacent = embryo(ii-1,jj) > 0.5;
                // }
                // if (ii != ROWS-1)
                // {
                //     none_adjacent = embryo(ii+1,jj) > 0.5;
                // }
                // if (jj != 0)
                // {
                //     none_adjacent = embryo(ii,jj-1) > 0.5;
                // }
                // if (jj != COLS-1)
                // {
                //     none_adjacent = embryo(ii,jj+1) > 0.5;
                // }

                // if (none_adjacent)
                // {
                //     continue;
                // }

                new_embryo(ii,jj) += identity_weights(ii,jj)*embryo(ii,jj);
                if (jj != 0)
                {
                    new_embryo(ii,jj) += left_weights(ii,jj-1)*embryo(ii,jj-1);
                }
                if (jj != COLS-1)
                {
                    new_embryo(ii,jj) += right_weights(ii,jj)*embryo(ii,jj+1);
                }
                if (ii != 0)
                {
                    new_embryo(ii,jj) += up_weights(ii-1,jj)*embryo(ii-1,jj);
                }
                if (ii != ROWS-1)
                {
                    new_embryo(ii,jj) += down_weights(ii,jj)*embryo(ii+1,jj);
                }
            }
        }

        embryo = new_embryo;
        round();
    }

    void generateImage(bitmap_image& img)
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

template<int M, int N>
void slice_mat(Eigen::Matrix<float,M,N>& out, const Eigen::Matrix<float,M,N>& in_a, const Eigen::Matrix<float,M,N>& in_b)
{
    std::random_device rd;  //Uniform random integer distribution
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib_row(0, M-1);
    std::uniform_int_distribution<size_t> distrib_col(0, N-1);
    size_t rowmax = distrib_row(gen);
    size_t colmax = distrib_col(gen);
    
    //Need to handle the case when rowmax or colmax is 0 because ::seq cant handle that

    if (rowmax != 0)
    {
        out(Eigen::seq(0,rowmax-1), Eigen::seq(0,N-1)) = in_a(Eigen::seq(0,rowmax-1), Eigen::seq(0,N-1));
    }
    
    if (colmax != 0)
    {
        out(rowmax, Eigen::seq(0,colmax-1)) = in_a(rowmax, Eigen::seq(0,colmax-1));
    }
    
    out(rowmax, Eigen::seq(colmax,N-1)) = in_b(rowmax, Eigen::seq(colmax,N-1));
    if (rowmax+1 != M)
    {
        out(Eigen::seq(rowmax+1,M-1), Eigen::seq(0,N-1)) = in_b(Eigen::seq(rowmax+1,M-1), Eigen::seq(0,N-1));
    }
}

template<int M, int N>
void splice(Embryo<M,N>& out, const Embryo<N,M>& em_a, const Embryo<N,M>& em_b)
{
    slice_mat(out.left_weights, em_a.left_weights, em_b.left_weights);
    slice_mat(out.right_weights, em_a.right_weights, em_b.right_weights);
    slice_mat(out.up_weights, em_a.up_weights, em_b.up_weights);
    slice_mat(out.down_weights, em_a.down_weights, em_b.down_weights);
    slice_mat(out.identity_weights, em_a.identity_weights, em_b.identity_weights);
}
