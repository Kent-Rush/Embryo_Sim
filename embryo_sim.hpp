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

// template<int M, int N>
// void splice(Embryo<M,N>& out, const Embryo<N,M>& em_a, const Embryo<N,M>& em_b)
// {
//     const int lr_max = M*(N-1);
//     const int ud_max = (M-1)*N;
//     const int i_max = M*N;

//     size_t lr_index = Eigen::Rand::UniformIntGen<size_t>(0,M*(N-1));
//     size_t ud_index = Eigen::Rand::UniformIntGen<size_t>(0,(M-1)*N);
//     size_t i_index  = Eigen::Rand::UniformIntGen<size_t>(0,M*N);

//     //I forgot these are not vectors q_q

//     out.left_weights(Eigen::seq(0,lr_index)) = em_a.left_weights(Eigen::seq(0,lr_index));
//     out.left_weights(Eigen::seq(lr_index+1, lr_max-1)) = em_b.left_weights(Eigen::seq(lr_index+1, lr_max-1));
//     out.right_weights(Eigen::seq(0,lr_index)) = em_a.right_weights(Eigen::seq(0,lr_index));
//     out.right_weights(Eigen::seq(lr_index+1,lr_max-1)) = em_b.right_weights(Eigen::seq(lr_index+1,lr_max-1));

//     out.up_weights(Eigen::seq(0,ud_index)) = em_a.up_weights(Eigen::seq(0,ud_index));
//     out.up_weights(Eigen::seq(ud_index+1,ud_max-1)) = em_b.up_weights(Eigen::seq(ud_index+1,ud_max-1));
//     out.down_weights(Eigen::seq(0,ud_index)) = em_a.down_weights(Eigen::seq(0,ud_index));
//     out.down_weights(Eigen::seq(ud_index+1,ud_max-1)) = em_b.down_weights(Eigen::seq(ud_index+1,ud_max-1));

//     out.identity_weights(Eigen::seq(0,i_index)) = em_a.identity_weights(Eigen::seq(0,i_index));
//     out.identity_weights(Eigen::seq(i_index,i_max-1)) = emb.identity_weights(Eigen::seq(i_index,i_max-1));
// }