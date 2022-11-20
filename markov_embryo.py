import numpy as np
import matplotlib.pyplot as plt
import pdb
import matplotlib.animation as animation

from embrio import Embryo


def unflatten(arr, n, m):
    output = np.zeros((n,m))

    for ii in range(n):
        output[ii,:] = arr[ii*m:(ii+1)*m]

    return output

def grow(arr, n, m, markov_mat):

    new_cells = []

    #local function to append the current point and adjacent points to a list
    #if they havent been already
    def append_maybe(index):

        row = np.floor(index/n)
        col = index%n

        #if we have already added this cell, dont do anything
        if row*m+col not in new_cells:
            new_cells.append(row*m+col)

        above = (row-1)*m + col
        below = (row+1)*m + col
        left  = row*m + col - 1
        right = row*m + col + 1

        if row == 0:
            if below not in new_cells:
                new_cells.append(below)
        elif row == n:
            if above not in new_cells:
                new_cells.append(above)
        else:
            if above not in new_cells:
                new_cells.append(above)
            if below not in new_cells:
                new_cells.append(below)

        if col == 0:
            if right not in new_cells:
                new_cells.append(right)
        elif col == m:
            if left not in new_cells:
                new_cells.append(left)
        else:
            if left not in new_cells:
                new_cells.append(left)
            if right not in new_cells:
                new_cells.append(right)

        
    arr2 = np.zeros(n*m)
    for ii in range(n*m):
        if arr[ii] == 1:
            append_maybe(ii)

            for cell in new_cells:
                arr2[int(cell)] = markov_mat[int(cell),:]@arr

    return arr2
                




def main():

    #grid size
    n = 100
    m = 100

    num_generations = 1

    num_growth_cycles = 20

    markov_matrix = np.random.rand(n*m,n*m)*2 - 1

    embryo = np.zeros(n*m)

    embryo[5050] = 1
    embryo[5051] = 1
    embryo[4950] = 1
    embryo[5150] = 1

    for ii in range(num_generations):

        animation_frames = np.zeros((num_growth_cycles, n*m))

        for jj in range(num_growth_cycles):


            print(jj)
            embryo = grow(embryo, n, m, markov_matrix)
            embryo = np.clip(embryo, 0, 1)
            embryo = np.round(embryo)

            animation_frames[jj] = 255.0*embryo

        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)

        def animate(index):
            ax1.clear()
            ax1.imshow(unflatten(animation_frames[index,:],n,m))
            ax1.set_title(index)

        ani = animation.FuncAnimation(fig, animate, interval=1000)
        plt.show()

            






if __name__ == "__main__":
    main()
