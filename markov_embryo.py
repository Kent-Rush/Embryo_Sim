import numpy as np
import matplotlib.pyplot as plt
import pdb
import matplotlib.animation as animation
from PIL import Image

from embrio import Embryo

#takes the embryo vector and turns it into a 2D grid for plotting
def unflatten(arr, n, m):
    output = np.zeros((n,m))

    for ii in range(n):
        output[ii,:] = arr[ii*m:(ii+1)*m]

    return output

#grows the embryo by multiplying the current state vector times the markov matrix,
#but only the squares adjacent to 'alive' squares so that it 'grows'. Also is more
#efficient. It is assumed that markov_mat is either 0s or 1s
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

        #if above > 99 or below > 99 or left > 99 or right > 99:
            #pdb.set_trace()

        if row == 0:
            if below not in new_cells:
                new_cells.append(below)
        elif row == n-1:
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
        elif col == m-1:
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

            #calculate the values for only the new cells
            for cell in new_cells:
                arr2[int(cell)] = markov_mat[int(cell),:]@arr

    return arr2


def main():

    np.random.seed(100)

    #grid size
    n = 10
    m = 10

    num_generations = 100
    num_growth_cycles = 20
    num_contestants = 10

    markov_matrix = np.random.rand(n*m,n*m)*2 - 1

    embryo0 = np.zeros(n*m)
    #inital embryo
    embryo0[0:10] = np.ones(10)

    #get a final goal body shape from a .bmp file I made in ms paint
    goal = np.array(Image.open('coral.bmp'), dtype = int).flatten()
    goal = np.ones(n*m) - goal

    errors = np.zeros(num_generations)
    for ii in range(num_generations):

        new_matrices = np.zeros((num_contestants, n*m, n*m))
        max_error = np.inf
        max_error_index = 0
        for kk in range(num_contestants):
            new_matrices[kk] = markov_matrix + (np.random.rand(n*m, n*m)*2 - 1)

            embryo = embryo0.copy()
            for jj in range(num_growth_cycles):
                print('generation {}, contestant {}, cycle {}, current max error {}'.format(ii, kk, jj, max_error))
                
                #iterate on the current embryo
                embryo = grow(embryo, n, m, new_matrices[kk])
                #make sure the values are all between 0 and 1
                embryo = np.clip(embryo, 0, 1)
                #make sure all value are 0 or 1
                embryo = np.round(embryo)

            error = np.linalg.norm(embryo-goal)
            print(error)
            if error < max_error:
                max_error = error
                max_error_index = kk

        markov_matrix = new_matrices[max_error_index]
        errors[ii] = max_error

    print(errors)
    animation_frames = np.zeros((num_growth_cycles, n*m))
    embryo = embryo0.copy()
    for jj in range(num_growth_cycles):

        #iterate on the current embryo
        embryo = grow(embryo, n, m, markov_matrix)
        #make sure the values are all between 0 and 1
        embryo = np.clip(embryo, 0, 1)
        #make sure all value are 0 or 1
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
