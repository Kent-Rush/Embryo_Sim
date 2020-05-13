import matplotlib.pyplot as plt
from numpy import *

class Embryo():

    def __init__(self, grid_size, switch_dna, action_dna, chemicals):

        if grid_size[0]+grid_size[1] != len(switch_dna):
            print('Grid Dimensions do not match switch DNA lenght.')

        if len(chemicals)*len(chemicals)*9 != len(action_dna):
            print("Action DNA length incorrect size for chemical combos")

        self.grid_size = grid_size
        self.switch_dna = switch_dna
        self.action_dna = action_dna
        self.chemicals = chemicals
        self.switches = {}
        self.cells = {}
        self.behaviour = {}
        self.colors = {}
        self.initialize_switches()
        self.initialize_cells()
        self.initialize_behaviour()
        self.initialize_chemcolors()


    def mutate(self):
        pass

    def initialize_chemcolors(self):
        for a in self.chemicals:
            for b in self.chemicals:
                self.colors[a+b] = random.randint(0,255,3, dtype = int)

    def get_image(self, show_chems = False):

        arr = zeros((self.grid_size[0], self.grid_size[1], 3), dtype = int)

        if show_chems:
            for x in range(self.grid_size[0]):
                for y in range(self.grid_size[1]):
                    chem = self.get_switch((x,y))
                    arr[x,y,:] = self.colors[chem]

        for pos in self.cells.keys():
            arr[pos[0], pos[1],:] = array([255,255,255], dtype = int)

        return arr


    def grow_step(self):
        new_cells = self.cells.copy()
        for cell_pos in self.cells.keys():
            chemical = self.get_switch(cell_pos)
            neighbors = self.get_neighbors(cell_pos)
            action = self.get_action(chemical, neighbors)
            self.perform_action(action, cell_pos, new_cells)

        self.cells = new_cells

    def perform_action(self, action, pos, cells):

        if action == 8:
            cells.pop(pos)
            return 0

        if action == 0:
            x = pos[0] + 0
            y = pos[1] + 1
        elif action == 1:
            x = pos[0] + 1
            y = pos[1] + 1
        elif action == 2:
            x = pos[0] + 1
            y = pos[1] + 0
        elif action == 3:
            x = pos[0] + 1
            y = pos[1] + -1
        elif action == 4:
            x = pos[0] + 0
            y = pos[1] + -1
        elif action == 5:
            x = pos[0] + -1
            y = pos[1] + -1
        elif action == 6:
            x = pos[0] + -1
            y = pos[1] + 0
        elif action == 7:
            x = pos[0] + -1
            y = pos[1] + 1

        if x > 0 and y > 0 and x < self.grid_size[0] and y < self.grid_size[1]:
            cells[(x,y)] = 1


    def get_action(self, chemical, neighbors):
        return self.behaviour[(chemical, neighbors)]

    def initialize_behaviour(self):
        dna_index = 0
        for i in chemicals:
            for j in chemicals:
                for k in range(9):
                    self.behaviour[(i+j,k)] = self.action_dna[dna_index]
                    dna_index += 1
    
    def initialize_cells(self):
        x = int(self.grid_size[0]/2)
        y = int(self.grid_size[1]/2)
        self.cells[(x,y)] = True
        
    def get_neighbors(self, pos):

        neighbors = 0
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                coord = (pos[0]+i, pos[1]+j)
                if (coord != pos) and (coord in self.cells.keys()):
                    neighbors+=1

        return neighbors

    def initialize_switches(self):
        for x in range(self.grid_size[0]):
            for y in range(self.grid_size[1]):
                self.switches[(x,y)] = self.switch_dna[x] + self.switch_dna[self.grid_size[0] + y]

    def get_switch(self, pos):
        return self.switches[pos]

    


if __name__ == '__main__':

    grid_size = (100,100)
    chemicals = ['A', 'B', 'C', 'D']

    switch_dna = []
    for i in range(grid_size[0] + grid_size[1]):
        switch_dna.append(random.choice(chemicals))

    action_dna = []
    for i in range(len(chemicals)*len(chemicals)*9):
        action_dna.append(random.choice([0,1,2,3,4,5,6,7,8]))

    run1 = Embryo(grid_size, switch_dna, action_dna, chemicals)

    fig = plt.figure()
    for i in range(100):
        plt.imshow(run1.get_image(show_chems = True))
        plt.title(i)
        plt.draw()
        plt.pause(.001)
        run1.grow_step()

