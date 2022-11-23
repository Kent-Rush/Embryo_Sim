COMPILER      = g++ 
OPTIONS       = -ansi -Wall -Wextra -g -std=c++17
INCLUDES      = -I ./Eigen/  -I ./Eigen/Dense
LINKER_OPT    = -L/usr/lib -lstdc++ -lm

all: embryo_sim

bitmap_test: bitmap_test.cpp bitmap_image.hpp
	$(COMPILER) $(OPTIONS) -o bitmap_test bitmap_test.cpp $(LINKER_OPT)

embryo_sim: embryo_sim.cpp bitmap_image.hpp
	$(COMPILER) $(OPTIONS) $(INCLUDES) -o embryo_sim embryo_sim.cpp $(LINKER_OPT)

valgrind_check:
	valgrind --leak-check=full --show-reachable=yes --track-origins=yes -v ./bitmap_test

clean:
	rm -f *.o
	sudo rm embryo_sim
	sudo rm bitmap_test

video:
	sudo ffmpeg -framerate 30 -i frames/fractal_%d.bmp gifs/output.mp4


