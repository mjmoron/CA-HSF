CC = g++
CFLAGS = -std=c++11 $(shell pkg-config --cflags opencv4) -lm -lgomp -lopencv_highgui -O3
LIBS = `pkg-config --libs opencv4`

all: 
	$(CC) -o CA_HSF  main_CA_HSF_2D_C_4stages_obj_dim.cpp CA_HSF_2D_C_4stages_obj_dim_func.cpp ca_infectious_process.cpp $(CFLAGS) $(LIBS)

clean:
	rm -rf CA_HSF   *.o

