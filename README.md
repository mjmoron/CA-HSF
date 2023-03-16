# CA-HSF
-__Authors__: Fernando Diaz-del-Rio, P. Sanchez-Cuevas, M. J. Moron-Fernández, José-Luis Guisado-Lizar, Senior Member,D. Cagigas-Muñiz, Pedro Real Jurado

-__Dptos__: ATC (www.atc.us.es) and MA1 (www.ma1.us.es) University of Seville. 2023 

-__Title__: Fully Parallel Cellular Automata for Topological Analysis of Color Digital Images 

-__Submitted to__: TRANSACTIONS ON IMAGE PROCESSING

-__Instructions__: Run make with the existent Makefile or modify it according to your machine. 

-__Requirements__: This project requires OpenCv4 and OpenMP libraries

-__Input__: image to be processed (Insert file name) that can be a real image (that is converted into a grayscale image) or a simple .txt file containing a synthetic image. 

-__Output__ : The output is the number of critical cells of each dimension and the minimum Time of a number of repetitions given by the constant MIN_NOF_REPETITIONS  

-__Options__ : Additionally, the user can modify the following constants:  

- #define MIN_NOF_REPETITIONS  1 // the minimun number of repetitions that the tests are done. Please increase this MIN_NOF_REPETITIONS   to obtain more accurate timing results. 

- #define MAX_NUM_THREADS 8. Please change this MAX_NUM_THREADS  according to the number of threads that you want to test

- In addition, the OpenMP SCHEDULING can be modified by these constants:  
  #define _OMP_SCHEDULING dynamic 
  #define CHUNK_SIZE 32
