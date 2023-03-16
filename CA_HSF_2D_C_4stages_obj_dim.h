// Copyright(c) 2023 Fernando Diaz-del-Rio , P. Sanchez-Cuevas, M. J. Moron-Fernández, José-Luis Guisado-Lizar, Senior Member,D. Cagigas-Muñiz, Pedro Real Jurado
// Submited to TRANSACTIONS ON IMAGE PROCESSING:
// TITLE: Fully Parallel Cellular Automata for Topological Analysis of Color Digital Images
// Dpto ATC: www.atc.us.es and Dpto MA1 www.ma1.us.es
// University of Seville.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :
// 
// *Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
// 
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and / or other materials provided with the distribution.
// 
// * Neither the name of CA-HSF nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

typedef  char  Orientation_4_adj_type;
#include <iostream>
#include <fstream>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <windows.h>


#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <omp.h>
using namespace std;


////////////////////
// for testing and debugging purposes. 
////////////////////
#define  dim_I 2
#define MIN_NOF_REPETITIONS  1// the minimun number of repetitions that the tests are done 

//#define _INTERNAL_TIMER_STEP4
#define  DEBUG_INNER_STATISTICS_no

#define M2C (1)  // a -1 is inserted in all indexations from 1 in MATLAB to 0 in C


// nof threads
extern int num_th;

////////////////////
////////////////////

#define MAX_NUM_THREADS 8
#define _OMP_SCHEDULING dynamic 
#define CHUNK_SIZE 32

#pragma once

////////////////////////////////  typedef_constants.h
typedef unsigned char Bool;

// binary images are allowed
typedef unsigned char ImageType;
typedef unsigned int LabelType;
typedef int JumpType;

// due to this eror message "Error	4	error C3016: 'c' : index variable in OpenMP 'for' statement must have signed integral type	"
// the index type has sign
typedef int RowColType;



//////////////////
// colors for :
#define FOREGROUND 1
#define BACKGROUND 0
#define FG FOREGROUND 
#define BG BACKGROUND 

#define INEXISTENT_ADJ_CELL (-88)
#define INEXISTENT_AFF_TREE (0)



#define row_ind(X) ( (X-1)%n_rows)
#define col_ind(X) ( (X-1)/n_rows)

#define tag_ind(row, col) ( (col)+(n_cols*row) )
#define tag_ind_MATLAB(row, col) ( (row)+(n_rows*col) )
// (row_new_tags_NE + 1 + (col_new_tags_NE * n_rows));

#define  ADD_ZERO_BORDER_MATRIX(X) \
for (RowColType c = 0; c < n_cols; c++) X[0][c] = 0;  \
for (RowColType c = 0; c < n_cols; c++) X[n_rows - 1][c] = 0;  \
for (RowColType r = 0; r < n_rows; r++) X[r][0] = 0; \
for (RowColType r = 0; r < n_rows; r++) X[r][n_cols - 1] = 0;





///////////  prototypes1.h
// init image functions

//just write negative numbers in the limits of previous loop to use another images 
//OJO!!!!void inic_synthetic_image_rand1(void);
void ca_infectious_process(ImageType **I, int	total_nof_cells_ext, int* nof_cells_per_nxel, int* nof_bound_neigh, int* nof_cobound_neigh, int** ext_acc, int** ext_acc_init, int** ext_acc_previous, Orientation_4_adj_type ** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** bound_position_numb, int*** cobound_position_numb, int*** cobound_row_increm, int*** cobound_col_increm,int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, int** ext_max_cobound, Orientation_4_adj_type * *ext_max_cobound_numb, Orientation_4_adj_type * *ext_obj_dim_belonging, Orientation_4_adj_type * *ext_aff_tree, Orientation_4_adj_type * *ext_s_k, Orientation_4_adj_type * *ext_t_k, Orientation_4_adj_type **ext_being_tail);


#if dim_I == 2
#define MAX_NOF_K_CELLS 2
// this is the max of the Tartaglia's Triangle for this dim_I , 
//  which can be computed as (dim_I  over (dim_I/2) ), combinations of dim_I  elem taken in (dim_I/2)  sets
#else
#if dim_I == 3
#define MAX_NOF_K_CELLS 3
#else
#if dim_I == 4
#define MAX_NOF_K_CELLS 6
#else
///
#endif
#endif
#endif

#define MAX_NOF_BOUNDS (2* dim_I)
// NOTE THAT MAX_NOF_COBOUNDS  is also (2* dim_I)
// the exact nof bounds for a k-cell is simply (2*k): 
// the exact nof cobounds for a k-cell is simply (2* (dim_I-k)): 

// the next tables are [4][2][3];  in R2
// and the next tables would be [6][3][4];  in R3



//////////////////
