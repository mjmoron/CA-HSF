// Copyright(c) 2021 Fernando Diaz del Rio and Pedro Real Jurado
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
// * Neither the name of CCLHSF nor the names of its
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

#include "CA_HSF_2D_C_4stages_obj_dim.h"

extern int n_rows;
extern int n_cols;

extern int n_rows_ext;
extern int n_cols_ext;




////////////////////
////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>




#include <iostream>
#include <fstream>

using namespace std;
using namespace cv;


int * primeraPosicionFuera;

int first_row_of_cells[dim_I + 1][MAX_NOF_K_CELLS] = {
1 - M2C, 1 - M2C,    // k = 0
1 - M2C, 2 - M2C,    // k = 1
2 - M2C, 2 - M2C,    // k = 2
};
int first_col_of_cells[dim_I + 1][MAX_NOF_K_CELLS] = {
	1 - M2C,1 - M2C,    // k = 0
	2 - M2C,1 - M2C,    // k = 1
	2 - M2C,2 - M2C,    // k = 2
};



void muestraMatrizResultado0(cv::Mat1i matrizBinaria) {
	printf("\nMatriz resultado:\n");
	for (int i = 0; i < matrizBinaria.rows; i++) {
		int* const fila = (int*)matrizBinaria.ptr<uint>(i);
		for (int j = 0; j < matrizBinaria.cols; j++) {
			printf("%d\t", fila[j]);
		}
		printf("\n");
	}
}

void muestraMatrizBinaria0(cv::Mat* matrizBinaria) {
	printf("\nMatriz binaria:\n");
	for (int i = 0; i < matrizBinaria->rows; i++) {
		uchar* const fila = (uchar*)matrizBinaria->ptr<uchar>(i);
		for (int j = 0; j < matrizBinaria->cols; j++) {
			printf("%d ", fila[j]);
		}
		printf("\n");
	}
}

////////////////////// --------------------

int f_nof_cells_per_nxel(int k, int dim) {

	if (k > dim || k < 0)
		exit(-1);

	//if (dim == 2)
	// 1, 2, 1
	int LUT_cells[3] = { 1, 2, 1 };
	/*else if (dim == 3)
	int LUT_cells = {1, 3, 3, 1};
	else if (dim == 4)
	int LUT_cells = {1, 4, 6, 4, 1};
	*/
	return LUT_cells[k + 1 - M2C];
}

int f_nof_boundary_neigh(int k, int dim) {

	if (k > dim || k < 0)
		exit(-1);

	//if (dim == 2)
	// 1, 2, 1
	int     LUT_neigh_cells[3] = { 0, 2, 4 };
	/*else if (dim == 3)
	else if (dim == 4)
	*/
	return LUT_neigh_cells[k + 1 - M2C];
}
int f_nof_coboundary_neigh(int k, int dim) {

	if (k > dim || k < 0)
		exit(-1);

	//if (dim == 2)
	// 1, 2, 1
	int     LUT_neigh_cells[3] = { 4,2,0 };
	/*else if (dim == 3)
	else if (dim == 4)
	*/
	return LUT_neigh_cells[k + 1 - M2C];
}

//////////////////////////////////////////////////////
void 	ca_preliminars(int* nof_cells_per_nxel, int* nof_bound_neigh, int* nof_cobound_neigh, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** bound_position_numb, int*** cobound_position_numb, int*** cobound_row_increm, int*** cobound_col_increm,int ***increm_boundary_neigh_row, int*** increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col) {

	//Preliminar file that builds several matrixes and structs('cells' in MATLAB)

		// dim_I = 2; // Number of array dimensions
	int size_I[2] = { n_rows, n_cols };


	int k = 1;
	int pp0[4][2] = { // (bound_numb, c, k + 1);
0, 1,  // bound 1  c = 1, c = 2
	0, -1,  // bound 2
	0, 0,  // bound 3 does not exist
	0, 0,  // bound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_boundary_neigh_row[i][j][k + 1 - M2C] = pp0[i][j];

	int pp1[4][2] = { // (bound_numb, c, k + 1);
	1, 0,  // bound 1
		-1, 0,  // bound 2
		0, 0,  // bound 3 does not exist
		0, 0,  // bound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_boundary_neigh_col[i][j][k + 1 - M2C] = pp1[i][j];

	k = 2;
	int pp2[4][2] = { // (bound_numb, c, k + 1);
	1, 0,  // bound 1  (only c = 1)
		0, 0,  // bound 2
		-1, 0,  // bound 3
		0, 0,  // bound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_boundary_neigh_row[i][j][k + 1 - M2C] = pp2[i][j];

	int pp3[4][2] = { // (bound_numb, c, k + 1);
	0, 0,  // bound 1 (only c = 1)
		1, 0,  // bound 2
		0, 0,  // bound 3
		-1, 0,  // bound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_boundary_neigh_col[i][j][k + 1 - M2C] = pp3[i][j];

	//////////////////////////////////////////////////////////////////////////////////////////
		////
//	int increm_coboundary_neigh_row[4][2][3];
	//int increm_coboundary_neigh_col[4][2][3];

	k = 0;
	int pp4[4][2] = { // (cobound_numb, c, k + 1);
	1, 0,  // cobound 1  (only c = 1)
		0, 0,  // cobound 2
		-1, 0,  // cobound 3
		0, 0,  // cobound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_coboundary_neigh_row[i][j][k + 1 - M2C] = pp4[i][j];

	int pp5[4][2] = { // (cobound_numb, c, k + 1);
	0, 0,  // cobound 1 (only c = 1)
		1, 0,  // cobound 2
		0, 0,  // cobound 3
		-1, 0,  // cobound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_coboundary_neigh_col[i][j][k + 1 - M2C] = pp5[i][j];

	k = 1;
	int pp6[4][2] = { // (cobound_numb, c, k + 1);
	1, 0,  // cobound 1  c = 1, c = 2
		-1, 0,  // cobound 2
		0, 0,  // cobound 3 does not exist
		0, 0,  // cobound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_coboundary_neigh_row[i][j][k + 1 - M2C] = pp6[i][j];

	int pp7[4][2] = { // (cobound_numb, c, k + 1);
	0, 1,  // cobound 1
		0, -1,  // cobound 2
		0, 0,  // cobound 3 does not exist
		0, 0,  // cobound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) increm_coboundary_neigh_col[i][j][k + 1 - M2C] = pp7[i][j];

	k = 2;
	// k = 2 does not have cobounds

//////////////////////////////////////////////////////////////////////////////////////////
		////
		// for each c, bound of a k - cell, this matrix returns the corresponding position c_bound of the k -1cell
	//int bound_position_numb[4][2][3]; // (bound_numb, c, k + 1)
	k = 0; //never ask for a bound

	k = 1;
	int pp8[4][2] = { // (cobound_numb, c, k + 1); cobound_numb @@
	1, 1,  // bound_numb1  c = 1, c = 2
	1, 1,  // bound_numb2  c = 1, c = 2
	0, 0,  // bound_numb 3 does  not exist
	0, 0,  // bound_numb 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) bound_position_numb[i][j][k + 1 - M2C] = pp8[i][j];

	k = 2;
	int pp9[4][2] = { // (cobound_numb, c, k + 1); cobound_numb @@
		1, 0,  // bound_numb1  c = 1, c = 2
		2, 0,  // bound_numb2  c = 1, c = 2
		1, 0,  // bound_numb 3 does  not exist
		2, 0,  // bound_numb 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) bound_position_numb[i][j][k + 1 - M2C] = pp9[i][j];
	////

	// for each cobound of a k -1 - cell, this matrix returns the corresponding position of the k cell
	//int cobound_position_numb[4][2][3] ; // (cobound_numb, c, k + 1)
	k = 0;
	int ppq[4][2] = { // (cobound_numb, c, k + 1); cobound_numb @@
	2, 0,  // cobound_numb1  c = 1,
	1, 0,  // cobound_numb2  c = 1,
	2, 0,  // cobound_numb 3
	1, 0,  // cobound_numb 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) cobound_position_numb[i][j][k + 1 - M2C] = ppq[i][j];

	k = 1;
	int ppw[4][2] = { // (cobound_numb, c, k + 1); cobound_numb @@
		1, 1,  // cobound_numb1  c = 1, c = 2
		1, 1,  // cobound_numb2  c = 1, c = 2
		0, 0,  // cobound_numb 3 does  not exist
		0, 0,  // cobound_numb 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) cobound_position_numb[i][j][k + 1 - M2C] = ppw[i][j];

	k = 2; //never ask for a cobound


	////////////////////////////////////////////////////////////////////////////////////////
		//// for each cobound of the k -1 cell, these matrices returns the corresponding row / col increm of the kcell

//	int cobound_row_increm[4][2][3]; // (cobound3_numb_grant, c_position_bound, k + 1 -1);
	//int cobound_col_increm[4][2][3]; // (cobound3_numb_grant, c_position_bound, k + 1 -1);
	k = 0;
	int ppe[4][2] = { // (cobound_numb, c, k + 1);
	1, 0,  // cobound 1  (only c = 1)
		0, 0,  // cobound 2
		-1, 0,  // cobound 3
		0, 0,  // cobound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) cobound_row_increm[i][j][k + 1 - M2C] = ppe[i][j];
	int ppr[4][2] = { // (cobound_numb, c, k + 1);
	0, 0,  // cobound 1  (only c = 1)
		1, 0,  // cobound 2
		0, 0,  // cobound 3
		-1, 0,  // cobound 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) cobound_col_increm[i][j][k + 1 - M2C] = ppr[i][j];

	k = 1;
	int ppt[4][2] = { // (cobound_numb, c, k + 1);
	1, 0,  // cobound_numb1  c = 1, c = 2
		-1, 0,  // cobound_numb2  c = 1, c = 2
		0, 0,  // cobound_numb 3 does  not exist
		0, 0,  // cobound_numb 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) cobound_row_increm[i][j][k + 1 - M2C] = ppt[i][j];

	int ppy[4][2] = { // (cobound_numb, c, k + 1);
	0, 1,  // cobound_numb1  c = 1, c = 2
		0, -1,  // cobound_numb2  c = 1, c = 2
		0, 0,  // cobound_numb 3 does  not exist
		0, 0,  // cobound_numb 4
	};
	for (int i = 0; i < 4; i++)  for (int j = 0; j < 2; j++) cobound_col_increm[i][j][k + 1 - M2C] = ppy[i][j];

	k = 2; //never ask for a cobound

	
	for (int i = 0; i < 3; i++)  for (int j = 0; j < 2; j++) ext_first_row_of_cells[i][j] = first_row_of_cells[i][j] + 2;
	//int ext_first_col_of_cells[3][2];
	for (int i = 0; i < 3; i++)  for (int j = 0; j < 2; j++) ext_first_col_of_cells[i][j] = first_col_of_cells[i][j] + 2;

	/////////////////////////////////////////////////////////////////////
	//// nof_cells_per_nxel : nof cells of each dim that a nxel has, and boundary and coboundary cells
	// for all dims

	for (int d = 1 - M2C; d <= dim_I; d++)
	{		// these little functions can be converted into LUT(Look Up Tables)
		nof_cells_per_nxel[d] = f_nof_cells_per_nxel(d, dim_I);
		nof_bound_neigh[d] = f_nof_boundary_neigh(d, dim_I);
		nof_cobound_neigh[d] = f_nof_coboundary_neigh(d, dim_I);
	}

	// the maximum for all dim :
	int max_nof_cells_per_nxel = 0;
	for (int d = 1 - M2C; d < dim_I; d++)
		max_nof_cells_per_nxel = max(max_nof_cells_per_nxel, nof_cells_per_nxel[d + 1 - M2C]);

	////////////////////////////////////////////////////////////////////
			//// ACC with new extended(double) borders to ease the CA computation.WE
			// need one additional borders for the primal vector searching around the real image,
			// and another additional borders for the dual vector searching around the previous additional border

// the biggest TMAX when using the infection process is :
#define TMAX (n_rows_ext  + n_cols_ext );//sum(size_acc); // max nof steps for the iterative CA

		//// in order to introduce the values(extended linear indexes) of the ext_acc
		// we initialize with :
	// linear idxs useful to mark cells :
/*	pp_ext_lin_idx = 1 : (size_ext_acc(1)*size_ext_acc(2));
	pp_ext_acc_init = reshape(pp_ext_lin_idx, size_ext_acc);
	*/

	//@ copying colors to the new extended borders
	// IT DOES NOT WORK : ext_acc_init(idx_all_cells{ :, : , : }) = acc_init;
	
	//end

}

/////////////////////////////////////////////////////////////////////////////////////////
//1st stage
#ifdef DEBUG_SYNTHETIC_IMAGE
#endif
#ifndef DEBUG_SYNTHETIC_IMAGE
#endif

#ifdef DEBUG_SYNTHETIC_IMAGE
//void Jonly_transports(const cv::Mat1b &I, cv::Mat1i &J_CRIT);
//void Jonly_transports_no_borders(const cv::Mat1b &I, int J_CRIT[n_rows][n_cols]);
#endif
#ifndef DEBUG_SYNTHETIC_IMAGE
//void Jonly_transports(const cv::Mat1b &I, cv::Mat1i &J_CRIT);
void Jonly_transports_no_borders(const cv::Mat1b &I, cv::Mat1i &J_CRIT);
#endif


