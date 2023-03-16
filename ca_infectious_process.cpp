// Copyright(c) 2023 Fernando Diaz-del-Rio , P. Sanchez-Cuevas, M. J. Moron-Fernández, José-Luis Guisado-Lizar, Senior Member,D. Cagigas-Muñiz, Pedro Real Jurado
// Submited to TRANSACTIONS ON IMAGE PROCESSING:
// Fully Parallel Cellular Automata for Topological Analysis of Color Digital Images
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

#include "CA_HSF_2D_C_4stages_obj_dim.h"
////////////////////


//#include "opencv/cv.h"
//#include "opencv/highgui.h"



using namespace std;
using namespace cv;

/////////////////////////////////////////////////////////////////////////////////////////

/**/ // ext_obj_dim_belonging must be global
//Orientation_4_adj_type ext_obj_dim_belonging[n_rows_ext][n_cols_ext];// = zeros(size_ext_acc);

extern int n_rows;
extern int n_cols;

extern int n_rows_ext;
extern int n_cols_ext;


// extended matrices :

extern int size_ext_acc[dim_I];



void init_n1_frontier_cells(int row, int col, int k, int c, int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int ***increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, ImageType** I, int nof_cobound_neigh[dim_I + 1], Orientation_4_adj_type** ext_obj_dim_belonging);
void init_obj_belonging_cells_color(int k, int c, int row, int col, int nof_cobound_neigh[dim_I + 1], int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, Orientation_4_adj_type** ext_obj_dim_belonging);

void stage0_k1_cells(int k, int total_nof_cells_ext,int* nof_cells_per_nxel, int** ext_acc,int** ext_first_row_of_cells, int** ext_first_col_of_cells,int nof_cobound_neigh[dim_I + 1], int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, int** ext_max_cobound, Orientation_4_adj_type** ext_max_cobound_numb, Orientation_4_adj_type** ext_being_tail, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree);
void stage1_k_cells(int k, int* nof_cells_per_nxel, int** ext_acc, int** ext_acc_init, Orientation_4_adj_type** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int nof_bound_neigh[dim_I + 1], int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col,int** reqmax_cobound, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_s_k, Orientation_4_adj_type** ext_t_k);
void stage2_k1_cells(int k, int* nof_cells_per_nxel, int* nof_cobound_neigh, int** ext_acc, int** ext_acc_init, Orientation_4_adj_type** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** cobound_position_numb,int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, int** ext_max_cobound, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_s_k, Orientation_4_adj_type** ext_t_k );
void stage3_k_cells(int k, int* nof_cells_per_nxel, int** ext_acc, Orientation_4_adj_type** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** bound_position_numb, int*** cobound_row_increm, int*** cobound_col_increm, int*** increm_boundary_neigh_row, int*** increm_boundary_neigh_col, int** ext_max_cobound, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_s_k, Orientation_4_adj_type** ext_t_k);


void ending_k2_cells_only_tail_cells(int k, int c, int row, int col, int* nof_cobound_neigh,int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int ***increm_coboundary_neigh_col, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_being_tail);


void ca_infectious_process(ImageType **I, int	total_nof_cells_ext, int* nof_cells_per_nxel, int* nof_bound_neigh, int* nof_cobound_neigh, int** ext_acc, int** ext_acc_init, int** ext_acc_previous,Orientation_4_adj_type** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** bound_position_numb, int*** cobound_position_numb, int*** cobound_row_increm, int*** cobound_col_increm,int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col,int** ext_max_cobound, Orientation_4_adj_type** ext_max_cobound_numb, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_s_k, Orientation_4_adj_type** ext_t_k,  Orientation_4_adj_type** ext_being_tail)
{
	//int step_counter[dim_I];//= zeros(dim_I, 1);
	int step_total = 0;


	//setInitValues(ext_obj_dim_belonging, ext_aff_tree, ext_s_k, ext_t_k);
	// dim_I - 1 CELLS
	int k = dim_I - 1;
	for (int c = 1 - M2C; c < nof_cells_per_nxel[k + 1 - M2C]; c++) { // for each position in a nxel(of k - cell)
		int first_row = (2 * (c == (1 - M2C))) + ext_first_row_of_cells[k + 1 - M2C][c];
		int first_col = (2 * (c == (2 - M2C))) + ext_first_col_of_cells[k + 1 - M2C][c];

		for (int col = first_col; col < (size_ext_acc[2 - M2C] - 2 - (2 * (c == (2 - M2C)))); col += 2) {
			for (int row = first_row; row < (size_ext_acc[1 - M2C] - 2 - (2 * (c == (1 - M2C)))); row += 2) {
			init_n1_frontier_cells(row, col, k, c, increm_boundary_neigh_row, increm_boundary_neigh_col, increm_coboundary_neigh_row, increm_coboundary_neigh_col,I, nof_cobound_neigh, ext_obj_dim_belonging);

			}
		}
	}// end of for c = 1 : nof_cells_per_nxel(k + 1 -1] // for each cell of dim k in a nxel


	//% rest of  k - CELLS
	//	% that is, 0 - cells in the case of R2, have(n - k) * 2 coborders.Depending on
	//	% them, they have a object dimension.

	
	k = dim_I - 2;
	for (int c = 1 - M2C; c < nof_cells_per_nxel[k + 1 - M2C]; c++) { // for each position in a nxel(of k - cell)
		int first_row = ext_first_row_of_cells[k + 1 - M2C][c];
		int first_col = ext_first_col_of_cells[k + 1 - M2C][c];

		for (int col = first_col; col < (size_ext_acc[2 - M2C] - 2); col += 2) {
			for (int row = first_row; row < (size_ext_acc[1 - M2C] - 2); row += 2) {
			init_obj_belonging_cells_color(k, c, row, col, nof_cobound_neigh, increm_boundary_neigh_row, increm_boundary_neigh_col,increm_coboundary_neigh_row, increm_coboundary_neigh_col,ext_obj_dim_belonging);

			}
		}
	}// end of for c = 1 : nof_cells_per_nxel(k + 1 -1] // for each cell of dim k in a nxel



	for (int k_k1_tree = dim_I; k_k1_tree >= 1; k_k1_tree--)
	{
		//// CA HSF building
		//	for k_k1_tree = dim_I : -1 : 1 // that is, k - (k - 1) tree

				////
				//     instead of the following lists it is better the coord o f the first
				//     k - cell in the position c.However this old list can be better for a
				//     GPU, and for n - dim images

		int step = 0;
		int TMAX = size_ext_acc[1 - M2C] * size_ext_acc[2 - M2C]; // 15; // @ provisional
		int ppsum = 1;

		while ((step <= TMAX) && (ppsum > 0))
			/*	while (step <= TMAX) && ((ppsum) > 0) // @a sufficient number of setps to guarantee the CA has converged(also, a while loop finishing when the CA has not changed)*/
		{
			step = step + 1;
			// ext_acc_previous = ext_acc;
#pragma omp parallel for default(none) shared (ext_acc, ext_acc_previous, ext_max_cobound_numb, ext_max_cobound, ext_requested_s_k,n_rows_ext,n_cols_ext)  					
			for (int i = 0; i < n_rows_ext; i++) {
				for (int j = 0; j < n_cols_ext; j++) {
					ext_acc_previous[i][j] = ext_acc[i][j];
					ext_max_cobound[i][j] = 0;
					ext_max_cobound_numb[i][j] = INEXISTENT_ADJ_CELL;
					ext_requested_s_k[i][j] = INEXISTENT_ADJ_CELL;
					/**/
					//ext_max_cobound = zeros(size_ext_acc);
					//ext_max_cobound_numb = zeros(size_ext_acc); //this is like the state of a k - 1 cell, the biggest label it has on its cobounds
				}
			}


			//// stage 0
					// k - 1 CELLS
			k = k_k1_tree - 1;
			
			stage0_k1_cells(k, total_nof_cells_ext,nof_cells_per_nxel, ext_acc,ext_first_row_of_cells, ext_first_col_of_cells, nof_cobound_neigh, increm_coboundary_neigh_row, increm_coboundary_neigh_col, ext_max_cobound, ext_max_cobound_numb, ext_being_tail, ext_obj_dim_belonging, ext_aff_tree);

			for (int i = 0; i < n_rows_ext; i++) {
				fill_n(ext_requested_s_k[i], n_cols_ext, INEXISTENT_ADJ_CELL);
			}
			  //// stage 1
			// k - CELLS
			k = k_k1_tree;
		
			stage1_k_cells(k, nof_cells_per_nxel, ext_acc, ext_acc_init,ext_requested_s_k, ext_first_row_of_cells, ext_first_col_of_cells, nof_bound_neigh, increm_boundary_neigh_row, increm_boundary_neigh_col, increm_coboundary_neigh_row, increm_coboundary_neigh_col,ext_max_cobound, ext_obj_dim_belonging, ext_aff_tree, ext_s_k, ext_t_k);

			 //// stage 2
			// k - 1 CELLS
			k = k_k1_tree - 1;

			
			stage2_k1_cells(k, nof_cells_per_nxel, nof_cobound_neigh, ext_acc, ext_acc_init,ext_requested_s_k, ext_first_row_of_cells, ext_first_col_of_cells, cobound_position_numb,increm_boundary_neigh_row, increm_boundary_neigh_col, increm_coboundary_neigh_row, increm_coboundary_neigh_col, ext_max_cobound, ext_obj_dim_belonging, ext_aff_tree, ext_s_k, ext_t_k);
			
			//// stage 3
			// k CELLS AGAIN
			k = k_k1_tree;
			
			stage3_k_cells(k, nof_cells_per_nxel, ext_acc,  ext_requested_s_k, ext_first_row_of_cells, ext_first_col_of_cells, bound_position_numb, cobound_row_increm, cobound_col_increm,increm_boundary_neigh_row, increm_boundary_neigh_col, ext_max_cobound, ext_obj_dim_belonging, ext_aff_tree, ext_s_k, ext_t_k);
				
			ppsum = 0;
#pragma omp parallel for default(none) shared (ext_acc_previous, ext_acc, n_rows_ext, n_cols_ext) reduction (+:ppsum) 
			for (int i = 0; i < n_rows_ext; i++) {
				for (int j = 0; j < n_cols_ext; j++) {
					ppsum += (ext_acc_previous[i][j] != ext_acc[i][j]);
				}
			}
		
			/*
			 //just a summatory of the changes in the CA to have an idea of the nof iterations of the CA :
			 ext_acc_differences = (ext_acc ~= ext_acc_previous);
			 ppsum = sum(ext_acc_differences);
			 for pp = 2:dim_I
				 ppsum = sum(ppsum, pp);
			 end
			 //// plotting
				 if dim_I == 2 && DEBUG_PLOTS == 9
					 // DEBUG
					 f_ca_plot2D_image(I, FIG_IMAGE_AND_HSF + step + step_total);
			 f_ca_plot2D_primal_vectors(ext_s_k, FIG_IMAGE_AND_HSF + step + step_total);
			 f_ca_bookkeep_dual_vectors();
			 f_ca_plot2D_dual_vectors(ext_u_k, FIG_IMAGE_AND_HSF + step + step_total);
			 //             f_ca_plot2D_critical_cells(ext_aff_tree, FIG_IMAGE_AND_HSF + step);
			 end
				 step
				 end // end of WHILE step
				 */


		}  //end of while

		printf(" (%d-%d) Tree in %d steps \n", k_k1_tree, k_k1_tree - 1, step);

		//// ENDINGS

		// NO ENDINGS for    // k-1 CELLS
		// % ONLY ENDINGS FOR TAIL CELLS ext_being_tail
			// % ONLY ENDINGS FOR TAIL CELLS ext_being_tail
		if (k_k1_tree > 1) {
			// k-2 CELLS
			k = k_k1_tree - 2;
			for (int c = 1 - M2C; c < nof_cells_per_nxel[k + 1 - M2C]; c++) // for each position in a nxel(of k - 1 cell)
			{
				int first_row = ext_first_row_of_cells[k + 1 - M2C][c];
				int first_col = ext_first_col_of_cells[k + 1 - M2C][c];
#pragma omp parallel for shared (k, first_row, first_col, size_ext_acc,  \
				cobound_row_increm, cobound_col_increm, nof_cobound_neigh, increm_coboundary_neigh_col, bound_position_numb,  \
				ext_s_k, ext_t_k, ext_requested_s_k, ext_aff_tree, ext_obj_dim_belonging, \
				ext_max_cobound_numb, ext_max_cobound, ext_being_tail)  					
				for (int row = first_row; row < (size_ext_acc[1 - M2C] - 2); row += 2) {
					for (int col = first_col; col < (size_ext_acc[2 - M2C] - 2); col += 2) {
						// 	%					ending_k2_cells;
//#include "ending_k2_cells_only_tail_cells.cpp" // it needs : row, col, k, c and global matrices
						ending_k2_cells_only_tail_cells(k, c, row, col, nof_cobound_neigh, increm_boundary_neigh_row, increm_boundary_neigh_col, increm_coboundary_neigh_row, increm_coboundary_neigh_col, ext_obj_dim_belonging, ext_aff_tree, ext_being_tail);
					}//end
				}//	end

			} //end // end of for c = 1 : nof_cells_per_nxel(k + 1] // for each cell of dim k in a nxel
		} //end of  if (k_k1_tree > 1) 

	} //end of for




}
                           
void init_n1_frontier_cells(int row, int col, int k, int c, int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col,ImageType** I, int nof_cobound_neigh[dim_I + 1], Orientation_4_adj_type** ext_obj_dim_belonging) {

	int nof_equal_colors = 1; // the first is idem to itself

	int cobound_numb = 1 - 1; // first cobound
	int row_increm = increm_coboundary_neigh_row[cobound_numb][c][k + 1 - 1];
	int col_increm = increm_coboundary_neigh_col[cobound_numb][c][k + 1 - 1];
	int row_cobound = row + row_increm;
	int col_cobound = col + col_increm;
	//int first_color = I[(row_cobound) / 2 - 1][(col_cobound) / 2 - 1]; // nof_regions  of the k - cell
	int first_color = I[(row_cobound + 1) / 2 - 1 - 1][(col_cobound + 1) / 2 - 1 - 1]; // nof_regions  of the k - cell

	// Any cobound cell of k - 1 cells is compared to find their diff
	for (int cobound_numb = 2 - 1; cobound_numb < nof_cobound_neigh[k + 1 - 1]; cobound_numb++)
	{
		int row_increm = increm_coboundary_neigh_row[cobound_numb][c][k + 1 - 1];
		int col_increm = increm_coboundary_neigh_col[cobound_numb][c][k + 1 - 1];
		int row_cobound = row + row_increm;
		int col_cobound = col + col_increm;

		//	int next_color = I[(row_cobound) / 2 - 1][(col_cobound) / 2 - 1]; // nof_regions  of the k - cell
		int next_color = I[(row_cobound + 1) / 2 - 1 - 1][(col_cobound + 1) / 2 - 1 - 1]; // nof_regions  of the k - cell

		if (first_color == next_color)
			// the cobound must be of the same color than the first n - xel
			nof_equal_colors = nof_equal_colors + 1;
	}//end //end of for cobound_neigh_numb

	if (nof_equal_colors == 2)
		ext_obj_dim_belonging[row][col] = dim_I;
	else
		ext_obj_dim_belonging[row][col] = dim_I - 1;


}


void init_obj_belonging_cells_color(int k, int c, int row, int col, int nof_cobound_neigh[dim_I + 1], int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col,Orientation_4_adj_type** ext_obj_dim_belonging) {
	
	int nof_k_dim_cob[dim_I]{}; // % to compute the sum of n-dim cobounds 
	for (int j = 0; j < dim_I; j++) {
		nof_k_dim_cob[j] = 0;
	}

	// Any cobound cell of k - 1 cells is compared to find their diff
	for (int cobound_numb = 1 - 1; cobound_numb < nof_cobound_neigh[k + 1 - 1]; cobound_numb++)
	{
		int row_increm = increm_coboundary_neigh_row[cobound_numb][c][k + 1 - 1];
		int col_increm = increm_coboundary_neigh_col[cobound_numb][c][k + 1 - 1];
		int row_cobound = row + row_increm;
		int col_cobound = col + col_increm;

		Orientation_4_adj_type obj_dim_cob = ext_obj_dim_belonging[row_cobound][col_cobound];
		nof_k_dim_cob[obj_dim_cob - 1] = nof_k_dim_cob[obj_dim_cob - 1] + 1;

	}//end //end of for cobound_neigh_numb

	ext_obj_dim_belonging[row][col] = floor(nof_k_dim_cob[k + 1 + 1 - 1] / (2));
}

void stage0_k1_cells(int k, int total_nof_cells_ext, int* nof_cells_per_nxel, int** ext_acc, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int nof_cobound_neigh[dim_I + 1], int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, int** ext_max_cobound, Orientation_4_adj_type** ext_max_cobound_numb, Orientation_4_adj_type** ext_being_tail, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree) {
	for (int c = 1 - M2C; c < nof_cells_per_nxel[k + 1 - M2C]; c++) // for each position in a nxel(of k - 1 cell)
	{
		int first_row = ext_first_row_of_cells[k + 1 - M2C][c];
		int first_col = ext_first_col_of_cells[k + 1 - M2C][c];
		
#pragma omp parallel for default(none) shared (c, k, first_row, first_col, size_ext_acc,  \
		nof_cobound_neigh, increm_coboundary_neigh_row, increm_coboundary_neigh_col, \
		total_nof_cells_ext,\
		ext_acc, ext_obj_dim_belonging, ext_max_cobound, ext_being_tail, ext_aff_tree,ext_max_cobound_numb)  
		// wwww quito estas var q no se usan aqui :   ext_s_k, ext_t_k,
		for (int row = first_row; row < (size_ext_acc[1 - M2C] - 2); row += 2) {
			for (int col = first_col; col < (size_ext_acc[2 - M2C] - 2); col += 2) {
				// if (ext_aff_tree(row, col) <= k_k1_tree) // $$ && v_color ~= 0)
			// NOT NECESSRAY for the k - 1 cells check the cell_belonging of this k - cell(it can belong to the superior tree) and that it is not an interface(or border) cell(v_color == 0)
	//#include "stage0_k1_cells.cpp" // it needs : row, col, k, c and global matrices

				int L_k1_ext = ext_acc[row][col]; // value(label) of the k - 1 - cell
				int max0_L_k = 0;
				// // // min0_L_k = inf;
				//// we ask for the coboundary k - cells(dual vectors)
				int biggest_L_k_cobound_numb = INEXISTENT_ADJ_CELL;
				int v0_obj_dim = ext_obj_dim_belonging[row][col];
				// Any cobound cell is compared to find the biggest values
				for (int cobound_numb = 1; cobound_numb <= nof_cobound_neigh[k + 1 - 1]; cobound_numb++)
				{
					// one of these cobound cells(of the previous cobound) is the
					// same k - cell whose maximum we are computing

					int row_increm = increm_coboundary_neigh_row[cobound_numb - M2C][c][k + 1 - M2C];
					int col_increm = increm_coboundary_neigh_col[cobound_numb - M2C][c][k + 1 - M2C];
					int row_cobound = row + row_increm;
					int col_cobound = col + col_increm;

					int v0_ext_acc_cobound = ext_acc[row_cobound][col_cobound]; // value of the k - cell

					int v0_cobound_obj_dim = ext_obj_dim_belonging[row_cobound][col_cobound];
					if ((v0_cobound_obj_dim == v0_obj_dim) && (ext_aff_tree[row_cobound][col_cobound] <= k + 1))
					{
						// the cobound must be of the same color than the cell
						//         and it must have the same color than the cell and not be affiliated

						// exploring the max values for all the cobounds
						if (max0_L_k < v0_ext_acc_cobound) {
							max0_L_k = v0_ext_acc_cobound;
							biggest_L_k_cobound_numb = cobound_numb;
						}//end
					}//end // end of if (v_nof_regions_cobound == v_nof_regions)
				}//end of for cobound_neigh_numb
				ext_max_cobound[row][col] = max0_L_k;
				ext_max_cobound_numb[row][col] = biggest_L_k_cobound_numb;
				// biggest possible to obbly the k cell to be coupled with it
				if (ext_being_tail[row][col] == 1) {
					ext_max_cobound[row][col] = L_k1_ext + total_nof_cells_ext;
					ext_being_tail[row][col] = 0; // the mark is deleted : the tail condition is to be applied now
				}//end;

			}//end
		}//	end

	} //end // end of for c = 1 : nof_cells_per_nxel(k + 1] // for each cell of dim k in a nxel


}


void stage1_k_cells(int k, int* nof_cells_per_nxel, int** ext_acc, int** ext_acc_init,Orientation_4_adj_type** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int nof_bound_neigh[dim_I + 1], int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col,int** ext_max_cobound, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_s_k, Orientation_4_adj_type** ext_t_k) {
	for (int c = 1 - M2C; c < nof_cells_per_nxel[k + 1 - M2C]; c++) // for each position in a nxel(of k - 1 cell)
	{
		//for c = 1:nof_cells_per_nxel[k + 1] // for each position in a nxel(of k - cell)
			//int first_row = ext_first_row_of_cells[k + 1][c];
		//int first_col = ext_first_col_of_cells[k + 1][c];


		int first_row = ext_first_row_of_cells[k + 1 - M2C][c];
		int first_col = ext_first_col_of_cells[k + 1 - M2C][c];
		//stage1_k_cells(k, c, first_row, first_col, ext_requested_s_k, nof_bound_neigh, increm_boundary_neigh_row, increm_boundary_neigh_col, ext_max_cobound, ext_obj_dim_belonging, ext_aff_tree, ext_s_k, ext_t_k, k_k1_tree);
#pragma omp parallel for default(none) shared (k, c,first_row, first_col, size_ext_acc,  \
		increm_coboundary_neigh_col, increm_boundary_neigh_row, increm_boundary_neigh_col, \
		ext_acc,ext_acc_init,ext_s_k, ext_t_k, ext_requested_s_k, ext_aff_tree, ext_obj_dim_belonging, \
		ext_max_cobound,nof_bound_neigh)  					
		for (int row = first_row; row < (size_ext_acc[1 - M2C] - 2); row += 2) {
			for (int col = first_col; col < (size_ext_acc[2 - M2C] - 2); col += 2) {
				if (ext_aff_tree[row][col] <= k)
				{
					ext_requested_s_k[row][col] = INEXISTENT_ADJ_CELL;
					// Initialize these max  with the cell label itself :
					int L_k_ext = ext_acc[row][col]; // value of the k - cell
					int v_k_aff_tree = ext_aff_tree[row][col]; // value of the k - cell
					// // first we ask for the boundary cells(primal vector)

					int nof_Ak1_all = 0;
					int max_L_k1_all = 0;
					int v1_obj_dim = ext_obj_dim_belonging[row][col];
					int selected_bound_cell_max_Ak1_all = 0;//Inicializado por MJ

					for (int bound_numb = 1; bound_numb <= nof_bound_neigh[k + 1 - M2C]; bound_numb++) // for each bound of this k - cell
					{
						int row_increm = increm_boundary_neigh_row[bound_numb - M2C][c][k + 1 - M2C];
						int col_increm = increm_boundary_neigh_col[bound_numb - M2C][c][k + 1 - M2C];
						int row_bound = row + row_increm;
						int col_bound = col + col_increm;

						int st_ext_acc_bound = ext_max_cobound[row_bound][col_bound]; // state of the k - 1 cell, that is, the max value of all the k cell surrounding the k - 1 cell

						// we proceed to look for the max of the bounds

						int v1_bound_obj_dim = ext_obj_dim_belonging[row_bound][col_bound];
						if (v1_bound_obj_dim == v1_obj_dim)
						{
							// if (v1_bound_being_frontier == 0)
							// if (v_nof_regions_bound == v_k_nof_regions) // the bound must be of the same color than the cell
							// not necessary to check if the bound is in another tree.It is
							// impossible because of the tree building(from k to 1)

							nof_Ak1_all = nof_Ak1_all + 1;
							if (max_L_k1_all <= st_ext_acc_bound)
							{
								max_L_k1_all = st_ext_acc_bound;
								selected_bound_cell_max_Ak1_all = bound_numb;
							}//end
						} //end //end of if (v_nof_regions_bound == v_nof_regions) ...
					} //end //end of for bound_neigh_numb2 = 1:nof_bound_neigh(k + 1)

					// //

					if (max_L_k1_all > L_k_ext)
					{
						// continuous updating : k cell always asks, and even steals, for a k - 1
						ext_requested_s_k[row][col] = selected_bound_cell_max_Ak1_all;

						// old couple is broken, because a new one is requested
						// even if the old couple is the same.
						if (ext_s_k[row][col] != INEXISTENT_ADJ_CELL) // actually this 'if' is not necessary
						{
							ext_s_k[row][col] = INEXISTENT_ADJ_CELL;
							ext_acc[row][col] = ext_acc_init[row][col];
							ext_aff_tree[row][col] = INEXISTENT_AFF_TREE;

						}//end
					}// end							

				}
			}//end
		}//	end

	} //end // end of for c = 1 : nof_cells_per_nxel(k + 1] // for each cell of dim k in a nxel
}



void stage2_k1_cells(int k, int* nof_cells_per_nxel, int* nof_cobound_neigh, int** ext_acc, int** ext_acc_init,Orientation_4_adj_type** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** cobound_position_numb,int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, int** ext_max_cobound, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_s_k, Orientation_4_adj_type** ext_t_k) {
	for (int c = 1 - M2C; c < nof_cells_per_nxel[k + 1 - M2C]; c++) // for each position in a nxel(of k - 1 cell)
	{
		int first_row = ext_first_row_of_cells[k + 1 - M2C][c];
		int first_col = ext_first_col_of_cells[k + 1 - M2C][c];
		
#pragma omp parallel for default(none) shared (k,c, first_row, first_col, size_ext_acc,  \
		nof_cobound_neigh, increm_coboundary_neigh_col, cobound_position_numb, increm_boundary_neigh_row, increm_boundary_neigh_col, \
		ext_acc,ext_acc_init,ext_s_k, ext_t_k, ext_requested_s_k, ext_aff_tree, ext_obj_dim_belonging, \
		ext_max_cobound,increm_coboundary_neigh_row)  					
		for (int row = first_row; row < (size_ext_acc[1 - M2C] - 2); row += 2) {
			for (int col = first_col; col < (size_ext_acc[2 - M2C] - 2); col += 2) {
				int L_k1_ext = ext_acc[row][col]; // value of the k - 1 - cell
				int v_k1_aff_tree = ext_aff_tree[row][col]; // value of the k - 1 - cell

				// // first we ask for the boundary cells(primal vector)

				int max_L_k_req_s_k = 0; //min label of cobounds(k cells) that requested a coupling

				int nof_req_s_k = 0; //cobounds are k cells
				int selected_cobound_cell_max = 0 + INEXISTENT_ADJ_CELL;

				// not necesary to detect cycles
				int v2_obj_dim = ext_obj_dim_belonging[row][col];

				//Then any cobound cell is compared to find the lowest values
				for (int cobound_numb = 1; cobound_numb <= nof_cobound_neigh[k + 1 - M2C]; cobound_numb++)
				{
					// one of these cobound cells(of the previous cobound) is the
					// same k - cell whose maximum we are computing

					int row_increm = increm_coboundary_neigh_row[cobound_numb - M2C][c][k + 1 - M2C];
					int col_increm = increm_coboundary_neigh_col[cobound_numb - M2C][c][k + 1 - M2C];
					int row_cobound = row + row_increm;
					int col_cobound = col + col_increm;

					//     v_nof_regions_cobound = ext_nof_regions[row_bound][col_bound]; // nof_regions  of the k - cell
					int v_ext_acc_cobound = ext_acc[row_cobound][col_cobound]; // value of the k - cell // &

					int v2_cobound_obj_dim = ext_obj_dim_belonging[row_cobound][col_cobound];
					if (v2_cobound_obj_dim == v2_obj_dim)
					{
						// the cobound must be of the same color than the cell

						// exploring the max values for all the cobounds that have requested a
						// couple
						int requested2_by_cobound = ext_requested_s_k[row_cobound][col_cobound];

						if (requested2_by_cobound != +INEXISTENT_ADJ_CELL) // this k cell requested a coupling
						{
							int c_pos2_cobound = cobound_position_numb[cobound_numb - M2C][c][k + 1 - M2C];

							int req_row_increm = increm_boundary_neigh_row[requested2_by_cobound - M2C][c_pos2_cobound - M2C][k + 1 + 1 - M2C];
							int req_col_increm = increm_boundary_neigh_col[requested2_by_cobound - M2C][c_pos2_cobound - M2C][k + 1 + 1 - M2C];
							int req_row_bound = row_cobound + req_row_increm;
							int req_col_bound = col_cobound + req_col_increm;

							if ((req_row_bound == row) && (req_col_bound == col))
							{
								// the requested cell by the cobound(k - cell)  is this same k - 1 cell
								nof_req_s_k = nof_req_s_k + 1;
								if (max_L_k_req_s_k < v_ext_acc_cobound)
								{
									max_L_k_req_s_k = v_ext_acc_cobound;
									selected_cobound_cell_max = cobound_numb;
								} //end
							} //end

						} //end // end of if (requested2_by_cobound ~= 0)

					} //end // end of if (v_nof_regions_cobound == v_nof_regions)

				} //end //end of for cobound_neigh_numb

				if (nof_req_s_k > 0)
				{
					ext_t_k[row][col] = selected_cobound_cell_max;
					ext_aff_tree[row][col] = k + 1;
					// and infecting the k - 1 cell with its previously computed  state :
					ext_acc[row][col] = ext_max_cobound[row][col];
					//% old couple(maybe from another tree) is broken, because a new one is requested
					//% even if the old couple is the same.
					if (ext_s_k[row][col] != 0)// % actually this 'if' is not necessary
						ext_s_k[row][col] = 0;

				}
				else
				{
					int prev_t_k = ext_t_k[row][col];
					if ((prev_t_k != +INEXISTENT_ADJ_CELL))
						//if ((prev_t_k != +INEXISTENT_ADJ_CELL) && (prev_t_k != 0))
					{
						// it has previously a couple
						int row_increm = increm_coboundary_neigh_row[prev_t_k - M2C][c][k + 1 - M2C];
						int col_increm = increm_coboundary_neigh_col[prev_t_k - M2C][c][k + 1 - M2C];
						int row_cobound = row + row_increm;
						int col_cobound = col + col_increm;
						int prev_L_k_couple = ext_acc[row_cobound][col_cobound]; // it was the previous couple_bound_numb

						//two cases :
						if (prev_L_k_couple == L_k1_ext)
						{		// the k cell continues to be its couple : nothing to do
						}
						else
						{
							//no k cell requested this cell and its old couple has requested anohter
							//k - 1 cell(it has s_k == 0), thus it has no couple now :
							ext_t_k[row][col] = INEXISTENT_ADJ_CELL;
							ext_aff_tree[row][col] = INEXISTENT_AFF_TREE;
							// and infecting the k - 1 cell with its previously computed  state :
							ext_acc[row][col] = ext_acc_init[row][col];

						} //end
					} //end

				} //end


			}//end
		}//	end
	}
}



void stage3_k_cells(int k,int* nof_cells_per_nxel, int** ext_acc, Orientation_4_adj_type ** ext_requested_s_k, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** bound_position_numb, int*** cobound_row_increm, int*** cobound_col_increm, int ***increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int** ext_max_cobound, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_s_k, Orientation_4_adj_type** ext_t_k) {

	for (int c = 1 - M2C; c < nof_cells_per_nxel[k + 1 - M2C]; c++) // for each position in a nxel(of k - 1 cell)
	{
		//for c = 1:nof_cells_per_nxel[k + 1] // for each position in a nxel(of k - cell)
			//int first_row = ext_first_row_of_cells[k + 1][c];
		//int first_col = ext_first_col_of_cells[k + 1][c];


		int first_row = ext_first_row_of_cells[k + 1 - M2C][c];
		int first_col = ext_first_col_of_cells[k + 1 - M2C][c];

		//stage3_k_cells(k, c, first_row, first_col, ext_requested_s_k, increm_boundary_neigh_row, increm_boundary_neigh_col, ext_max_cobound, ext_obj_dim_belonging, ext_aff_tree, ext_s_k, ext_t_k, k_k1_tree);
		
#pragma omp parallel for default(none) shared (k, c, first_row, first_col, size_ext_acc,  \
		cobound_row_increm, cobound_col_increm,  bound_position_numb, \
		ext_acc,ext_s_k, ext_t_k, ext_requested_s_k, ext_aff_tree, ext_obj_dim_belonging, \
		ext_max_cobound,increm_boundary_neigh_row,increm_boundary_neigh_col)  					
		for (int row = first_row; row < (size_ext_acc[1 - M2C] - 2); row += 2) {
			for (int col = first_col; col < (size_ext_acc[2 - M2C] - 2); col += 2) {
				//for col = first_col:2 : (size_ext_acc[2] - 2)
				//for row = first_row : 2 : (size_ext_acc[1] - 2)

				if (ext_aff_tree[row][col] <= k)
				{
		
					int L3_k_ext = ext_acc[row][col]; // value of the k - cell
					int v3_k_aff_tree = ext_aff_tree[row][col]; // value of the k - cell
					//cout << k << ";" << c << ";" << row << ";" << col << endl;
					// they are inititalized here to prevent errors.

					int bound3_numb_req = ext_requested_s_k[row][col]; //for the previously requested

					if (bound3_numb_req == INEXISTENT_ADJ_CELL) {
						// the k cell does not request : if its couple was stolen, it must reset its state
						int s3_k_ext = ext_s_k[row][col]; // bound of previous couple
						//Condición modificada por MJMF
						if (s3_k_ext != (INEXISTENT_ADJ_CELL) && (s3_k_ext != 0))
						{
							int row_increm_prev = increm_boundary_neigh_row[s3_k_ext - M2C][c][k + 1 - M2C];
							int col_increm_prev = increm_boundary_neigh_col[s3_k_ext - M2C][c][k + 1 - M2C];


							int row_bound_prev = row + row_increm_prev;
							int col_bound_prev = col + col_increm_prev;

							// the grant may has changed.The current one is :
							int cobound3_numb_current = ext_t_k[row_bound_prev][col_bound_prev];

							int c_position_bound_current = bound_position_numb[s3_k_ext - M2C][c][k + 1 - M2C];

							int inc_row_k_grant = cobound_row_increm[cobound3_numb_current - M2C][c_position_bound_current - M2C][k + 1 - 1 - M2C];
							int inc_col_k_grant = cobound_col_increm[cobound3_numb_current - M2C][c_position_bound_current - M2C][k + 1 - 1 - M2C];
							int row_cobound_grant_current = row_bound_prev + inc_row_k_grant;
							int col_cobound_grant_current = col_bound_prev + inc_col_k_grant;

							if ((row_cobound_grant_current == row) && (col_cobound_grant_current == col))
							{
								int kk = 9999; // nothing to do.THis cell still is the granted one
							}
							else 		// THis cell couple was stolen : reseting
							{
								ext_s_k[row][col] = INEXISTENT_ADJ_CELL;
								ext_aff_tree[row][col] = INEXISTENT_AFF_TREE;
								//             ext_acc[row][col] = ext_acc_init[row][col]; // @if this value
								//             is reseted, the CA enters in resonance !!!@
							}
						}//end
					}//end

					else {
						// the k cell request sth

						//// first we ask if the request was granted by the boundary cell

						int row_increm = increm_boundary_neigh_row[bound3_numb_req - M2C][c][k + 1 - M2C];
						int col_increm = increm_boundary_neigh_col[bound3_numb_req - M2C][c][k + 1 - M2C];
						int row_bound_req = row + row_increm;
						int col_bound_req = col + col_increm;

						int cobound3_numb_grant = ext_t_k[row_bound_req][col_bound_req]; //the answer of the previously requested bound cell
						
						int nof3_valid_bounds = 0;
						// we proceed to look for the max of the bounds, and to detect cycles

						//// cycle CHECKING IS NOT NECESSARY

						//// infection and accepting the couple

						// it is infected by its couple($k - 1$ - cell)
						// label in case it was bigger.
						if (cobound3_numb_grant != INEXISTENT_ADJ_CELL) {
							int c_position_bound = bound_position_numb[bound3_numb_req - M2C][c][k + 1 - M2C];
							int inc_row_k_grant = cobound_row_increm[cobound3_numb_grant - M2C][c_position_bound - M2C][k + 1 - 1 - M2C];
							int inc_col_k_grant = cobound_col_increm[cobound3_numb_grant - M2C][c_position_bound - M2C][k + 1 - 1 - M2C];

							// detecting if old couple remains :
							if ((row_bound_req + inc_row_k_grant == row) && (col_bound_req + inc_col_k_grant == col)) {
								// note that RESETING THE REQUESTs is done at the beginning of each step

								ext_acc[row][col] = ext_max_cobound[row_bound_req][col_bound_req]; //@9
								if (v3_k_aff_tree == INEXISTENT_AFF_TREE) {
									// Each $k$ - cell that was granted by some $k - 1$ - cell accepts the
									// provisional couple if this does not introduce a cycle.
									ext_s_k[row][col] = bound3_numb_req;
									ext_aff_tree[row][col] = k;
								}
								else {
									//If this $k$ - cell was previously coupled(with a different
									//$k - 1$ - cell), it breaks the old couple and points to the new one.
									//This last situation can happen for an upgrade or a promotion of stage 1.
									ext_s_k[row][col] = bound3_numb_req;
									// ext_aff_tree[row][col] = k;  // it was already set

								}// end

							}// end
						}// end

					}// end


				}
			}//end
		}//	end


	} //end // end of for c = 1 : nof_cells_per_nxel(k + 1] // for each cell of dim k in a nxel

}


void ending_k2_cells_only_tail_cells(int k, int c, int row, int col, int* nof_cobound_neigh,int *** increm_boundary_neigh_row, int ***increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col, Orientation_4_adj_type** ext_obj_dim_belonging, Orientation_4_adj_type** ext_aff_tree, Orientation_4_adj_type** ext_being_tail)
{
	//#include "ending_k2_cells_only_tail_cells.cpp" // it needs : row, col, k, c and global matrices
	int ppacum_nof_valid_cobounds = 0;
	// Any cobound cell is compared to find the biggest values
	for (int cobound_numb = 1; cobound_numb <= nof_cobound_neigh[k + 1 - M2C]; cobound_numb++) {
		// one of these cobound cells(of the previous cobound) is the
		// same k - cell whose maximum we are computing

		int row_increm = increm_coboundary_neigh_row[cobound_numb - M2C][c][k + 1 - M2C];
		int col_increm = increm_coboundary_neigh_col[cobound_numb - M2C][c][k + 1 - M2C];
		int row_cobound = row + row_increm;
		int col_cobound = col + col_increm;

		int pp_ext_cobound_obj_dim = ext_obj_dim_belonging[row_cobound][col_cobound];

		// we need for ext_being_tail :
		if ((ext_aff_tree[row_cobound][col_cobound] <= k + 1) && (pp_ext_cobound_obj_dim != (k + 1))) {
			ppacum_nof_valid_cobounds = ppacum_nof_valid_cobounds + 1;
		}//end

	} //end// end of for cobound_neigh_numb

	if ((ppacum_nof_valid_cobounds == 1) && (ext_obj_dim_belonging[row][col] == (k + 2))) {
		ext_being_tail[row][col] = 1;
	}//end	
}



