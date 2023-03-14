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
int num_th;
////////////////////
////////////////////

#include "CA_HSF_2D_C_4stages_obj_dim.h"



using namespace std;
using namespace cv;


//#include "opencv/cv.h"
//#include "opencv/highgui.h"
#include <direct.h>


int n_rows = 6;
int n_cols = 6;

int n_rows_ext;
int n_cols_ext;
///////////
//int size_acc[dim_I] = { n_rows * 2 + 1, n_cols * 2 + 1 }; // = (size(I) * 2 + 1);
int size_acc [dim_I];

int	total_nof_cells_ext;


int ***increm_boundary_neigh_row;
int ***increm_boundary_neigh_col;

int ***increm_coboundary_neigh_row;
int ***increm_coboundary_neigh_col;

int ***bound_position_numb; // (bound_numb, c, k + 1)
int ***cobound_position_numb; // (cobound_numb, c, k + 1)

int ***cobound_row_increm; // (cobound3_numb_grant, c_position_bound, k + 1 -1);
int ***cobound_col_increm; // (cobound3_numb_grant, c_position_bound, k + 1 -1);

int* nof_cells_per_nxel;// [dim_I + 1] ;
int* nof_bound_neigh;// [dim_I + 1] ;
int* nof_cobound_neigh;// [dim_I + 1] ;


int** ext_first_row_of_cells; //[dim_I + 1] [MAX_NOF_K_CELLS] ;
int** ext_first_col_of_cells;// [dim_I + 1] [MAX_NOF_K_CELLS] ;


#define DIM_0CELLS 0
#define DIM_1CELLS 1
#define DIM_2CELLS 2

#define NOF_GHOST_COLS_ACC 2
#define NOF_GHOST_ROWS_ACC 2  // THERE ARE TWO ghost cols and rows
//PP_MAX_NOF_CRIT_CELLS (EXT_N_ROWS* EXT_N_COLS)
//#define PP_MAX_NOF_CRIT_CELLS 128  // THERE ARE TWO ghost cols and rows


int size_ext_acc[dim_I]; // = size(ext_acc_init);

int** ext_acc;
int** ext_acc_init; // size_acc + 4
Orientation_4_adj_type** ext_being_tail;// = ZEROS(size_ext_acc);
int** ext_acc_previous;// [n_rows_ext] [n_cols_ext] ;
int** ext_max_cobound;// [n_rows_ext] [n_cols_ext] ;
Orientation_4_adj_type** ext_max_cobound_numb;// [n_rows_ext] [n_cols_ext] ;
Orientation_4_adj_type** ext_requested_s_k;// [n_rows_ext] [n_cols_ext] ;


// This matrix will hold which tree each cell belongs to(or which affiliation has) :
Orientation_4_adj_type** ext_aff_tree;// [n_rows_ext] [n_cols_ext] ;
Orientation_4_adj_type** ext_s_k;// [n_rows_ext] [n_cols_ext] ; //bound of the k cell for a primal vector
Orientation_4_adj_type** ext_t_k;// [n_rows_ext] [n_cols_ext] ; //cobound of the k -1 cell for a primal vector

Orientation_4_adj_type** ext_obj_dim_belonging;//[n_rows_ext][n_cols_ext]; /**/
int max_nof_crit_cells = 128;
//int crit_cell_list[dim_I + 1][PP_MAX_NOF_CRIT_CELLS];  //this must be done with std vector . Or even with a dynamic matrix for the 3 types of cells: e.g. crit_cell_list[dim_I+1][dynamic_length]
//they contain the labels of the crit cells
int** crit_cell_list; //[dim_I + 1] [PP_MAX_NOF_CRIT_CELLS] ;
int** crit_cell_nof_cob; //[dim_I] [PP_MAX_NOF_CRIT_CELLS] ;  //this must be done with std vector . Or even with a dynamic matrix for the 3 types of cells: e.g. crit_cell_nof_cob[dim_I][dynamic_length]
//they contain the nof  valid cobounds of each crit cells
// 2 cells have no cobound: thus the first size is [dim_I] instead of [dim_I+1]
/* In order to mark crit cells that represent a hole, the nof cobounds are added this value: 2 * (dim_I)
This is to codify the hole condition in this crit_cell_nof_cob instead of defining another mattrix sth like "crit_cell_hole_condition"
*/

int*** crit_cell_cob_list;// [dim_I] [PP_MAX_NOF_CRIT_CELLS] [2 * (dim_I - DIM_0CELLS)] ;  //this must be done with std vector . Or even with a dynamic matrix for the 3 types of cells: e.g. crit_cell_cob_list[dim_I][dynamic_length][2 * (dim_I - 0)] ; third size because the max nof cobounds of any cell
//Note that 2 cells have no cobound
// these are the labels of the cobounds (of the crit cells)

ImageType** I;



void createBidimArrays();
/////////////////////
////////////////////
// Get binary image given a image's FileName; 
bool getBinaryImage(const string FileName, cv::Mat1b& binaryMat) {

	// Image load
	cv::Mat image;

	//CV_LOAD_IMAGE_GRAYSCALE
	image = cv::imread(FileName, cv::IMREAD_GRAYSCALE);   // Read the file

	// Check if image exist
	if (image.empty())
		return false;

	// Adjust the threshold to actually make it binary
	//CV_THRESH_BINARY
	threshold(image, binaryMat, 100, 1, cv::THRESH_BINARY);

	return true;
}

bool get_color_img(const string& filename, Mat& color_mat)
{
	// Image load
	color_mat = cv::imread(filename, IMREAD_COLOR);   // Read the file
												  // Check if image exist
	if (color_mat.empty()) {
		return false;
	}

	return true;
}
bool get_grayscale_img(const string& filename, Mat& color_mat)
{
	// Image load
	color_mat = cv::imread(filename, IMREAD_GRAYSCALE);   // Read the file
												  // Check if image exist
	if (color_mat.empty()) {
		return false;
	}

	return true;
}

bool get_and_convert_into_grayscale_img(const string& filename, Mat1b& image_gray)
{
	// Image load
	Mat img_color;
	if (!get_color_img(filename, img_color)) {
		cout << "Unable to check on '" + filename + "', file does not exist" << endl;
		exit(0);
	}
	image_gray = Mat1b(img_color.size());
	cvtColor(img_color, image_gray, cv::COLOR_BGR2GRAY);

	return true;
}
////////////////////////////////////////////
void createBidimArrays()
{


	n_rows_ext = n_rows * 2 + 1 + 4;
	n_cols_ext = n_cols * 2 + 1 + 4;

	size_ext_acc[0] = n_rows_ext;
	size_ext_acc[1] = n_cols_ext;
	size_acc[0] = n_rows * 2 + 1;
	size_acc[1] = n_cols * 2 + 1;

	ext_acc = new int* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ;
// extended matrices :
	ext_acc_init = new int* [n_rows_ext]; //[n_rows_ext] [n_cols_ext] ; // size_acc + 4
	

	ext_being_tail = new Orientation_4_adj_type* [n_rows_ext];// = ZEROS(size_ext_acc);
	ext_acc_previous = new int* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ;
	ext_max_cobound = new int* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ;
	ext_max_cobound_numb = new Orientation_4_adj_type* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ;
	ext_requested_s_k = new Orientation_4_adj_type* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ;


	// This matrix will hold which tree each cell belongs to(or which affiliation has) :
	ext_aff_tree = new Orientation_4_adj_type* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ;
	ext_s_k = new Orientation_4_adj_type* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ; //bound of the k cell for a primal vector
	ext_t_k = new Orientation_4_adj_type* [n_rows_ext];// [n_rows_ext] [n_cols_ext] ; //cobound of the k -1 cell for a primal vector

	ext_obj_dim_belonging = new Orientation_4_adj_type* [n_rows_ext];//[n_rows_ext][n_cols_ext];
	

	for (int i = 0; i < n_rows_ext; i++) {
		ext_acc[i] = new int[n_cols_ext]();// [n_rows_ext] [n_cols_ext] ;
		ext_acc_init[i] = new int[n_cols_ext]();
		
		ext_being_tail[i] = new Orientation_4_adj_type[n_cols_ext]();
		ext_acc_previous[i] = new int[n_cols_ext]();// [n_rows_ext] [n_cols_ext] ;
		ext_max_cobound[i] = new int[n_cols_ext]();// [n_rows_ext] [n_cols_ext] ;
		ext_max_cobound_numb[i] = new Orientation_4_adj_type [n_cols_ext]();// [n_rows_ext] [n_cols_ext] ;
		ext_requested_s_k[i] = new Orientation_4_adj_type [n_cols_ext]();// [n_rows_ext] [n_cols_ext] ;
		ext_aff_tree[i] = new Orientation_4_adj_type [n_cols_ext]();// [n_rows_ext] [n_cols_ext] ;
		ext_s_k[i] = new Orientation_4_adj_type [n_cols_ext]();// [n_rows_ext] [n_cols_ext] ; //bound of the k cell for a primal vector
		ext_t_k[i] = new Orientation_4_adj_type [n_cols_ext]();// [n_rows_ext] [n_cols_ext] ; //cobound of the k -1 cell for a primal vector
		ext_obj_dim_belonging[i] = new Orientation_4_adj_type[n_cols_ext]();

	}
	
	nof_cells_per_nxel= new int [dim_I + 1]();
	nof_bound_neigh= new int[dim_I + 1]();
	nof_cobound_neigh =new int[dim_I + 1]();

	ext_first_row_of_cells = new int* [dim_I + 1];  //[dim_I + 1] [MAX_NOF_K_CELLS] ;
	ext_first_col_of_cells = new int* [dim_I + 1];
	for (int i = 0; i < (dim_I + 1); i++) {
		ext_first_row_of_cells [i] = new int[MAX_NOF_K_CELLS]();
		ext_first_col_of_cells [i] = new int[MAX_NOF_K_CELLS]();
	}
	increm_boundary_neigh_row=new int **[MAX_NOF_BOUNDS];
	increm_boundary_neigh_col = new int**[MAX_NOF_BOUNDS];
	increm_coboundary_neigh_row = new int** [MAX_NOF_BOUNDS];
	increm_coboundary_neigh_col = new int**[MAX_NOF_BOUNDS];

	bound_position_numb = new int** [MAX_NOF_BOUNDS];
	cobound_position_numb = new int** [MAX_NOF_BOUNDS];
	cobound_row_increm = new int** [MAX_NOF_BOUNDS];
	cobound_col_increm = new int** [MAX_NOF_BOUNDS];
	for (int i = 0; i < MAX_NOF_BOUNDS; i++) {
		increm_boundary_neigh_row[i] = new int* [MAX_NOF_K_CELLS];
		increm_boundary_neigh_col[i]= new int* [MAX_NOF_K_CELLS];
		increm_coboundary_neigh_row[i] = new int* [MAX_NOF_K_CELLS];
		increm_coboundary_neigh_col[i] = new int* [MAX_NOF_K_CELLS];
		bound_position_numb[i] = new int* [MAX_NOF_K_CELLS];
		cobound_position_numb[i] = new int* [MAX_NOF_K_CELLS];
		cobound_row_increm[i] = new int* [MAX_NOF_K_CELLS];
		cobound_col_increm[i] = new int* [MAX_NOF_K_CELLS];
		for (int j = 0; j < MAX_NOF_K_CELLS; j++) {
			increm_boundary_neigh_row[i][j] = new int [dim_I + 1]();
			increm_boundary_neigh_col[i][j] = new int[dim_I + 1]();
			increm_coboundary_neigh_row[i][j] = new int[dim_I + 1]();
			increm_coboundary_neigh_col[i][j] = new int[dim_I + 1]();
			bound_position_numb[i][j] = new int [dim_I + 1]();
			cobound_position_numb[i][j] = new int [dim_I + 1]();
			cobound_row_increm[i][j] = new int [dim_I + 1]();
			cobound_col_increm[i][j] = new int [dim_I + 1]();
		}
	}
}
void createCrit_Cell_Arrays()
{
	max_nof_crit_cells = n_rows_ext * n_cols_ext;
	crit_cell_list = new int* [dim_I + 1]; //[dim_I + 1] [PP_MAX_NOF_CRIT_CELLS] ;
	crit_cell_nof_cob = new int* [dim_I + 1]; //[dim_I] [PP_MAX_NOF_CRIT_CELLS] ;  //this must be done with std vector . Or even with a dynamic matrix for the 3 types of cells: e.g. crit_cell_nof_cob[dim_I][dynamic_length]
	crit_cell_cob_list = new int** [dim_I + 1];// [dim_I] [PP_MAX_NOF_CRIT_CELLS] [2 * (dim_I - DIM_0CELLS)] ;  //this must be done with std vector . Or even with a dynamic matrix for the 3 types of cells: e.g. crit_cell_cob_list[dim_I][dynamic_length][2 * (dim_I - 0)] ; third size because the max nof cobounds of any cell
	for (int i = 0; i < (dim_I + 1); i++) {
		crit_cell_list[i] = new int[max_nof_crit_cells]();
		crit_cell_nof_cob[i] = new int[max_nof_crit_cells]();
	}
	for (int i = 0; i < (dim_I + 1); i++) {
		crit_cell_cob_list[i] = new int* [max_nof_crit_cells];
		for (int j = 0; j < max_nof_crit_cells; j++)
		{
			crit_cell_cob_list[i][j] = new int[2 * (dim_I - DIM_0CELLS)]();
		}
	}
}
void deleteCrit_Cells_Arrays()
{
	for (int i = 0; i < (dim_I + 1); i++) {
		delete[] crit_cell_list[i];
		delete[] crit_cell_nof_cob[i];
	}
	delete[] crit_cell_list;
	delete[] crit_cell_nof_cob;
	for (int i = 0; i < (dim_I + 1); i++) {

		for (int j = 0; j < max_nof_crit_cells; j++)
		{
			delete[] crit_cell_cob_list[i][j];
		}
		delete[] crit_cell_cob_list[i];
	}
	delete[] crit_cell_cob_list;
}
void setInitValues()
{

	for (int i = ext_first_row_of_cells[dim_I + 1 - M2C][1 - M2C]; i < n_rows_ext - 2; i += 2)
		for (int j = ext_first_col_of_cells[dim_I + 1 - M2C][1 - M2C]; j < n_cols_ext - 2; j += 2)
			ext_obj_dim_belonging[i][j] = dim_I;


	//// 1 - cell in the frontier with external dummy object are 1 - object obviously
	int c = 1;  // //for each position in a pixel(of 1 - cell)

	int ppfr = ext_first_row_of_cells[dim_I + 1 - 1 - M2C][c - M2C];
	int ppfc = ext_first_col_of_cells[dim_I + 1 - 1 - M2C][c - M2C];
	//ext_obj_dim_belonging(ppfr, ppfc: 2 : end - ppfc + 1) = dim_I - 1;
	int i = ppfr;
	for (int j = ppfc; j < n_cols_ext - (ppfc); j += 2)
		ext_obj_dim_belonging[i][j] = dim_I - 1;
	//ext_obj_dim_belonging(end - ppfr + 1, ppfc: 2 : end - ppfc + 1) = dim_I - 1;
	i = (n_rows_ext - M2C) - (ppfr)+1 - M2C;
	for (int j = ppfc; j < n_cols_ext - (ppfc); j += 2)
		ext_obj_dim_belonging[i][j] = dim_I - 1;


	c = 2;  // //for each position in a pixel(of 1 - cell)

	ppfr = ext_first_row_of_cells[dim_I + 1 - 1 - M2C][c - M2C];
	ppfc = ext_first_col_of_cells[dim_I + 1 - 1 - M2C][c - M2C];

	//ext_obj_dim_belonging(ppfr: 2 : end - ppfr + 1, ppfc) = dim_I - 1;
	int j = ppfc;
	for (int i = ppfr; i < n_rows_ext - (ppfr)+1; i += 2) {
		ext_obj_dim_belonging[i][j] = dim_I - 1;
	}

	//ext_obj_dim_belonging(ppfr: 2 : end - ppfr + 1, end - ppfc + 1) = dim_I - 1;
	j = (n_cols_ext - M2C) - (ppfc)+1 - M2C;
	for (int i = ppfr; i < n_rows_ext - (ppfr)+1; i += 2) {
		ext_obj_dim_belonging[i][j] = dim_I - 1;
	}

	////external dummy object is monochrome :
	//ext_obj_dim_belonging(1:2, 1 : end) = dim_I;
	for (int i = 1 - M2C; i <= 2 - M2C; i++) {
		fill_n(ext_obj_dim_belonging[i], n_cols_ext, dim_I);
	}
	//ext_obj_dim_belonging(end - 1:end, 1 : end) = dim_I;
	for (int i = n_rows_ext - 2; i < n_rows_ext; i++) {
		fill_n(ext_obj_dim_belonging[i], n_cols_ext, dim_I);
	}
	//ext_obj_dim_belonging(1:end, 1 : 2) = dim_I;
	for (int i = 0; i < n_rows_ext; i++) {
		for (int j = 1 - M2C; j < 2; j++) {
			ext_obj_dim_belonging[i][j] = dim_I;
		}
	}
	//ext_obj_dim_belonging(1:end, end - 1 : end) = dim_I;
	for (int i = 0; i < n_rows_ext; i++) {
		for (int j = n_cols_ext - 2; j < n_cols_ext; j++) {
			ext_obj_dim_belonging[i][j] = dim_I;
		}
	}

	//init: matrices set to 0:
	for (int i = 0; i < n_rows_ext; i++) {
		for (int j = 0; j < n_cols_ext; j++) {
			ext_aff_tree[i][j] = INEXISTENT_AFF_TREE;
			ext_s_k[i][j] = INEXISTENT_ADJ_CELL; //bound of the k cell for a primal vector
			ext_t_k[i][j] = 0 + INEXISTENT_ADJ_CELL; //cobound of the k -1 cell for a primal vector
		}
	}

	//it determines if it is a tail cell. 1 means that it is; bigger than 1 means that it is not  nxels are never tails
			//init: matrices set to 0:
	for (int i = 0; i < n_rows_ext; i++) {
		fill_n(ext_being_tail[i], n_cols_ext, 0);
	}
	if (dim_I == 2)
	{
		//	ext_acc_init(3:end - 2, 3 : end - 2) = pp_ext_acc_init(3:end - 2, 3 : end - 2); // @ only 2D
			// // //     ext_acc_color(3:end - 2, 3 : end - 2) = pp_colors; //&
		//for (int i = 0; i < size_ext_acc[0] * size_ext_acc[1]; i++)
		//	*((int *)ext_acc_init + i) = 0;
		//memset(ext_acc_init, 0, n_rows_ext*n_cols_ext);
		for (int j = 3 - M2C; j < size_ext_acc[1] - 2; j++)
			for (int i = 3 - M2C; i < size_ext_acc[0] - 2; i++)
				ext_acc_init[i][j] = (j)*size_ext_acc[0] + i + 1;
	}
	else
	{
		exit(0);
		//input('*** NOT done yet. Press Intro to finish');
	}
	//ext_acc = ext_acc_init;
	for (int i = 0; i < n_rows_ext; i++)  for (int j = 0; j < n_cols_ext; j++)
		ext_acc[i][j] = ext_acc_init[i][j];



	total_nof_cells_ext = 1;
	for (int d = 1 - M2C; d < dim_I; d++)
		//	for d = 1:dim_I
		total_nof_cells_ext = total_nof_cells_ext * size_ext_acc[d];
}

void DeleteArrasys() {
	for (int i = 0; i < n_rows_ext; i++) {
		delete[] ext_acc [i];
		delete[] ext_acc_init[i];
		delete[] ext_being_tail[i];
		delete[] ext_acc_previous[i];
		delete[] ext_max_cobound [i];
		delete[] ext_max_cobound_numb[i];
		delete[] ext_requested_s_k[i];
		delete[] ext_aff_tree[i];
		delete[] ext_s_k[i];
		delete[] ext_t_k[i];
		delete[] ext_obj_dim_belonging[i];
		
	}
	for (int i = 0; i < n_rows; i++) {
		delete[] I[i];
	}

	delete[] ext_acc;
	delete[] ext_acc_init;
	delete[] ext_being_tail;
	delete[] ext_acc_previous;
	delete[] ext_max_cobound;
	delete[] ext_max_cobound_numb;
	delete[] ext_requested_s_k;
	delete[] ext_aff_tree;
	delete[] ext_s_k;
	delete[] ext_t_k;
	delete[] ext_obj_dim_belonging;
	delete[] I;
	
	
	delete[] nof_cells_per_nxel;
	delete[] nof_bound_neigh;
	delete[] nof_cobound_neigh;

	for (int i = 0; i < (dim_I + 1); i++) {
		delete[] ext_first_row_of_cells[i];
		delete[] ext_first_col_of_cells[i];
	}
	delete[] ext_first_row_of_cells;
	delete[] ext_first_col_of_cells;
	for (int i = 0; i < MAX_NOF_BOUNDS; i++) {
		
		for (int j = 0; j < MAX_NOF_K_CELLS; j++) {
			delete [] increm_boundary_neigh_row[i][j];
			delete [] increm_boundary_neigh_col[i][j];
			delete [] increm_coboundary_neigh_row[i][j];
			delete [] increm_coboundary_neigh_col[i][j];
			delete[] bound_position_numb[i][j];
			delete[] cobound_position_numb[i][j];
			delete[] cobound_row_increm[i][j];
			delete[] cobound_col_increm[i][j];
		}
		delete [] increm_boundary_neigh_row[i];
		delete [] increm_boundary_neigh_col[i];
		delete [] increm_coboundary_neigh_row[i];
		delete [] increm_coboundary_neigh_col[i];
		delete[] bound_position_numb[i];
		delete[] cobound_position_numb[i];
		delete[] cobound_row_increm[i];
		delete[] cobound_col_increm[i];
		}
	delete [] increm_boundary_neigh_row;
	delete[] increm_boundary_neigh_col;
	delete [] increm_coboundary_neigh_row;
	delete [] increm_coboundary_neigh_col;
	delete[] bound_position_numb;
	delete[] cobound_position_numb;
	delete[] cobound_row_increm;
	delete[] cobound_col_increm;

}

ImageType* procesarLinea(string linea, int ncs) {
	string tmp;
	stringstream sstream(linea);
	ImageType*  ims = new ImageType[ncs]();
	for (int c = 0; c < ncs;c++) {
		getline(sstream, tmp, ',');
		ims[c] = stoi(tmp);
	}
	return ims;
}

ImageType **loadSynthIm(string path) {
	
		string fila,tmp;
		int nr, r=0;
		int nc;
		ifstream fe;
		
		

		ImageType **im;
		fe.open(".\\"+path);
		if (fe.is_open())
		{

			getline(fe, fila);
			stringstream sstream(fila);
			getline(sstream, tmp, ';');
			nr = stoi(tmp);
			getline(sstream, tmp);
			nc = stoi(tmp);

			im = new ImageType *[nr];
			getline(fe, fila);
			im[r++]=procesarLinea(fila, nc);
			while (!fe.eof() && (r < nr)) {
				getline(fe, fila);
				im[r++] = procesarLinea(fila, nc);
			}
			fe.close();
		}
		else
		{
			std::cout << "File not Found!" << endl;
			exit(-1);
		}
		n_rows = nr;
		n_cols = nc;

		return im;

}
int main(int argc, char** argv) {

	string fname;
	boolean synthetic = false;
	int num_repet = MIN_NOF_REPETITIONS;
	
	if (argc >= 2)
	{
		fname =  string(argv[1]);
		num_repet = atoi(argv[2]);
	}
	else
	{
		std::cout << "Insert file name: ";
		std::cin >> fname;
	}
	if (fname.find(".txt")!= std::string::npos)
	{
		synthetic = true;
		cout << "Loading synthetic image from file " << endl;
	}
	else
		cout << "Loading real image from file " << endl;


	void 	ca_preliminars(int* nof_cells_per_nxel, int* nof_bound_neigh, int* nof_cobound_neigh, int** ext_first_row_of_cells, int** ext_first_col_of_cells, int*** bound_position_numb, int*** cobound_position_numb, int*** cobound_row_increm, int*** cobound_col_increm,int*** increm_boundary_neigh_row, int*** increm_boundary_neigh_col, int*** increm_coboundary_neigh_row, int*** increm_coboundary_neigh_col);

	void printing_machine_info(void);
	printing_machine_info();

	
	Mat1b  image_gray;
	//cv::Mat1b binaryImg;



#define NOF_COLORS 2

	int kk = 111;
	//scanf("%d", &kk);
	// filling  .data with rand() numbers to build a synthetic image 
	srand(kk);  // this guarantees the same random matrix
	/*for (RowColType r = 0; r < n_rows; r++) {
		for (RowColType c = 0; c < n_cols; c++) {
			I[r][c] = 1 + rand() % NOF_COLORS;
		}
	}*/
	if (synthetic)
	{
		I = loadSynthIm(fname);
		cout << "--Matrix size r, c =  " << n_rows << "," << n_cols << endl;

	}
	else
	{

		Mat img_color;
		char* cwd = _getcwd(NULL, 0);
		fname = "\\" + fname;
		
		
		if (!get_and_convert_into_grayscale_img(cwd + fname, image_gray)) {
			std::cout << "Unable to check on '" + fname + "', file does not exist" << endl;
			free(cwd);
			exit(0);
		}	
		
		free(cwd);
		//std::cout << "  -- Matrix size r,c = " << image_gray.rows << ", " << image_gray.cols << endl;


		n_rows = image_gray.rows;
		n_cols = image_gray.cols;
		cv::imshow("Image",image_gray);
		waitKey(0);
		destroyAllWindows();

		I = new  ImageType * [n_rows];
	

		for (int r = 0; r < n_rows; r++) {
			I[r] = new  ImageType[n_cols]();
			for (int c = 0; c < n_cols; c++) {
				I[r][c] = image_gray(r, c);
			}
		}

	}
	std::cout << "------------Testing Image : ------------ " << fname.c_str() << endl;
	std::cout << "ROWS: " << n_rows << endl;
	std::cout << "COLS:" << n_cols << endl;
	
	createBidimArrays();
	
	
	for (num_th = 1; num_th <= MAX_NUM_THREADS; num_th++) {
		std::cout << " ------------ Testing " << num_th << "threads. ------------ "  << endl;
		omp_set_num_threads(num_th);

		double t_inc_min_complete = 1.0e+35, t_inc_complete, t0_complete;

		

		// a group of repetitions to compute several times the execution time and print the minimum
		//   (with more threads, more repetitions are done)


		for (int rep = 0; rep < num_repet; rep++) {
			
			/// MEASUREMENTS BEGIN
			// FUNCTIONS TO BE TESTED :
			
			t0_complete = omp_get_wtime();
			ca_preliminars( nof_cells_per_nxel, nof_bound_neigh, nof_cobound_neigh, ext_first_row_of_cells, ext_first_col_of_cells, bound_position_numb, cobound_position_numb, cobound_row_increm, cobound_col_increm,increm_boundary_neigh_row, increm_boundary_neigh_col, increm_coboundary_neigh_row, increm_coboundary_neigh_col);
			setInitValues();

			ca_infectious_process(I,total_nof_cells_ext,nof_cells_per_nxel, nof_bound_neigh, nof_cobound_neigh, ext_acc, ext_acc_init,ext_acc_previous,ext_requested_s_k, ext_first_row_of_cells, ext_first_col_of_cells, bound_position_numb,  cobound_position_numb, cobound_row_increm, cobound_col_increm,increm_boundary_neigh_row, increm_boundary_neigh_col, increm_coboundary_neigh_row, increm_coboundary_neigh_col,ext_max_cobound, ext_max_cobound_numb, ext_obj_dim_belonging, ext_aff_tree, ext_s_k, ext_t_k, ext_being_tail);
				
			
			t_inc_complete = omp_get_wtime() - t0_complete;
				
			cout << rep << ": "<< t_inc_complete << endl;
			// choosing minimum time:
			if (t_inc_min_complete > t_inc_complete) t_inc_min_complete = t_inc_complete;

		}  //end of for (int rep = 0; rep < num_repet; rep++)

	void hgig_counting_crit_cells(int*);
	int nof_crit_cells[dim_I + 1];
	hgig_counting_crit_cells(nof_crit_cells);

	printf(" nof_crit_cells 0,1,2: %u, %u, %u \n", nof_crit_cells[0], nof_crit_cells[1], nof_crit_cells[2]);

	printf(" * Total minimum Time (%3d threads):  %lf \n", num_th, t_inc_min_complete);

	}
	
	DeleteArrasys();
	getchar();
	return 0;
}



////////////////


//////////////////////
void  hgig_counting_crit_cells(int* nof_crit) {
	//must return: int nof_crit_cells[dim_I + 1] ;
// the rest of matrixes should be embedded into a class, so that each object would contain the info of a crit cell
	nof_crit[DIM_0CELLS] = 0;
	nof_crit[DIM_1CELLS] = 0;
	nof_crit[DIM_2CELLS] = 0; //this must be done with std vector . Or even with a dynamic matrix for the 3 types of cells: e.g. nof_crit_cells[dynamic_length][dim_I+1]
	//they contain the number of crit cells
	// 
	createCrit_Cell_Arrays();
	
	//obtained from:  int ext_aff_tree[n_rows_ext][n_cols_ext];
		// de momento no le pongo omp porque no esta dentro de la medicion de tieimpo s xxxxx
	for (int col = NOF_GHOST_COLS_ACC; col < n_cols_ext - NOF_GHOST_COLS_ACC; col++) { // The loop order "col outer, row inner" matches that of matlab; better for debugging 
		for (int row = NOF_GHOST_COLS_ACC; row < n_rows_ext - NOF_GHOST_COLS_ACC; row++) { //beginning end ending in 2 or -2 because of the 2 ghost rows
				//int crit = (ext_aff_tree[row][col] == 0);
			if (ext_aff_tree[row][col] == 0) {
				//the dim of the first cell is 0 (upper left corner)
				// if row is even we are in 0 or 1 cells
				// if row is odd we are in 1 or 2 cells
				// Thus, these  &1 computes the cell dim:
				int cell_dim = (row & 1) + (col & 1);
				int prev_nof_crit_cells = nof_crit[cell_dim];
				crit_cell_list[cell_dim][prev_nof_crit_cells] = ext_acc[row][col];
				if (cell_dim < dim_I) { //note that dim_I-cells have no cobounds
					// exploring and logging cobounds of crit cells
					int cell_pos = (row & 1);
					/* a little trick valid for planar images and for 0,1 cells: only for odd rows, the cell position (of 1 cells is 1)
					For 0 cells, always row and col are even; thus cell_pos = 0 */
					for (int cobound_numb = 1; cobound_numb <= nof_cobound_neigh[cell_dim + 1 - M2C]; cobound_numb++) {
						//exploring the cobounds  of the recent crit cell
						int row_increm = increm_coboundary_neigh_row[cobound_numb - M2C][cell_pos][cell_dim + 1 - M2C];
						int col_increm = increm_coboundary_neigh_col[cobound_numb - M2C][cell_pos][cell_dim + 1 - M2C];
						int row_cobound = row + row_increm;
						int col_cobound = col + col_increm;

						if ((cell_dim + 1) == ext_obj_dim_belonging[row_cobound][col_cobound]) {
							// if the ext_obj_dim_belonging of the boundary were superior to (ext_obj_dim_belonging of the cell plus 1, this logging would not be useful
							//e.g if a 0 cell had in its cobound a region (1 cell with ext_obj_dim_belonging=2), the logging is not useful
							//e.g a 0 cell can have in its cobound frontiers (1 cell with ext_obj_dim_belonging=1), these are valid  loggings 
							
							int prev_nof_cob = crit_cell_nof_cob[cell_dim][prev_nof_crit_cells];
							crit_cell_cob_list[cell_dim][prev_nof_crit_cells][prev_nof_cob] = ext_acc[row_cobound][col_cobound];
							crit_cell_nof_cob[cell_dim][prev_nof_crit_cells]++;
							//  xxxxxx
						}
					}
				}  //end of if (cell_dim < dim_I) 
				nof_crit[cell_dim] ++;
			}
		}
	}
	//now that all crit cells have been logged, we detect those that represent a hole (and not a cross or a frontier between regions)
	/* In order to mark crit cells that represent a hole, the nof cobounds are added this value: 2 * (dim_I)
This is to codify the hole condition in this crit_cell_nof_cob instead of defining another mattrix sth like "crit_cell_hole_condition"
*/
	for (int d = 0; d < dim_I; d++) { // The loop order "col outer, row inner" matches that of matlab; better for debugging 
		for (int cc = 0; cc < nof_crit[d]; cc++) { //beginning end ending in 2 or -2 because of the 2 ghost rows
			if (crit_cell_nof_cob[d][cc] == 2 &&
				crit_cell_cob_list[d][cc][0] == crit_cell_cob_list[d][cc][1]) {
				// a simple way to detect hole in R2: only 2 cobounds and they are the same
				crit_cell_nof_cob[d][cc] += dim_I * 2;
			}
		}
	}
	// now we can obtain el graph HGIG: for each k cell of crit_cell_list[k][number], I get its cobounds from crit_cell_cob_list[k].  
	//Then, a binary search in the (k+1)-list of the crit_cell_list[k+1], which are ordered in increasing order.
	deleteCrit_Cells_Arrays();
};


/////////////////////////////////

void printing_machine_info(void) {

	int CPUInfo[4] = { -1 };
	unsigned   nExIds, i = 0;
	char CPUBrandString[0x40];
	// Get the information associated with each extended ID.
	__cpuid(CPUInfo, 0x80000000);
	nExIds = CPUInfo[0];
	for (i = 0x80000000; i <= nExIds; ++i)
	{
		__cpuid(CPUInfo, i);
		// Interpret CPU brand string
		if (i == 0x80000002)
			memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
		else if (i == 0x80000003)
			memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
		else if (i == 0x80000004)
			memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
	}
	//string includes manufacturer, model and clockspeed
	std::cout << "CPU Type: " << CPUBrandString << endl;

	SYSTEM_INFO sysInfo;
	GetSystemInfo(&sysInfo);
	std::cout << "Number of Cores: " << sysInfo.dwNumberOfProcessors << endl;

	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);
	std::cout << "Total System Memory: " << (statex.ullTotalPhys / 1024) / 1024 << "MB" << endl;
}

///////////////////////////
