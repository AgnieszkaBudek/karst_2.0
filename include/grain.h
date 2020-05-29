/**
* @file grain.h
* @brief This file contains class grains.
* @author Agnieszka Budek
* @date 25/09/2019
*/

#ifndef GRAIN_H
#define GRAIN_H 

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>

#include "pore.h"
#include "node.h"
#include "constants.h"



using namespace std;



class Node;
class Pore;
class Network;

/**
 * @file
 * @author  Agnieszka Budek
 * @date 25/09/2019
 *
 * @class Grain
 *
 * @section DESCRIPTION
 *
 * Class grain represents one grain - the piece of material surrounded by pores.
 * Grain can consists of the materials A and E. The first one can be dissolved by acid in dissolution reaction while the second
 * will precipitate in the second reaction.
 *
 */

/// represents one grains/// represents one pore/// represents one pore
class Grain{

	public:

		const int bm;    ///< number of solid species building a grain
		double *V;		 ///< table with volume of particular species
		double *V_old;   ///< table with volume of particular species in previous time step

		double tmp;		 ///< temporal information (name or type of grain)///

		int a;           ///< name of a grain
		int bN;          ///< numbers of nodes in vicinity (used for non-triangular network)
		int bP;          ///< numbers of pores in vicinity (used for non-triangular network)

		Node **n;		///< node list
		Pore **p;		///< pore list


	public:

		Grain (float tmp,  Node* nn0, Node* nn1, Node* nn2,int bm);                              //constructor for triangular network
		Grain (float tmp,  Node* nn0, Node* nn1, Node* nn2, Node* nn3,int bm);                   //constructor for square network
		Grain (float tmp=0,  int bbP=3, int bbN=3, Node** nn0=NULL, Pore** pp0=NULL, int bm);    //constructor for a general network
		Grain (float name, double V_a_tmp, double V_e_tmp, double V_f_tmp,int bb_N, int bb_P, int bm);          //constructor for a network read from file
		Grain (Grain &g);
		Grain (int bw);
		~Grain ()  {delete[] n; n=NULL; delete[] p; p=NULL; delete[] V; V=NULL; delete[] V_old; V_old=NULL;}

		Grain& operator = (Grain &g);

		void   calculate_initial_volume (Network *S);
		double calculate_maximal_volume (Network *S);
		void   change_nodes(Node *n_old, Node *n_new);
		void   change_pores(Pore *p_old,Pore * p_new);
		void   set_effective_d_and_l(Pore *p,Network *S);
		void   remove_Pore(Pore *p);
		void   add_Pore   (Pore *p);
		void   remove_Node(Node *n);
		void   add_Node   (Node *n);
		bool   to_be_merge();
		bool   do_contain_Node(Node *n_tmp);
		bool   do_contain_Pore(Pore *p_tmp);
		bool   is_pathological ();
};


ofstream_txt & operator << (ofstream_txt & stream, Grain &g);
ostream      & operator << (ostream      & stream, Grain &g);

Pore* findPore(Node*, Node*);


inline void swap_grains(Grain **g1, Grain **g2){

	Grain *g_tmp = *g1;
	*g1=*g2;
	*g2=g_tmp;
	double tmptmp = (*g2)->tmp;
	(*g2)->tmp = (*g1)->tmp;
	(*g1)->tmp = tmptmp;

}

#endif


