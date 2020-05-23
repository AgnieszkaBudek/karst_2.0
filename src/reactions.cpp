/*
 * reactions.cpp
 *
 *  Created on: May 17, 2020
 *      Author: agnieszka
 */

#include "../include/reactions.h"

Reactions::Reactions(string model_name) {

	if      (model_name == "pure_dissolution") 					Set_reactions_to_pure_dissolution();
	else if (model_name == "dissolution_and_precipitation_A") 	Set_reactions_to_dissolution_and_precipitation_A();
	else if (model_name == "dissolution_and_precipitation_F1") 	Set_reactions_to_dissolution_and_precipitation_F1();
	else    {cerr<<"ERROR: Wrong reaction model!!!"<<endl; exit(666);}

}

Reactions::~Reactions() {

	delete [] k;  k = NULL;
	delete [] D;  D = NULL;
	delete [] DD; DD= NULL;
	delete [] C;  C = NULL;
	delete [] V;  V = NULL;
	delete [] gamma; gamma= NULL;

}

void Reactions::Set_reactions_to_pure_dissolution(){

	bw = 1;    ///< number of soluble species
	bm = 2;    ///< one material reacting and one immune to reactions
	br = 1;    ///< one reaction

	w = new Soluble* [bw];    		///< table with soluble materials
	m = new Solid*   [bm];      		///< table with solid materials

	C = new double [bw];    		///< table with soluble materials concentration
	V = new double [bm];      		///< table with total amount of volume for all species

	k = new double [br];      		///< table with reactions rate (not used now)
	D = new double [bw];      		///< table with diffusion coefficient (not used now)
	DD= new double [bw];      		///< table with transversal diffusion coefficient (not used now)
	gamma = new double[bm];         ///< table with capacity number for all solids

	C[0]     = 1;                   ///< by definition in this system inlet concentration of acid is 1
	gamma[0] = 1;					///< by definition in this system gamma is 1

	w[0] = new Soluble_B_Agnieszkas_model();
	m[0] = new Solid_A_Agnieszkas_model();
	m[1] = new Solid_non_reacting();

}

void Reactions::Set_reactions_to_dissolution_and_precipitation_A(){


	bw = 2;    ///< number of soluble species
	bm = 3;    ///< one material reacting and one immune to reactions
	br = 2;    ///< one reaction

	w = new Soluble* [bw];    		///< table with soluble materials
	m = new Solid*   [bm];      		///< table with solid materials

	C = new double [bw];    		///< table with soluble materials concentration
	V = new double [bm];      		///< table with total amount of volume for all species

	k = new double [br];      		///< table with reactions rate (not used now)
	D = new double [bw];      		///< table with diffusion coefficient (not used now)
	DD= new double [bw];      		///< table with transversal diffusion coefficient (not used now)
	gamma = new double[bm];         ///< table with capacity number for all solids

	C[0]     = 1;                   ///< can be change later in network creator (while reading setup file)
	C[1]     = 0;                   ///< can be change later in network creator (while reading setup file)
	gamma[0] = 1;					///< can be change later in network creator (while reading setup file)
	gamma[1] = 1;					///< can be change later in network creator (while reading setup file)

	w[0] = new Soluble_B_Agnieszkas_model();
	w[1] = new Soluble_C_Agnieszkas_model();
	m[0] = new Solid_A_Agnieszkas_model();
	m[1] = new Solid_E_Agnieszkas_model();
	m[2] = new Solid_non_reacting();
}

void Reactions::Set_reactions_to_dissolution_and_precipitation_F1(){

	//to be implemented

	cerr<<"ERROR: Florians model has not been implemented yet."<<endl;
	exit(666);

}
