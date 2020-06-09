/*
 * reactions.cpp
 *
 *  Created on: May 17, 2020
 *      Author: agnieszka
 */

#include "reactions.h"

Single_Reaction::Single_Reaction(){

	sw = NULL;
	sm = NULL;
	pw = NULL;
	pm = NULL;
	k  = -1;

}

Single_Reaction::Single_Reaction(Soluble * ssw, Solid *ssm, Soluble *ppw, Solid *ppm){

	sw = ssw;
	sm = ssm;
	pw = ppw;
	pm = ppm;
	k  = -1;

}


Reactions::Reactions(string model_name) {

	if      (model_name == "pure_dissolution") 					Set_reactions_to_pure_dissolution();
	else if (model_name == "dissolution_and_precipitation_A") 	Set_reactions_to_dissolution_and_precipitation_A();
	else if (model_name == "dissolution_and_precipitation_F1") 	Set_reactions_to_dissolution_and_precipitation_F1();
	else    {cerr<<"ERROR: Wrong reaction model!!!"<<endl; exit(666);}

}

Reactions::~Reactions() {

	delete [] w;
	delete [] m;
	delete [] r;

}

void Reactions::Set_reactions_to_pure_dissolution(){

	bw = 1;    ///< number of soluble species
	bm = 2;    ///< one material reacting and one immune to reactions
	br = 1;    ///< one reaction

	w = new Soluble* [bw];    		///< table with soluble materials
	m = new Solid*   [bm];          ///< table with solid materials
	r = new Singla_Reaction* [br];  ///< list of reactions

	w[0] = new Soluble_B_Agnieszkas_model();
	m[0] = new Solid_A_Agnieszkas_model();
	m[1] = new Solid_non_reacting();

	r[0] = new Single_Reaction(w[0],m[0],NULL,NULL);

////to be deleted
////giving info about reactions to Soluble and Solids
//	w[0]->br   = 1;
//	w[0]->r    = new Single_Reaction* [1];
//	w[0]->r[0] = r[0];
//
//	m[0]->br   = 1;
//	m[0]->r    = new Single_Reaction* [1];
//	m[0]->r[0] = r[0];
//
//	m[1]->br   = 1;
//	m[1]->r    = new Single_Reaction* [1];
//	m[1]->r[0] = r[0];
}

void Reactions::Set_reactions_to_dissolution_and_precipitation_A(){


	bw = 2;    ///< number of soluble species
	bm = 3;    ///< one material reacting and one immune to reactions
	br = 2;    ///< one reaction

	w = new Soluble* [bw];    		///< table with soluble materials
	m = new Solid*   [bm];          ///< table with solid materials
	r = new Singla_Reaction* [br];  ///< list of reactions


	w[0] = new Soluble_B_Agnieszkas_model();
	w[1] = new Soluble_C_Agnieszkas_model();
	m[0] = new Solid_A_Agnieszkas_model();
	m[0] = new Solid_E_Agnieszkas_model();
	m[1] = new Solid_non_reacting();

	r[0] = new Single_Reaction(w[0],m[0],w[1],NULL);
	r[1] = new Single_Reaction(w[1],NULL,NULL,m[1]);


//to be deleted
////giving info about reactions to Soluble and Solids
//	w[0]->br   = 1;
//	w[0]->r    = new Single_Reaction* [1];
//	w[0]->r[0] = r[0];
//
//	m[0]->br   = 1;
//	m[0]->r    = new Single_Reaction* [1];
//	m[0]->r[0] = r[0];
//
//	w[1]->br   = 2;
//	w[1]->r    = new Single_Reaction* [2];
//	w[1]->r[0] = r[0];  w[1]->r[1] = r[1];
//
//
//	m[1]->br   = 1;
//	m[1]->r    = new Single_Reaction* [1];
//	m[1]->r[0] = r[1];



}

void Reactions::Set_reactions_to_dissolution_and_precipitation_F1(){

	//to be implemented

	cerr<<"ERROR: Florians model has not been implemented yet."<<endl;
	exit(666);

}


double Reactions::default_diameter_change(Pore* p, Network *S){

	double dd = 0;
	for(int im=0; im<bm; im++) dd+=m[im]->dd(p,S);

	return dd;
}



void Reactions::final_geometry_change(Pore* p0, Network *S){

//checking if reaction occures

	bool if_reaction = false;
    for(int iw=0;iw<bw;iw++) if(w[iw]->is_reacting(p0,S)) if_reaction = true;
    if(!if_reaction) return;

    double d_old      = p0->d;
	double dd         = 0;  //total change in diameter due to all reactions
	double dd_plus    = 0;  //all dissolution effect
	double dd_minus   = 0;  //all precipitation effect
	double dd_minus_f = 0;  //final amount of precipitation effect that has enough space to precipitate

	for(int im=0; im<bm; im++){
		double d_tmp = m[im]->dd(p0,S);
		if(dd>0) dd_plus  += d_tmp;
		else     dd_minus += d_tmp;
		dd += d_tmp;
	}

	//Checking if there is enough space for full dissolution
	if(p0->d + dd < S->d_min){		//there is not enough space for all precipitating material
		p0->d=S->d_min;
		dd_minus_f = p0->d + dd_plus - S->d_min;
	}
	else{											//there is enough space for all precipitating material
		p0->d += dd;
		dd_minus_f = dd_minus;
	}

	double d_V_tot = 0;
	//updating grains volume
	if (S->if_track_grains) for(int im=0;im<bm;im++){
		double dd_tmp = m[im]->dd(p0,S);
		double d_V_tmp = -(M_PI*d_old*dd_tmp*(S->d0/2)*p0->l);
		if(d_V_tmp>0) d_V_tmp = d_V_tmp*dd_minus_f/dd_minus;
		d_V_tot += d_V_tmp;
		int bG_tmp=0;
		for(int s=0; s<p0->bG;s++) if(p0->g[s]->V[im]>0 || dd_tmp<0) bG_tmp++;
		for(int s=0; s<p0->bG;s++) if(p0->g[s]->V[im]>0 || dd_tmp<0) p0->g[s]->V_old[im]+=d_V_tmp/bG_tmp;
		}

	//check if to adapt time step
	if(S->if_adaptive_dt)      S->set_adaptive_dt(dd/p0->d, d_V_tot);
}


