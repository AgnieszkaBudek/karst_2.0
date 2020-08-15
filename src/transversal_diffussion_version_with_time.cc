#include "network.h"
#include "printing.h"
#include "constants.h"
#include "algorithms_cc.h"

/**
* This function calculates the concentration field for species B in the entire network.
* First linear equation for concentration in each node are created and
* then there are solved using function Network::solve_matrix().
*This version takes into account transversal diffusion through the pores.
*
* @author Agnieszka Budek
* @date 12/06/2020
*/
// evolution with transversal diffusion
void Network::calculate_concentrations_b_diff_T(){

	//cerr<<"Calculating concentrations for species B taking into account diffusion and time..."<<endl;


	//calculating new cb for each node
	for(int i=0;i<NN;i++){

		Node *nn = n[i];

		if(nn->t==1){      nn->cb = Cb_0;    continue;}
		if(tot_steps==0){if(nn->xy.y < 0) nn->cb = Cb_init; continue;}

		double J_in   = 0;    //amount of species B flowing into the node due to convection
		double J_diff = 0;    //amount of species B flowing in/out of the node due to diffusion
		double Q_out  = 0;    //total flow through the node

		for(int s=0; s<nn->b; s++){		//eqs for normal and outlet nodes

			Pore *pp = findPore(nn,nn->n[s]);
			if (pp->d == 0 || pp->l<=l_min)  continue;    //no reaction in tiny grain or in pore with no flow
			double qq;
			if(nn==pp->n[0]) qq = -pp->q;
			else             qq =  pp->q;
			double l_tmp = point_distance(nn->xy,nn->n[s]->xy)/2.;

			if(qq>0) { J_in += (qq * pp->cb_old); Q_out+=qq;}
			J_diff += (nn->p[s]->cb_old - nn->cb_old) * M_PI*pow(pp->d,2)/4.*D1/l_tmp;
		}

		nn->cb += dt * (J_in + J_diff - Q_out*nn->cb_old) / nn->V;
		if(nn->cb<0)  nn->cb = 0;

	}


	//cerr<<"Calculating concentrations for species B in pores taking into account diffusion and time..."<<endl;


	//calculating new cb for each pore
	for(int i=0;i<NP;i++){

		Pore *pp = p[i];
		if (pp->d == 0 || pp->l<=l_min)     continue;    //no reaction in tiny grain or in pore with no flow
		if(tot_steps==0){ if(pp->n[0]->xy.y < 0.) pp->cb = Cb_init; continue;}

		double J_in   = 0;    //amount of species B flowing into the pore in one time step due to convection
		double J_diff = 0;    //amount of species B flowing in/out of the pore due to diffusion

		Node* nn = NULL;
		if(pp->n[0]->u > pp->n[1]->u) nn = pp->n[0];
		else             			  nn = pp->n[1];

		double l_tmp = point_distance(pp->n[0]->xy,pp->n[1]->xy)/2.;

		J_in   =  nn->cb_old * fabs(pp->q);
		J_diff = (pp->n[0]->cb_old + pp->n[1]->cb_old - 2*pp->cb_old) * M_PI*pow(pp->d,2)/4.*D1/l_tmp;

		pp->cb +=  dt*( J_in + J_diff - fabs(pp->q)*pp->cb_old) / pp->volume();
		if(pp->cb<0)  pp->cb = 0;
	}

	//additional printing for debugging
	print_network_for_debugging ("After calculating concentration B field ","acid concentration", "concentration B");

}

void Network::calculate_concentrations_c_diff_T(){

	//cerr<<"Calculating concentrations for species C taking into account diffusion and time..."<<endl;


	//calculating new cc for each node
	for(int i=0;i<NN;i++){

		Node *nn = n[i];

		if(nn->t==-1){     nn->cc = Cc_0;     continue;}
		if(tot_steps==0){ if(nn->xy.y >= 0.) nn->cc = Cc_init;  continue;}

		double J_in   = 0;    //amount of species C flowing into the node due to convection
		double J_diff = 0;    //amount of species C flowing in/out of the node due to diffusion
		double Q_out  = 0;    //total flow through the node

		for(int s=0; s<nn->b; s++){		//eqs for normal and outlet nodes

			Pore *pp = findPore(nn,nn->n[s]);
			if (pp->d == 0 || pp->l<=l_min)  continue;    //no reaction in tiny grain or in pore with no flow
			double qq;
			if(nn==pp->n[0]) qq = -pp->q;
			else             qq =  pp->q;
			double l_tmp = point_distance(nn->xy,nn->n[s]->xy)/2.;

			if(qq>0) { J_in += (qq * pp->cc_old); Q_out+=qq;}
			J_diff += (nn->p[s]->cc_old - nn->cc_old) * M_PI*pow(pp->d,2)/4.*D2/l_tmp;
		}

		nn->cc += dt * (J_in + J_diff - Q_out*nn->cc_old) / nn->V;
		if(nn->cc<0)  nn->cc = 0;

	}


	//cerr<<"Calculating concentrations for species B in pores taking into account diffusion and time..."<<endl;


	//calculating new cc for each node
	for(int i=0;i<NP;i++){

		Pore *pp = p[i];
		if (pp->d == 0 || pp->l<=l_min)     continue;    //no reaction in tiny grain or in pore with no flow
		if(tot_steps==0){if(pp->n[0]->xy.y >= 0) pp->cc = Cc_init; continue;}

		double J_in   = 0;    //amount of species C flowing into the pore in one time step due to convection
		double J_diff = 0;    //amount of species C flowing in/out of the pore due to diffusion

		Node* nn = NULL;
		if(pp->n[0]->u > pp->n[1]->u) nn = pp->n[0];
		else             			  nn = pp->n[1];
		double l_tmp = point_distance(pp->n[0]->xy,pp->n[1]->xy)/2.;

		J_in   =  nn->cc_old * fabs(pp->q);
		J_diff = (pp->n[0]->cc_old + pp->n[1]->cc_old - 2*pp->cc_old) * M_PI*pow(pp->d,2)/4.*D2/l_tmp;


		pp->cc += dt*(J_in + J_diff - fabs(pp->q)*pp->cc_old) / pp->volume();
		if(pp->cc<0)  pp->cc = 0;
	}

	//additional printing for debugging
	print_network_for_debugging ("After calculating concentration B field ","acid concentration", "concentration B");


}


void Network::calculate_concentrations_f_diff_T(){

	//cerr<<"Calculating concentrations for species C taking into account diffusion and time..."<<endl;


	//calculating new cf for each node
	for(int i=0;i<NN;i++){

		Node *nn = n[i];
		if(tot_steps==0)  {nn->cf = Cf_init; continue;}

		double J_in   = 0;    //amount of species C flowing into the node due to convection
		double J_diff = 0;    //amount of species C flowing in/out of the node due to diffusion
		double Q_out  = 0;    //total flow through the node

		for(int s=0; s<nn->b; s++){		//eqs for normal and outlet nodes

			Pore *pp = findPore(nn,nn->n[s]);
			if (pp->d == 0 || pp->l<=l_min)  continue;    //no reaction in tiny grain or in pore with no flow
			double qq;
			if(nn==pp->n[0]) qq = -pp->q;
			else             qq =  pp->q;
			double l_tmp = point_distance(nn->xy,nn->n[s]->xy)/2.;

			if(qq>0) { J_in += (qq * pp->cf_old); Q_out+=qq;}
			J_diff += (nn->p[s]->cf_old - nn->cf_old) * M_PI*pow(pp->d,2)/4.*D3/l_tmp;
		}

		nn->cf += dt * (J_in + J_diff - Q_out*nn->cf_old) / nn->V;
		if(nn->cf<0)  nn->cf = 0;

	}



	//cerr<<"Calculating concentrations for species B in pores taking into account diffusion and time..."<<endl;


	//calculating new cf for each node
	for(int i=0;i<NP;i++){

		Pore *pp = p[i];
		if (pp->d == 0 || pp->l<=l_min)     continue;    //no reaction in tiny grain or in pore with no flow
		if(tot_steps==0){ pp->cf = Cf_init; continue;}

		double J_in   = 0;    //amount of species C flowing into the pore in one time step due to convection
		double J_diff = 0;    //amount of species C flowing in/out of the pore due to diffusion

		Node* nn = NULL;
		if(pp->n[0]->u > pp->n[1]->u) nn = pp->n[0];
		else             			  nn = pp->n[1];
		double l_tmp = point_distance(pp->n[0]->xy,pp->n[1]->xy)/2.;

		J_in   =  nn->cf_old * fabs(pp->q);
		J_diff = (pp->n[0]->cf_old + pp->n[1]->cf_old - 2*pp->cf_old) * M_PI*pow(pp->d,2)/4.*D3/l_tmp;


		pp->cf += dt*(J_in + J_diff - fabs(pp->q)*pp->cf_old) / pp->volume();
		if(pp->cf<0)  pp->cf = 0;
	}

	//additional printing for debugging
	print_network_for_debugging ("After calculating concentration B field ","acid concentration", "concentration B");


}



void Network::precipitate(){


	cerr<<"Precipitating in zeolite..."<<endl;

	//for updating grains volume;
	if(if_track_grains) for (int i=0;i<NG;i++) {g[i]->tmp=0; g[i]->tmp2=0; }

	//avoid negative concentrations
	for (int i=0; i<NN; i++) {
		if(n[i]->cb<0) n[i]->cb=0;
		if(n[i]->cc<0) n[i]->cc=0;
		if(n[i]->cf<0) n[i]->cf=0;
	}
	//avoid negative concentrations
	for (int i=0; i<NP; i++) {
		if(p[i]->cb<0) p[i]->cb=0;
		if(p[i]->cc<0) p[i]->cc=0;
		if(p[i]->cf<0) p[i]->cf=0;
	}


	//checking dt
	for (int i=0;i<NN;i++){
		Node* nn = n[i];
		if(nn->cb>0) set_adaptive_dt_for_cT(fabs((nn->cb - nn->cb_old)/nn->cb));
		if(nn->cc>0) set_adaptive_dt_for_cT(fabs((nn->cc - nn->cc_old)/nn->cc));
		if(nn->cf>0) set_adaptive_dt_for_cT(fabs((nn->cf - nn->cf_old)/nn->cf));
	}

	for (int i=0;i<NP;i++){
		Pore* pp = p[i];
		if(pp->cb>0) set_adaptive_dt_for_cT(fabs((pp->cb - pp->cb_old)/pp->cb));
		if(pp->cc>0) set_adaptive_dt_for_cT(fabs((pp->cc - pp->cc_old)/pp->cc));
		if(pp->cf>0) set_adaptive_dt_for_cT(fabs((pp->cf - pp->cf_old)/pp->cf));
	}


	//final reaction


	for(int i=0;i<NP;i++){ //for each pore...


		Pore* p0 = p[i];
		if (p0->d == 0 || p0->l<=l_min)  continue;    //no reaction in tiny grain or in pore with no flow

		//update cb, cc and cf due to reactions:
		double R_1_tmp =  p0->R_1(this);
		double C1_tmp  =  R_1_tmp / p0->volume();
		p0->cb -=  C1_tmp;
		p0->cc -=  C1_tmp;
		p0->cf +=  C1_tmp;
		double R_2_tmp =  p0->R_2(this);
		p0->cf -=  R_2_tmp / p0->volume();


		//avoiding negative concentration
		if(p0->cb<0) p0->cb=0;
		if(p0->cc<0) p0->cc=0;
		if(p0->cf<0) p0->cf=0;


		double d_old    = p0->d;
		double dd_plus  = 0;
		double dd_minus = R_2_tmp/ (M_PI*p0->l*p0->d/2.);;

		//update geometry
		//p0->d += (dd_plus - dd_minus);
		if(p0->d<0) {Ve_tot+= M_PI*pow(p0->d,2) * p0->l; p0->d = 0;}



		//updating Va and Ve volumes
		int bG_tmp_A=0;
		double d_V_A = (M_PI*(d_old)*(dd_plus *d0)/2*p0->l);
		double d_V_E = (M_PI*(d_old)*(dd_minus*d0)/2*p0->l);
		for(int s=0; s<p0->bG;s++) if(p0->g[s]->Va >0) bG_tmp_A++;
		for(int s=0; s<p0->bG;s++) {
			if(p0->g[s]->Va >0) p0->g[s]->tmp -=d_V_A/bG_tmp_A;
			if (true)           p0->g[s]->tmp2+=d_V_E/p0->bG;
		}

		if(if_adaptive_dt)      set_adaptive_dt((dd_plus - dd_minus)*d0/p0->d, d_V_A + d_V_E);

		//temporal line: I don;t want temporally to change the geometry
		if(p0->d<d_min) p0->d = d_min*0.999;
	}

	//updating Va and Vc (must be done after main dissolution for c_out to be calculated correctly)
	if(if_track_grains){
		for (int i=0;i<NG;i++) {g[i]->Va+=g[i]->tmp;  if(g[i]->Va<0) {Va_tot-=g[i]->Va; g[i]->Va = 0;}}
		for (int i=0;i<NG;i++)  g[i]->Ve+=g[i]->tmp2;
	}

	//update concentration_old
	for (int i=0; i<NN; i++) {n[i]->cb_old = n[i]->cb; n[i]->cc_old = n[i]->cc; n[i]->cf_old = n[i]->cf;}
	for (int i=0; i<NP; i++) {p[i]->cb_old = p[i]->cb; p[i]->cc_old = p[i]->cc; p[i]->cf_old = p[i]->cf;}


	//updating pore lengths
	if(if_dynamical_length) for(int i=0; i<NP; i++) p[i]->calculate_actual_length(this);
	//updating node volume
	if(if_track_grains)     for(int i=0; i<NN; i++) n[i]->calculate_volume(this);

	//additional printing for debugging
	print_network_for_debugging ("After dissolution and precipitation","pressure", "diameter","volume A");

}


void Network::set_adaptive_dt_for_cT(double dc){

	if (set_new_dt == 1) return;
	if (d_d_max > 0 && dc > d_d_max)   {set_new_dt = 1; return;}
	if (set_new_dt == -1 && d_d_min > 0 && dc > d_d_min) set_new_dt = 0;
}
