#include "network.h"
#include "printing.h"

//functions related to evolution


/**
* This function simulates the evolution of the system.
*
* @param T maximal no of time steps to be simulated
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::evolution(long int T){

	cerr<<"\n\n\nEvolution:"<<endl;

	if (T==0) T = T_max;
	double 	t  = 0; //total time
	int 	s  = 0; //nr of steps

	if(if_save_txt)     print_net_txt();
	save_all_data       (true);
	analyze_diss_pattern(true);

	while(s<T){
		cerr<<endl;
		cerr<<s<<". step of evolution (t = "<<t<<")"<<endl;

		if(if_leapfrog)  do_one_leapfrog_step();	//do frog leap version (not implemented yet)
		else			 do_one_euler_step();		//do normal Euler version

		s++; t+=dt;
		tot_steps++; tot_time+=dt;

		save_all_data();
		analyze_diss_pattern();

		if(if_system_dissolved){
			save_all_data       (true);
			analyze_diss_pattern(true);
			break;
		}

	}

}

/**
* This function do one euler time step.
* The order of calculations is as follows:
*
* 1. Calculate pressure field
*
* 2. Calculate flow field
*
* 3. Calculate the concentration profiles.
*
* 4. Change pore sizes according to the dissolution (and precipitation).
*
* 5. Check the flow and mass balance.
*
* 6. Do merging (optionally).
*
* This is the simplest method for integrating of differential equations describing the system dynamics.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::do_one_euler_step(){

//Calculate pressure and flow field.
	calculate_pressure_and_flow_filed();

//Calculate concentration fields for all species
	calculate_concentration_field();

//Calculate new values of diameters and lengths
	update_geometry();

//Check mass balance
	check_mass_balance();

//Topology changes
	do_merging();

//
	if(if_adaptive_dt)      adapt_dt();           //adapt dt
	if(if_full_dissolution) check_if_dissolved(); //check if the system is dissolved

}


/**
* This function will do one leap-frog step.
* Will be implemented later.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::do_one_leapfrog_step(){

//trzeba będzie też znak pradu dyfuzji lepiej liczyć w leaf frogu!!!
//bo tam biore znak z wczesniejszego kroku (mozna by to iteracyjnie, ale moim zdaniem bez sensu strata czasu)

}





void Network::calculate_pressure_and_flow_filed(){

//calculating pressures
	if(if_smarter_calculation_of_pressure)	{
		calculate_pressures_and_flows_smarter(d_max_for_u);}

	else {
		calculate_pressures();
		calculate_flows();}
}



void Network::calculate_concentration_field(){

	if(Pe==-1){  //without transversal diffusion

		if(if_precipitation){
			calculate_concentrations_b();
			calculate_concentrations_c();
			}
		else{
			if(if_streamtube_mixing) 	calculate_concentrations_streamtube_mixing();
			else						calculate_concentrations_b();
		}
	}

	else {   //version with transversal diffusion
		calculate_concentrations_b_diff();
		if(if_precipitation) calculate_concentrations_c_diff();
	}

}


void Network::update_geometry(){

	//updating diameters and lengths

	if(if_precipitation)	dissolve_and_precipitate();
	else                    dissolve();

}

void Network::check_mass_balance(){


	check_flow_balance();
	if(if_track_grains) check_acid_balance();
	if(if_precipitation&&if_track_grains) check_precipitating_balance();
	//check_GMash_connections();  //additional checking for G_Mash network

}
