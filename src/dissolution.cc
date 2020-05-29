#include "network.h"
#include "printing.h"
#include "algorithms_cc.h"

//functions related to evolution

/**
* This function calculates the pressure field for entilre network.
* First linear equation for pressure in each node are created and
* then there are solved using function Network::solve_matrix().
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_pressures(){

	cerr<<"Calculating pressures..."<<endl;
	
	for(int i=0;i<NN;i++) n[i]->tmp=i;

	int R_no = 0;
	for(int i=0;i<NN;i++)
		if(n[i]->t==0) R_no+=n[i]->b+1;
		else           R_no++;

	int R_m  = NN;				    //rank of the matrix to be solved

	int* ww_r = new int [R_no];		//raw indexes
	int* ww_c = new int [R_no];		//column indexes
	double* B = new double [R_no];	//non zero elements
	double* y = new double [NN];	//rhs

	//filling rhs (pressures in nodes)
	for(int i=0;i<NN;i++) if(n[i]->t==1) y[i]=P_in; else y[i]=0;
	
	//filing matrix of eqs for each node
	int r_no=0; double S=0;
	for(int i=0;i<NN;i++){
		S=0;
		if(abs(n[i]->t)==1) S=-1;			     //if node is an I/O node
		else for(int s=0; s<n[i]->b;s++){	 //eqs for normal nodes
			ww_r[r_no] 	= i;
			ww_c[r_no] 	= n[i]->n[s]->tmp;
			B[r_no]		= n[i]->p[s]->perm(mu_0);
			S+=B[r_no];
			r_no++;}
		ww_r[r_no] 		= i;
		ww_c[r_no] 		= i;
		B[r_no]			=-S;
		r_no++;	
	}
	if(r_no!=R_no) {cerr<<"Problem with filling linear equations for pressures!"<<endl; exit(666);}

	int M_out = solve_matrix(R_m, R_no, ww_r, ww_c, B, y);
	if(M_out!=0) cerr<<"Problem with solving linear equations; M_out = "<<M_out<<endl;
	
	for(int i=0;i<NN;i++) n[i]->u = y[i];  //filling the solution
	
	//checking if the matrix has been solved correctly
	double S_tmp=0;
	for(int i=0;i<NN;i++) if(y[i]!=P_in) S_tmp+= y[i];
	if(S_tmp==0) {
		cerr<<"Pressure has not been calculated properly."<<endl;
		ofstream_txt B_err;
		B_err     .open("B_err.out",	      ios_base::out | ios_base::trunc );
		for(int i=0;i<R_no;i++) B_err<<ww_r[i]<<"\t"<<ww_c[i]<<"\t"<<B[i]<<endl;
		B_err.close();
		exit(123);
	}

	//additional printing for debugging
	print_network_for_debugging ("After calculating pressure field","pressure", "diameter");
	
	delete[] ww_r;
	delete[] ww_c;
	delete[] B;
	delete[] y;
}

/**
* This function calculates the pressure and flow field in two step method.
* First the system is divided into two sets: small and large pore  using definition of dissolution pattern (with threshold equal d_max).
* Then the pressure field for small pores are calculated assuming no pressure drop in large pores.
* Next the flow field is calculated for this region.
* In the and the pressure and flow field is calculated for the rest of the system.
* @param d_max threshold for definition of dissolution pattern
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_pressures_and_flows_smarter(double d_max){


	check_diss_pattern              (d_max);	//select nodes connected to the inlet by thick pores

	//additional printing for debugging
	print_network_for_debugging ("After calculating pattern","x", "flow");

	calculate_pressures_for_small_d (d_max);
	calculate_flows ();
	calculate_pressures_for_large_d (d_max);
	calculate_flows_for_large_d     (d_max);

}

/**
* This function calculates pressure field for the small pores region.
* Is is used only in Network::calculate_pressures_and_flows_smarter() function.
*
* @param d_max threshold for definition of dissolution pattern
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_pressures_for_small_d(double d_max){

	cerr<<"Calculating pressures smarter for small pores..."<<endl;

	for(int i=0;i<NN;i++) n[i]->tmp=i;


	int R_no = 0;
	for(int i=0;i<NN;i++)
		if(n[i]->t==0 && n[i]->x<2) R_no+=n[i]->b+1;
		else           R_no++;

	int R_m  = NN;						//rank of the matrix to be solved

	int* ww_r = new int [R_no];		//raw indexes
	int* ww_c = new int [R_no];		//column indexes
	double* B = new double [R_no];	//non zero elements
	double* y = new double [NN];	//rhs

	//filling rhs (pressures in nodes)
	for(int i=0;i<NN;i++) if(n[i]->t==1 || n[i]->x==2) y[i]=P_in; else y[i]=0;

	//filing matrix of eqs for each node
	int r_no=0; double S=0;
	for(int i=0;i<NN;i++){
		S=0;
		if(abs(n[i]->t)==1 || n[i]->x==2) S=-1;			     //if node is an I/O node
		else for(int s=0; s<n[i]->b;s++){	 //eqs for normal nodes
			ww_r[r_no] 	= i;
			ww_c[r_no] 	= n[i]->n[s]->tmp;
			B[r_no]		= n[i]->p[s]->perm(mu_0);
			S+=B[r_no];
			r_no++;}
		ww_r[r_no] 		= i;
		ww_c[r_no] 		= i;
		B[r_no]			=-S;
		r_no++;
	}
	if(r_no!=R_no) {cerr<<"Problem with filling linear equations for pressures!"<<endl; exit(666);}

	int M_out = solve_matrix(R_m, R_no, ww_r, ww_c, B, y);
	if(M_out!=0) cerr<<"Problem with solving linear equations; M_out = "<<M_out<<endl;

	for(int i=0;i<NN;i++) n[i]->u = y[i];  //filling the solution

	//additional printing for debugging
	print_network_for_debugging ("After calculating pressure smarter for small d ","pressure", "flow");

	delete[] ww_r;
	delete[] ww_c;
	delete[] B;
	delete[] y;

}


/**
* This function calculates pressure field for the large pores region.
* Is is used only in Network::calculate_pressures_and_flows_smarter() function.
*
* @param d_max threshold for definition of dissolution pattern
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_pressures_for_large_d(double d_max){

	cerr<<"Calculating pressures smarter for large pores..."<<endl;

	for(int i=0;i<NN;i++) n[i]->tmp=i;

	int R_no = 0;       //nr of nonzero elements of the matrix
	int R_m  = 0;		//nr of nodes to solve
	for(int i=0;i<NN;i++)
		if(n[i]->x==2) {
			n[i]->tmp = R_m++;
			if(n[i]->t==0) for(int s=0; s<n[i]->b;s++) if(n[i]->n[s]->x==2) R_no++;
			R_no++;
		}

	if (R_no == R_m) {cerr<<"Pressure for large d is not neccesery." <<endl; return;}


	int* ww_r = new int [R_no];		//raw indexes
	int* ww_c = new int [R_no];		//column indexes
	double* B = new double [R_no];	//non zero elements
	double* y = new double [R_m];	//rhs

	//filling rhs (pressures in nodes)
	for(int i=0;i<NN;i++){
		if(n[i]->x==2 && n[i]->t!=0)  y[i]=n[i]->u;
	    if(n[i]->x==2 && n[i]->t==0)  y[i]=0;}


	//filing matrix of eqs for each node
	int r_no=0; double S=0;
	for(int i=0;i<NN;i++){
		if(n[i]->x!=2) continue;

		S=0;
		if(abs(n[i]->t)==1) S=-1;			 //if node is an I/O node
		else for(int s=0; s<n[i]->b;s++){	 //eqs for normal nodes
			if(n[i]->n[s]->x==2){
				ww_r[r_no] 	= i;
				ww_c[r_no] 	= n[i]->n[s]->tmp;
				B[r_no]		= n[i]->p[s]->perm(mu_0);
				S+=B[r_no];
				r_no++;}
			else y[i] +=fabs(n[i]->p[s]->q);
		}
		ww_r[r_no] 		= i;
		ww_c[r_no] 		= i;
		B[r_no]			=-S;
		r_no++;
	}
	if(r_no!=R_no) {cerr<<"Problem with filling linear equations for pressures!"<<endl; exit(666);}

	int M_out = solve_matrix(R_m, R_no, ww_r, ww_c, B, y);
	if(M_out!=0) cerr<<"Problem with solving linear equations; M_out = "<<M_out<<endl;

	for(int i=0;i<NN;i++) if(n[i]->x==2) n[i]->u = y[(int)n[i]->tmp];  //filling the solution

	//additional printing for debugging
	print_network_for_debugging ("After calculating pressure smarter for large d ","pressure", "flow");

	delete[] ww_r;
	delete[] ww_c;
	delete[] B;
	delete[] y;
}


/**
* This function calculates flow field for the entilre network.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_flows(){ //version without pore merging

	cerr<<"Calculating flows..."<<endl;
	
	for(int i=0;i<NP;i++) p[i]->q = (p[i]->n[0]->u - p[i]->n[1]->u)*p[i]->perm(mu_0);
	for(int i=0;i<NP;i++) if(not isfinite(p[i]->q)) {
		cerr<<"Problem with flow through the system. Check if the system is not clogged."<<endl;
		s_save_data = 1;
		save_all_data();
		exit(3);}
	

	// additionally we can neglect extremely small q as they can cause numerical problems (equivalent to neglect to small d in precipitation)
	// if (q_min>0) for(int i=0;i<NP;i++) if (fabs(p[i]->q) < q_min) p[i]->q = 0;

	if(Q_tot!=0) recalculate_flows_to_keep_Q_tot("outlet");

	double Q_tot_tmp=0;
	for(int i=0;i<N_wi;i++) for(int s=0;s<wi[i]->b;s++) Q_tot_tmp+= fabs(wi[i]->p[s]->q);
	cerr<<"Q_tot = "<<Q_tot_tmp<< "    Pressure drop = " << wi[0]->u<<endl;
	//additional printing for debugging
	print_network_for_debugging ("After calculating flows ","pressure", "flow");
}

/**
* This function calculates flow field for the large pores region.
* Is is used only in Network::calculate_pressures_and_flows_smarter() function.
*
* @param d_max threshold for definition of dissolution pattern
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_flows_for_large_d(double d_max){
	cerr<<"Calculating flows..."<<endl;

	for(int i=0;i<NP;i++) if(p[i]->n[0]->x==2 && p[i]->n[1]->x==2) p[i]->q = (p[i]->n[0]->u - p[i]->n[1]->u)*p[i]->perm(mu_0);
	for(int i=0;i<NP;i++) if(!isfinite(p[i]->q)) {
			cerr<<"Problem with flow through the system. Check if the system is not clogged."<<endl;
			s_save_data = 1;
			save_all_data();
			exit(3);}

	// additionally we can neglect extremely small q as they can cause numerical problems (equivalent to neglect to small d in precipitation)
	// if (q_min>0) for(int i=0;i<NP;i++) if (fabs(p[i]->q) < q_min) p[i]->q = 0;


	//additional printing for debugging
	print_network_for_debugging ("After calculating flows for large d","pressure", "flow");
}

/**
* This function rescale the flow field to obtain the demanded total flow through the system.
* Is is used only in Network::calculate_pressures_and_flows_smarter() function.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::recalculate_flows_to_keep_Q_tot(string type){

	double Q_tmp = 0;


	if(type == "inlet"){  //to be used for large d only !!!
		//cerr<<"Recalculating flows for large d only!!!"<<endl;
		for(int i=0;i<N_wi;i++){
			Node* n_tmp = wi[i];
			for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) Q_tmp+=fabs(n_tmp->p[j]->q);
			}
		double factor = Q_tot/Q_tmp;
		for(int i=0;i<NP;i++) if(p[i]->n[0]->x==2 && p[i]->n[1]->x==2)  p[i]->q *=factor;   //rescaling all flows
		for(int i=0;i<NN;i++) if(n[i]->x==2)                            n[i]->u *=factor;   //rescaling all pressures
	}
	else if(type == "outlet"){ //to be used for small d or in general case
		for(int i=0;i<N_wo;i++){
			Node* n_tmp = wo[i];
			for (int j=0; j<n_tmp->b;j++) if(n_tmp->p[j]->d>0) Q_tmp+=fabs(n_tmp->p[j]->q);
			}
		if (fabs(Q_tmp)<1e-10) {
			cerr<<"Problem with flow through the system. Check if the system is not clogged."<<endl;
			s_save_data = 1;
			save_all_data();
			exit(3);}

		double factor = Q_tot/Q_tmp;
		for(int i=0;i<NP;i++) p[i]->q *=factor;   //rescaling all flows
		for(int i=0;i<NN;i++) n[i]->u *=factor;   //rescaling all pressures
	}
	else {
		cerr<<"ERROR: wrong value of type in function \
				Network::recalculate_flows_to_keep_Q_tot(string type), use \"inlet\" or \"outlet\"."<<endl;
		exit(777);
	}

}

/**
* This function returns the drop of the concentration of species B in a given pore as a result of dissolution reaction.
* The real outlet concentration is the inlet concentration times the factor returned by this function.
*
* @param pore given pore
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::outlet_c_b (Pore *p0){
	
	if(p0->q==0 || p0->d == 0) 						return 0;    //pores with no flow
	if(p0->l==l_min)           						return 1;    //no reaction in tiny pore
	if(if_track_grains && !(p0->is_Va_left()))      return 1;    //no dissolution if there is no A material

	double f = p0->local_Da_eff(this);      //effective reaction rate (taking into account both reaction and transversal diffusion)

	return exp(-f);
}   


/**
* This function returns the rise of the concentration of species C in a given pore as a result to dissolution reaction.
* The real value of the C_c_outlet  = C_c_inlet*outlet_c_c_2 (p0) + C_b_inlet*outlet_c_c_1 (p0)
* WARNING: The d_min and d_max needed to determined if the system will clogged
*  are calculated based on old concentration field.
*  There is no way to synchronized it. One can only use relatively small time step.
* @param pore given pore
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::outlet_c_c_1 (Pore *p0){

	if(p0->q==0 || p0->d ==0) return 0;   //pore with no flow
	if(!(p0->is_Va_left()))   return 0;   //no contribution if there is no A material
	if(p0->l == l_min)        return 0;   //no reaction in tiny grain


	double f1 = p0->local_Da_eff   (this);  //effective reaction rate (taking into account both reaction and transversal diffusion)
	double f2 = p0->local_Da_eff_2 (this);  //effective reaction rate (taking into account both reaction and transversal diffusion)

	double dd_plus  = p0->default_dd_plus (this);
	double dd_minus = p0->default_dd_minus(this);

	double q_tmp = fabs(p0->q);
	double c_tmp_in;
	if(if_streamtube_mixing) c_tmp_in = p0->c_in;
	else                     c_tmp_in = p0->calculate_inlet_cb();



	double x = 0;   //value to be returned
	if(p0->d + (dd_plus - dd_minus)*d0 < d_min)   { //if there is no space for full precipitation
		dd_minus=(p0->d/d0 + dd_plus - d_min/d0);
		x = c_tmp_in*q_tmp*(1-exp(-f1)) - (M_PI*(p0->d)*(dd_minus*d0)/2*p0->l)/(gamma*dt/dt_unit);
		if(x<0) {
			if (if_verbose) cerr<< "WARNING: STH wrong in outlet_c_c_1 with calculating x!!!"<<endl;
			x=0;}
		}
	else{	//if there is no problem with a space for precipitation
		if (f1!=f2) x = c_tmp_in*q_tmp*    (exp(-f1) - exp(-f2))*(f1)/(f2-f1);
		else        x = c_tmp_in*q_tmp*    f2*exp(-f2);
		}

	return x;
}


/**
* This function returns the drop of the concentration of species C in a given pore as a result to precipitation reaction.
* The real value of the C_c_outlet  = C_c_inlet*outlet_c_c_2 (p0) + C_b_inlet*outlet_c_c_1 (p0)
* @param pore given pore
* @author Agnieszka Budek
* @date 25/09/2019
*/
double Network::outlet_c_c_2 (Pore *p0){

	if(p0->q==0 || p0->d ==0) 				   return 0;      //pore with no flow
	if(p0->l==l_min)          				   return 1;	  //no reactions in tiny grain
	if(p0->d<=d_min && (!(p0->is_Va_left())))  return 1;


	double f2       = p0->local_Da_eff_2  (this);
	double dd_plus  = p0->default_dd_plus (this);
	double dd_minus = p0->default_dd_minus(this);


	//Checking if there is enough space for full dissolution
	if(p0->d + (dd_plus - dd_minus) *d0<d_min){
		if(p0->is_Va_left()) return 1;            //if there is no space for full precipitation I use ugly but working formula form outlet_c_c_2
		else                 return 1 - (p0->d-d_min)/d0/dd_minus*(1-exp(-f2));
	}
	else 		return  exp(-f2);	  //if there is enough space for precipitation

}




/**
* This function calculates the concentration field for all soluble species in the entire network.
*
* @author Agnieszka Budek
* @date 25/05/2019
*/
void Network::calculate_all_concentrations(){

	for(int i=0;i<R->bw;i++)
		calculate_concentrations(i);

}



/**
* This function calculates the concentration field for  species i.
*
* @author Agnieszka Budek
* @date 25/05/2019
*/
void Network::calculate_concentrations(int iw){

	cerr<<"Calculating concentrations for species "<<iw<<" ..."<<endl;

	for(int i=0;i<NN;i++) n[i]->tmp=i;

	//calculating nr of non zero elements of linear equations
	int R_no = 0;
	for(int i=0;i<NN;i++)
		if(n[i]->t==0) R_no+=n[i]->b+1;
		else           R_no++;


	int R_m  = NN;						    //rank of the matrix to be solved

	//for matrix containing linear equations for concentration to be solved
	int* ww_r = new int [R_no];		//raw indexes
	int* ww_c = new int [R_no];		//column indexes
	double* B = new double [R_no];	//non zero elements
	double* y = new double [NN];	//rhs


	int r_no=0;
	//filing matrix of eqs for each node
	for(int i=0;i<NN;i++){
		double S_q=0;
		if(n[i]->t==1) 	y[i] = R->C0[iw];
		else 			y[i] = 0;

		if(n[i]->t==1) S_q=-1;					//if node is an inlet node
		else for(int s=0; s<n[i]->b; s++){		//eqs for normal and outlet nodes
			Pore *pp = findPore(n[i],n[i]->n[s]);
			double qq;
			if(n[i]==pp->n[0]) qq = -pp->q;
			else               qq =  pp->q;

			if(qq > 0){
				ww_r[r_no] 	= i;
				ww_c[r_no] 	= n[i]->n[s]->tmp;
				B[r_no]		= R->w[iw]->c1(pp,this);
				S_q+=qq;
				y[i]-= R->w[iw]->c0(pp,this);
				r_no++;}}
		ww_r[r_no] 		= i;
		ww_c[r_no] 		= i;
		B[r_no]			= -S_q;

		if(S_q==0) B[r_no]=1;		//nothing is flowing into the pore
		r_no++;
	}

	//if(r_no!=R_no) {cerr<<"Problem with filling linear equations for concentration! R_no = "<<R_no<<" r_no = "<<r_no<<endl; exit(666);}
	cerr<<"Calculating concentrations: solving matrix..."<<endl;
	int M_out = solve_matrix(R_m, r_no, ww_r, ww_c, B, y);
	if(M_out!=0) cerr<<"Problem with solving linear equations; M_out = "<<M_out<<endl;

	cerr<<"Filling solution..."<<endl;
	for(int i=0;i<NN;i++) n[i]->c[iw] = y[i];     //filling the solution


	//additional printing for debugging
	print_network_for_debugging ("After calculating concentration C field ","concentration C", "flow");

	delete[] ww_r;
	delete[] ww_c;
	delete[] B;
	delete[] y;


}



/**
* This function calculates the concentration field for species B in the entire network.
* First linear equation for concentration in each node are created and
* then there are solved using function Network::solve_matrix().
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_concentrations(){

	cerr<<"Calculating concentrations for species B..."<<endl;

	//calculating no of non zero elements of linear equations
	for(int i=0;i<NN;i++) n[i]->tmp=i;

	int R_no = 0;
	for(int i=0;i<NN;i++)
		if(n[i]->t==0) R_no+=n[i]->b+1;
		else           R_no++;


	int R_m  = NN;						    //rank of the matrix to be solved

	//for matrix containing linear equations for concentration to be solved
	int* ww_r = new int [R_no];		//raw indexes
	int* ww_c = new int [R_no];		//column indexes
	double* B = new double [R_no];	//non zero elements
	double* y = new double [NN];	//rhs


	int r_no=0;
	//filing matrix with eqs for each node
	for(int i=0;i<NN;i++){
		double S_q=0;
		if(n[i]->t==1) 	y[i] = Cb_0;
		else 			y[i] = 0;

		if(n[i]->t==1) S_q=-1;					//if node is an inlet one
		else for(int s=0; s<n[i]->b; s++){		//eqs for normal and outlet nodes
			Pore *pp = findPore(n[i],n[i]->n[s]);
			double qq;
			if(n[i]==pp->n[0]) qq = -pp->q;
			else               qq =  pp->q;

			if(qq > 0){
				ww_r[r_no] 	= i;
				ww_c[r_no] 	= n[i]->n[s]->tmp;
				B[r_no]		= qq*outlet_c_b(pp);
				S_q+=qq;
				r_no++;}}
		ww_r[r_no] 		= i;
		ww_c[r_no] 		= i;
		B[r_no]			= -S_q;

		if(S_q==0) B[r_no]=1;		//nothing is flowing into the pore
		r_no++;
	}
	//if(r_no!=R_no) {cerr<<"Problem with filling linear equations for concentration! R_no = "<<R_no<<" r_no = "<<r_no<<endl; exit(666);}
	cerr<<"Calculating concentrations: solving matrix..."<<endl;
	int M_out = solve_matrix(R_m, r_no, ww_r, ww_c, B, y);
	if(M_out!=0) cerr<<"Problem with solving linear equations; M_out = "<<M_out<<endl;

	cerr<<"Filling solution..."<<endl;
	for(int i=0;i<NN;i++) n[i]->cb = y[i];     //filling the solution

	//additional printing for debugging
	print_network_for_debugging ("After calculating concentration B field ","acid concentration", "flow");

	delete[] ww_r;
	delete[] ww_c;
	delete[] B;
	delete[] y;

}

/**
* This function calculates the concentration field for species C in the entire network.
* First linear equation for concentration in each node are created and
* then there are solved using function Network::solve_matrix().
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::calculate_concentrations_c(){

	cerr<<"Calculating concentrations for species C..."<<endl;
	for(int i=0;i<NN;i++) n[i]->tmp=i;

	//calculating nr of nono zero elements of linear equations
	int R_no = 0;
	for(int i=0;i<NN;i++)
		if(n[i]->t==0) R_no+=n[i]->b+1;
		else           R_no++;


	int R_m  = NN;						    //rank of the matrix to be solved

	//for matrix containing linear equations for concentration to be solved
	int* ww_r = new int [R_no];		//raw indexes
	int* ww_c = new int [R_no];		//column indexes
	double* B = new double [R_no];	//non zero elements
	double* y = new double [NN];	//rhs

	
	int r_no=0;  
	//filing matrix of eqs for each node
	for(int i=0;i<NN;i++){ 
		double S_q=0;
		if(n[i]->t==1) 	y[i] = Cc_0;
		else 			y[i] = 0;
		
		if(n[i]->t==1) S_q=-1;					//if node is an inlet node
		else for(int s=0; s<n[i]->b; s++){		//eqs for normal and outlet nodes
			Pore *pp = findPore(n[i],n[i]->n[s]);
			double qq;
			if(n[i]==pp->n[0]) qq = -pp->q;
			else               qq =  pp->q;

			if(qq > 0){			
				ww_r[r_no] 	= i;
				ww_c[r_no] 	= n[i]->n[s]->tmp;
				B[r_no]		= qq*outlet_c_c_2(pp);
				S_q+=qq;
				y[i]-= outlet_c_c_1(pp);
				r_no++;}}
		ww_r[r_no] 		= i;
		ww_c[r_no] 		= i;
		B[r_no]			= -S_q;
		
		if(S_q==0) B[r_no]=1;		//nothing is flowing into the pore
		r_no++;	
	}

	//if(r_no!=R_no) {cerr<<"Problem with filling linear equations for concentration! R_no = "<<R_no<<" r_no = "<<r_no<<endl; exit(666);}
	cerr<<"Calculating concentrations: solving matrix..."<<endl;
	int M_out = solve_matrix(R_m, r_no, ww_r, ww_c, B, y);
	if(M_out!=0) cerr<<"Problem with solving linear equations; M_out = "<<M_out<<endl;

	cerr<<"Filling solution..."<<endl;
	for(int i=0;i<NN;i++) n[i]->cc = y[i];     //filling the solution


	//additional printing for debugging
	print_network_for_debugging ("After calculating concentration C field ","concentration C", "flow");
	
	delete[] ww_r;
	delete[] ww_c;
	delete[] B;
	delete[] y;

}


/**
* This function update the diameters of all pores due to the dissolution process.
* WARNING: No precipitation is done in this version.
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::dissolve(){

	cerr<<"Dissolving..."<<endl;

	//for updating grains volume;
	if(if_track_grains) for (int i=0;i<NG;i++) g[i]->tmp=0;

	for(int i=0;i<NP;i++){ //for each pore...

		Pore* p0 = p[i];
		if (p0->q == 0 || p0->d == 0) 				continue;
		if(if_track_grains && !p0->is_Va_left()) 	continue;

		double dd  = p0->default_dd_plus (this);

		//calculation updates for grain volumes
		int bG_tmp=0;
		for(int s=0; s<p0->bG;s++) if(p0->g[s]->Va >0) bG_tmp++;
		if(if_track_grains && bG_tmp>0) for(int s=0; s<p0->bG;s++) p0->g[s]->tmp-=(M_PI*(p0->d)*(dd*d0)/2*p0->l)/p0->bG;

		//updating diameter
		p0->d += (dd*d0);   //increasing pore diameter

		if(if_adaptive_dt)      set_adaptive_dt(dd*d0/p0->d,0);   //check if adapting dt is in needed
	}

	//updating Va and Ve (must be done after main dissolution for c_out to be calculated correctly)
	if(if_track_grains)     for(int i=0;i<NG;i++)   {g[i]->Va+=g[i]->tmp; if(g[i]->Va<0) {Va_tot-=g[i]->Va; g[i]->Va = 0;}}
	if(if_dynamical_length) for(int i=0; i<NP; i++) p[i]->calculate_actual_length(this);


	//additional printing for debugging
	print_network_for_debugging ("After dissolution ","pressure", "diameter","volume A");

}


/**
* This function update the diameters of all pores due to both dissolution and precipitation reactions.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::dissolve_and_precipitate(){

	cerr<<"Dissolving and precipitating..."<<endl;


	for(int i=0;i<NP;i++){ //for each pore...
		Pore* p0 = p[i];
		if (p0->q == 0 || p0->d == 0 || p0->l<=l_min)  continue;    //no reaction in tiny grain or in pore with no flow
		R->final_geometry_change(p0,this);
	}

	//final volumes update
	if(if_track_grains)     for(int i=0; i<NG; i++) for(int im; im<R->bm; im++) g[i]->V[im] = g[i]->V_old[im];

	//final length u[date
	if(if_dynamical_length) for(int i=0; i<NP; i++) p[i]->calculate_actual_length(this);


	//additional printing for debugging
	print_network_for_debugging ("After dissolution and precipitation","pressure", "diameter","volume A");

}


/**
* This function update the diameters of all pores due to both dissolution and precipitation reactions.
*
* @author Agnieszka Budek
* @date 25/09/2019
*/
void Network::dissolve_and_precipitate_old(){

	cerr<<"Dissolving and precipitating..."<<endl;

	//for updating grains volume;
	if(if_track_grains) for (int i=0;i<NG;i++) {g[i]->tmp=0; g[i]->tmp2=0; }

	for(int i=0;i<NP;i++){ //for each pore...


		Pore* p0 = p[i];
		if (p0->q == 0 || p0->d == 0 || p0->l==l_min)  continue;    //no reaction in tiny grain or in pore with no flow
		if (p0->d<=d_min && (!(p0->is_Va_left())))     continue;    //no reactions at all in this pore

		double d_old    = p0->d;
		double dd_plus  = p0->default_dd_plus (this);
		double dd_minus = p0->default_dd_minus(this);

		//Checking if there is enough space for full dissolution
		if(p0->d + (dd_plus - dd_minus) *d0<d_min){		//there is not enough space for all precipitating material
			dd_minus = p0->d/d0 + dd_plus - d_min/d0;
			p0->d=d_min;

		}
		else{											//there is enough space for all precipitating material
			p0->d += (dd_plus - dd_minus) *d0;
		}

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
	}


	//updating Va and Vc (must be done after main dissolution for c_out to be calculated correctly)
	if(if_track_grains){
		for (int i=0;i<NG;i++) {g[i]->Va+=g[i]->tmp;  if(g[i]->Va<0) {Va_tot-=g[i]->Va; g[i]->Va = 0;}}
		for (int i=0;i<NG;i++)  g[i]->Ve+=g[i]->tmp2;
	}

	if(if_dynamical_length) for(int i=0; i<NP; i++) p[i]->calculate_actual_length(this);


	//additional printing for debugging
	print_network_for_debugging ("After dissolution and precipitation","pressure", "diameter","volume A");

}


