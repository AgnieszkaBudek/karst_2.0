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
void Network::calculate_concentrations_b_diff(){

	cerr<<"Calculating concentrations for species B taking into account diffusion..."<<endl;

	//calculating no of non zero elements of linear equations
	for(int i=0;i<NN;i++) n[i]->tmp=i;

	int R_no = 0;
	for(int i=0;i<NN;i++)
		if(n[i]->t!=1) R_no+=n[i]->b+1;
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

		if(n[i]->t==1) S_q=1;					//if node is an inlet one
		else for(int s=0; s<n[i]->b; s++){		//eqs for normal and outlet nodes

			Pore *pp = findPore(n[i],n[i]->n[s]);
			double qq;
			if(n[i]==pp->n[0]) qq = -pp->q;
			else               qq =  pp->q;
			if(fabs(qq)<=1e-10)  continue;  //epsilon for numerical stability, one has to check the best value

			ww_r[r_no] 	= i;
			ww_c[r_no] 	= n[i]->n[s]->tmp;
			if (n[i]->t==-1) {
				//y[i]=0.5;                                        //for debugging
				//B[r_no++] = qq*outlet_c_b(pp);  S_q += -qq;      //first attempt, but mass not conserved
				B[r_no++] = qq*outlet_c_b_1_d(pp,_sign(qq)); S_q += qq*outlet_c_b_0_d(pp,_sign(qq))-qq;}
			else           {
				B[r_no++] = qq*outlet_c_b_1_d(pp,_sign(qq)); S_q += qq*outlet_c_b_0_d(pp,_sign(qq));}
			}

		ww_r[r_no] 		= i;
		ww_c[r_no] 		= i;
		B[r_no]			= S_q;

		if(S_q==0) B[r_no]=1;		//nothing is flowing into the pore
		r_no++;
	}

// Printing entire matrix for debugging purpose
//	double M[NN][NN];
//	cerr<<endl<<endl;
//	for (int i=0;i<NN;i++) for(int j=0;j<NN;j++)
//		M[i][j]=0;
//	for (int i=0;i<r_no;i++) M[ww_r[i]][ww_c[i]] = B[i];
//	for (int i=0;i<NN;i++){
//		for(int j=0;j<NN;j++) cerr<<setw(7)<<setprecision(4)<<M[i][j];
//		cerr<<endl;
//	}



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


void Network::calculate_concentrations_c_diff(){

	cerr<<"Has to be implemented!!!"<<endl;
	exit(123);

}



double Network::outlet_c_b_1_d   (Pore* p, int s) {

	if(p->d==0) return 0;
	if(p->q==0) return 0;

	//pore connected to the outlet does not track diffusion
	//if(p->n[0]->t==-1 || p->n[1]->t==-1)	return 1; //problems when dc=0 in the last pore (with small q)

	double pe = p->local_Pe    (this);
	double da = p->local_Da_eff(this);

	double a = sqrt(pe*(4*da+pe));
	double b = s*pe/2.;

	return - (a*exp(a/2.+b)) / (2.*b*(1-exp(a)));
}



double Network::outlet_c_b_0_d   (Pore* p, int s) {

	if(p->d==0) return 0;
	if(p->q==0) return 0;


	//pore connected to the outlet does not track diffusion
	//if(p->n[0]->t==-1 || p->n[1]->t==-1)	return 1; //problems when dc=0 in the last pore (with small q)


	double pe = p->local_Pe    (this);
	double da = p->local_Da_eff(this);

	double a = sqrt(pe*(4*da+pe));
	double b = s*pe/2.;

	return  (a + 2*b + (a - 2*b)*exp(a)) / (4.*b*(1-exp(a)));

}


double Network::outlet_c_c_1_d   (Pore* p,int s) {return 0;}   //to be implemented
double Network::outlet_c_c_0_d   (Pore* p,int s) {return 0;}   //to be implemented
