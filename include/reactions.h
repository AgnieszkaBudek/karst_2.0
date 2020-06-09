/*
 * reactions.h
 *
 *  Created on: May 17, 2020
 *      Author: agnieszka
 */

#ifndef INCLUDE_REACTIONS_H_
#define INCLUDE_REACTIONS_H_

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>

#include "soluble.h"
#include "solid.h"
#include "network.h"

class Soluble;
class Solid;
class Pore;
class Network;



class Single_Reaction{

public:
	Single_Reaction();
	Single_Reaction(Soluble *, Solid *, Soluble *, Solid *);

public:

	Soluble   *sw;  ///<soluble substrate  (there may be more soluble substrates in general case, but it is not implemented yet)
	Solid     *sm;	///<solid substrate
	Soluble   *pw;  ///<soluble product
	Solid     *pm;  ///<solid product

	double      k;  //<reaction rate

};



class Reactions {

public:

	int bw;    ///< number of soluble species
	int bm;    ///< number of solid species
	int br;    ///< number of independent reactions

	Soluble** w;    	  ///< table with soluble materials
	Solid  ** m;   		  ///< table with solid materials
	Single_Reaction** r;  ///< table with all reactions in the system


public:
	// Different models that are implemented:
	void Set_reactions_to_pure_dissolution();  						///< Model with pure dissolution due to linear reaction with acid injected to the sysytem
	void Set_reactions_to_dissolution_and_precipitation_A();        ///< Model with both dissolution and precipitation (my model form PhD studies)
	void Set_reactions_to_dissolution_and_precipitation_F1();       ///< First approach to simulate Florians model

	Reactions(string model_name);
	virtual ~Reactions();

	//Functions connected with network evolution
	double default_diameter_change (Pore *, Network*);
	void   final_geometry_change   (Pore *, Network*);
	void   update_volumes          (Grain*, Network*);




};

#endif /* INCLUDE_REACTIONS_H_ */
