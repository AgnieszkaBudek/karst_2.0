/*
 * solid.h
 *
 *  Created on: May 18, 2020
 *      Author: agnieszka
 */

#ifndef INCLUDE_SOLID_H_
#define INCLUDE_SOLID_H_

class Pore;
class Network;

class Solid {

public:
	virtual bool    is_reacting(Pore* p, Network* S)=0;  //if true the species is reacting in a pore
	virtual double  dd(Pore *p, Network *S)=0;           //the change in diameter due to reaction

public:
	Solid();
	virtual ~Solid();
};

class Solid_A_Agnieszkas_model : public Solid {

	public:
	virtual bool    is_reacting(Pore* p, Network* S);  //if true the species is reacting in a pore
	virtual double  dd(Pore *p, Network *S);           //the change in diameter due to reaction

};


class Solid_E_Agnieszkas_model : public Solid {

	public:
	virtual bool    is_reacting(Pore* p, Network* S);  //if true the species is reacting in a pore
	virtual double  dd(Pore *p, Network *S);           //the change in diameter due to reaction

};


class Solid_non_reacting : public Solid {

	public:
	virtual bool    is_reacting(Pore* p, Network* S)  {return false;}  //if true the species is reacting in a pore
	virtual double  dd(Pore *p, Network *S)           {return 0;}      //the change in diameter due to reaction

};


#endif /* INCLUDE_SOLID_H_ */
