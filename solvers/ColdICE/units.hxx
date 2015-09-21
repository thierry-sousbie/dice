#ifndef __COLDICE_UNITS_HXX__
#define __COLDICE_UNITS_HXX__

#include <dice/cosmo/cosmology.hxx>

struct Units {
  typedef dice::cosmo::Cosmology Cosmology;
  
  double G;
  double H; // hubble param for cosmo

  double length;
  double velocity;
  //double time;
  double mass;

  int useCosmo;
  Cosmology cosmology;

  //double Gconv;
  //double Hconv;

  void update() {/*updateTime();*/updateG();updateH();}
  //void updateTime() {time=length/velocity;}
  void updateG() {G*=pow(length,-3)*mass;}
  void updateH() {};//H*=length/velocity;}

  //void updateG() {Gconv=pow(length,-3)*mass*pow(time,2);}
  //void updateH() {Hconv=length/velocity;}
};

#endif
