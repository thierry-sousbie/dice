#ifndef __DICE_COSMOLOGY_HXX__
#define __DICE_COSMOLOGY_HXX__

#include "./cosmoParameters.hxx"

#include "../internal/namespace.header" 

namespace cosmo {

  class Cosmology
  {
  public:
    typedef CosmoParameters Params;
    static std::string classHeader() {return "cosmology";}
    static float classVersion() {return 0.10;}
    static float compatibleSinceClassVersion() {return 0.10;}

    Cosmology()
    {}

    ~Cosmology()
    {}

    Cosmology(const CosmoParameters &cosmoParams)
    {
      initialize(cosmoParams);
    }

    void initialize(const CosmoParameters &cosmoParams)
    {
      params = cosmoParams;
      initPrivate();
    }

    template <class R, class PM>
    void initializeFromManager(R *reader, PM &manager,Params p=Params())
    {
      //Params p;
      p.parse(reader,manager);
      initialize(p);
    }

    template <class PP>
    void initializeFromParser(PP &parser, Params p=Params())
    {
      //Params p;
      p.parse(parser);
      initialize(p);
    }

    double dtau_da(double a)
    {
      return f_dtau_da(a);
    }
    
    double a_of_tau(double tau, 
		    double a_start, 
		    double tolerance=1.E-8, 
		    double a_end=2.0) const
    {            
      double tol2 = tolerance / 10;
      double aMin = a_start;
      double aMax = a_end;
      int level=0;      
      //printf("Computing a : tau=%g a_start=%g \n",tau,a_start);
      /*
      double deltaTau;      
      if (!cache.check(a_start,tol2))
	{	  
	  //const F_dtau_da dtau_da(params);
	  deltaTau=rombergsIntegration(a_start, 1.0, tol2, f_dtau_da);
	  cache.setDeltaTau(deltaTau);
	}
      else deltaTau = cache.getDeltaTau();
      */
      // printf("  --> deltaTau = %g\n",deltaTau);

      std::pair<bool,double> cached = cache.getRange(tau,aMin,aMax,tolerance);
      if (cached.first) 
	{
	  //printf("  --> CACHED => %g\n",cached.second);
	  return cached.second;
	}

      /*
      if ((a_start == cache.a_start)&&(tolerance==cache.tolerance))
	deltaTau=cache.deltaTau;
      else
	{
	  cache.clear();
	  const F_dtau_da dtau_da(params);
	  cache.a_start=a_start;
	  cache.tolerance=tolerance;
	  cache.deltaTau=deltaTau=rombergsIntegration(a_start, 1.0, tol2, dtau_da);
	}
      */
      double aCur = (aMax+aMin)/2.0L;    
      double t = tau_of_a(aCur,tol2);
      cache.insert(t,aCur);

      //printf("  (0)--> %g< aCur=%g <%g, a_start=%g => t=%g \n",aMin,aCur,aMax,a_start,t);

      while (fabs(t-tau) > tolerance)
	{
	  if (t>tau) aMax = aCur;
	  else aMin = aCur;
	 
	  aCur = (aMax+aMin)/2.0L;
	  t = tau_of_a(aCur,tol2);
	  if (level<=10) cache.insert(t,aCur);
	  level++;
	  //printf("  (%d)--> %g< aCur=%g <%g, a_start=%g => t=%g (%g>tol=%g) \n",level,aMin,aCur,aMax,a_start,t,fabs(t-tau),tolerance);
	};
      //printf("LEVEL = %d\n",level);
      return aCur;
    }

    double tau_of_a(double a, double tolerance = 1.E-9) const
    { 
      return rombergsIntegration(1.0,a,tolerance,f_dtau_da);     
    }

    const Params& getParams() {return params;}

    double d1_plus(double a) const
    {
      double omegaA=omega(a);
      double lambdaA=lambda(a);
      
      double tmp=pow(omegaA,4.0/7.0)-lambdaA+(1.0+0.5*omegaA)*(1.0+lambdaA/70.0);
      return 2.5*a*omegaA/tmp;
    }

    double d2_plus(double a) const
    {
      double d1p2=pow(d1_plus(a),2);
      double omegaA=omega(a);

      return -3.0*d1p2*g(omegaA)/7.0;
    }

    double f1(double a) const
    {
      return pow(omega(a),5.0/9.0);      
    }

    double f2(double a) const
    {
      return 2.0*pow(omega(a),6.0/11.0);
    }

    double h_over_h0(double a) const
    {
      double a2=a*a;
      double a3=a2*a;

      return sqrt(params.oM/a3+params.oK/a2+params.oL);
    }

  private:
    
    static double g(double omega)
    {
      return pow(omega,-1.0/143.0);
    }
    
    double lambda(double a) const
    {
      double a3=a*a*a;
      return (a3*params.oL)/(a+params.oM*(1.0-a)+params.oL*(a3-a));
    }

    double omega(double a) const
    {
      double a3=a*a*a;
      return params.oM/(a+params.oM*(1.0-a)+params.oL*(a3-a));      
    }

    void initPrivate()
    {
      f_dtau_da.initialize(params);
    }
   
    ////////////////////////////////////////////////////////////////////////////////
    //  double Rombergs_Integration_Method( double a, double h, double tolerance, //
    //                             int max_cols, double (*f)(double), int *err ); //
    //                                                                            //
    //  Description:                                                              //
    //    If T(f,h,a,b) is the result of applying the trapezoidal rule to approx- //
    //    imating the integral of f(x) on [a,b] using subintervals of length h,   //
    //    then if I(f,a,b) is the integral of f(x) on [a,b], then                 //
    //                           I(f,a,b) = lim T(f,h,a,b)                        //
    //    where the limit is taken as h approaches 0.                             //
    //    The classical Romberg method applies Richardson Extrapolation to the    //
    //    limit of the sequence T(f,h,a,b), T(f,h/2,a,b), T(f,h/4,a,b), ... ,     //
    //    in which the limit is approached by successively deleting error terms   //
    //    in the Euler-MacLaurin summation formula.                               //
    //                                                                            //
    //  Arguments:                                                                //
    //     double a          The lower limit of the integration interval.         //
    //     double h          The length of the interval of integration, h > 0.    //
    //                       The upper limit of integration is a + h.             //
    //     double tolerance  The acceptable error estimate of the integral.       //
    //                       Iteration stops when the magnitude of the change of  //
    //                       the extrapolated estimate falls below the tolerance. //
    //     int    max_cols   The maximum number of columns to be used in the      //
    //                       Romberg method.  This corresponds to a minimum       //
    //                       integration subinterval of length 1/2^max_cols * h.  //
    //     double *f         Pointer to the integrand, a function of a single     //
    //                       variable of type double.                             //
    //     int    *err       0 if the extrapolated error estimate falls below the //
    //                       tolerance; -1 if the extrapolated error estimate is  //
    //                       greater than the tolerance and the number of columns //
    //                       is max_cols.                                         //
    //                                                                            //
    //  Return Values:                                                            //
    //     The integral of f(x) from a to a +  h.                                 //
    //                                                                            //
    ////////////////////////////////////////////////////////////////////////////////
    //static const int rombergs_max_columns = 27;
    
    template <class F>
    static double rombergsIntegration( double a, double b, double tolerance, 
				       const F &f,//F *f, //double (*f)(double)
				       int &err , int max_cols=27) 
    { 
      //printf("Integrating from %g to %g\n",a,b);

      static const double rich[27] = {  
	3.333333333333333333e-01, 6.666666666666666667e-02, 1.587301587301587302e-02,
	3.921568627450980392e-03, 9.775171065493646139e-04, 2.442002442002442002e-04,
	6.103888176768601599e-05, 1.525902189669642176e-05, 3.814711817595739730e-06,
	9.536752259018191355e-07, 2.384186359449949133e-07, 5.960464832810451556e-08,
	1.490116141589226448e-08, 3.725290312339701922e-09, 9.313225754828402544e-10,
	2.328306437080797376e-10, 5.820766091685553902e-11, 1.455191522857861004e-11,
	3.637978807104947841e-12, 9.094947017737554185e-13, 2.273736754432837583e-13,
	5.684341886081124604e-14, 1.421085471520220567e-14, 3.552713678800513551e-15,
	8.881784197001260212e-16, 2.220446049250313574e-16
      };
     
      double sign=1.0;
      if (b<a)
	{
	  std::swap(a,b);
	  sign=-1.0;
	}

      double h = b-a;
      double upper_limit = a + h;     // upper limit of integration      
      double dt[27];         // dt[i] is the last element in column i.
      double integral = 0.5 * ( f(a) + f(a+h) );
      double x, old_h, delta;
      int j,k;

      // Initialize err and the first column, dt[0], to the numerical estimate //
      // of the integral using the trapezoidal rule with a step size of h.     //
 
      err = 0;
      dt[0] = 0.5 * h *  ( f(a) + f(a+h) );

      // For each possible succeeding column, halve the step size, calculate  //
      // the composite trapezoidal rule using the new step size, and up date  //
      // preceeding columns using Richardson extrapolation.                   //

      max_cols = std::min(std::max(max_cols,0),27);
      for (k = 1; k < max_cols; k++) {
	old_h = h;

	// Calculate T(f,h/2,a,b) using T(f,h,a,b) //
 
	h *= 0.5;
	integral = 0.0;
	for (x = a + h; x < upper_limit; x += old_h) integral +=  f(x);
	integral = h * integral + 0.5 * dt[0];

	//  Calculate the Richardson Extrapolation to the limit //

	for (j = 0; j < k; j++) {
	  delta =  integral - dt[j];
	  dt[j] = integral;
	  integral += rich[j] * delta;
	} 

	//  If the magnitude of the change in the extrapolated estimate //
	//  for the integral is less than the preassigned tolerance,    //
	//  return the estimate with err = 0.                           //

	if ( fabs( delta ) < tolerance ) {
	  return integral*sign;
	}
      
	//  Store the current esimate in the kth column. //

	dt[k] = integral;
      }

      // The process didn't converge within the preassigned tolerance //
      // using the maximum number of columns designated.              //
      // Return the current estimate of integral and set err = -1.    //
   
      err = -1;
      return integral*sign;
    }
    template <class F>
    static double rombergsIntegration( double a, double b, double tolerance, const F &f)
    {
      int err;
      return rombergsIntegration(a,b,tolerance,f,err);
    }
    /*
    double E_z(double z)
    {
      return 
	sqrt(params.oM*pow(1.0L+z,3)+params.oK*pow(1.0L+z,2.)+params.oL*pow(1+z,3*(1.0L+params.w)));
    }

    double E_a(double a)
    {
      return 
	sqrt(params.oM*pow(a,-3.0L)+params.oK*pow(a,-2.0L)+params.oL*pow(a,-3.0L*(1.0L+params.w)));
    }    
    */
    class F_dtau_da
    {
    public:
      
      F_dtau_da()
      {}

      explicit F_dtau_da(const CosmoParameters &p)
      {initialize(p);}

      void initialize(const CosmoParameters &p)
      {
	params=p;
      }

      double operator()(double a) const
      {
	return pow((a*a*a)*E_a(a),-1.0L);//2.0L/3.0L*
      }
    private:
      double E_a(double a) const
      {
	return 
	  sqrt(params.oM*pow(a,-3.0L)+params.oK*pow(a,-2.0L)+
	       params.oL*pow(a,-3.0L*(1.0L+params.w)));
      }    
      
      CosmoParameters params;
    };
   
    class Cache
    {
    public:
      static const int MAX_ELEMENTS=1000000;

      Cache():
	a_start(-1),
	tol2(0)
      {}

      std::pair<bool,double> getRange(double val, double &min, double &max, double tolerance)
      {
	std::pair<bool,double> result(false,double(0));
	if (m.size()==0) return result;
	//printf("SIZE : %ld\n  [%e .. %e]",m.size(),min,max);
	Map::iterator it = m.lower_bound(val);
	
	if (it!=m.end())
	  {
	    max=(*it).second;
	    if (fabs((*it).first - val)<tolerance)
	      {
		result.first=true;
		result.second=max;
	      }
	  }

	if (it!=m.begin())
	  {
	    --it;
	    min=(*it).second;
	    if (fabs((*it).first - val)<tolerance)
	      {
		result.first=true;
		result.second=min;
	      }
	  }
	//printf("=> [%e .. %e]\n",min,max);
	return result;
      }
      
      void insert(double key, double val)
      {
	if (m.size()>MAX_ELEMENTS) return;
	m.insert(std::make_pair(key,val));
      }

      void clear() {m.clear();}

      bool check(double a_start_p, double tol2_p)
      {
	if ((a_start_p!=a_start)||
	    (tol2_p!=tol2))
	  {
	    clear();
	    a_start=a_start_p;
	    tol2=tol2_p;
	    return false;
	  }
	return true;
      }

      void setDeltaTau(double dt) {deltaTau = dt;}
      double getDeltaTau() {return deltaTau;}

    private: 
      double a_start;
      double tol2;
      double deltaTau;
    
      typedef std::map<double,double> Map;
      Map m;
    };
    
    mutable Cache cache;

    F_dtau_da f_dtau_da;
    CosmoParameters params;    
  };

} // namespace cosmo

#include "../internal/namespace.footer"
#endif
