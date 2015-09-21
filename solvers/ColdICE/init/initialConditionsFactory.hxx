#ifndef __INITIAL_CONDITIONS_FACTORY_HXX__
#define __INITIAL_CONDITIONS_FACTORY_HXX__

#include <string> 
#include <iostream> 
#include <sstream> 

#include <dice/tools/helpers/helpers.hxx>
#include <dice/tools/MPI/mpiCommunication.hxx>

template <class ICL>
class InitialConditionsFactoryT
{
public:
  typedef typename ICL::Type::Interface IC;

  template <class MT, class PM, class R>
  static IC* create(const std::string &which, PM &manager, R *reader, 
		    dice::MpiCommunication *com=0) // nullptr
  {
    typedef typename dice::hlp::ConstantValue< (ICL::SIZE>0) > Continue;
    return createHelper<ICL,MT>(which,manager,reader,com,Continue());
  }

  static void getList(std::vector<std::string> &result)
  {
    typedef typename dice::hlp::ConstantValue< (ICL::SIZE>0) > Continue;

    result.clear();
    getListHelper<ICL>(result,Continue());
  }

  static std::string getList(const char separator[]=", ")
  {
    std::vector<std::string> result;
    std::stringstream ss;
    getList(result);

    auto it = result.begin();
    if (it==result.end()) return std::string("");

    ss << *(it++);
    while (it!=result.end())
      ss << separator << *(it++);

    return ss.str();
  }

private: 
  template <class TL, class MT, class PM, class R>
  static IC* createHelper(const std::string &which, PM &manager, R *reader, 
			  dice::MpiCommunication *com,
			  dice::hlp::ConstantValue<false>)
  {
    return static_cast<IC*>(NULL);
  }

  template <class TL, class MT, class PM, class R>
  static IC* createHelper(const std::string &which, PM &manager, R *reader, 
			  dice::MpiCommunication *com,
			  dice::hlp::ConstantValue<true>)
  {
    typedef typename dice::hlp::ConstantValue< (TL::INDEX>0) > Continue;
    typedef typename dice::hlp::ConstantType< MT > MeshTraits;
    typedef typename TL::Next Next;
    typedef typename TL::Type Type;
    
     if (which == Type::name())
      {
	Type *ic = new Type(manager,reader,com);
	return static_cast<IC*>(ic);
      }

     return createHelper<Next,MT>(which,manager,reader,com,Continue());
  }  

  template <class TL>
  static void getListHelper(std::vector<std::string> &result, 
			    dice::hlp::ConstantValue<false>)
  {}

  template <class TL>
  static void getListHelper(std::vector<std::string> &result, 
			    dice::hlp::ConstantValue<true>)
  {
    typedef typename dice::hlp::ConstantValue< (TL::INDEX>0) > Continue;
    
    typedef typename TL::Next Next;
    typedef typename TL::Type Type;
    
    result.push_back(Type::name());

    getListHelper<Next>(result,Continue());    
  }  
};

#endif
