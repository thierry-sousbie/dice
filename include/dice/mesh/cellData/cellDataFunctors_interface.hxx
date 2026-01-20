#ifndef __CELL_DATA_FUNCTORS_INTERFACE_HXX__
#define __CELL_DATA_FUNCTORS_INTERFACE_HXX__

#include "../../internal/namespace.header"

namespace cellDataFunctors
{

  enum Flags
  {
    F_NO_FLAG = (0),
    F_SKIP_ON_FILE_DUMP = (1 << 0)
  };

  template <class M, class CT>
  class InterfaceT
  {
  public:
    typedef CT Cell;
    InterfaceT(const M *m, int flags = F_NO_FLAG) : mesh(m), flags_(flags) {}
    virtual ~InterfaceT() {}
    virtual double get(const Cell *s) const
    {
      return get(s, 0);
    }
    virtual double get(const Cell *s, int index) const = 0;
    virtual int get(const Cell *s, double *result) const
    {
      for (int i = 0; i < getSize(); ++i)
        result[i] = get(s, i);
      return getSize();
    }
    virtual std::string getName() const = 0;
    double operator()(const Cell *s) const { return this->get(s); }
    double operator()(const Cell *s, int index) const { return this->get(s, index); }
    int operator()(const Cell *s, double *result) const { return this->get(s, result); }
    virtual int getSize() const = 0;

    int getFlags() const { return flags_; }
    void setFlags(int flags) { flags_ = flags; }

  protected:
    const M *mesh;
    int flags_;
  };

}

#include "../../internal/namespace.footer"
#endif
