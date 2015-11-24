#ifndef GATHERED_GRID_SUBSET_HXX__
#define GATHERED_GRID_SUBSET_HXX__

#include "../dice_globals.hxx"

/**
 * @file 
 * @brief A regular grid shared among MPI processes
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"

/** \addtogroup GRID
 *   \{
 */

/**
 * \class GatheredGridSubsetT
 * \brief
 */

template <class CC, class G>
class GatheredGridSubsetT {
public:
  typedef typename G::Data Data;
  typedef typename CC::Coord Coord;
  //friend class G;

private:
  enum Mode {indexed=0, full=1};
  
  explicit GatheredGridSubsetT(long count, int mode[], int nReceive[])
  {
    this->count=count;
    indexedStart.assign(count+1,0);
    fullStart.assign(count+1,0);
    
    for (int i=1;i<count+1;++i)
      {
	if (mode[i-1]==indexed)
	  indexedStart[i]=nReceive[i-1];
	else
	  fullStart[i]=nReceive[i-1];

	fullStart[i]+=fullStart[i-1];
	indexedStart[i]+=indexedStart[i-1];
      }
    
    indexedData.resize(indexedStart.back());
    fullData.resize(fullStart.back());
  }
  
  long count;
  std::vector<Data> indexedData;
  std::vector<Data> fullData;

  std::vector<long> indexedStart;
  std::vector<long> fullStart;
};

/** \}*/
#include "../internal/namespace.footer"
#endif
