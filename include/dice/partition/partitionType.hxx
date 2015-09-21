#ifndef __PARTITION_TYPE_HXX__
#define __PARTITION_TYPE_HXX__

#include "../tools/types/typeSelect.hxx"

/**
 * @file 
 * @brief  Definition of mesh partition and partition refinement methods 
 * @author Thierry Sousbie
 */


#include "../internal/namespace.header"

struct PartitionTypeV {
  enum Type {PH=0,   //!< Peano-hilbert (position space)
	     KWAY=1, //!< K-Way topological partitionning
	     PH_KWAY=2, //!< hybrid method (see parMetis)
	     UNDEFINED=100
  };
};

typedef PartitionTypeV::Type PartitionType;

struct PartitionTypeSelect : public TypeSelectT<PartitionTypeV> {
  PartitionTypeSelect()
  {
    this->insert("PH",PartitionTypeV::PH);
    this->insert("KWAY",PartitionTypeV::KWAY);
    this->insert("PH_KWAY",PartitionTypeV::PH_KWAY);
  }
  std::string name() {return "partitionType";}
};

struct RefinePartitionTypeV {
  enum Type {KWAY=0, //!< K-Way topological partitionning (full)
	     REFINE_KWAY=1, //!< Improve an already well balanced partionning using K-WAY
	     ADAPTIVE=2, //!< Repartition to compensate imbalance due to successive refining / coarsening
	     PH=3, //!< Peano-hilbert (position space)
	     UNDEFINED=100};
};

typedef RefinePartitionTypeV::Type RefinePartitionType;

struct RefinePartitionTypeSelect : public TypeSelectT<RefinePartitionTypeV> {
  RefinePartitionTypeSelect()
  {    
    this->insert("KWAY",RefinePartitionTypeV::KWAY);
    this->insert("REFINE_KWAY",RefinePartitionTypeV::REFINE_KWAY);
    this->insert("ADAPTIVE",RefinePartitionTypeV::ADAPTIVE);
    this->insert("PH",RefinePartitionTypeV::PH);
  }
  std::string name() {return "refinePartitionType";}
};

#include "../internal/namespace.footer"
#endif
