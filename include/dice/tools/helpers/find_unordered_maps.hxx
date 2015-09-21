#ifndef _FIND_UNORDERED_MAPS_HEADER_
#define _FIND_UNORDERED_MAPS_HEADER_

/**
 * @file 
 * @brief  Generic support for unordered maps and sparsehash : import hash table support
 * from in C++11, TR1, boost or regular map in that order depending on availability. Sparse
 * hash wrapper is also defined, using regular hash map if unavailable.
 * @author Thierry Sousbie
 */

#ifdef HAVE_CPP11
#include <unordered_map>

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 * \{
 */
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef std::unordered_map<KeyType, MappedType> type;
};
/** \}*/
#include "../../internal/namespace.footer"

#elif HAVE_TR1 // No CPP11 support

#ifdef HAVE_TR1_HEADER_PREFIX
#include <tr1/unordered_map>
#else //HAVE_TR1_HEADER_PREFIX
#include <unordered_map>
#endif //HAVE_TR1_HEADER_PREFIX


#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 * \{
 */
#ifdef HAVE_TR1_NAMESPACE_PREFIX
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef std::tr1::unordered_map<KeyType, MappedType> type;
};
#else //HAVE_TR1_NAMESPACE_PREFIX
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef std::unordered_map<KeyType, MappedType> type;
};
#endif//HAVE_TR1_NAMESPACE_PREFIX
/** \}*/
#include "../../internal/namespace.footer"

#elif defined(HAVE_BOOST) // No CPP11 or TR1 support -> use boost

#include <boost/unordered_map.hpp>

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 * \{
 */
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef boost::unordered_map<KeyType, MappedType> type;
};
/** \}*/
#include "../../internal/namespace.footer"

#else // No CPP11, No TR1 and No boost !!! -> replace with a standard map

#include <map>

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 * \{
 */
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef std::map<KeyType, MappedType> type;
};
/** \}*/
#include "../../internal/namespace.footer"

#endif // HAVE_CPP11


#ifdef HAVE_SPARSEHASH 
// specialized hash maps, for high speed (dense) or low memory (sparse)

// REMEMBER : iterators are invalidated on insert, but safe on delete
// One has to call set_hash_empty_key() exactly once before anything is inserted.
// One has to call set_hash_deleted_key exactly once before anything is deleted (only if anything is ever deleted)

#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/dense_hash_set>
#include <sparsehash/sparse_hash_set>

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 * \{
 */
template<typename KeyType, typename MappedType>
struct my_dense_hash
{
  typedef typename google::dense_hash_map<KeyType,MappedType> type;  
};

template<typename KeyType, typename MappedType>
struct my_sparse_hash
{
  typedef typename google::sparse_hash_map<KeyType,MappedType> type;  
};

template<typename KeyType>
struct my_dense_set
{
  typedef typename google::dense_hash_set<KeyType> type;  
};

template<typename KeyType>
struct my_sparse_set
{
  typedef typename google::sparse_hash_set<KeyType> type;  
};

template <class T>
void set_set_empty_key(T &set, typename T::key_type key)
{
  set.set_empty_key(key);
}

template <class T>
void set_hash_empty_key(T &hash, typename T::key_type key)
{
  hash.set_empty_key(key);
}

template <class T>
void set_hash_deleted_key(T &hash, typename T::key_type key)
{
  hash.set_deleted_key(key);
}
/** \}*/
#include "../../internal/namespace.footer"
#else // no sparsehash, emulate it ...

#include "../../internal/namespace.header"
/** \addtogroup TOOLS
 * \{
 */
template<typename KeyType, typename MappedType>
struct my_dense_hash
{
  typedef typename my_unordered_map<KeyType,MappedType>::type type;  
};

template<typename KeyType, typename MappedType>
struct my_sparse_hash
{
  typedef typename my_unordered_map<KeyType,MappedType>::type type;  
};

template <class T>
void set_hash_empty_key(T &hash, typename T::key_type key)
{}

template <class T>
void set_hash_deleted_key(T &hash, typename T::key_type key)
{}
/** \}*/
#include "../../internal/namespace.footer"
#endif // HAVE_SPARSEHASH

#endif
