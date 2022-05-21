// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_STORAGE_HASH_MAP_H_
#define BCL_STORAGE_HASH_MAP_H_

// include the namespace header
#include "bcl_storage.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_storage_pair.h"

// external includes - sorted alphabetically
#if defined (__GNUC__)
  #if (__GNUC__ > 3 && __GNUC_MINOR__ > 2)
    #include <backward/hash_map>
  #else
    #include <ext/hash_map>
  #endif
  using namespace __gnu_cxx;
#elif defined (_MSC_VER)
  #include <hash_map>
  using namespace std;
  using namespace stdext;
#endif

#if defined (__MINGW64__)
#include <hash_fun.h>

namespace __gnu_cxx
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //! @class hash
  //! @brief a specialization of the hash necessary for the mingw compiler
  //! @author karakam
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<>
  struct hash< size_t>
  {
    size_t
    operator()( size_t __x) const
    { return __x;}
  };
} // namespace __gnu_cxx
#endif

namespace __gnu_cxx
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //! @class hash
  //! @brief specialization of the hash function for std::string
  //! @author fischea
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<>
  struct hash< std::string>
  {
    size_t operator()( const std::string &ARG) const
    {
      return hash < const char *> ()( ARG.c_str());
    }
  };

} // namespace __gnu_cxx

namespace bcl
{
  namespace storage
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HashMap
    //! @brief This is the HashMap Associative Container.
    //! @details There are several caveats in order to overcome differences
    //! between windows and Linux implementations. There are differences in the template parameters. Windows takes
    //! four: t_KeyType, t_DataType, t_Hasher, t_Alloc. GNU takes five: t_KeyType, t_DataType, t_Hasher, t_KeyComp,
    //! t_Alloc. Windows only takes four template parameters because it combines the Hash functor and KeyComp functor
    //! into one object. So the bcl HashMap will only have two template parameters in order to ensure compatibility:
    //! t_KeyType, t_DataType. The HashMap allows the key types of size_t and char *.
    //!
    //! attributes:
    //! - stores pairs of key and value
    //! - representations of the key are stored rather than the objects themselves
    //! - faster access to values than a regular Map
    //! - only single objects can be hashed into a single key.
    //!
    //! @tparam t_KeyType indicates the type of key which will be used to access data
    //! @tparam t_DataType indicates the type of data that will be accessed by the t_KeyType
    //!
    //! @see @link example_storage_hash_map.cpp @endlink
    //! @author karakam, fischea
    //! @date 04.09.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_KeyType, typename t_DataType>
    class HashMap :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    protected:

      hash_map< t_KeyType, t_DataType> m_Data;

    public:

      //! typedef for iterator
      typedef typename hash_map< t_KeyType, t_DataType>::iterator               iterator;

      //! typedef for const_iterator
      typedef typename hash_map< t_KeyType, t_DataType>::const_iterator         const_iterator;

      //! typedef for key type
      typedef typename hash_map< t_KeyType, t_DataType>::key_type               key_type;

      //! typedef for mapped type
      typedef typename hash_map< t_KeyType, t_DataType>::mapped_type            mapped_type;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      HashMap< t_KeyType, t_DataType>() :
        m_Data()
      {
      }

      //! @brief construct HashMap with two iterators [begin, end)
      //! @param ITR_BEGIN is the iterator denoting what will be the first element of the HashMap
      //! @param ITR_END is the iterator denoting the end of elements to be put into the HashMap (not inclusive)
      template< typename t_InputIterator>
      HashMap< t_KeyType, t_DataType>( t_InputIterator ITR_BEGIN, t_InputIterator ITR_END) :
        m_Data( ITR_BEGIN, ITR_END)
      {
      }

      //! @brief copy constructor
      HashMap< t_KeyType, t_DataType>
      (
        const HashMap< t_KeyType, t_DataType> &HASH_MAP
      ) :
        m_Data( HASH_MAP.m_Data)
      {
      }

      //! @brief clone copy constructor
      HashMap< t_KeyType, t_DataType> *Clone() const
      {
        return new HashMap< t_KeyType, t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief operator [] allows access to the data stored for a given key.
      //! If the key does not exist it will be used to create the data
      //! @param KEY t_KeyType for which access or creation is desired
      //! @return returns t_DataType which is the data accessed by KEY
      virtual t_DataType &operator []( const t_KeyType &KEY)
      {
        return m_Data[ KEY];
      }

      //! @brief Erase removes all the elements which have the specified key
      //! @param KEY denotes the key which will be removed
      //! @return returns the number of keys which were removed
      virtual size_t Erase( const t_KeyType &KEY)
      {
        return m_Data.erase( KEY);
      }

      //! @brief RemoveElement removes the element which is pointed to by the given iterator
      //! @param ITR iterator which denotes the element to be deleted
      virtual void RemoveElement( iterator ITR)
      {
        m_Data.erase( ITR);
      }

      //! @brief Erase removes ell elements in range denoted by two iterators
      //! @param ITR_A iterator denoting the start of elements to be deleted
      //! @param ITR_B iterator denoting the end position of elements to be deleted (non-inclusive)
      virtual void Erase( iterator ITR_A, iterator ITR_B)
      {
        m_Data.erase( ITR_A, ITR_B);
      }

      //! @brief Clear erases all elements of the associative container
      virtual void Reset()
      {
        m_Data.clear();
      }

      //! @brief Find takes a key and returns an iterator to the element denoted by key, or an off the end iterator
      //! @param KEY the key for which an iterator is desired
      //! @return returns an iterator to the element
      virtual iterator Find( const t_KeyType &KEY)
      {
        return m_Data.find( KEY);
      }

      //! @brief Find takes a key and returns a const_iterator to the element denoted by key, or off the end iterator
      //! @param KEY the key for which a const_iterator is desired
      //! @return returns an iterator to the element
      virtual const_iterator Find( const t_KeyType &KEY) const
      {
        return m_Data.find( KEY);
      }

      //! @brief Count takes a key and gives the number of times it appears in the associative container
      //! @param KEY the key for which a count is desired
      //! @return returns the number of instances of KEY
      virtual size_t Count( const t_KeyType &KEY) const
      {
        return m_Data.count( KEY);
      }

      //! @brief EqualRange gives a pair of iterators that denotes the range of the first and last instance of a key
      //! @param KEY for which the range is desired
      //! @return returns a pair of iterators denoting the first and last instance of a key
      virtual std::pair< iterator, iterator>
      EqualRange( const t_KeyType &KEY)
      {
        return m_Data.equal_range( KEY);
      }

      //! @brief returns size of the container
      //! @return size, i.e. number of elements stored
      virtual size_t GetSize() const
      {
        return m_Data.size();
      }

      //! @brief return the maximal size of the container
      //! @return maximum size, i.e. maximal number of elements to store
      virtual size_t MaxSize() const
      {
        return m_Data.max_size();
      }

      //! @brief returns a const reference to the used stl container
      //! @return const reference to the internal stl container
      virtual hash_map< t_KeyType, t_DataType> const &InternalData() const
      {
        return m_Data;
      }

      //! @brief returns a changeable reference to the used stl container
      //! @return changeable reference to the internal stl container
      virtual hash_map< t_KeyType, t_DataType> &InternalData()
      {
        return m_Data;
      }

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      virtual iterator Begin()
      {
        return m_Data.begin();
      }

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      virtual const_iterator Begin() const
      {
        return m_Data.begin();
      }

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      virtual iterator End()
      {
        return m_Data.end();
      }

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      virtual const_iterator End() const
      {
        return m_Data.end();
      }

      //! @brief Empty tests whether the container is empty or not
      //! @return boolean value indicating if the container is empty or not
      virtual bool IsEmpty() const
      {
        return m_Data.empty();
      }

      //! @brief Swap exchanges the contents of two associative containers
      //! @param SWAPPER is the container with which the contents of m_Data will be swapped
      virtual void Swap( HashMap< t_KeyType, t_DataType> &SWAPPER)
      {
        m_Data.swap( SWAPPER.m_Data);
      }

      //! @brief Insert takes a standard pair and inserts it into the map
      //! @param STD_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      virtual
      std::pair< iterator, bool>
      Insert( const std::pair< t_KeyType, t_DataType> &STD_PAIR)
      {
        return m_Data.insert( STD_PAIR);
      }

      //! @brief Insert takes a standard pair and inserts it into the map
      //! @param BCL_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      virtual
      std::pair< iterator, bool>
      Insert( const Pair< t_KeyType, t_DataType> &BCL_PAIR)
      {
        return m_Data.insert
        (
          std::pair< t_KeyType, t_DataType>( BCL_PAIR.First(), BCL_PAIR.Second())
        );
      }

      //! @brief InsertElement takes a standard pair and inserts it into the map
      //! @param STD_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      void
      InsertElement( const std::pair< t_KeyType, t_DataType> &STD_PAIR)
      {
        m_Data.insert( STD_PAIR);
      }

      //! @brief Insert takes a BCL pair and inserts it into the map
      //! @param BCL_PAIR is the pair that is desired to be inserted into the map
      //! @return a std pair with an iterator pointing to the place of insertion and a bool indicating success or not
      //! if the element already existed then the iterator points to the previous existance and the bool is false
      void
      InsertElement( const Pair< t_KeyType, t_DataType> &BCL_PAIR)
      {
        m_Data.insert
        (
          std::pair< t_KeyType, t_DataType>( BCL_PAIR.First(), BCL_PAIR.Second())
        );
      }

      //! @brief InsertElements inserts a number of elements into the container [INCLUSIVE, NOT_INCLUSIVE)
      //! @param ITR_A iterator denoting the first element to be inserted
      //! @param ITR_B iterator denoting the end of the elements to be inserted (not included in the insertions)
      template< typename t_Iterator>
      void InsertElements( const t_Iterator ITR_A, const t_Iterator ITR_B)
      {
        m_Data.insert( ITR_A, ITR_B);
      }

      //! @brief InsertElements insert all elements of argument CONTAINER
      //! @param MAP container in which objects are inserted
      virtual void
      InsertElements( const HashMap< t_KeyType, t_DataType> &MAP)
      {
        m_Data.insert
        (
          MAP.Begin(), MAP.End()
        );
      }

      //! @brief InsertElement takes a standard pair and inserts it into the UniqueSortedAssociativeContainer
      //! with the passed iterator as a suggestion on where to start looking for the insertion
      //! @tparam t_Iterator is the type of iterator that suggests where the insert should start looking
      //! @param ITR is the iterator denoting the position to start looking for the correct insertion point
      //! @param STD_PAIR is the pair that is desired to be inserted into the UniqueSortedAssociativeContainer
      //! @return an iterator pointing to the place of insertion
      iterator
      InsertElement( iterator ITR, const std::pair< t_KeyType, t_DataType> &STD_PAIR)
      {
        return m_Data.insert( ITR, STD_PAIR);
      }

      //! @brief InsertElement takes a BCL pair and inserts it into the UniqueAssociativeContainer
      //! with the passed iterator as a suggestion on where to start looking for the insertion
      //! @tparam t_Iterator is the type of iterator that suggests where the insert should start looking
      //! @param ITR is the iterator denoting the position to start looking for the correct insertion point
      //! @param BCL_PAIR is the pair that is desired to be inserted into the UniqueSortedAssociativeContainer
      //! @return an iterator pointing to the place of insertion
      iterator
      InsertElement( iterator ITR, const Pair< t_KeyType, t_DataType> &BCL_PAIR)
      {
        return m_Data.insert
        (
          ITR,
          std::pair< t_KeyType, t_DataType>( BCL_PAIR.First(), BCL_PAIR.Second())
        );
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write size of container
        io::Serialize::Write( GetSize(), OSTREAM, INDENT) << '\n';

        // write each element of container to new line
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          io::Serialize::Write( std::pair< t_KeyType, t_DataType>( *itr), OSTREAM, INDENT) << '\n';
        }

        // end
        return OSTREAM;
      }

      //! @brief read container from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // empty the container
        Reset();

        // create int "size" for holding the size of the container
        size_t size;

        // read in size of container
        io::Serialize::Read( size, ISTREAM);

        // read in each element from stream
        for( size_t count( 0); count < size; ++count)
        {
          std::pair< t_KeyType, t_DataType> current;

          BCL_Assert( io::Serialize::Read( current, ISTREAM), "Error reading element!");

          // check that insert was successful
          BCL_Assert( Insert( current).second, "error inserting");
        }

        // return
        return ISTREAM;
      }

    }; // template class HashMap

    // instantiate s_Instance
    template< typename t_KeyType, typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> HashMap< t_KeyType, t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new HashMap< t_KeyType, t_DataType>())
    );

    //! @brief operator == checks if two containers are the same
    //! @param HASHMAP_A first container
    //! @param HASHMAP_B second container
    //! @return returns boolean value
    template< typename t_KeyType, typename t_DataType>
    inline bool operator ==
    (
      const HashMap< t_KeyType, t_DataType> &HASHMAP_A,
      const HashMap< t_KeyType, t_DataType> &HASHMAP_B
    )
    {
      return HASHMAP_A.InternalData() == HASHMAP_B.InternalData();
    }

  } // namespace storage
} // namespace bcl

#endif //BCL_STORAGE_HASH_MAP_H_
