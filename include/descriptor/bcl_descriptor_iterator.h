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

#ifndef BCL_DESCRIPTOR_ITERATOR_H_
#define BCL_DESCRIPTOR_ITERATOR_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_type.h"
#include "function/bcl_function_member.h"
#include "iterate/bcl_iterate_generic.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Iterator
    //! @brief provides the interface for generating descriptions for a given argument
    //! @details Iterator class provides the interface for any type of descriptor class to be derived from in order to
    //! be used in frameworks developed. It requires overwriting of following functions: operator() and GetLength()
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_iterator.cpp @endlink
    //! @author mendenjl
    //! @date Dec 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Iterator :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! elements of the internal iterator
      storage::Vector< iterate::Generic< const t_DataType> > m_Elements;

      //! iterator/descriptor type
      Type m_Type;

      //! Function used to change m_Elements when operator++ is called, based on Type
      typename function::Member< Iterator, void>::MemberFunctionPtr m_Increment;

      //! Function used to change m_Elements when operator-- is called, based on Type
      typename function::Member< Iterator, void>::MemberFunctionPtr m_Decrement;

      //! Current position
      size_t m_Position;

      //! Size = total # of t_DataType tuples that will be iterated over
      size_t m_Size;

    public:

      typedef typename storage::Vector< iterate::Generic< const t_DataType> >::const_iterator         const_iterator;
      typedef typename storage::Vector< iterate::Generic< const t_DataType> >::const_reverse_iterator const_reverse_iterator;
      typedef typename storage::Vector< iterate::Generic< const t_DataType> >::iterator               iterator;
      typedef typename storage::Vector< iterate::Generic< const t_DataType> >::reverse_iterator       reverse_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Iterator();

      //! @brief constructor from type
      Iterator( const Type &TYPE);

      //! @brief constructor from a type and sequence
      Iterator( const Type &TYPE, const SequenceInterface< t_DataType> &SEQUENCE);

      //! @brief constructor for elementwise iterator given a simple generic iterator
      Iterator( const iterate::Generic< const t_DataType> &ITR);

      //! @brief constructor for elementwise iterator given a simple reflecting iterator
      Iterator( const iterate::Reflecting< const t_DataType> &ITR);

      //! @brief constructor for pairwise iterator given two generic iterators and a type
      Iterator
      (
        const Type &TYPE,
        const iterate::Generic< const t_DataType> &ITR_A,
        const iterate::Generic< const t_DataType> &ITR_B
      );

      //! @brief constructor from a type and vector of iterators
      Iterator
      (
        const Type &TYPE,
        const storage::Vector< iterate::Generic< const t_DataType> > &ITERATORS
      );

      //! @brief Clone function
      //! @return pointer to new Iterator
      Iterator< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the type of descriptor / iterator
      //! @return the type of this iterator / descriptor
      const Type &GetType() const;

      //! @brief get the size, e.g. the number of valid element tuples, given the type of iteration
      //! @return the number of features
      size_t GetSize() const;

      //! @brief get the position, e.g. the number of prior sequence element tuples that have been iterated over
      //! @return the number of prior sequence element tuples that have been iterated over
      size_t GetPosition() const;

      //! @brief goto a particular position (slow unless dimension is 0 or 1)
      //! @param POSITION the position of interest
      void GotoPosition( const size_t &POSITION);

      //! @brief Set the object that this sequence iterator is iterating over
      //! @param SEQUENCE the sequence to iterate over
      void SetObject( const SequenceInterface< t_DataType> &SEQUENCE);

      //! @brief access the beginning of this sequence iterator : first of the internal sequence iterators
      //! @return iterator to the beginning of the iterator vector
      const_iterator Begin() const;

      //! @brief access the end of this sequence iterator array
      //! @return iterator to the end of the iterator vector
      const_iterator End() const;

      //! @brief determine whether the sequence iterator is at its logical end for the current sequence
      //! @return true if this iterator is not at the end
      bool NotAtEnd() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief dereference to gain access to the underlying iterator vector
      //! @return the underlying iterator vector
      const storage::Vector< iterate::Generic< const t_DataType> > &operator *() const;

      //! @brief test for less than
      bool operator <( const Iterator &ITR) const;

      //! @brief Access by index
      //! @param POS the position to access
      //! @return the underlying iterator
      const iterate::Generic< const t_DataType> &operator()( const size_t &POS) const;

      //! @brief increment; move to the next tuple in the sequence based on the type
      //! @return this object
      Iterator< t_DataType> &operator++();

      //! @brief decrement; move to the previous tuple in the sequence based on the type
      //! @return this object
      Iterator< t_DataType> &operator--();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief determine the correct incrementation function
      void SetupIncrementer();

      //! @brief increment/decrement the first (and only) iterator
      void IncrementFirst();
      void DecrementFirst();

      //! @brief increment/decrement over all pairs (permutations with replacement)
      void IncrementAllPairsAllElements();
      void DecrementAllPairsAllElements();

      //! @brief increment over all pairs with unique elements (permutations without replacement)
      void IncrementAllPairsUniqueElements();
      void DecrementAllPairsUniqueElements();

      //! @brief increment over unique pairs (combinations with replacement)
      void IncrementUniquePairsAllElements();
      void DecrementUniquePairsAllElements();

      //! @brief increment over unique pairs with unique elements (combinations without replacement)
      void IncrementUniquePairsUniqueElements();
      void DecrementUniquePairsUniqueElements();

      //! @brief increment/decrement over all triplets (permutations with replacement)
      void IncrementAllTripletsAllElements();
      void DecrementAllTripletsAllElements();

      //! @brief increment over all triplets with unique elements (permutations without replacement)
      void IncrementAllTripletsUniqueElements();

      //! @brief increment over unique triplets (combinations with replacement)
      void IncrementUniqueTripletsAllElements();

      //! @brief increment over unique triplets with unique elements (combinations without replacement)
      void IncrementUniqueTripletsUniqueElements();

      //! @brief > 3 dimensional iteration for permutations with replacement
      void IncrementAllToAllElements();

      //! @brief > 3 dimensional iteration for permutations without replacement
      void IncrementAllToUniqueElements();

      //! @brief > 3 dimensional iteration for combinations with replacement
      void IncrementUniqueToAllElements();

      //! @brief > 3 dimensional iteration for combinations without replacement
      void IncrementUniqueToUniqueElements();

      //! @brief a trivial function that warns that a particular function is not available, then exits
      void NotImplemented();

    }; // template class Iterator

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Iterator< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Iterator< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Iterator< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Iterator< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ITERATOR_H_
