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

#ifndef BCL_DESCRIPTOR_ITERATOR_HPP_
#define BCL_DESCRIPTOR_ITERATOR_HPP_

// include the header of this class
#include "bcl_descriptor_iterator.h"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_sequence_interface.h"
#include "iterate/bcl_iterate_reflecting.h"

// external includes - sorted alphabetically
#include <set>

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Iterator< t_DataType>::Iterator() :
      m_Elements(),
      m_Type(),
      m_Increment( &Iterator< t_DataType>::NotImplemented),
      m_Decrement( &Iterator< t_DataType>::NotImplemented),
      m_Position( 0),
      m_Size( 0)
    {
    }

    //! @brief constructor from type
    template< typename t_DataType>
    Iterator< t_DataType>::Iterator( const Type &TYPE) :
      m_Elements( TYPE.GetDimension()),
      m_Type( TYPE),
      m_Increment( &Iterator< t_DataType>::IncrementFirst),
      m_Decrement( &Iterator< t_DataType>::NotImplemented),
      m_Position( 0),
      m_Size( 0)
    {
      SetupIncrementer();
    }

    //! @brief constructor from type and sequence interface
    template< typename t_DataType>
    Iterator< t_DataType>::Iterator( const Type &TYPE, const SequenceInterface< t_DataType> &SEQUENCE) :
      m_Elements( TYPE.GetDimension()),
      m_Type( TYPE),
      m_Increment( &Iterator< t_DataType>::IncrementFirst),
      m_Decrement( &Iterator< t_DataType>::NotImplemented),
      m_Position( 0),
      m_Size( 0)
    {
      SetupIncrementer();
      SetObject( SEQUENCE);
    }

    //! @brief constructor for elementwise iterator given a simple generic iterator
    template< typename t_DataType>
    Iterator< t_DataType>::Iterator( const iterate::Generic< const t_DataType> &ITR) :
      m_Elements( size_t( 1), ITR),
      m_Type( size_t( 1), false, Type::e_Symmetric),
      m_Increment( &Iterator< t_DataType>::IncrementFirst),
      m_Decrement( &Iterator< t_DataType>::DecrementFirst),
      m_Position( ITR.GetPosition()),
      m_Size( ITR.GetSize())
    {
    }

    //! @brief constructor for elementwise iterator given a simple reflecting iterator
    template< typename t_DataType>
    Iterator< t_DataType>::Iterator( const iterate::Reflecting< const t_DataType> &ITR) :
      m_Elements( size_t( 1), iterate::Generic< const t_DataType>( ITR)),
      m_Type( size_t( 1), false, Type::e_Symmetric),
      m_Increment( &Iterator< t_DataType>::IncrementFirst),
      m_Decrement( &Iterator< t_DataType>::DecrementFirst),
      m_Position( ITR.GetPosition()),
      m_Size( ITR.GetSize())
    {
    }

    //! @brief constructor for pairwise iterator given two generic iterators and a type
    template< typename t_DataType>
    Iterator< t_DataType>::Iterator
    (
      const Type &TYPE,
      const iterate::Generic< const t_DataType> &ITR_A,
      const iterate::Generic< const t_DataType> &ITR_B
    ) :
      m_Elements( size_t( 2), ITR_A),
      m_Type( size_t( 2), TYPE.ConsiderRepeatedObjects(), TYPE.GetSymmetry()),
      m_Increment( &Iterator< t_DataType>::IncrementFirst),
      m_Decrement( &Iterator< t_DataType>::NotImplemented),
      m_Position( 0),
      m_Size( m_Type.GetNumberFeatures( ITR_A.GetSize()))
    {
      m_Elements( 1) = ITR_B;
      SetupIncrementer();
      m_Position =
        m_Type.GetPosition
        (
          storage::Vector< size_t>::Create( ITR_A.GetPosition(), ITR_B.GetPosition()),
          ITR_A.GetSize()
        );
    }

    //! @brief constructor from a type and vector of iterators
    template< typename t_DataType>
    Iterator< t_DataType>::Iterator
    (
      const Type &TYPE,
      const storage::Vector< iterate::Generic< const t_DataType> > &ITERATORS
    ) :
      m_Elements( ITERATORS),
      m_Type( TYPE),
      m_Increment( &Iterator< t_DataType>::IncrementFirst),
      m_Decrement( &Iterator< t_DataType>::NotImplemented),
      m_Position( 0),
      m_Size( 1)
    {
      BCL_Assert
      (
        m_Type.GetDimension() == ITERATORS.GetSize(),
        "Iterator constructor was given different # of element iterators from given type"
      );
      SetupIncrementer();
      if( ITERATORS.GetSize())
      {
        const size_t sequence_size( ITERATORS.FirstElement().GetSize());
        // need to determine the size and position for the iterators
        m_Size = m_Type.GetNumberFeatures( sequence_size);
        storage::Vector< size_t> positions( m_Elements.GetSize());
        for( size_t i( 0); i < TYPE.GetDimension(); ++i)
        {
          positions( i) = ITERATORS( i).GetPosition();
        }
        m_Position = m_Type.GetPosition( positions, sequence_size);
      }
    }

    //! @brief Clone function
    //! @return pointer to new Iterator
    template< typename t_DataType>
    Iterator< t_DataType> *Iterator< t_DataType>::Clone() const
    {
      return new Iterator( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Iterator< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the type of descriptor / iterator
    //! @return the type of this iterator / descriptor
    template< typename t_DataType>
    const Type &Iterator< t_DataType>::GetType() const
    {
      return m_Type;
    }

    //! @brief get the size, e.g. the number of valid element tuples, given the type of iteration
    //! @return the number of features
    template< typename t_DataType>
    size_t Iterator< t_DataType>::GetSize() const
    {
      return m_Size;
    }

    //! @brief get the position, e.g. the number of prior sequence element tuples that have been iterated over
    //! @return the number of prior sequence element tuples that have been iterated over
    template< typename t_DataType>
    size_t Iterator< t_DataType>::GetPosition() const
    {
      return m_Position;
    }

    //! @brief goto a particular position (slow unless dimension is 0 or 1)
    //! @param POSITION the position of interest
    template< typename t_DataType>
    void Iterator< t_DataType>::GotoPosition( const size_t &POSITION)
    {
      if( POSITION != m_Position)
      {
        BCL_Assert( POSITION < m_Size, "Cannot goto position outside of range!");
        if( m_Type.GetDimension() == size_t( 1))
        {
          m_Elements( 0).GotoPosition( POSITION);
          m_Position = POSITION;
        }
        else
        {
          if( POSITION < m_Position)
          {
            for( iterator itr( m_Elements.Begin()), itr_end( m_Elements.End()); itr != itr_end; ++itr)
            {
              itr->Restart();
            }
            m_Position = 0;
          }
          while( m_Position < POSITION)
          {
            operator++();
          }
        }
      }
    }

    //! @brief Set the object that this sequence iterator is iterating over
    //! @param SEQUENCE the sequence to iterate over
    template< typename t_DataType>
    void Iterator< t_DataType>::SetObject( const SequenceInterface< t_DataType> &SEQUENCE)
    {
      m_Position = 0;
      if( m_Type.GetDimension())
      {
        // initialize the elements appropriately
        if( m_Type.ConsiderRepeatedObjects())
        {
          std::fill( m_Elements.Begin(), m_Elements.End(), SEQUENCE.GetIterator());
        }
        else
        {
          m_Elements.FirstElement() = SEQUENCE.GetIterator();
          for( size_t i( 1), size( m_Type.GetDimension()); i < size; ++i)
          {
            m_Elements( i) = m_Elements( i - 1);
            ++m_Elements( i);
          }
        }
        m_Size = m_Type.GetNumberFeatures( m_Elements.FirstElement().GetSize());
      }
      else
      {
        // nothing to be done for empty sequence
        m_Size = size_t( 1);
      }
    }

    //! @brief access the beginning of this sequence iterator : first of the internal sequence iterators
    //! @return iterator to the beginning of the iterator vector
    template< typename t_DataType>
    typename Iterator< t_DataType>::const_iterator Iterator< t_DataType>::Begin() const
    {
      return m_Elements.Begin();
    }

    //! @brief access the end of this sequence iterator array
    //! @return iterator to the end of the iterator vector
    template< typename t_DataType>
    typename Iterator< t_DataType>::const_iterator Iterator< t_DataType>::End() const
    {
      return m_Elements.End();
    }

    //! @brief determine whether the sequence iterator is at its logical end for the current sequence
    //! @return true if this iterator is not at the end
    template< typename t_DataType>
    bool Iterator< t_DataType>::NotAtEnd() const
    {
      return m_Position < m_Size;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief dereference to gain access to the underlying iterator vector
    //! @return the underlying iterator vector
    template< typename t_DataType>
    const storage::Vector< iterate::Generic< const t_DataType> > &Iterator< t_DataType>::operator *() const
    {
      return m_Elements;
    }

    //! @brief test for less than
    template< typename t_DataType>
    bool Iterator< t_DataType>::operator <( const Iterator &ITR) const
    {
      return m_Position < ITR.m_Position;
    }

    //! @brief Access by index
    //! @param POS the position to access
    //! @return the underlying iterator
    template< typename t_DataType>
    const iterate::Generic< const t_DataType> &Iterator< t_DataType>::operator()( const size_t &POS) const
    {
      return m_Elements( POS);
    }

    //! @brief increment; move to the next tuple in the sequence based on the type
    //! @return this object
    template< typename t_DataType>
    Iterator< t_DataType> &Iterator< t_DataType>::operator++()
    {
      if( NotAtEnd() && ++m_Position != m_Size)
      {
        ( this->*m_Increment)();
      }
      return *this;
    }

    //! @brief decrement; move to the previous tuple in the sequence based on the type
    //! @return this object
    template< typename t_DataType>
    Iterator< t_DataType> &Iterator< t_DataType>::operator--()
    {
      if( m_Position)
      {
        --m_Position;
        ( this->*m_Decrement)();
      }
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Iterator< t_DataType>::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Elements, ISTREAM);
      io::Serialize::Read( m_Type, ISTREAM);
      *this = Iterator( m_Type);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_DataType>
    std::ostream &Iterator< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Elements, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Type, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determine the correct incrementation function
    template< typename t_DataType>
    void Iterator< t_DataType>::SetupIncrementer()
    {
      // choose the function that will increment this iterator
      if( m_Type.GetDimension() > size_t( 1))
      {
        if( m_Type.GetDimension() == size_t( 2))
        {
          // pairwise descriptors
          if( m_Type.ConsiderRepeatedObjects())
          {
            // all elements
            if( m_Type.GetSymmetry() == Type::e_Symmetric)
            {
              // symmetric (combinations) with replacement
              m_Increment = &Iterator< t_DataType>::IncrementUniquePairsAllElements;
              m_Decrement = &Iterator< t_DataType>::DecrementUniquePairsAllElements;
            }
            else
            {
              // permutations with replacement
              m_Increment = &Iterator< t_DataType>::IncrementAllPairsAllElements;
              m_Decrement = &Iterator< t_DataType>::DecrementAllPairsAllElements;
            }
          }
          else if( m_Type.GetSymmetry() == Type::e_Symmetric)
          {
            // combinations without replacement
            m_Increment = &Iterator< t_DataType>::IncrementUniquePairsUniqueElements;
            m_Decrement = &Iterator< t_DataType>::DecrementUniquePairsUniqueElements;
          }
          else
          {
            // permutations without replacement
            m_Increment = &Iterator< t_DataType>::IncrementAllPairsUniqueElements;
            m_Decrement = &Iterator< t_DataType>::DecrementAllPairsUniqueElements;
          }
        }
        else if( m_Type.GetDimension() == size_t( 3))
        {
          // triplet descriptors
          if( m_Type.ConsiderRepeatedObjects())
          {
            // all elements
            if( m_Type.GetSymmetry() == Type::e_Symmetric)
            {
              // symmetric (combinations) with replacement
              m_Increment = &Iterator< t_DataType>::IncrementUniqueTripletsAllElements;
            }
            else
            {
              // permutations with replacement
              m_Increment = &Iterator< t_DataType>::IncrementAllTripletsAllElements;
              m_Decrement = &Iterator< t_DataType>::DecrementAllTripletsAllElements;
            }
          }
          else if( m_Type.GetSymmetry() == Type::e_Symmetric)
          {
            // combinations without replacement
            m_Increment = &Iterator< t_DataType>::IncrementUniqueTripletsUniqueElements;
          }
          else
          {
            // permutations without replacement
            m_Increment = &Iterator< t_DataType>::IncrementAllTripletsUniqueElements;
          }
        }
        else
        {
          // higher dimensional descriptors
          if( m_Type.ConsiderRepeatedObjects())
          {
            // all elements
            if( m_Type.GetSymmetry() == Type::e_Symmetric)
            {
              // symmetric (combinations) with replacement
              m_Increment = &Iterator< t_DataType>::IncrementUniqueToAllElements;
            }
            else
            {
              // permutations with replacement
              m_Increment = &Iterator< t_DataType>::IncrementAllToAllElements;
            }
          }
          else if( m_Type.GetSymmetry() == Type::e_Symmetric)
          {
            // combinations without replacement
            m_Increment = &Iterator< t_DataType>::IncrementUniqueToUniqueElements;
          }
          else
          {
            // permutations without replacement
            m_Increment = &Iterator< t_DataType>::IncrementAllToUniqueElements;
          }
        }
      }
      else
      {
        // dimension = 0-1, scalar or elementwise descriptor
        m_Increment = &Iterator< t_DataType>::IncrementFirst;
        m_Decrement = &Iterator< t_DataType>::DecrementFirst;
      }
    }

    //! @brief increment the first (and only) iterator
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementFirst()
    {
      iterate::Generic< const t_DataType> &itr( m_Elements.FirstElement());
      ++itr;
    }

    //! @brief decrement the first (and only) iterator
    template< typename t_DataType>
    void Iterator< t_DataType>::DecrementFirst()
    {
      iterate::Generic< const t_DataType> &itr( m_Elements.FirstElement());
      --itr;
    }

    //! @brief increment over all pairs (permutations with replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementAllPairsAllElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      ++second_itr;
      if( !second_itr.NotAtEnd())
      {
        iterate::Generic< const t_DataType> &first_itr( m_Elements( 0));
        ++first_itr;
        second_itr.Restart();
      }
    }

    //! @brief increment over all pairs (permutations with replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::DecrementAllPairsAllElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      if( second_itr.AtBegin())
      {
        second_itr.GotoEnd();
        --m_Elements( 0);
      }
      --second_itr;
    }

    //! @brief increment over all pairs with unique elements (permutations without replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementAllPairsUniqueElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      iterate::Generic< const t_DataType> &first_itr( m_Elements( 0));
      ++second_itr;
      if( second_itr.NotAtEnd())
      {
        if( second_itr == first_itr)
        {
          ++second_itr;
        }
      }
      else
      {
        ++first_itr;
        second_itr.Restart();
      }
    }

    //! @brief decrement over all pairs with unique elements (permutations without replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::DecrementAllPairsUniqueElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      if( second_itr.AtBegin())
      {
        --m_Elements( 0);
        second_itr.GotoEnd();
        --second_itr;
      }
      else if( --second_itr == m_Elements( 0))
      {
        --second_itr;
      }
    }

    //! @brief increment over unique pairs (combinations with replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementUniquePairsAllElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      ++second_itr;
      if( !second_itr.NotAtEnd())
      {
        iterate::Generic< const t_DataType> &first_itr( m_Elements( 0));
        ++first_itr;
        second_itr = first_itr;
      }
    }

    //! @brief decrement over unique pairs (combinations with replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::DecrementUniquePairsAllElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      iterate::Generic< const t_DataType> &first_itr( m_Elements( 0));
      if( second_itr == first_itr)
      {
        second_itr.GotoEnd();
        --first_itr;
      }
      --second_itr;
    }

    //! @brief increment over unique pairs with unique elements (combinations without replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementUniquePairsUniqueElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      ++second_itr;
      if( !second_itr.NotAtEnd())
      {
        iterate::Generic< const t_DataType> &first_itr( m_Elements( 0));
        ++first_itr;
        second_itr = first_itr;
        ++second_itr;
      }
    }

    //! @brief increment over unique pairs with unique elements (combinations without replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::DecrementUniquePairsUniqueElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      iterate::Generic< const t_DataType> &first_itr( m_Elements( 0));
      if( --second_itr == first_itr)
      {
        second_itr.GotoEnd();
        --first_itr;
        --second_itr;
      }
    }

    //! @brief increment over all triplets (permutations with replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementAllTripletsAllElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &third_itr( m_Elements( 2));
      ++third_itr;
      if( !third_itr.NotAtEnd())
      {
        third_itr.Restart();
        IncrementAllPairsAllElements();
      }
    }

    //! @brief increment over all triplets (permutations with replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::DecrementAllTripletsAllElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &third_itr( m_Elements( 2));
      if( third_itr.AtBegin())
      {
        DecrementAllPairsAllElements();
        third_itr.GotoEnd();
      }
      --third_itr;
    }

    //! @brief increment over all triplets with unique elements (permutations without replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementAllTripletsUniqueElements()
    {
      // for 4 elements, this should produce 24 elements, e.g.:
      // 0 1 2; 0 1 3; 0 2 1; 0 2 3; 0 3 1; 0 3 2
      // 1 0 2; 1 0 3; 1 2 0; 1 2 3; 1 3 0; 1 3 2;
      // 2 0 1; 2 0 3; 2 1 0; 2 1 3; 2 3 0; 2 3 1;
      // 3 0 1; 3 0 2; 3 1 0; 3 1 2; 3 2 0; 3 2 1;

      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &third_itr( m_Elements( 2));
      ++third_itr;
      iterate::Generic< const t_DataType> &second_itr( m_Elements( 1));
      iterate::Generic< const t_DataType> &first_itr( m_Elements( 0));
      if( !third_itr.NotAtEnd())
      {
        // reached the end of the third, reset it and increment the second
        third_itr.Restart();
        ++second_itr;

        // it is not possible that 2nd hits the end here too, since that would imply that the second and third
        // elements were equal originally

        // because 2nd was just incremented, and 3rd is at 0, it is impossible that third == second
        if( second_itr == first_itr)
        {
          // 1 0 3 -> 1 2 0
          ++second_itr;
        }
        else if( third_itr == first_itr)
        {
          // 0 1 3 -> 0 2 1
          ++third_itr;
        }
      }
      else if( third_itr == first_itr)
      {
        ++third_itr;
        if( !third_itr.NotAtEnd())
        {
          // 3 1 2 -> 3 2 0
          // 1st and 3rd are right before end, only thing that needs to happen is incrementing 2nd and restarting 3rd
          third_itr.Restart();
          ++second_itr;
        }
        else if( third_itr == second_itr)
        {
          // 1 2 0 -> 1 2 3
          ++third_itr;
          if( !third_itr.NotAtEnd())
          {
            // 1 2 0 -> 2 0 1 (given only 3 elements)
            ++first_itr;
            third_itr.Restart();
            second_itr.Restart();
            ++third_itr;
          }
        }
      }
      else if( third_itr == second_itr)
      {
        // 0 2 1 -> 0 2 3
        ++third_itr;
        if( !third_itr.NotAtEnd())
        {
          // 0 2 1 -> 1 0 1
          ++first_itr;
          third_itr.Restart();
          second_itr.Restart();
          ++third_itr;
        }
        if( third_itr == first_itr)
        {
          // 1 0 1 -> 1 0 2
          // 2 1 0 -> 2 1 3
          ++third_itr;
        }
      }
    }

    //! @brief increment over unique triplets (combinations with replacement)
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementUniqueTripletsAllElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &third_itr( m_Elements( 2));
      ++third_itr;
      if( !third_itr.NotAtEnd())
      {
        IncrementUniquePairsAllElements();
        third_itr = m_Elements( 1);
      }
    }

    //! @brief increment over unique triplets with unique elements (combinations without replacement)
    //! @return true if *this can be incremented further
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementUniqueTripletsUniqueElements()
    {
      // common case, unrolled for performance
      iterate::Generic< const t_DataType> &third_itr( m_Elements( 2));
      ++third_itr;
      if( !third_itr.NotAtEnd())
      {
        IncrementUniquePairsUniqueElements();
        third_itr = m_Elements( 1);
        ++third_itr;
      }
    }

    //! @brief > 3 dimensional iteration for permutations with replacement
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementAllToAllElements()
    {
      for( reverse_iterator itr( m_Elements.ReverseBegin()), itr_end( m_Elements.ReverseEnd()); itr != itr_end; ++itr)
      {
        ++*itr;
        if( itr->NotAtEnd())
        {
          break;
        }
        itr->Restart();
      }
    }

    //! @brief > 3 dimensional iteration for permutations without replacement
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementAllToUniqueElements()
    {
      // lookup next_k_permutation code on google if this functionality is ever desired
      // the below algorithm is recursive and insanely slow in that it actually iterates over all n^k permutations
      bool found_repeat( true);
      while( found_repeat)
      {
        found_repeat = false;
        long pos( m_Elements.GetSize());
        for
        (
          reverse_iterator itr( m_Elements.ReverseBegin()), itr_end( m_Elements.ReverseEnd());
          itr != itr_end;
          ++itr, --pos
        )
        {
          ++*itr;
          if( itr->NotAtEnd())
          {
            break;
          }
          else
          {
            itr->Restart();
          }
        }

        // continue walking down the sequence, looking for the first time that incrementing an iterator does
        // not make it equal to the next iterator.  This is guaranteed to happen, so there is no loop conditional
        ++pos;
        for( const long size( m_Elements.GetSize()); pos < size; ++pos)
        {
          m_Elements( pos) = m_Elements( pos - 1);
          ++m_Elements( pos);
        }

        std::set< size_t> positions;
        for( iterator itr( m_Elements.Begin()), itr_end( m_Elements.End()); itr != itr_end; ++itr)
        {
          if( !positions.insert( itr->GetPosition()).second)
          {
            ++itr;
            // set all later elements to just before the end point to avoid needless repetitions of the above process
            while( itr != itr_end)
            {
              itr->GotoEnd();
              --*itr;
              ++itr;
            }
            found_repeat = true;
            break;
          }
        }
      }
    }

    //! @brief > 3 dimensional iteration for combinations with replacement
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementUniqueToAllElements()
    {
      // all to all
      long pos( m_Elements.GetSize() - 1);
      for
      (
        reverse_iterator itr( m_Elements.ReverseBegin()), itr_end( m_Elements.ReverseEnd());
        itr != itr_end;
        ++itr, --pos
      )
      {
        ++*itr;
        if( itr->NotAtEnd())
        {
          std::fill( m_Elements.Begin() + pos + 1, m_Elements.End(), *itr);
          break;
        }
      }
    }

    //! @brief > 3 dimensional iteration for combinations without replacement
    template< typename t_DataType>
    void Iterator< t_DataType>::IncrementUniqueToUniqueElements()
    {
      // try incrementing the last element
      reverse_iterator itr( m_Elements.ReverseBegin());
      ++*itr;
      if( itr->NotAtEnd())
      {
        // return; no problems
        return;
      }
      // return the iterator to its previous position
      --*itr;

      reverse_iterator itr_next( itr);
      ++itr_next;

      // continue walking down the sequence, looking for the first time that incrementing an iterator does
      // not make it equal to the next iterator.  This is guaranteed to happen, so there is no loop conditional
      for( long pos( m_Elements.GetSize() - 1);; --pos)
      {
        ++*itr_next;
        if( *itr_next == *itr)
        {
          --*itr_next;
          ++itr_next;
          ++itr;
          continue;
        }
        for( const long dim( m_Elements.GetSize()); pos < dim; ++pos)
        {
          m_Elements( pos) = m_Elements( pos - 1);
          ++m_Elements( pos);
        }
        break;
      }
    }

    //! @brief a trivial function that warns that a particular function is not available, then exits
    template< typename t_DataType>
    void Iterator< t_DataType>::NotImplemented()
    {
      BCL_Exit( "operator-- has not yet been implemented for type : " + util::Format()( m_Type), -1);
    }

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ITERATOR_HPP_
