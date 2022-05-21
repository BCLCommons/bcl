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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_type.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Type::Type() :
      m_Dimension( 0),
      m_ConsiderRepeats( false),
      m_Symmetry( Type::e_Symmetric)
    {
    }

    //! @brief constructor from members: dimension, whether to consider repeats, and symmetry
    Type::Type( const size_t DIMENSION, const bool &REPEATS, const Symmetry &SYMMETRY) :
      m_Dimension( DIMENSION),
      m_ConsiderRepeats( REPEATS),
      m_Symmetry( SYMMETRY)
    {
      // if the dimension is 0 or 1, m_ConsiderRepeats and m_Symmetry are irrelevant, so set them to the most general
      // value
      if( m_Dimension < size_t( 2))
      {
        m_ConsiderRepeats = false;
        m_Symmetry = Type::e_Symmetric;
      }
    }

    //! @brief Clone function
    //! @return pointer to new Type
    Type *Type::Clone() const
    {
      return new Type( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Type::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the symmetry of the descriptor
    //! @return the symmetry of the descriptor
    const Type::Symmetry &Type::GetSymmetry() const
    {
      return m_Symmetry;
    }

    //! @brief get whether to consider repeated elements for training examples
    //! @detail for example, whether there should be a feature for elements B, B
    //! @return true if there should be features for repeated elements
    const bool &Type::ConsiderRepeatedObjects() const
    {
      return m_ConsiderRepeats;
    }

    //! @brief get the dimension of the descriptor
    //! @return the # of sub objects that go into calculating the descriptor
    const size_t &Type::GetDimension() const
    {
      return m_Dimension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the number of features for the object, given the number of sub objects
    //! @param SIZE the number of sub objects for the object (e.g. atoms in the molecule, AAs in the sequences)
    //! @return the number of features for the object given the dimension and symmetry
    size_t Type::GetNumberFeatures( const size_t &SIZE) const
    {
      if( m_Dimension == size_t( 0))
      {
        return size_t( 1);
      }
      else if( m_Dimension == size_t( 1))
      {
        return SIZE;
      }
      else if( SIZE < m_Dimension)
      {
        return 0;
      }

      // 2+ objects, have to consider symmetry and whether or not repeats are allowed
      if( m_ConsiderRepeats)
      {
        if( m_Symmetry == e_Symmetric)
        {
          // number of combinations with replacement
          return math::BinomialCoefficient( SIZE + m_Dimension - 1, m_Dimension);
        }
        else
        {
          // asymmetric, so SIZE ^ m_Dimension different ways of choosing things
          return math::Pow( SIZE, m_Dimension);
        }
      }

      // not considering repeats
      if( m_Symmetry == e_Symmetric)
      {
        // number of combinations without replacement
        return math::BinomialCoefficient( SIZE, m_Dimension);
      }

      // asymmetric, no repeats == SIZE! / (SIZE - DIM)! == B(SIZE,DIM) * DIM!
      return math::BinomialCoefficient( SIZE, m_Dimension) * math::Factorial( m_Dimension);
    }

    //! @brief given a vector of element positions, return the overall position for the iterator of this type
    //! @param POSITIONS vector of size_t positions for each element
    //! @param SIZE the number of sub objects for the object (e.g. atoms in the molecule, AAs in the sequences)
    //! @return the overall position for the iterator of this type
    size_t Type::GetPosition( const storage::Vector< size_t> &POSITIONS, const size_t &SIZE) const
    {
      if( m_Dimension == size_t( 0))
      {
        return size_t( 0);
      }
      else if( m_Dimension == size_t( 1))
      {
        return POSITIONS( 0);
      }
      else if( SIZE < m_Dimension)
      {
        return 0;
      }

      size_t position( 0);

      // handle permutations first, as they are relatively easy

      if( m_Symmetry == e_Asymmetric)
      {
        // 2+ objects, assume <= relationship for combinations
        // have to consider symmetry and whether or not repeats are allowed
        if( m_ConsiderRepeats)
        {
          // asymmetric, so SIZE ^ m_Dimension different ways of choosing things
          for
          (
            storage::Vector< size_t>::const_iterator itr( POSITIONS.Begin()), itr_end( POSITIONS.End());
            itr != itr_end;
            ++itr
          )
          {
            // Horner's algorithm, with x = SIZE, coefficients given by POSITIONS
            position = position * SIZE + *itr;
          }
        }
        else
        // asymmetric, no repeats == SIZE! / (SIZE - DIM)! == B(SIZE,DIM) * DIM!
        {
          size_t n_choices( SIZE - m_Dimension);
          size_t multiplier( 1);

          for
          (
            storage::Vector< size_t>::const_reverse_iterator
              itr( POSITIONS.ReverseBegin()), itr_end( POSITIONS.ReverseEnd());
            itr != itr_end;
            ++itr, ++n_choices, multiplier *= n_choices
          )
          {
            // because replacements are not allowed, ids are reduced by the number of previous choices
            // in the permutation that are less than the current value
            size_t effective_id( *itr);
            for( storage::Vector< size_t>::const_reverse_iterator itr_b( itr); itr_b != itr_end; ++itr_b)
            {
              if( *itr > *itr_b)
              {
                --effective_id;
              }
            }
            position += effective_id * multiplier;
          }
        }
      }
      else
      {
        // combinations
        // these algorithms could definitely be optimized better by using table lookups or more deterministic equations
        // if they come into common use
        const size_t one_less_than_dimension( m_Dimension - 1);
        size_t denominator( one_less_than_dimension);
        if( !m_ConsiderRepeats)
        {
          // no replacement
          for( size_t curr_val( SIZE - POSITIONS( 0)); curr_val < SIZE; ++curr_val)
          {
            position += math::BinomialCoefficient( curr_val, denominator);
          }
          --denominator;
          for( size_t dimension( 1); dimension < one_less_than_dimension; ++dimension, --denominator)
          {
            const size_t max_val( POSITIONS( dimension));
            for( size_t curr_val( POSITIONS( dimension - 1) + 1); curr_val < max_val; ++curr_val)
            {
              position += math::BinomialCoefficient( curr_val, denominator);
            }
          }
          position += POSITIONS( one_less_than_dimension) - POSITIONS( one_less_than_dimension - 1) - 1;
        }
        else
        {
          // combinations with replacement; size == 2
          const size_t max_val( SIZE + one_less_than_dimension);
          for( size_t curr_val( SIZE - POSITIONS( 0) + one_less_than_dimension); curr_val < max_val; ++curr_val)
          {
            position += math::BinomialCoefficient( curr_val, denominator);
          }
          --denominator;
          for( size_t dimension( 1); dimension < one_less_than_dimension; ++dimension, --denominator)
          {
            const size_t diluted_size( SIZE + denominator - POSITIONS( dimension - 1));
            const size_t diff( POSITIONS( dimension) - POSITIONS( dimension - 1));
            for( size_t curr_val( diluted_size - diff); curr_val < diluted_size; ++curr_val)
            {
              position += math::BinomialCoefficient( curr_val, denominator);
            }
          }
          position += POSITIONS( one_less_than_dimension) - POSITIONS( one_less_than_dimension - 1);
        }
      }
      return position;
    }

    //! @brief generalize this type to allow it to generate descriptors appropriate for another type
    //! @details for example become asymmetric if the other descriptor is asymmetric
    //! @param OTHER the type to generalize for
    //! @return reference to this
    Type &Type::GeneralizeToHandle( const Type &OTHER)
    {
      if( OTHER.m_Symmetry == e_Asymmetric)
      {
        m_Symmetry = e_Asymmetric;
      }
      if( OTHER.m_ConsiderRepeats)
      {
        m_ConsiderRepeats = true;
      }
      if( OTHER.m_Dimension > m_Dimension)
      {
        m_Dimension = OTHER.m_Dimension;
      }
      return *this;
    }

    //! @brief equality test
    //! @param OTHER other type to test for equality
    bool Type::operator ==( const Type &OTHER) const
    {
      return m_Dimension == OTHER.m_Dimension
             && m_ConsiderRepeats == OTHER.m_ConsiderRepeats
             && m_Symmetry == OTHER.m_Symmetry;
    }

    //! @brief inequality test
    //! @param OTHER other type to test for inequality
    bool Type::operator !=( const Type &OTHER) const
    {
      return m_Dimension != OTHER.m_Dimension
             || m_ConsiderRepeats != OTHER.m_ConsiderRepeats
             || m_Symmetry != OTHER.m_Symmetry;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Type::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Dimension, ISTREAM);
      bool symmetric;
      io::Serialize::Read( symmetric, ISTREAM);
      m_Symmetry = symmetric ? e_Symmetric : e_Asymmetric;
      io::Serialize::Read( m_ConsiderRepeats, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Type::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Dimension, OSTREAM, INDENT) << '\t' << bool( m_Symmetry) << '\t' << m_ConsiderRepeats;
      return OSTREAM;
    }

  } // namespace descriptor
} // namespace bcl
