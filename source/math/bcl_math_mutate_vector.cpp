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
#include "math/bcl_math_mutate_vector.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_mutate_result.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically
#include <set>

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateVector::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateVector())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateVector::MutateVector() :
      m_Constant( util::GetUndefinedSize_t())
    {
    }

    //! @brief constructor from mutate range, number params mutated and a boolean for using percentages
    //! @param NUMBER_PARAMS Total number of parameters
    //! @param MUTATE_RANGE the range allowed for all parameters
    //! @param NUMBER_PARAMS_MUTATED the number of params mutated at each call to the operator function
    //! @param USE_PERCENTAGES boolean to indicate whether the mutate ranges are percentages
    //! @param KEEP_POSITIVE keep the values as positive numbers
    //! @param CONSTANT index to keep constant. By default, no such index is defined
    MutateVector::MutateVector
    (
      const size_t NUMBER_PARAMS,
      const double MUTATE_RANGE,
      const size_t NUMBER_PARAMS_MUTATED,
      const bool USE_PERCENTAGE,
      const bool KEEP_POSITIVE,
      const size_t &CONSTANT
    ) :
      m_Range( NUMBER_PARAMS, MUTATE_RANGE),
      m_NumberParamsMutated( NUMBER_PARAMS_MUTATED),
      m_UsePercentage( USE_PERCENTAGE),
      m_KeepPositive( KEEP_POSITIVE),
      m_Constant( CONSTANT)
    {
      BCL_Assert
      (
        NUMBER_PARAMS >= NUMBER_PARAMS_MUTATED,
        "Request to mutate " + util::Format()( NUMBER_PARAMS_MUTATED) + " parameters but only " +
        util::Format()( NUMBER_PARAMS) + " parameters exist"
      );

      // make sure the given range is non-zero
      BCL_Assert( MUTATE_RANGE != double( 0.0), "given range is 0.0");
    }

    //! @brief constructor from mutate range, number params mutated and a boolean for using percentages
    //! @param MUTATE_RANGES a vector of acceptable ranges for each parameter
    //! @param NUMBER_PARAMS_MUTATED the number of params mutated at each call to the operator function
    //! @param USE_PERCENTAGE boolean to indicate whether the mutate ranges are percentages
    //! @param KEEP_POSITIVE keep the values as positive numbers
    //! @param CONSTANT index to keep constant. By default, no such index is defined
    MutateVector::MutateVector
    (
      const linal::Vector< double> &MUTATE_RANGES,
      const size_t NUMBER_PARAMS_MUTATED,
      const bool USE_PERCENTAGE,
      const bool KEEP_POSITIVE,
      const size_t &CONSTANT
    ) :
      m_Range( MUTATE_RANGES),
      m_NumberParamsMutated( NUMBER_PARAMS_MUTATED),
      m_UsePercentage( USE_PERCENTAGE),
      m_KeepPositive( KEEP_POSITIVE),
      m_Constant( CONSTANT)
    {
      BCL_Assert
      (
        MUTATE_RANGES.GetSize() >= NUMBER_PARAMS_MUTATED,
        "Request to mutate " + util::Format()( NUMBER_PARAMS_MUTATED) + " parameters but only " +
        util::Format()( MUTATE_RANGES.GetSize()) + " parameters exist"
      );

      size_t mutates_not_zero( 0);
      for( const double *ptr( m_Range.Begin()), *ptr_end( m_Range.End()); ptr != ptr_end; ++ptr)
      {
        if( *ptr != double( 0.0))
        {
          ++mutates_not_zero;
        }
      }

      // check that there are larger equal ranges not zero than parameters to be mutated
      BCL_Assert
      (
        mutates_not_zero >= NUMBER_PARAMS_MUTATED,
        "there are less ranges different from zero that parameters to be mutated: " +
        util::Format()( mutates_not_zero) + " < " + util::Format()( NUMBER_PARAMS_MUTATED)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateVector
    MutateVector *MutateVector::Clone() const
    {
      return new MutateVector( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateVector::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking an VECTOR and returning a mutated Vector
    //! @param VECTOR Vector of parameters to be mutated
    //! @return MutateResult containing the mutated Vector
    MutateResult< linal::Vector< double> > MutateVector::operator()( const linal::Vector< double> &VECTOR) const
    {
      // indices that indicate which parameters to mutate
      std::set< size_t> indices;

      const size_t have_constant( m_Constant < VECTOR.GetSize() ? 1 : 0);
      if( m_NumberParamsMutated >= m_Range.GetSize() - have_constant)
      {
        // mutate all parameters being optimized
        for( size_t i( 0); i < m_Range.GetSize(); ++i)
        {
          // only insert index if param is actually to be mutated given by range
          if( m_Range( i) != double( 0.0) && i != m_Constant)
          {
            indices.insert( i);
          }
        }
      }
      else
      {
        // generate random indices to select which parameters to mutate
        while( indices.size() < m_NumberParamsMutated)
        {
          // select a value in the desired range; if there is a constant, the range is effectively one smaller
          // if the value ends up being >= the constant, we add 1 then to bring it back to the original range of
          // interest
          size_t chosen_value
          (
            random::GetGlobalRandom().Random< size_t>( 0, m_Range.GetSize() - 1 - have_constant)
          );
          if( chosen_value >= m_Constant)
          {
            ++chosen_value;
          }
          indices.insert( chosen_value);
        }
      }

      // create the vector of mutated parameters to return
      // initially a copy of original parameters - we overwrite the ones we mutate
      util::ShPtr< linal::Vector< double> > sp_mutated_vector( VECTOR.Clone());

      // loop over and mutate all parameters indicated by indices
      for
      (
        std::set< size_t>::const_iterator index_itr( indices.begin()), index_itr_end( indices.end());
        index_itr != index_itr_end;
        ++index_itr
      )
      {
        // get the index and look up the associated mutate range
        const size_t current_index( *index_itr);
        const double current_range( m_Range( current_index));

        if( m_KeepPositive)
        {
          // using random number generator and the mutate range, calculate the change to be applied
          double change( random::GetGlobalRandom().Random< double>( 0.0, 1.0));
          if( change < 0.45)
          {
            sp_mutated_vector->operator()( current_index) *= 0.75;
          }
          else
          {
            sp_mutated_vector->operator()( current_index) *= 1.333333;
          }
        }
        else
        {
          // using random number generator and the mutate range, calculate the change to be applied
          double change
          (
            random::GetGlobalRandom().Random< double>( -current_range, current_range)
          );

          // if percentages is being used then multiply the value with the percentage
          if( m_UsePercentage)
          {
            change *= sp_mutated_vector->operator()( current_index);
          }
          // change the value at the current_index correspondingly
          sp_mutated_vector->operator()( current_index) += change;
        }
      }

      // return the result
      return MutateResult< linal::Vector< double> >( sp_mutated_vector, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateVector::Read( std::istream &ISTREAM)
    {
      // read member variables
      io::Serialize::Read( m_Range              , ISTREAM);
      io::Serialize::Read( m_NumberParamsMutated, ISTREAM);
      io::Serialize::Read( m_UsePercentage      , ISTREAM);
      io::Serialize::Read( m_KeepPositive       , ISTREAM);
      io::Serialize::Read( m_Constant           , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateVector::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member variables
      io::Serialize::Write( m_Range              , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberParamsMutated, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UsePercentage      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_KeepPositive       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Constant           , OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
