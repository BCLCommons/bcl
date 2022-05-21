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

#ifndef BCL_MATH_RUNNING_SUM_H_
#define BCL_MATH_RUNNING_SUM_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_limits.h"
#include "bcl_math_running_average.h"
#include "io/bcl_io_serialize.h"
#include "type/bcl_type_enable_if.h"
#include "type/bcl_type_is_sequence.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RunningSum
    //! @brief class maintains a running sum of values that are passed to it
    //! This class allows summation over values or vectors of values without any need for vector interconversion
    //!
    //! @tparam t_SumType the type in which to store the average
    //!
    //! @see @link example_math_running_sum.cpp @endlink
    //! @author mendenjl
    //! @date Feb 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_SumType>
    class RunningSum :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_SumType m_Sum;           //!< Current sum
      bool          m_IsInitialized; //!< True if the sum has been initialized

    public:
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RunningSum() :
        m_Sum(),
        m_IsInitialized( false)
      {
        Reset();
      }

      //! @brief constructor from initial value
      RunningSum( const t_SumType &AVERAGE) :
        m_Sum( AVERAGE),
        m_IsInitialized( true)
      {
      }

      //! virtual copy constructor
      RunningSum *Clone() const
      {
        return new RunningSum( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns an alias for the data member returned by the given pointer
      //! @param FUNCTION the function for which to get the alias interest
      //! @return the alias as const ref std::string
      const std::string &GetAliasFromFunctionPtr
      (
        const t_SumType &( RunningSum< t_SumType>::*FUNCTION)() const
      ) const
      {
        static const std::string s_alias( "Sum");
        return s_alias;
      }

      //! @brief get the current min
      //! @return the current min
      const t_SumType &GetSum() const
      {
        return m_Sum;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset the running average
      void Reset()
      {
        RunningAverage< t_SumType>::Zero( m_Sum);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief add an observation from the running average
      //! @param VALUE the value to add to the running average
      template< typename t_ArgumentType>
      const RunningSum &operator +=( const t_ArgumentType &VALUE)
      {
        if( !m_IsInitialized)
        {
          m_Sum = VALUE;
          m_IsInitialized = true;
          return *this;
        }
        Add( m_Sum, VALUE);
        return *this;
      }

      //! @brief remove an observation from the running average
      //! @param VALUE the value to add to the running average
      template< typename t_ArgumentType>
      const RunningSum &operator -=( const t_ArgumentType &VALUE)
      {
        if( !m_IsInitialized)
        {
          m_Sum = VALUE;
          m_IsInitialized = true;
          Reset();
        }
        Subtract( m_Sum, VALUE);
        return *this;
      }

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param VALUE the value to add to the running average
      //! @param WEIGHT the weight to give the observation (may be negative to remove an old observation)
      //! @param the weight to add to
      template< typename t_ArgumentType>
      const RunningSum &AddWeightedObservation( const t_ArgumentType &VALUE, double WEIGHT)
      {
        // handle uninitialized sums
        if( !m_IsInitialized)
        {
          m_Sum = VALUE; // use operator =; this will set m_Sum to have the same dimensions as VALUE
          Reset();
          m_IsInitialized = true;
        }
        if( !WEIGHT)
        {
          // no weight added
          return *this;
        }

        // get the ratio of this observations weight to the new weight = WEIGHT / (d + WEIGHT)
        // update the average, using the weight ratio
        AddProduct( m_Sum, VALUE, WEIGHT);

        return *this;
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        io::Serialize::Read( m_Sum, ISTREAM);
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::Write( m_Sum, OSTREAM, INDENT);
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief Add to the given sum with the new value
      //! @param SCALAR_SUM the scalar sum to update
      //! @param VALUE the new value to update the scalar sum with
      //! @details overloaded version for t_DataTypes without iterators
      template< typename t_DataType, typename t_OtherDataType>
      static
      typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
      Add( t_DataType &SCALAR_SUM, const t_OtherDataType &VALUE)
      {
        SCALAR_SUM += VALUE;
      }

      //! @brief Add to the given sum with the new value
      //! @param SCALAR_SUM the scalar sum to update
      //! @param VALUE the new value to update the scalar sum with
      //! @details overloaded version for t_DataTypes without iterators
      template< typename t_DataType, typename t_OtherDataType>
      static
      typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
      AddProduct( t_DataType &SCALAR_SUM, const t_OtherDataType &VALUE, const double &WEIGHT)
      {
        SCALAR_SUM += VALUE * WEIGHT;
      }

      //! @brief Subtract from the given sum with the new value
      //! @param SCALAR_SUM the scalar sum to update
      //! @param VALUE value to subtract from the scalar sum
      //! @details overloaded version for t_DataTypes without iterators
      template< typename t_DataType, typename t_OtherDataType>
      static typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
      Subtract( t_DataType &SCALAR_SUM, const t_OtherDataType &VALUE)
      {
        SCALAR_SUM -= VALUE;
      }

      //! @brief Add to the given sum with the new values
      //! @param VECTOR_SUM the vector sum to update
      //! @param SEQUENCE the new sequence of values to add to the vector sum
      //! @details overloaded version for t_DataTypes with iterators
      template< typename t_DataType, typename t_OtherDataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
      Add( t_DataType &VECTOR_SUM, const t_OtherDataType &SEQUENCE)
      {
        // call add on each iterated-over object
        typename t_OtherDataType::const_iterator itr_seq( SEQUENCE.Begin()), itr_seq_end( SEQUENCE.End());
        for( typename t_DataType::iterator itr( VECTOR_SUM.Begin()); itr_seq != itr_seq_end; ++itr, ++itr_seq)
        {
          RunningSum< t_SumType>::Add( *itr, *itr_seq);
        }
      }

      //! @brief Add to the given sum with the new values
      //! @param VECTOR_SUM the vector sum to update
      //! @param SEQUENCE the new sequence of values to add to the vector sum
      //! @details overloaded version for t_DataTypes with iterators
      template< typename t_DataType, typename t_OtherDataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
      AddProduct( t_DataType &VECTOR_SUM, const t_OtherDataType &SEQUENCE, const double &WEIGHT)
      {
        // call add on each iterated-over object
        typename t_OtherDataType::const_iterator itr_seq( SEQUENCE.Begin()), itr_seq_end( SEQUENCE.End());
        for( typename t_DataType::iterator itr( VECTOR_SUM.Begin()); itr_seq != itr_seq_end; ++itr, ++itr_seq)
        {
          RunningSum< t_SumType>::AddProduct( *itr, *itr_seq, WEIGHT);
        }
      }

      //! @brief Subtract from the given sum with the new values
      //! @param VECTOR_SUM the vector sum to update
      //! @param SEQUENCE the new sequence of values to subtract from the vector sum
      //! @details overloaded version for t_DataTypes with iterators
      template< typename t_DataType, typename t_OtherDataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
      Subtract( t_DataType &VECTOR_SUM, const t_OtherDataType &SEQUENCE)
      {
        // call add on each iterated-over object
        typename t_OtherDataType::const_iterator itr_seq( SEQUENCE.Begin()), itr_seq_end( SEQUENCE.End());
        for( typename t_DataType::iterator itr( VECTOR_SUM.Begin()); itr_seq != itr_seq_end; ++itr, ++itr_seq)
        {
          RunningSum< t_SumType>::Subtract( *itr, *itr_seq);
        }
      }

    }; // template class RunningSum

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RUNNING_SUM_H_
