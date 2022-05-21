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

#ifndef BCL_MATH_RUNNING_MIN_MAX_H_
#define BCL_MATH_RUNNING_MIN_MAX_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_limits.h"
#include "bcl_math_running_average.h"
#include "bcl_math_running_sum.h"
#include "bcl_math_statistics.h"
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
    //! @class RunningMinMax
    //! @brief class maintains a running min and max of observations given to it via +=
    //! RunningMinMax does not store observations, so it is ideal for
    //! computing the min / max of vectors or large sets of numbers that are generated if only the extrema are needed
    //!
    //! RunningMinMax works for vector/list/etc. and scalar types; for container types, the running min/max is position
    //! specific, so RunningMinMax< Vector< int> > x( Vector< double>( 2, 0));
    //!              x += MakeVector< int>(1, 3); x.GetMin()( 0) == 0, x.GetMax()( 0) == 1
    //! Like the other running stats classes, the RunningMinMax should be constructed with an example object of the
    //! proper size for container types
    //!
    //! @tparam t_AverageType the type in which to store the average
    //!
    //! @see @link example_math_running_min_max.cpp @endlink
    //! @author mendenjl
    //! @date Aug 17, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_AverageType>
    class RunningMinMax :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_AverageType m_Min; //!< Running minimum
      t_AverageType m_Max; //!< Running maximum
      mutable t_AverageType m_Range; //!< m_Max - m_Min; computed only when GetRange is called

    public:
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RunningMinMax() :
        m_Min(),
        m_Max()
      {
        Reset();
      }

      //! @brief constructor from an initial min or max, needed so that m_Average has the right size if t_ArgumentType is a container
      //! @param INITIAL_AVERAGE initial average
      //! If t_ArgumentType is a vector, INITIAL_AVERAGE must be the same size as all vectors that will be averaged
      explicit RunningMinMax( const t_AverageType &INITIAL) :
        m_Min( INITIAL),
        m_Max( INITIAL)
      {
      }

      //! virtual copy constructor
      RunningMinMax *Clone() const
      {
        return new RunningMinMax( *this);
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
        const t_AverageType &( RunningMinMax< t_AverageType>::*FUNCTION)() const
      ) const
      {
        static const std::string s_min( "Min"), s_max( "Max"), s_range( "Range"), s_unk( "Unknown");
        return FUNCTION == &RunningMinMax< t_AverageType>::GetMin ? s_min
               : FUNCTION == &RunningMinMax< t_AverageType>::GetMax ? s_max
               : FUNCTION == &RunningMinMax< t_AverageType>::GetRange ? s_range
               : s_unk;
      }

      //! @brief get the current min
      //! @return the current min
      const t_AverageType &GetMin() const
      {
        return m_Min;
      }

      //! @brief get the current max
      //! @return the current max
      const t_AverageType &GetMax() const
      {
        return m_Max;
      }

      //! @brief get the current range (max-min)
      //! @return the current range (max-min)
      const t_AverageType &GetRange() const
      {
        m_Range = m_Max;
        RunningSum< t_AverageType>::Subtract( m_Range, m_Min);
        return m_Range;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset the running average
      void Reset()
      {
        Set( m_Min, false);
        Set( m_Max, true);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief add an observation from the running average
      //! @param VALUE the value to add to the running average
      template< typename t_ArgumentType>
      const RunningMinMax &operator +=( const t_ArgumentType &VALUE)
      {
        UpdateMinMax( m_Min, m_Max, VALUE);
        return *this;
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        io::Serialize::Read( m_Min, ISTREAM);
        io::Serialize::Read( m_Max, ISTREAM);
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::Write( m_Min, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Max, OSTREAM, INDENT);
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief equals operator for built-in types, but sets value to min or max depending on PARITY
      //! @param SCALAR the scalar to = with
      //! @param PARITY true for min value, false for max value
      template< typename t_DataType>
      static typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
      Set( t_DataType &SCALAR, const bool &PARITY)
      {
        SCALAR = PARITY ? GetLowestUnboundedValue< t_DataType>() : GetHighestUnboundedValue< t_DataType>();
      }

      //! @brief equals operator for vector types, but sets value to min or max depending on PARITY
      //! @param VECTOR the vector to perform element-wise setting on
      //! @param PARITY true for min value, false for max value
      template< typename t_DataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
        Set( t_DataType &VECTOR, const bool &PARITY)
      {
        for( typename t_DataType::iterator itr( VECTOR.Begin()), itr_end( VECTOR.End()); itr != itr_end; ++itr)
        {
          RunningMinMax< t_AverageType>::Set( *itr, PARITY);
        }
      }

      //! @brief update min and max with a new value
      //! @param MIN the current min value
      //! @param MAX the current max value
      //! @param TEST the current value to consider
      template< typename t_DataType>
      static typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
        UpdateMinMax( t_DataType &MIN, t_DataType &MAX, const t_DataType &TEST)
      {
        if( TEST < MIN)
        {
          MIN = TEST;
        }
        if( TEST > MAX)
        {
          MAX = TEST;
        }
      }

      //! @brief take the square root of a vector value
      //! @param VECTOR the vector to perform element-wise sqrt on
      template< typename t_DataType, typename t_DataTypeB>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
        UpdateMinMax( t_DataType &MIN, t_DataType &MAX, const t_DataTypeB &VECTOR)
      {
        if( MIN.Begin() == MIN.End())
        {
          MIN = MAX = VECTOR;
          return;
        }
        BCL_Assert
        (
          std::distance( VECTOR.Begin(), VECTOR.End()) == std::distance( MIN.Begin(), MIN.End()),
          "Incompatible vector sizes for Min/Max calculation!"
        );
        typename t_DataType::iterator itr_min( MIN.Begin()), itr_max( MAX.Begin());
        for
        (
          typename t_DataTypeB::const_iterator itr( VECTOR.Begin()), itr_end( VECTOR.End());
          itr != itr_end;
          ++itr, ++itr_min, ++itr_max
        )
        {
          RunningMinMax< t_AverageType>::UpdateMinMax( *itr_min, *itr_max, *itr);
        }
      }

    }; // template class RunningMinMax

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RUNNING_MIN_MAX_H_
