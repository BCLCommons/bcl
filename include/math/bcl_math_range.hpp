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

#ifndef BCL_MATH_RANGE_HPP_
#define BCL_MATH_RANGE_HPP_

// include the header of this class
#include "bcl_math_range.h"

// includes from bcl - sorted alphabetically
#include "bcl_math_comparisons.h"
#include "io/bcl_io_serialize.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Range< t_DataType>::s_Instance( GetObjectInstances().AddInstance( new Range< t_DataType>()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor - initialize min and max to undefined and the borders to a closed interval
    template< typename t_DataType>
    Range< t_DataType>::Range() :
      m_Borders(),
      m_Min( util::GetUndefined< t_DataType>()),
      m_Max( util::GetUndefined< t_DataType>())
    {
    }

    //! @brief construct from min and max and two border conditions
    //! @param LEFT_BORDER_CONDITION border condition to the left
    //! @param MIN lower (left-hand) bound of closed interval
    //! @param MAX upper (right-hand) bound of closed interval
    //! @param RIGHT_BORDER_CONDITION border condition to the right
    template< typename t_DataType>
    Range< t_DataType>::Range
    (
      const RangeBorders::ConditionLeft LEFT_BORDER_CONDITION,
      const t_DataType &MIN,
      const t_DataType &MAX,
      const RangeBorders::ConditionRight RIGHT_BORDER_CONDITION
    ) :
      m_Borders( LEFT_BORDER_CONDITION, RIGHT_BORDER_CONDITION),
      m_Min( MIN),
      m_Max( MAX)
    {
      // check that range is valid
      BCL_Assert
      (
        m_Min <= m_Max,                              // min < max for one side open interval
        "invalid range: " + GetString()
      );
    }

    //! @brief construct from min and max [min,max] in a closed interval
    //! @param MIN lower (left-hand) bound of closed interval
    //! @param MAX upper (right-hand) bound of closed interval
    template< typename t_DataType>
    Range< t_DataType>::Range
    (
      const t_DataType &MIN,
      const t_DataType &MAX
    ) :
      m_Borders( RangeBorders::e_LeftClosed, RangeBorders::e_RightClosed),
      m_Min( MIN),
      m_Max( MAX)
    {
      // check that min is smaller than max
      BCL_Assert( m_Min <= m_Max, "invalid range: " + GetString());
    }

    //! @brief Clone function
    //! @return pointer to new Range
    template< typename t_DataType>
    Range< t_DataType> *Range< t_DataType>::Clone() const
    {
      return new Range< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Range< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to min of range
    //! @return min
    template< typename t_DataType>
    const t_DataType &Range< t_DataType>::GetMin() const
    {
      return m_Min;
    }

    //! @brief access to max of range
    //! @return max
    template< typename t_DataType>
    const t_DataType &Range< t_DataType>::GetMax() const
    {
      return m_Max;
    }

    //! @brief access to left border condition
    //! @return ConditionLeft open or closed
    template< typename t_DataType>
    RangeBorders::ConditionLeft Range< t_DataType>::GetLeftCondition() const
    {
      return m_Borders.GetLeftCondition();
    }

    //! @brief access to right border condition
    //! @return ConditionRight open or closed
    template< typename t_DataType>
    RangeBorders::ConditionRight Range< t_DataType>::GetRightCondition() const
    {
      return m_Borders.GetRightCondition();
    }

    //! @brief return the width
    //! @return max-min
    template< typename t_DataType>
    t_DataType Range< t_DataType>::GetWidth() const
    {
      return m_Max - m_Min;
    }

    //! @brief return the middle value of the interval
    //! @return the middle of the interval
    template< typename t_DataType>
    t_DataType Range< t_DataType>::GetMiddle() const
    {
      return t_DataType( m_Min + GetWidth() * 0.5);
    }

    //! @brief convert to a human readable string with formatted numbers
    //! @return string like "[2.5,20.1]"
    template< typename t_DataType>
    std::string Range< t_DataType>::GetString() const
    {
      return this->GetString( util::Format());
    }

    //! @brief convert to a human readable string with formatted numbers
    //! @param FORMAT format for the numbers
    //! @return string like "[2.5,20.1]"
    template< typename t_DataType>
    std::string Range< t_DataType>::GetString( const util::Format &FORMAT) const
    {
      std::string range_string;
      range_string += m_Borders.GetConditionLeftChar();
      range_string += FORMAT( m_Min) + ",";
      range_string += FORMAT( m_Max);
      range_string += m_Borders.GetConditionRightChar();

      // end
      return range_string;
    }

    //! @brief get the label containing only initialization parameters
    //! @param WITH_DATA whether to include any data members, else, only include initialization members
    template< typename t_DataType>
    util::ObjectDataLabel Range< t_DataType>::GetLabel( const bool &WITH_DATA) const
    {
      return util::ObjectDataLabel( "", GetString());
    }

    //! @brief change the min by VALUE
    //! @param MIN the new minimum value
    //! @return true, if range could be changed (range is valid after operation)
    template< typename t_DataType>
    bool Range< t_DataType>::SetMin( const t_DataType &MIN)
    {
      if( MIN <= m_Min || IsWithin( MIN))
      {
        m_Min = MIN;
        return true;
      }

      // end
      return false;
    }

    //! @brief change the max by VALUE
    //! @param MAX the new maximum value
    //! @return true, if range could be changed (range is valid after operation)
    template< typename t_DataType>
    bool Range< t_DataType>::SetMax( const t_DataType &MAX)
    {
      if( MAX >= m_Max || IsWithin( MAX))
      {
        m_Max = MAX;
        return true;
      }

      // end
      return false;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief rescale a value from a source range to this range
    //! @param VALUE the value to be rescaled to this range
    //! @param SOURCE_RANGE the range the VALUE comes from
    //! @return rescaled value
    template< typename t_DataType>
    t_DataType Range< t_DataType>::Rescale( const t_DataType &VALUE, const Range< t_DataType> &SOURCE_RANGE) const
    {
      return t_DataType( ( VALUE - SOURCE_RANGE.m_Min) * ( GetWidth() / double( SOURCE_RANGE.GetWidth())) + m_Min);
    }

    //! @brief check value for being in this range
    //! @param VALUE
    //! @return true, if within range [min,max]
    template< typename t_DataType>
    bool Range< t_DataType>::IsWithin( const t_DataType &VALUE) const
    {
      // check left border of interval (min)
      if( m_Borders.IsLeftBorderClosed())
      {
        if( VALUE < m_Min) return false;
      }
      else if( VALUE <= m_Min)
      {
        return false;
      }

      // check right border of interval (max)
      if( m_Borders.IsRightBorderClosed())
      {
        if( VALUE > m_Max) return false;
      }
      else if( VALUE >= m_Max)
      {
        return false;
      }

      // did not return false yet
      return true;
    }

    //! @brief check whether the range contains any values
    //! @return true if there are no values in the range, or it is undefined
    template< typename t_DataType>
    bool Range< t_DataType>::IsEmpty() const
    {
      // min cannot be larger than max, since this is asserted in the constructor

      // if min is undefined, the range is empty
      // note that max may be undefined without the range being empty,
      // e.g. Range< size_t>( size_t( 0), std::numeric_limits< size_t>())
      if( util::IsNaN( m_Min))
      {
        return true;
      }

      // if min and max are equal, the range is empty if either border is open
      if( m_Max == m_Min)
      {
        return !m_Borders.IsClosed();
      }

      // m_Min < m_Max, so the only way that the range could be empty is if both borders are open
      if( m_Borders.IsOpen())
      {
        // in this case, the range is not empty if the middle is within the range
        return std::numeric_limits< t_DataType>::is_integer && t_DataType( m_Min + 1) == m_Max;
      }

      // the range is not empty, because at least one side is closed and m_Max > m_Min
      return false;
    }

    //! @brief check if one range falls within another range
    //! @param VALUE
    //! @return true, if there is overlap within a range
    template< typename t_DataType>
    bool Range< t_DataType>::DoesOverlap( const Range< t_DataType> &VALUE) const
    {
      // determine comparison of border conditions
      typename Comparisons< t_DataType>::Comparison valueleft_thisleft
      (
        VALUE.m_Borders.IsLeftBorderClosed() && m_Borders.IsLeftBorderClosed() ?
          Comparisons< t_DataType>::GetEnums().e_GreaterEqual : Comparisons< t_DataType>::GetEnums().e_Greater
      );
      typename Comparisons< t_DataType>::Comparison valueleft_thisright
      (
        VALUE.m_Borders.IsLeftBorderClosed() && m_Borders.IsRightBorderClosed() ?
          Comparisons< t_DataType>::GetEnums().e_LessEqual : Comparisons< t_DataType>::GetEnums().e_Less
      );
      typename Comparisons< t_DataType>::Comparison valueright_thisleft
      (
        VALUE.m_Borders.IsRightBorderClosed() && m_Borders.IsLeftBorderClosed() ?
          Comparisons< t_DataType>::GetEnums().e_GreaterEqual : Comparisons< t_DataType>::GetEnums().e_Greater
      );
      typename Comparisons< t_DataType>::Comparison valueright_thisright
      (
        VALUE.m_Borders.IsRightBorderClosed() && m_Borders.IsRightBorderClosed() ?
          Comparisons< t_DataType>::GetEnums().e_LessEqual : Comparisons< t_DataType>::GetEnums().e_Less
      );
      typename Comparisons< t_DataType>::Comparison thisleft_valueright
      (
        VALUE.m_Borders.IsRightBorderClosed() && m_Borders.IsLeftBorderClosed() ?
          Comparisons< t_DataType>::GetEnums().e_LessEqual : Comparisons< t_DataType>::GetEnums().e_Less
      );

      typename Comparisons< t_DataType>::Comparison thisright_valueleft
      (
        VALUE.m_Borders.IsLeftBorderClosed() && m_Borders.IsRightBorderClosed() ?
          Comparisons< t_DataType>::GetEnums().e_GreaterEqual : Comparisons< t_DataType>::GetEnums().e_Greater
      );

      bool does_overlap( false);

      // true if min_val falls within the range of *this Range
      does_overlap |= ( *valueleft_thisleft)->operator()( VALUE.m_Min, m_Min) &&
        ( *valueleft_thisright)->operator()( VALUE.m_Min, m_Max);

      // true if m_Min value falls in the range of "VALUE"
      does_overlap |= ( *valueleft_thisleft)->operator()( m_Min, VALUE.m_Min) &&
        ( *thisleft_valueright)->operator()( m_Min, VALUE.m_Max);

      // true if max_val falls in the range of *this Range
      does_overlap |= ( *valueright_thisleft)->operator()( VALUE.m_Max, m_Min) &&
        ( *valueright_thisright)->operator()( VALUE.m_Max, m_Max);

      // true if m_Max falls in the range of "VALUE"
      does_overlap |= ( *thisright_valueleft)->operator()( m_Max, VALUE.m_Min) &&
        ( *valueright_thisright)->operator()( m_Max, VALUE.m_Max);

      return does_overlap;
    }

    //! @brief combine overlapping ranges into a single range
    //! @param RANGE_1 the first range to combine
    //! @param RANGE_2 the second range to combine
    //! @return the overlapping ranges
    template< typename t_DataType>
    Range< t_DataType> Range< t_DataType>::CombineRanges
    (
      const Range< t_DataType> &RANGE_1,
      const Range< t_DataType> &RANGE_2
    )
    {
      // make sure the ranges overlapped
      BCL_Assert( RANGE_1.DoesOverlap( RANGE_2), "Cannot combine non-overlapping ranges!");

      // make a new range to hold the information
      Range< t_DataType> new_range;

      // determine the new minimum of the range
      if( RANGE_1.m_Min < RANGE_2.m_Min)
      {
        new_range.m_Min = RANGE_1.m_Min;
        new_range.m_Borders.Set( RANGE_1.GetLeftCondition());
      }
      else if( RANGE_2.m_Min < RANGE_1.m_Min)
      {
        new_range.m_Min = RANGE_2.m_Min;
        new_range.m_Borders.Set( RANGE_2.GetLeftCondition());
      }
      else // RANGE_2.m_Min == RANGE_1.m_Min
      {
        new_range.m_Min = RANGE_1.m_Min;

        // was the first range closed ?
        if( RANGE_1.m_Borders.IsLeftBorderClosed())
        {
          // yep, so the new range is closed
          new_range.m_Borders.Set( RANGE_1.GetLeftCondition());
        }
        else
        {
          // nope, so use the condition on the second range
          new_range.m_Borders.Set( RANGE_2.GetLeftCondition());
        }
      }

      // determine the new maximum of the range
      if( RANGE_1.m_Max > RANGE_2.m_Max)
      {
        new_range.m_Max = RANGE_1.m_Max;
        new_range.m_Borders.Set( RANGE_1.GetRightCondition());
      }
      else if( RANGE_2.m_Max > RANGE_1.m_Max)
      {
        new_range.m_Max = RANGE_2.m_Max;
        new_range.m_Borders.Set( RANGE_2.GetRightCondition());
      }
      else // RANGE_2.m_Max == RANGE_1.m_Max
      {
        new_range.m_Max = RANGE_1.m_Max;

        // was the first range closed ?
        if( RANGE_1.m_Borders.IsRightBorderClosed())
        {
          // yep, so the new range is closed
          new_range.m_Borders.Set( RANGE_1.GetRightCondition());
        }
        else
        {
          // nope, so use the other range's condition
          new_range.m_Borders.Set( RANGE_2.GetRightCondition());
        }
      }

      return new_range;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @param ERROR_STREAM stream to write errors to
    //! @return true on success
    template< typename t_DataType>
    bool Range< t_DataType>::FromStream( std::istream &ISTREAM, std::ostream &ERROR_STREAM)
    {
      // read member
      char current_character;
      std::string::size_type border_cond;

      // read left border
      ISTREAM >> current_character;
      border_cond = RangeBorders::GetConditionLeftChars().find( current_character);
      if( border_cond == std::string::npos)
      {
        ERROR_STREAM << "expected either character: \"" << RangeBorders::GetConditionLeftChars()
                     << "\" as left-border of range but found '" << current_character << "'";
        return false;
      }
      m_Borders.Set( RangeBorders::ConditionLeft( border_cond));

      // read min
      io::Serialize::Read( m_Min, ISTREAM);

      // read comma separator, or the right fence, if max = min
      ISTREAM >> current_character;
      border_cond = RangeBorders::GetConditionRightChars().find( current_character);
      if( border_cond != std::string::npos)
      {
        m_Max = m_Min;
      }
      else // read comma separator, then max
      {
        if( current_character != ',')
        {
          ERROR_STREAM << "range min and max should be separated by ',' but this was found '"
                       << current_character << "'";
          return false;
        }

        // read max
        io::Serialize::Read( m_Max, ISTREAM);

        // read right border
        ISTREAM >> current_character;
        border_cond = RangeBorders::GetConditionRightChars().find( current_character);

        if( border_cond == std::string::npos)
        {
          ERROR_STREAM << "expected either character: \"" << RangeBorders::GetConditionRightChars()
                       << "\" as right-border of range but found '" << current_character << "' actual range: "
                       << GetString();
          return false;
        }
      }

      m_Borders.Set( ConditionRight( border_cond));

      return true;
    }

    //! @brief set the value of the corresponding member based on the label
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
    template< typename t_DataType>
    bool Range< t_DataType>::TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      std::istringstream stream( LABEL.GetValue());
      return FromStream( stream, ERROR_STREAM);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief rescale a value from a source range to this range
    //! @param VALUE the value to be rescaled to this range
    //! @param SOURCE_RANGE the range the VALUE comes from
    //! @see Rescale
    //! @return rescaled value
    template< typename t_DataType>
    t_DataType Range< t_DataType>::operator()( const t_DataType &VALUE, const Range< t_DataType> &SOURCE_RANGE) const
    {
      return Rescale( VALUE, SOURCE_RANGE);
    }

    //! @brief equal operator for two ranges
    //! @param RANGE_LHS left hand side range
    //! @param RANGE_RHS right hand side range
    //! @return true if min, max and border conditions are equal
    template< typename t_DataType>
    bool Range< t_DataType>::operator ==( const Range< t_DataType> &RANGE_RHS) const
    {
      return
        m_Borders == RANGE_RHS.m_Borders &&
        m_Min     == RANGE_RHS.m_Min &&
        m_Max     == RANGE_RHS.m_Max;
    }

    //! @brief less than operator for two ranges.  Returns true only when one range is strictly less than other (no overlap)
    //! @param RANGE_RHS right hand side range
    //! @return true if the ranges have no overlap and the max LHS is < the other
    template< typename t_DataType>
    bool Range< t_DataType>::operator <( const Range< t_DataType> &RANGE_RHS) const
    {
      return m_Max < RANGE_RHS.m_Min
             ||
             (
               m_Max == RANGE_RHS.m_Min
               && ( m_Borders.IsRightBorderOpen() || RANGE_RHS.m_Borders.IsLeftBorderOpen())
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Range< t_DataType>::Read( std::istream &ISTREAM)
    {
      std::stringstream error_stream;
      // assert that there were no errors
      BCL_Assert( FromStream( ISTREAM, error_stream), "Error while reading range: " + error_stream.str());

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &Range< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // insert indent
      io::Serialize::InsertIndent( OSTREAM, INDENT);

      // write member
      OSTREAM << m_Borders.GetConditionLeftChar() << ' ';
      io::Serialize::Write( m_Min, OSTREAM, 0) << " , ";
      io::Serialize::Write( m_Max, OSTREAM, 0) << ' ';
      OSTREAM << m_Borders.GetConditionRightChar();

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes the help for the label
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent
    //! @return the given stream to which the help was written to
    template< typename t_DataType>
    std::ostream &Range< t_DataType>::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      OSTREAM << "A range of values, [] indicates closed lhs/rhs borders, () denotes open borders\n";
      io::Serialize::InsertIndent( OSTREAM, INDENT)
              << "',' separates low/high end of the border when a true range is specified\n";
      io::Serialize::InsertIndent( OSTREAM, INDENT)
              << "A single value may also be given e.g. [7]";
      return OSTREAM;
    }

    //! @brief standardize range for integral types to a left closed and right open range
    //! @return a range [new_min,new_max)
    template< typename t_DataType>
    Range< t_DataType> Range< t_DataType>::StandardizeRange() const
    {
      // ensure that the range is for integer types
      BCL_Assert
      (
        std::numeric_limits< t_DataType>::is_integer,
        "cannot standardize non-integral range"
      );

      // ensure that the max value will not wrap around
      BCL_Assert
      (
        m_Max != std::numeric_limits< t_DataType>::max() || m_Borders.IsRightBorderOpen(),
        "cannot standardize closed range with max == max value for " + GetStaticClassName< t_DataType>()
      );

      return
        Range< t_DataType>
        (
          RangeBorders::e_LeftClosed,
          m_Min + t_DataType( m_Borders.IsLeftBorderOpen()),
          m_Max + t_DataType( m_Borders.IsRightBorderClosed()),
          RangeBorders::e_RightOpen
        );
    }

    //! @brief close any open borders, while maintaining the same set of values for which IsWithin returns true
    //! @return a new range with closed borders
    //! @note this method should be preferred over StandardizeRange because standardize range does not work for
    //! @note floating point types, and in any event can cause m_Max to wrap around, e.g. if the range included
    //! @note the entire range of size_t's, then standardizing the range would make the range empty
    template< typename t_DataType>
    Range< t_DataType> Range< t_DataType>::CloseBorders() const
    {
      // handle the trivial case that both borders are already closed
      if( m_Borders.IsClosed())
      {
        return *this;
      }

      // handle the trivial case of an empty range
      if( IsEmpty())
      {
        return Range< t_DataType>( t_DataType( 0.0), t_DataType( 0.0));
      }

      // copy the current range, with the borders closed
      Range< t_DataType> new_range( m_Min, m_Max);

      // handle integer types
      if( std::numeric_limits< t_DataType>::is_integer)
      {
        if( m_Borders.IsLeftBorderOpen())
        {
          new_range.m_Min += t_DataType( 1);
        }
        if( m_Borders.IsRightBorderOpen())
        {
          new_range.m_Max -= t_DataType( 1);
        }
      }
      else // handle floating point types
      {
        // get the machine epsilon
        const t_DataType epsilon( std::numeric_limits< t_DataType>::epsilon());

        // get the smallest absolute value representable for t_DataType
        const t_DataType abs_min( std::numeric_limits< t_DataType>::min());

        // handle the left hand side
        if( m_Borders.IsLeftBorderOpen() && util::IsDefined( m_Min))
        {
          // handle values > 0.0
          if( m_Min > t_DataType( 0.0)) // positive values
          {
            // the first x for which x > m_Min is (1 + epsilon) * m_Min
            // note: for numeric reasons, it is NOT okay to use += here, e.g. epsilon * abs_min is too small to represent
            new_range.m_Min = ( 1.0 + epsilon) * m_Min;
          }
          else if( m_Min < t_DataType( 0.0)) // negative values
          {
            // note: for numeric reasons, it is NOT okay to use -= here, e.g. epsilon * abs_min is too small to represent
            new_range.m_Min = m_Min / ( 1.0 + epsilon);
          }
          else // if m_Min == 0.0
          {
            // abs_min is the first number past 0 that compares greater than 0
            new_range.m_Min = abs_min;
          }
        }

        // handle the right hand side
        if( m_Borders.IsRightBorderOpen() && util::IsDefined( m_Max))
        {
          // handle values > 0.0
          if( m_Max > t_DataType( 0.0)) // positive values
          {
            // the first x for which x > m_Max is (1 + epsilon) * m_Max
            // note: for numeric reasons, it is NOT okay to use += here, e.g. epsilon * abs_min is too small to represent
            new_range.m_Max = ( 1.0 + epsilon) * m_Max;
          }
          else if( m_Max < t_DataType( 0.0)) // negative values
          {
            // note: for numeric reasons, it is NOT okay to use -= here, e.g. epsilon * abs_min is too small to represent
            new_range.m_Max = m_Max / ( 1.0 + epsilon);
          }
          else // if m_Min == 0.0
          {
            // abs_min is the first number past 0 that compares greater than 0
            new_range.m_Max = abs_min;
          }
        }
      }

      return new_range;
    }

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RANGE_HPP_
