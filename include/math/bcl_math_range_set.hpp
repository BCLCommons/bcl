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

#ifndef BCL_MATH_RANGE_SET_HPP_
#define BCL_MATH_RANGE_SET_HPP_

// include the header of this class
#include "bcl_math_range_set.h"

// includes from bcl - sorted alphabetically
#include "bcl_math_limits.h"
#include "io/bcl_io_serialize.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_instances.h"
#include "util/bcl_util_string_functions.h"

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
    const util::SiPtr< const util::ObjectInterface> RangeSet< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new RangeSet< t_DataType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor - initialize min and max to undefined and the borders to a closed interval
    template< typename t_DataType>
    RangeSet< t_DataType>::RangeSet() :
      m_Ranges()
    {
    }

    //! @brief construct from a single range
    //! @param RANGE range to use to initialize the range set
    template< typename t_DataType>
    RangeSet< t_DataType>::RangeSet( const Range< t_DataType> &RANGE)
    {
      if( !RANGE.IsEmpty())
      {
        // add the range to the range set
        m_Ranges.Insert( RANGE);
      }
    }

    //! @brief Clone function
    //! @return pointer to new RangeSet
    template< typename t_DataType>
    RangeSet< t_DataType> *RangeSet< t_DataType>::Clone() const
    {
      return new RangeSet< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &RangeSet< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the label containing only initialization parameters
    //! @param WITH_DATA whether to include any data members, else, only include initialization members
    template< typename t_DataType>
    util::ObjectDataLabel RangeSet< t_DataType>::GetLabel( const bool &WITH_DATA) const
    {
      return util::ObjectDataLabel( "", AsString());
    }

    //! @brief get the minimum value for which IsWithin will return true
    //! @return the minimum value of the range set for which IsWithin will return true
    template< typename t_DataType>
    t_DataType RangeSet< t_DataType>::GetMin() const
    {
      if( IsEmpty())
      {
        return util::GetUndefined< t_DataType>();
      }
      return m_Ranges.Begin()->GetMin();
    }

    //! @brief get the maximum value for which IsWithin will return true
    //! @return the maximum value of the range set for which IsWithin will return true
    template< typename t_DataType>
    t_DataType RangeSet< t_DataType>::GetMax() const
    {
      if( IsEmpty())
      {
        return util::GetUndefined< t_DataType>();
      }
      return m_Ranges.ReverseBegin()->GetMax();
    }

    //! @brief check value for being in this RangeSet
    //! @param VALUE
    //! @return true, if within RangeSet [min,max]
    template< typename t_DataType>
    bool RangeSet< t_DataType>::IsWithin( const t_DataType &VALUE) const
    {
      // make a copy of the new range, which will be expanded if it overlaps any other ranges in the set
      return util::IsDefined( VALUE) && m_Ranges.Contains( Range< t_DataType>( VALUE, VALUE));
    }

    //! @brief check whether the RangeSet contains any ranges
    //! @return true if there are no values in the RangeSet, or it is undefined
    template< typename t_DataType>
    bool RangeSet< t_DataType>::IsEmpty() const
    {
      return m_Ranges.IsEmpty();
    }

    //! @brief read from a string
    //! @param STRING input string
    //! @param ERROR_STREAM stream to print errors to
    template< typename t_DataType>
    bool RangeSet< t_DataType>::FromString( const std::string &STRING, std::ostream &ERROR_STREAM)
    {
      // reset any existing ranges
      m_Ranges.Reset();

      // put the trimmed string into an input string stream
      std::string trimmed_string( util::TrimString( STRING));

      // handle - as the first character, which means the user wants the full range of the numeric type minus the specified range
      if( trimmed_string.size() > 1 && trimmed_string[ 0] == '-')
      {
        m_Ranges.Insert
        (
          Range< t_DataType>
          (
            RangeBorders::e_LeftClosed,
            GetLowestUnboundedValue< t_DataType>(),
            GetHighestUnboundedValue< t_DataType>(),
            RangeBorders::e_RightClosed
          )
        );
      }
      // '+' is prepended to keep the while loop below simple, because then there is always a sign before each range
      // only prepend the '+' if the string does not already begin with + and if there really is a range
      else if( trimmed_string.size() > 1 && trimmed_string[ 0] != '+')
      {
        trimmed_string = "+" + trimmed_string;
      }

      std::istringstream range_string_stream( trimmed_string);

      // determine the length of the string stream
      const std::istringstream::streampos string_length( trimmed_string.size());

      // while there are characters left in the stream, attempt to read another range
      while( range_string_stream.tellg() < string_length)
      {
        // load the + or - character
        char sign;
        range_string_stream >> sign;
        if( sign != '+' && sign != '-')
        {
          // bad sign character
          ERROR_STREAM << "ranges must be added/subtracted with +/-, but encountered " << sign;
          return false;
        }

        // initialize a new range
        Range< t_DataType> temp;

        // read the range from the stream
        if( !temp.FromStream( range_string_stream, ERROR_STREAM))
        {
          return false;
        }

        if( sign == '+')
        {
          // add the range to the set
          *this += temp;
        }
        else // if sign == '-'
        {
          // remove this range from the set
          *this -= temp;
        }
      }
      return true;
    }

    //! @brief Get the mapped subset of values chosen from this rangeset
    //! @param SELECTION the selected range of this rangeset; min must be >= 0, max <= total width of internal ranges
    //! @return SELECTION of this RangeSet
    template< typename t_DataType>
    RangeSet< t_DataType> RangeSet< t_DataType>::GetMappedSubset( const Range< t_DataType> &SELECTION)
    {
      BCL_Assert( SELECTION.GetMin() >= t_DataType( 0), "Selection must be positive for the mapped subset");
      BCL_Assert( std::numeric_limits< t_DataType>::is_integer, "MappedSubset is currently only implemented for integral types");
      storage::Vector< Range< t_DataType> > selected_ranges;
      selected_ranges.AllocateMemory( m_Ranges.GetSize());
      Range< t_DataType> closed_range( SELECTION.CloseBorders());

      t_DataType amount_till_start( closed_range.GetMin());

      // for integral types, the total width is really 1 more than GetWidth returns
      t_DataType remaining_width( closed_range.GetWidth() + 1);
      for
      (
        typename storage::Set< Range< t_DataType> >::const_iterator itr( m_Ranges.Begin()), itr_end( m_Ranges.End());
        itr != itr_end && remaining_width;
        ++itr
      )
      {
        t_DataType effective_width( itr->GetWidth() + 1);
        if( amount_till_start >= effective_width)
        {
          amount_till_start -= effective_width;
          continue;
        }
        if( amount_till_start)
        {
          effective_width -= amount_till_start;
        }
        effective_width = std::min( remaining_width, effective_width);
        selected_ranges.PushBack
        (
          Range< t_DataType>
          (
            RangeBorders::e_LeftClosed,
            itr->GetMin() + amount_till_start,
            itr->GetMin() + amount_till_start + effective_width,
            RangeBorders::e_RightOpen
          )
        );
        amount_till_start = 0;
        remaining_width -= effective_width;
      }
      return RangeSet< t_DataType>( selected_ranges.Begin(), selected_ranges.End());
    }

    //! @brief add a range to the range set
    //! @param RANGE a range to be added to the range set
    //! @return this
    template< typename t_DataType>
    RangeSet< t_DataType> &RangeSet< t_DataType>::operator +=( const Range< t_DataType> &RANGE)
    {
      // make a copy of the new range, which will be expanded if it overlaps any other ranges in the set
      Range< t_DataType> temp( RANGE.CloseBorders());

      // an empty range does not change the range set, so only add the range if it is non-empty
      if( !temp.IsEmpty())
      {
        // check whether this range overlaps with any other range in the set
        typename storage::Set< Range< t_DataType> >::const_iterator
          itr( m_Ranges.Find( temp));

        // until this range does not overlap with any other range in the set
        while( itr != m_Ranges.End())
        {
          // combine this range with the overlapping range in the set
          temp = Range< t_DataType>::CombineRanges( temp, *itr);

          // remove the old range from the set
          m_Ranges.RemoveElement( itr);

          // look to see if this range overlaps with any others in the set
          itr = m_Ranges.Find( temp);
        }

        m_Ranges.Insert( temp);
      }

      return *this;
    }

    //! @brief remove a range from the set
    //! @param RANGE a range of values which should be excluded from the set
    //! @return this
    template< typename t_DataType>
    RangeSet< t_DataType> &RangeSet< t_DataType>::operator -=( const Range< t_DataType> &RANGE)
    {
      // an empty range does not change the range set, so only remove the range if it is non-empty
      if( !RANGE.IsEmpty())
      {
        // look to see if this range overlaps with any other range in the set
        typename storage::Set< Range< t_DataType> >::const_iterator
          itr( m_Ranges.Find( RANGE));

        storage::Vector< Range< t_DataType> > newly_split_ranges;
        // until this range does not overlap with any other range in the set
        while( itr != m_Ranges.End())
        {
          // remove RANGE from itr
          storage::Vector< Range< t_DataType> > split_ranges( SplitRanges( *itr, RANGE));

          // remove the old range from the set
          m_Ranges.RemoveElement( itr);

          // add the ranges back into the set, no need to go through operator += since
          // the ranges returned will not overlap or be empty
          newly_split_ranges.Append( split_ranges);

          // look to see if this range overlaps with any others in the set
          itr = m_Ranges.Find( RANGE);
        }
        m_Ranges.InsertElements( newly_split_ranges.Begin(), newly_split_ranges.End());
      }
      return *this;
    }

    //! @brief compare range sets
    //! @param B a range set of values to compare this range set against
    //! @return true if all internal ranges are the same
    template< typename t_DataType>
    bool RangeSet< t_DataType>::operator ==( const RangeSet< t_DataType> &B) const
    {
      return m_Ranges == B.m_Ranges;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get a range set encompassing the entire range of the data type
    //! @return a range set encompassing the entire range of the data type
    template< typename t_DataType>
    RangeSet< t_DataType> RangeSet< t_DataType>::GetCompleteRange()
    {
      // for integral types, numeric limits min/max give the minimum and maximal values, respectively
      // for floating point values, the min is the negative of the max value, the value returned by min
      // is the smallest positive value the type can represent
      return
        RangeSet< t_DataType>
        (
          Range< t_DataType>
          (
            GetLowestBoundedValue< t_DataType>(),
            GetHighestBoundedValue< t_DataType>()
          )
        );
    }

    //! @brief writes the help for the label
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent
    //! @return the given stream to which the help was written to
    template< typename t_DataType>
    std::ostream &RangeSet< t_DataType>::WriteHelp
    (
      std::ostream &OSTREAM,
      const size_t INDENT
    ) const
    {
      OSTREAM << "A set of ranges, which can be written using addition (+) and subtraction (-) of range from one\n";
      io::Serialize::InsertIndent( OSTREAM, INDENT) << "another, e.g. [0,10) - [5,8) or [0,5) + [8,10)";
      return OSTREAM;
    }

    //! @brief write this range set as a string
    //! @return this range set as a string
    template< typename t_DataType>
    std::string RangeSet< t_DataType>::AsString() const
    {
      if( m_Ranges.IsEmpty())
      {
        return std::string();
      }
      std::string ranges_str( m_Ranges.Begin()->GetString());
      for
      (
        typename storage::Set< Range< t_DataType> >::const_iterator
          itr( ++m_Ranges.Begin()), itr_end( m_Ranges.End());
        itr != itr_end;
        ++itr
      )
      {
        ranges_str += "+" + itr->GetString();
      }
      return ranges_str;
    }

    //! @brief set the value of the corresponding member based on the label
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
    template< typename t_DataType>
    bool RangeSet< t_DataType>::TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      return FromString( LABEL.GetValue(), ERROR_STREAM);
    }

    //! @brief helper function to remove any overlap of one range from another
    //! @param ORIGINAL_RANGE the original range
    //! @param RANGE_TO_EXCLUDE the part of the original range to exclude
    //! @return a vector, containing 0-2 ranges that are the segments of the ORIGINAL_RANGE
    //!         that were not in the RANGE_TO_EXCLUDE
    template< typename t_DataType>
    storage::Vector< Range< t_DataType> > RangeSet< t_DataType>::SplitRanges
    (
      const Range< t_DataType> &ORIGINAL_RANGE,
      const Range< t_DataType> &RANGE_TO_EXCLUDE
    )
    {
      storage::Vector< Range< t_DataType> > split_range;

      // make versions of the ranges with closed borders
      Range< t_DataType> original_range( ORIGINAL_RANGE.CloseBorders());
      Range< t_DataType> range_to_exclude( RANGE_TO_EXCLUDE.CloseBorders());

      // handle the trivial case in which EXCLUDED_RANGE contains no part of the original range or is empty
      if
      (
        range_to_exclude.GetMax() < original_range.GetMin()
        || range_to_exclude.GetMin() > original_range.GetMax()
        || range_to_exclude.IsEmpty()
      )
      {
        // the split range is just the original range
        split_range.PushBack( original_range);
        // as it is used in this class, this branch will never be called, because only overlapping ranges are
        // removed
      }
      else // at least part of range_to_exclude overlaps the original range
      {
        // add any part of the original range before the excluded range
        if( range_to_exclude.GetMin() > original_range.GetMin())
        {
          split_range.PushBack
          (
            Range< t_DataType>
            (
              RangeBorders::e_LeftClosed,
              original_range.GetMin(),
              range_to_exclude.GetMin(),
              RangeBorders::e_RightOpen
            ).CloseBorders()
          );
        }

        // add any part of the original range that extends beyond the max of the excluded range
        if( range_to_exclude.GetMax() < original_range.GetMax())
        {
          split_range.PushBack
          (
            Range< t_DataType>
            (
              RangeBorders::e_LeftOpen,
              range_to_exclude.GetMax(),
              original_range.GetMax(),
              RangeBorders::e_RightClosed
            ).CloseBorders()
          );
        }
      }

      // return the range(s) that remain
      return split_range;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &RangeSet< t_DataType>::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Ranges, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &RangeSet< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Ranges, OSTREAM, INDENT);
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RANGE_SET_HPP_
