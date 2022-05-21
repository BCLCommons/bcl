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

#ifndef BCL_MATH_RANGE_SET_H_
#define BCL_MATH_RANGE_SET_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_range.h"
#include "io/bcl_io_serialization_interface.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_data_type.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RangeSet
    //! @brief a set of ranges, which provides much of the same behavior as Range
    //!
    //! @tparam t_DataType the type of data of that RangeSet. It requires comparison <=, <, >, >= operators, subtraction,
    //! scalar multiplication, summation and division for the t_DataType.
    //!
    //! @see @link example_math_range_set.cpp @endlink
    //! @author mendenjl
    //! @date Jan 31, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class RangeSet :
      public io::SerializationInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Ranges included in this range set
      storage::Set< Range< t_DataType> > m_Ranges;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor - no ranges
      RangeSet();

      //! @brief construct from a single range
      //! @param RANGE range to use to initialize the range set
      RangeSet( const Range< t_DataType> &RANGE);

      //! @brief construct from an iterator over a range of ranges
      template< typename t_InputIterator>
      RangeSet( const t_InputIterator &BEGIN, const t_InputIterator &END)
      {
        for( t_InputIterator itr( BEGIN); itr != END; ++itr)
        {
          *this += *itr;
        }
      }

      //! @brief Clone function
      //! @return pointer to new RangeSet
      RangeSet< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the label containing only initialization parameters
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the ranges in this range set
      //! @return the ranges in this range set
      const storage::Set< Range< t_DataType> > &GetRanges() const
      {
        return m_Ranges;
      }

      //! @brief get the minimum value for which IsWithin will return true
      //! @return the minimum value of the range set for which IsWithin will return true
      t_DataType GetMin() const;

      //! @brief get the maximum value for which IsWithin will return true
      //! @return the maximum value of the range set for which IsWithin will return true
      t_DataType GetMax() const;

      //! @brief check value for being in this RangeSet
      //! @param VALUE
      //! @return true, if within RangeSet [min,max]
      bool IsWithin( const t_DataType &VALUE) const;

      //! @brief check whether the RangeSet contains any ranges
      //! @return true if there are no values in the RangeSet, or it is undefined
      bool IsEmpty() const;

      //! @brief read from a string
      //! @param STRING input string
      //! @param ERROR_STREAM stream to print errors to
      bool FromString( const std::string &STRING, std::ostream &ERROR_STREAM);

      //! @brief Get the mapped subset of values chosen from this rangeset
      //! @param SELECTION the selected range of this rangeset; min must be >= 0, max <= total width of internal ranges
      //! @return SELECTION of this RangeSet
      RangeSet< t_DataType> GetMappedSubset( const Range< t_DataType> &SELECTION);

      //! @brief add a range to the range set
      //! @param RANGE a range to be added to the range set
      //! @return this
      RangeSet< t_DataType> &operator +=( const Range< t_DataType> &RANGE);

      //! @brief remove a range from the set
      //! @param RANGE a range of values which should be excluded from the set
      //! @return this
      RangeSet< t_DataType> &operator -=( const Range< t_DataType> &RANGE);

      //! @brief compare range sets
      //! @param B a range set of values to compare this range set against
      //! @return true if all internal ranges are the same
      bool operator ==( const RangeSet< t_DataType> &B) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a range set encompassing the entire range of the data type
      //! @return a range set encompassing the entire range of the data type
      static RangeSet< t_DataType> GetCompleteRange();

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

      //! @brief write this range set as a string
      //! @return this range set as a string
      std::string AsString() const;

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    protected:

      //! @brief helper function to remove any overlap of one range from another
      //! @param ORIGINAL_RANGE the original range
      //! @param RANGE_TO_EXCLUDE the part of the original range to exclude
      //! @return a vector, containing 0-2 ranges that are the segments of the ORIGINAL_RANGE
      //!         that were not in the RANGE_TO_EXCLUDE
      static storage::Vector< Range< t_DataType> > SplitRanges
      (
        const Range< t_DataType> &ORIGINAL_RANGE,
        const Range< t_DataType> &RANGE_TO_EXCLUDE
      );

    }; // template class RangeSet

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    // when adding types to this list, they also need to be added to math::Comparisons
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< bool>;
    BCL_EXPIMP_TEMPLATE template class BCL_API RangeSet< char>;

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RANGE_SET_H_
