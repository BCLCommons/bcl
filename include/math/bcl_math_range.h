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

#ifndef BCL_MATH_RANGE_H_
#define BCL_MATH_RANGE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization_interface.h"
#include "util/bcl_util_binary_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RangeBorders
    //! @brief is the base class for Ranges, containing the border condition for open and closed intervals
    //!
    //! @see @link example_math_range.cpp @endlink
    //! @author woetzen
    //! @date Jul 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RangeBorders
    {
    public:

    ///////////
    // types //
    ///////////

      //! enum that defines if the interval is open or closed on the left side (includes or excludes the min)
      enum ConditionLeft
      {
        e_LeftClosed, //!< closed to the left [MIN, MAX .. includes MIN
        e_LeftOpen    //!< open   to the left (MIN, MAX .. excludes MIN
      };

      //! enum that defines if the interval is open or closed on the right side (includes or excludes the max)
      enum ConditionRight
      {
        e_RightClosed, //!< closed to the left .. MIN, MAX] includes MAX
        e_RightOpen    //!< open   to the left .. MIN, MAX) excludes MAX
      };

      static const std::string &GetConditionLeftChars();  //!< character used for left border condition
      static const std::string &GetConditionRightChars(); //!< character used for right border condition

    private:

    //////////
    // data //
    //////////

      ConditionLeft  m_ConditionLeft;  //!< left border condition
      ConditionRight m_ConditionRight; //!< right border condition

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, sets both borders to closed
      //! @param LEFT_BORDER_CONDITION border condition to the left
      //! @param RIGHT_BORDER_CONDITION border condition to the right
      RangeBorders() :
        m_ConditionLeft( e_LeftClosed),
        m_ConditionRight( e_RightClosed)
      {
      }

      //! @brief construct from two border conditions
      //! @param LEFT_BORDER_CONDITION border condition to the left
      //! @param RIGHT_BORDER_CONDITION border condition to the right
      RangeBorders
      (
        const ConditionLeft LEFT_BORDER_CONDITION,
        const ConditionRight RIGHT_BORDER_CONDITION
      ) :
        m_ConditionLeft( LEFT_BORDER_CONDITION),
        m_ConditionRight( RIGHT_BORDER_CONDITION)
      {
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

      //! @brief access to left border condition
      //! @return ConditionLeft open or closed
      ConditionLeft GetLeftCondition() const
      {
        return m_ConditionLeft;
      }

      //! @brief access to right border condition
      //! @return ConditionRight open or closed
      ConditionRight GetRightCondition() const
      {
        return m_ConditionRight;
      }

      //! @brief check if left border is closed
      //! @return true if the left border is closed
      bool IsLeftBorderClosed() const
      {
        return m_ConditionLeft == e_LeftClosed;
      }

      //! @brief check if left border is open
      //! @return true if the left border is open
      bool IsLeftBorderOpen() const
      {
        return m_ConditionLeft == e_LeftOpen;
      }

      //! @brief set the left border
      //! @param LEFT the new border condition
      void Set( const ConditionLeft &LEFT)
      {
        m_ConditionLeft = LEFT;
      }

      //! @brief set the right border
      //! @param RIGHT the new border condition
      void Set( const ConditionRight &RIGHT)
      {
        m_ConditionRight = RIGHT;
      }

      //! @brief check if right border is closed
      //! @return true if the right border is closed
      bool IsRightBorderClosed() const
      {
        return m_ConditionRight == e_RightClosed;
      }

      //! @brief check if right border is open
      //! @return true if the right border is open
      bool IsRightBorderOpen() const
      {
        return m_ConditionRight == e_RightOpen;
      }

      //! @brief check whether both borders are closed
      //! @return true if both borders are closed
      bool IsClosed() const
      {
        return m_ConditionRight == e_RightClosed && m_ConditionLeft == e_LeftClosed;
      }

      //! @brief check whether both borders are open
      //! @return true if both borders are open
      bool IsOpen() const
      {
        return m_ConditionRight == e_RightOpen && m_ConditionLeft == e_LeftOpen;
      }

      //! @brief get the character for the left border condition
      const char &GetConditionLeftChar() const
      {
        return GetConditionLeftChars()[ m_ConditionLeft];
      }

      //! @brief get the character for the right border condition
      const char &GetConditionRightChar() const
      {
        return GetConditionRightChars()[ m_ConditionRight];
      }

      //! @brief comparison
      bool operator ==( const RangeBorders &OTHER) const
      {
        return m_ConditionLeft == OTHER.m_ConditionLeft && m_ConditionRight == OTHER.m_ConditionRight;
      }

    }; // class RangeBorders

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Range
    //! @brief represents an interval [MIN, MAX], (MIN, MAX], [MIN, MAX) or (MIN, MAX)
    //! @details it provides functions to check if a number is within an interval, maps a number from one interval to another
    //! interval
    //!
    //! @tparam t_DataType the type of data of that range. It requires comparison <=, <, >, >= operators, subtraction,
    //! scalar multiplication, summation and division for the t_DataType.
    //!
    //! @see @link example_math_range.cpp @endlink
    //! @author woetzen, mueller
    //! @date Apr 8, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Range :
      virtual public util::BinaryFunctionInterface< t_DataType, Range< t_DataType>, t_DataType>,
      virtual public io::SerializationInterface
    {

    /////////////
    // friends //
    /////////////

      //! @brief Range class is friend
      template< typename t_OtherDataType> friend class Range;

    private:

    //////////
    // data //
    //////////

      RangeBorders  m_Borders;  //!< Borders object
      t_DataType    m_Min;      //!< min of range
      t_DataType    m_Max;      //!< max of range

    public:

    //////////
    // data //
    //////////

      // import range border enums
      typedef RangeBorders::ConditionLeft ConditionLeft;
      typedef RangeBorders::ConditionRight ConditionRight;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor - initialize min and max to undefined and the borders to a closed interval
      Range();

      //! @brief construct from min and max and two border conditions
      //! @param LEFT_BORDER_CONDITION border condition to the left
      //! @param MIN lower (left-hand) bound of closed interval
      //! @param MAX upper (right-hand) bound of closed interval
      //! @param RIGHT_BORDER_CONDITION border condition to the right
      Range
      (
        const ConditionLeft LEFT_BORDER_CONDITION,
        const t_DataType    &MIN,
        const t_DataType    &MAX,
        const ConditionRight RIGHT_BORDER_CONDITION
      );

      //! @brief construct from min and max [min,max] in a closed interval
      //! @param MIN lower (left-hand) bound of closed interval
      //! @param MAX upper (right-hand) bound of closed interval
      Range
      (
        const t_DataType &MIN,
        const t_DataType &MAX
      );

      //! @brief Clone function
      //! @return pointer to new Range
      Range< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief access to min of range
      //! @return min
      const t_DataType &GetMin() const;

      //! @brief access to max of range
      //! @return max
      const t_DataType &GetMax() const;

      //! @brief access to left border condition
      //! @return ConditionLeft open or closed
      ConditionLeft GetLeftCondition() const;

      //! @brief access to right border condition
      //! @return ConditionRight open or closed
      ConditionRight GetRightCondition() const;

      //! @brief return the width
      //! @return max-min
      t_DataType GetWidth() const;

      //! @brief return the middle value of the interval
      //! @return the middle of the interval
      t_DataType GetMiddle() const;

      //! @brief convert to a human readable string with formatted numbers
      //! @return string like "[2.5,20.1]"
      std::string GetString() const;

      //! @brief convert to a human readable string with formatted numbers
      //! @param FORMAT format for the numbers
      //! @return string like "[2.5,20.1]"
      std::string GetString( const util::Format &FORMAT) const;

      //! @brief get the label containing only initialization parameters
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const;

      //! @brief change the min by VALUE
      //! @param MIN the new minimum value
      //! @return true, if range could be changed (range is valid after operation)
      bool SetMin( const t_DataType &MIN);

      //! @brief change the max by VALUE
      //! @param MAX the new maximum value
      //! @return true, if range could be changed (range is valid after operation)
      bool SetMax( const t_DataType &MAX);

    ////////////////
    // operations //
    ////////////////

      //! @brief rescale a value from a source range to this range
      //! @param VALUE the value to be rescaled to this range
      //! @param SOURCE_RANGE the range the VALUE comes from
      //! @return rescaled value
      t_DataType Rescale( const t_DataType &VALUE, const Range< t_DataType> &SOURCE_RANGE) const;

      //! @brief check value for being in this range
      //! @param VALUE
      //! @return true, if within range [min,max]
      bool IsWithin( const t_DataType &VALUE) const;

      //! @brief check whether the range contains any values
      //! @return true if there are no values in the range, or it is undefined
      bool IsEmpty() const;

      //! @brief check if one range falls within another range
      //! @param VALUE
      //! @return true, if there is overlap within a range
      bool DoesOverlap( const Range< t_DataType> &VALUE) const;

      //! @brief combine overlapping ranges into a single range
      //! @param RANGE_1 the first range to combine
      //! @param RANGE_2 the second range to combine
      //! @return the overlapping ranges
      static Range< t_DataType> CombineRanges
      (
        const Range< t_DataType> &RANGE_1,
        const Range< t_DataType> &RANGE_2
      );

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @param ERROR_STREAM stream to write errors to
      //! @return true on success
      bool FromStream( std::istream &ISTREAM, std::ostream &ERROR_STREAM);

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    ///////////////
    // operators //
    ///////////////

      //! @brief rescale a value from a source range to this range
      //! @param VALUE the value to be rescaled to this range
      //! @param SOURCE_RANGE the range the VALUE comes from
      //! @see Rescale
      //! @return rescaled value
      t_DataType operator()( const t_DataType &VALUE, const Range< t_DataType> &SOURCE_RANGE) const;

      //! @brief equal operator for two ranges
      //! @param RANGE_RHS right hand side range
      //! @return true if min, max and border conditions are equal
      bool operator ==( const Range< t_DataType> &RANGE_RHS) const;

      //! @brief less than operator for two ranges.  Returns true only when one range is strictly less than other (no overlap)
      //! @param RANGE_RHS right hand side range
      //! @return true if the ranges have no overlap and the max LHS is <= the other
      bool operator <( const Range< t_DataType> &RANGE_RHS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

      //! @brief convert range of one datatype to a range of another datatype
      //! @tparam t_OtherDataType datatype of new range to convert to
      //! @return range for other datatype, with casted min and max
      template< typename t_OtherDataType>
      Range< t_OtherDataType> Convert() const
      {
        return Range< t_OtherDataType>
        (
          GetLeftCondition(),
          t_OtherDataType( m_Min),
          t_OtherDataType( m_Max),
          GetRightCondition()
        );
      }

      //! @brief standardize range for integral types to a left closed and right open range
      //! @return a range [new_min,new_max)
      Range< t_DataType> StandardizeRange() const;

      //! @brief close any open borders, while maintaining the same set of values for which IsWithin returns true
      //! @return a new range with closed borders
      //! @note this method should be preferred over StandardizeRange because this function works for floating point types
      //! @note and always returns a valid range (StandardizeRange does not work if the range included
      //! @note the entire range of size_t's, then standardizing the range would make the range empty)
      Range< t_DataType> CloseBorders() const;

    }; // template class Range

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    // when adding types to this list, they also need to be added to math::Comparisons
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< bool>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Range< char>;

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RANGE_H_
