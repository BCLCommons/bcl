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

#ifndef BCL_MATH_RUNNING_AVERAGE_H_
#define BCL_MATH_RUNNING_AVERAGE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "type/bcl_type_compare.h"
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
    //! @class RunningAverage
    //! @brief class maintains a running average of observations
    //! Running average does not store observations, so it is ideal for
    //! averaging vectors or large sets of numbers that are generated if only the average is needed.
    //! RunningAverage is also more numerically stable than computing the sum and then dividing by the # of observations
    //!
    //! @tparam t_AverageType the type in which to store the average
    //! t_AverageType for non-integral types should be the same as t_ArgumentType
    //! t_AverageType for integer types should be double, to avoid roundoff
    //!
    //! @see @link example_math_running_average.cpp @endlink
    //! @author mendenjl
    //! @date Dec 3, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_AverageType>
    class RunningAverage :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_AverageType m_Average; //!< The current average
      double        m_Weight;  //!< Sum of observation weights, = the # of observations if observations are unweighted

    public:
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RunningAverage() :
        m_Average(),
        m_Weight( 0.0)
      {
        Zero( m_Average);
      }

      //! @brief constructor from an initial average, needed so that m_Average has the right size if t_ArgumentType is a container
      //! @param INITIAL_AVERAGE initial average
      //! If t_ArgumentType is a vector, INITIAL_AVERAGE must be the same size as all vectors that will be averaged
      explicit RunningAverage( const t_AverageType &INITIAL_AVERAGE) :
        m_Average( INITIAL_AVERAGE),
        m_Weight( 0.0)
      {
      }

      //! virtual copy constructor
      RunningAverage *Clone() const
      {
        return new RunningAverage( *this);
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
        const t_AverageType &( RunningAverage< t_AverageType>::*FUNCTION)() const
      ) const
      {
        static const std::string s_alias( "Mean");
        return s_alias;
      }

      //! @brief get the current weight
      //! @return the current weight
      const double &GetWeight() const
      {
        return m_Weight;
      }

      //! @brief change the current weight
      //! @param WEIGHT the new weight to use
      void SetWeight( const double &WEIGHT)
      {
        m_Weight = WEIGHT;
      }

      //! @brief get the current weight
      //! @return the current weight
      const t_AverageType &GetAverage() const
      {
        return m_Average;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset the running average
      void Reset()
      {
        m_Weight = 0.0;
        Zero( m_Average);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief convert the running average to the average type
      //! @return the current average
      operator const t_AverageType &() const
      {
        return m_Average;
      }

      //! @brief add an observation from the running average
      //! @param VALUE the value to add to the running average
      const RunningAverage &operator +=( const RunningAverage &VALUE)
      {
        return AddWeightedObservation( VALUE.GetAverage(), VALUE.GetWeight());
      }

      //! @brief add an observation from the running average
      //! @param VALUE the value to add to the running average
      template< typename t_ArgumentType>
      const RunningAverage &operator +=( const t_ArgumentType &VALUE)
      {
        return AddWeightedObservation( VALUE, 1.0);
      }

      //! @brief remove an observation from the running average
      //! @param VALUE the value to add to the running average
      template< typename t_ArgumentType>
      const RunningAverage &operator -=( const t_ArgumentType &VALUE)
      {
        return AddWeightedObservation( VALUE, -1.0);
      }

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param VALUE the value to add to the running average
      //! @param WEIGHT the weight to give the observation (may be negative to remove an old observation)
      //! @param the weight to add to
      template< typename t_ArgumentType>
      const RunningAverage &AddWeightedObservation( const t_ArgumentType &VALUE, double WEIGHT)
      {
        if( WEIGHT == -m_Weight && WEIGHT)
        {
          // if the weight is now zero, just make the average 0 too and return
          Reset();
          return *this;
        }
        else if( m_Weight == 0.0)
        {
          // if the current weight is 0.0, then just set m_Weight to WEIGHT and the average is just VALUE
          m_Average = VALUE;
          m_Weight = WEIGHT;
          return *this;
        }

        // add the new weight to m_WEIGHT
        m_Weight += WEIGHT;

        // get the ratio of this observations weight to the new weight = WEIGHT / (d + WEIGHT)
        // update the average, using the weight ratio
        UpdateAverage( m_Average, VALUE, WEIGHT / m_Weight);

        return *this;
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        io::Serialize::Read( m_Average, ISTREAM);
        io::Serialize::Read( m_Weight, ISTREAM);
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::Write( m_Average, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Weight, OSTREAM, INDENT);
        return OSTREAM;
      }

      //! @brief Update the given average with the new value according to the given weight ratio
      //! @param SCALAR_AVE the scalar average to update
      //! @param VALUE the new value to update the scalar average with
      //! @param WEIGHT_RATIO the ratio of weights (weight for NEW_VALUE) / (total new weight)
      //! @details overloaded version for t_DataTypes without iterators
      template< typename t_DataType, typename t_OtherDataType>
      static
      typename type::EnableIf
      <
        !type::IsSequence< t_DataType>::value && !type::Compare< t_OtherDataType, t_DataType>::e_Same,
        void
      >::Type
      UpdateAverage( t_DataType &SCALAR_AVE, const t_OtherDataType &VALUE, const double &WEIGHT_RATIO)
      {
        // a weighted average of a vector is always just
        // x_ave = sum(x[i]*weight_x[i]) / sum(weight_x[i]) for 0 <= i < x.GetSize()
        // abbreviating the numerator with n and the denominator with d:
        // x_ave = n / d
        // adding the new observation VALUE with weight WEIGHT gives:
        // x_ave_new = (sum(x[i]*weight_x[i]) + VALUE * WEIGHT) / ( sum(weight_x[i])) + WEIGHT)
        //           = (n + VALUE * WEIGHT) / (d + WEIGHT)
        // taking x_ave_new - x_ave gives us
        // x_ave_new - x_ave = (n + VALUE * WEIGHT) / (d + WEIGHT) - n / d
        //                   = ((n + VALUE * WEIGHT) * d - n * (d+WEIGHT)) / (d * (d + WEIGHT))
        // simplifying:
        // x_ave_new - x_ave = (VALUE - n / d) * WEIGHT / (d + WEIGHT)
        //                   = (VALUE - x_ave) * WEIGHT / (d + WEIGHT)
        // WEIGHT_RATIO = WEIGHT / (d + WEIGHT)
        // thus,
        // x_ave_new = x_ave + (VALUE - x_ave) * WEIGHT_RATIO
        // So updating is just x_ave += (VALUE - x_ave) * WEIGHT_RATIO

        SCALAR_AVE += WEIGHT_RATIO * ( t_DataType( VALUE) - SCALAR_AVE);
      }

      //! @brief Update the given average with the new value according to the given weight ratio
      //! @param SCALAR_AVE the scalar average to update
      //! @param VALUE the new value to update the scalar average with
      //! @param WEIGHT_RATIO the ratio of weights (weight for NEW_VALUE) / (total new weight)
      //! @details overloaded version for t_DataTypes without iterators and no type conversion
      template< typename t_DataType, typename t_OtherDataType>
      static
      typename type::EnableIf
      <
        !type::IsSequence< t_DataType>::value && type::Compare< t_OtherDataType, t_DataType>::e_Same,
        void
      >::Type
      UpdateAverage( t_DataType &SCALAR_AVE, const t_OtherDataType &VALUE, const double &WEIGHT_RATIO)
      {
        // in this version, t_OtherDataType and t_DataType are identical, so the cast needed above is unnecessary
        SCALAR_AVE += WEIGHT_RATIO * ( t_DataType( VALUE) - SCALAR_AVE);
      }

      //! @brief Update the given average with the new value according to the given weight ratio
      //! @param VECTOR_AVE the vector average to update
      //! @param SEQUENCE the new sequence of values to update the vector average with
      //! @param WEIGHT_RATIO the ratio of weights (weight for NEW_VALUE) / (total new weight)
      //! @details overloaded version for t_DataTypes with iterators
      template< typename t_DataType, typename t_OtherDataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
      UpdateAverage( t_DataType &VECTOR_AVE, const t_OtherDataType &SEQUENCE, const double &WEIGHT_RATIO)
      {
        // call the update average on each iterated-over object
        typename t_OtherDataType::const_iterator itr_seq( SEQUENCE.Begin()), itr_seq_end( SEQUENCE.End());
        for( typename t_DataType::iterator itr( VECTOR_AVE.Begin()); itr_seq != itr_seq_end; ++itr, ++itr_seq)
        {
          RunningAverage< t_AverageType>::UpdateAverage( *itr, *itr_seq, WEIGHT_RATIO);
        }
      }

    public:

      //! @brief = 0 for scalar types
      //! @param SCALAR the scalar to = with
      template< typename t_DataType>
      static typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
      Zero( t_DataType &SCALAR)
      {
        SCALAR = t_DataType( 0);
      }

      //! @brief = 0 for vector types
      //! @param VECTOR the vector to perform element-wise zeroing of
      template< typename t_DataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
        Zero( t_DataType &VECTOR)
      {
        for( typename t_DataType::iterator itr( VECTOR.Begin()), itr_end( VECTOR.End()); itr != itr_end; ++itr)
        {
          RunningAverage< t_AverageType>::Zero( *itr);
        }
      }
    }; // template class RunningAverage

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RUNNING_AVERAGE_H_
