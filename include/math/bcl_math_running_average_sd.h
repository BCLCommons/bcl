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

#ifndef BCL_MATH_RUNNING_AVERAGE_SD_H_
#define BCL_MATH_RUNNING_AVERAGE_SD_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RunningAverageSD
    //! @brief class maintains a running average and running variance of observations
    //! RunningAverageSD does not store observations, so it is ideal for
    //! computing the average / variance / std of vectors or large sets of numbers that are generated if only the
    //! average/variance/std is needed.
    //! RunningAverageSD works for vector and scalar types
    //!
    //! @tparam t_AverageType the type in which to store the average
    //!
    //! @see @link example_math_running_average_sd.cpp @endlink
    //! @author mendenjl
    //! @date Apr 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_AverageType>
    class RunningAverageSD :
      public util::ObjectInterface
    {

    public:

    //////////
    // enum //
    //////////

      //! The type of information currently held in the temporary variable (StdDev)
      enum StdDevType
      {
        e_None,     //!< Any information in the temporary variable is out of date
        e_Biased,   //!< Population standard deviation = sqrt(Sum((x(i)-ave_x)^2)/N)
        e_Unbiased, //!< Sample standard deviation = sqrt(Sum((x(i)-ave_x)^2)/(N-1))
        s_Unknown
      };

    private:

    //////////
    // data //
    //////////

      t_AverageType m_AveVariance;    //!< Running average of the variance
      t_AverageType m_Average;        //!< The current average
      double        m_Weight;         //!< Weight

      mutable t_AverageType m_StdDev;     //!< The last computed standard deviation
      mutable StdDevType    m_StdDevType; //!< Type of standard deviation information

    public:
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RunningAverageSD() :
        m_AveVariance(),
        m_Average(),
        m_Weight( 0),
        m_StdDev(),
        m_StdDevType( e_None)
      {
        Reset();
      }

      //! virtual copy constructor
      RunningAverageSD *Clone() const
      {
        return new RunningAverageSD( *this);
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
        const t_AverageType &( RunningAverageSD< t_AverageType>::*FUNCTION)() const
      ) const
      {
        static const std::string s_mean( "Mean"), s_var( "Variance"), s_std( "StandardDeviation"), s_unk( "Unknown");
        return
          FUNCTION == &RunningAverageSD< t_AverageType>::GetAverage ? s_mean
          : FUNCTION == &RunningAverageSD< t_AverageType>::GetVariance ? s_var
          : FUNCTION == &RunningAverageSD< t_AverageType>::GetStandardDeviation ? s_std
          : s_unk;
      }

      //! @brief get the current weight
      //! @return the current weight
      const double &GetWeight() const
      {
        return m_Weight;
      }

      //! @brief get the current average
      //! @return the current average
      const t_AverageType &GetAverage() const
      {
        return m_Average;
      }

      //! @brief get the current variance
      //! @return the current variance
      const t_AverageType &GetVariance() const
      {
        return m_AveVariance;
      }

      //! @brief et the current standard deviation
      //! @return the current standard deviation
      const t_AverageType &GetStandardDeviation() const
      {
        if( m_StdDevType != e_Biased)
        {
          // store the variance in Temp
          m_StdDev = m_AveVariance;

          // take the element-wise square root of temp
          ComputeSqrt( m_StdDev);

          m_StdDevType = e_Biased;
        }
        return m_StdDev;
      }

      //! @brief et the current standard deviation
      //! @return the current standard deviation
      const t_AverageType &GetSampleStandardDeviation() const
      {
        if( m_StdDevType != e_Unbiased)
        {
          // store the variance in Temp
          m_StdDev = m_AveVariance;

          // adjust for sample standard deviation using the bessel correction
          Multiply( m_StdDev, m_Weight / ( std::max( m_Weight, double( 1.5)) - 1.0));

          // take the element-wise square root of temp
          ComputeSqrt( m_StdDev);

          m_StdDevType = e_Unbiased;
        }
        return m_StdDev;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset the running average
      void Reset()
      {
        m_Weight = 0.0;
        RunningAverage< t_AverageType>::Zero( m_Average);
        RunningAverage< t_AverageType>::Zero( m_AveVariance);
        m_StdDevType = e_None;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief add an observation from the running average
      //! @param VALUE the value to add to the running average
      template< typename t_ArgumentType>
      const RunningAverageSD &operator +=( const t_ArgumentType &VALUE)
      {
        return AddWeightedObservation( VALUE, 1.0);
      }

      //! @brief remove an observation from the running average
      //! @param VALUE the value to add to the running average
      template< typename t_ArgumentType>
      const RunningAverageSD &operator -=( const t_ArgumentType &VALUE)
      {
        return AddWeightedObservation( VALUE, -1.0);
      }

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param VALUE the value to add to the running average
      //! @param WEIGHT the weight to give the observation (may be negative to remove an old observation)
      //! @param the weight to add to
      template< typename t_ArgumentType>
      const RunningAverageSD &AddWeightedObservation( const t_ArgumentType &VALUE, double WEIGHT)
      {
        m_StdDevType = e_None;
        if( WEIGHT == -m_Weight)
        {
          // if the weight is now zero, just make the average 0 too and return
          Reset();
          return *this;
        }
        else if( m_Weight == 0.0)
        {
          // if the current weight is 0.0, then just set m_Weight to WEIGHT and the average is just VALUE
          m_Average = VALUE;
          m_AveVariance = VALUE;
          RunningAverage< t_AverageType>::Zero( m_AveVariance);
          m_Weight = WEIGHT;
          return *this;
        }

        // see http://en.wikipedia.org/wiki/Standard_deviation#Weighted_calculation
        // given variables X = the new value
        //                 w = the weight for the new value
        //                 A = the prior average
        //                 W = the prior weight
        //                 V = the prior average variance
        // we can compute the updated quantities A', W', and V' as follows:
        // W' = W + w
        // dA = w / W' * (X - A)
        // A' = A + dA
        // ave variance = ave variance * W/W' + dA * (X - A')

        const double old_weight( m_Weight);
        m_Weight += WEIGHT;
        Update( m_Average, m_AveVariance, VALUE, WEIGHT / m_Weight, old_weight / m_Weight);
        return *this;
      }

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param VALUE the value to add to the running average
      //! @param WEIGHT the weight to give the observation (may be negative to remove an old observation)
      //! @param the weight to add to
      template< typename t_ArgumentType>
      const RunningAverageSD &AddEachElement( const t_ArgumentType &CONTAINER, double WEIGHT = 1.0)
      {
        for
        (
          typename t_ArgumentType::const_iterator itr( CONTAINER.Begin()), itr_end( CONTAINER.End());
          itr != itr_end;
          ++itr
        )
        {
          AddWeightedObservation( *itr, WEIGHT);
        }
        return *this;
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        io::Serialize::Read( m_Average, ISTREAM);
        io::Serialize::Read( m_AveVariance, ISTREAM);
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
        io::Serialize::Write( m_AveVariance, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Weight, OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief take the square root of a scalar value
      //! @param SCALAR the scalar to compute the square root of
      template< typename t_DataType>
      static typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
        ComputeSqrt( t_DataType &SCALAR)
      {
        SCALAR = Sqrt( std::max( SCALAR, t_DataType( 0.0)));
      }

      //! @brief multiply a vector by a value
      //! @param VECTOR the vector to perform element-wise *= on
      //! @param WEIGHT amount to multiply the vector by
      template< typename t_DataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
        Multiply( t_DataType &VECTOR, const double &WEIGHT)
      {
        for( typename t_DataType::iterator itr( VECTOR.Begin()), itr_end( VECTOR.End()); itr != itr_end; ++itr)
        {
          RunningAverageSD< t_AverageType>::Multiply( *itr, WEIGHT);
        }
      }
      //! @brief multiply a scalar by a value
      //! @param SCALAR the scalar to perform *= with
      //! @param WEIGHT amount to multiply the scalar by
      template< typename t_DataType>
      static typename type::EnableIf< !type::IsSequence< t_DataType>::value, void>::Type
        Multiply( t_DataType &SCALAR, const double &WEIGHT)
      {
        SCALAR *= WEIGHT;
      }

      //! @brief take the square root of a vector value
      //! @param VECTOR the vector to perform element-wise sqrt on
      template< typename t_DataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
        ComputeSqrt( t_DataType &VECTOR)
      {
        for( typename t_DataType::iterator itr( VECTOR.Begin()), itr_end( VECTOR.End()); itr != itr_end; ++itr)
        {
          RunningAverageSD< t_AverageType>::ComputeSqrt( *itr);
        }
      }

      //! @brief Update an average and average variance given a new value, weight ratio, and ratio of the prior weights
      //! @param SCALAR_AVE the scalar average to update
      //! @param SCALAR_AVE_VARIANCE the scalar variance to update
      //! @param VALUE the value to update the prior two variables with
      //! @param WEIGHT_RATIO the ratio of weights w / W'
      //! @param OLD_WEIGHT_RATIO the ratio of weights W / W'
      //! @details overloaded version for t_DataTypes without iterators
      template< typename t_DataType, typename t_OtherDataType>
      static
      typename type::EnableIf
      <
        !type::IsSequence< t_DataType>::value && !type::Compare< t_OtherDataType, t_DataType>::e_Same,
        void
      >::Type
      Update
      (
        t_DataType &SCALAR_AVE,
        t_DataType &SCALAR_AVE_VARIANCE,
        const t_OtherDataType &VALUE,
        const double &WEIGHT_RATIO,
        const double &OLD_WEIGHT_RATIO
      )
      {
        // see http://en.wikipedia.org/wiki/Standard_deviation#Weighted_calculation
        // given variables X = the new value
        //                 w = the weight for the new value
        //                 A = the prior average
        //                 W = the prior weight
        //                 V = the prior average variance
        // we can compute the updated quantities A', W', and V' as follows:
        // W' = W + w
        // dA = w / W' * (X - A)
        // A' = A + dA
        // ave variance = ave variance * W/W' + dA * (X - A')

        t_DataType value( VALUE);
        t_DataType delta_ave( WEIGHT_RATIO * ( value - SCALAR_AVE));
        SCALAR_AVE += delta_ave;
        SCALAR_AVE_VARIANCE *= OLD_WEIGHT_RATIO;
        SCALAR_AVE_VARIANCE += delta_ave * ( value - SCALAR_AVE);
      }

      //! @brief Update an average and average variance given a new value, weight ratio, and ratio of the prior weights
      //! @param SCALAR_AVE the scalar average to update
      //! @param SCALAR_AVE_VARIANCE the scalar variance to update
      //! @param VALUE the value to update the prior two variables with
      //! @param WEIGHT_RATIO the ratio of weights w / W'
      //! @param OLD_WEIGHT_RATIO the ratio of weights W / W'
      //! @details overloaded version for t_DataTypes without iterators and no type conversion
      template< typename t_DataType, typename t_OtherDataType>
      static
      typename type::EnableIf
      <
        !type::IsSequence< t_DataType>::value && type::Compare< t_OtherDataType, t_DataType>::e_Same,
        void
      >::Type
      Update
      (
        t_DataType &SCALAR_AVE,
        t_DataType &SCALAR_AVE_VARIANCE,
        const t_OtherDataType &VALUE,
        const double &WEIGHT_RATIO,
        const double &OLD_WEIGHT_RATIO
      )
      {
        // in this version, t_OtherDataType and t_DataType are identical, so the cast needed above is unnecessary
        t_DataType delta_ave( WEIGHT_RATIO * ( VALUE - SCALAR_AVE));
        SCALAR_AVE += delta_ave;
        SCALAR_AVE_VARIANCE *= OLD_WEIGHT_RATIO;
        SCALAR_AVE_VARIANCE += delta_ave * ( VALUE - SCALAR_AVE);
      }

      //! @brief Update an average and average variance given a new value, weight ratio, and ratio of the prior weights
      //! @param SCALAR_AVE the scalar average to update
      //! @param SCALAR_AVE_VARIANCE the scalar variance to update
      //! @param VALUE the value to update the prior two variables with
      //! @param WEIGHT_RATIO the ratio of weights w / W'
      //! @param OLD_WEIGHT_RATIO the ratio of weights W / W'
      //! @details overloaded version for t_DataTypes with iterators
      template< typename t_DataType, typename t_OtherDataType>
      static typename type::EnableIf< type::IsSequence< t_DataType>::value, void>::Type
      Update
      (
        t_DataType &VECTOR_AVE,
        t_DataType &VECTOR_AVE_VARIANCE,
        const t_OtherDataType &SEQUENCE,
        const double &WEIGHT_RATIO,
        const double &OLD_WEIGHT_RATIO
      )
      {
        // call the update average on each iterated-over object
        typename t_OtherDataType::const_iterator itr_seq( SEQUENCE.Begin()), itr_seq_end( SEQUENCE.End());
        for
        (
          typename t_DataType::iterator itr( VECTOR_AVE.Begin()), itr_var( VECTOR_AVE_VARIANCE.Begin());
          itr_seq != itr_seq_end;
          ++itr, ++itr_seq, ++itr_var
        )
        {
          RunningAverageSD< t_AverageType>::Update( *itr, *itr_var, *itr_seq, WEIGHT_RATIO, OLD_WEIGHT_RATIO);
        }
      }

    }; // template class RunningAverageSD

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LessThanAbsoluteMean
    //! @brief binary operator function struct for comparing the absolute value of the mean, or std's if means are equal
    //! @author alexanns
    //! @date Apr 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    struct BCL_API LessThanAbsoluteMean
    {
      template< typename t_AverageType>
      bool operator()( const RunningAverageSD< t_AverageType> &A, const RunningAverageSD< t_AverageType> &B) const
      {
        const t_AverageType lhs_abs_mean( Absolute( A.GetAverage()));
        const t_AverageType rhs_abs_mean( Absolute( B.GetAverage()));

        // true if rhs abs mean is undefined
        if( !util::IsDefined( rhs_abs_mean))
        {
          // return lhs is not less than rhs
          return false;
        }

        // true if lhs abs mean is undefined
        if( !util::IsDefined( lhs_abs_mean))
        {
          // return lhs is less than rhs
          return true;
        }

        // true if means are equal
        if( lhs_abs_mean == rhs_abs_mean)
        {
          // true if standard deviations are equal
          if( A.GetStandardDeviation() == B.GetStandardDeviation())
          {
            // return comparison of counts
            return A.GetWeight() < B.GetWeight();
          }

          // standard deviations are not equal return evaluation by standard deviation
          return A.GetStandardDeviation() < B.GetStandardDeviation();
        }

        // means are not equal, return evaluation by mean
        return lhs_abs_mean < rhs_abs_mean;
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LessThanCounts
    //! @brief binary operator function struct for comparing the counts/weights of the running averages
    //! @author alexanns
    //! @date Apr 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    struct BCL_API LessThanCounts
    {
      template< typename t_AverageType>
      bool operator()( const RunningAverageSD< t_AverageType> &A, const RunningAverageSD< t_AverageType> &B) const
      {
        return A.GetWeight() < B.GetWeight();
      }
    };

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_RUNNING_AVERAGE_SD_H_
