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
#include "math/bcl_math_roc_curve.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_binary_function_bind_second.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_linear_function.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_spline_border_type.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_logger_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from members
    //! @param FP # of false positives
    //! @param TP # of true positives
    //! @param CUTOFF actual cutoff value
    ROCCurve::Point::Point( const size_t &FP, const size_t &TP, const double &CUTOFF) :
      m_FalsePositives( FP),
      m_TruePositives( TP),
      m_Cutoff( CUTOFF)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ROCCurve
    ROCCurve::Point *ROCCurve::Point::Clone() const
    {
      return new ROCCurve::Point( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ROCCurve::Point::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get number of true positives
    //! @return number of true positives = TP
    size_t ROCCurve::Point::GetNumberTruePositives() const
    {
      return m_TruePositives;
    }

    //! @brief get number of false positives
    //! @return number of false positives = FP
    size_t ROCCurve::Point::GetNumberFalsePositives() const
    {
      return m_FalsePositives;
    }

    //! @brief number predicted positives
    //! @return sum of TP and FP = P'
    size_t ROCCurve::Point::GetNumberPredictedPositives() const
    {
      return m_TruePositives + m_FalsePositives;
    }

    //! @brief get the cutoff
    //! @return the cutoff
    double ROCCurve::Point::GetCutoff() const
    {
      return m_Cutoff;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Construct a contingency matrix, given the final point in the ROC
    //! @param FINAL_POINT the final point in the ROC
    //! @return a contingency matrix for this object
    ContingencyMatrix ROCCurve::Point::GetContingencyMatrix( const Point &FINAL_POINT) const
    {
      return
        ContingencyMatrix
        (
          m_TruePositives,
          m_FalsePositives,
          FINAL_POINT.m_TruePositives - m_TruePositives,
          FINAL_POINT.m_FalsePositives - m_FalsePositives
        );
    }

    //! @brief add true positives to the point
    //! @param RESULT true for true positive, false for false positive
    //! @param CUTOFF the new cutoff
    void ROCCurve::Point::Update( const bool &RESULT, const double &CUTOFF)
    {
      // update the counts
      RESULT ? ++m_TruePositives : ++m_FalsePositives;
      // update the threshold
      m_Cutoff = CUTOFF;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ROCCurve::Point::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_FalsePositives, ISTREAM);
      io::Serialize::Read( m_TruePositives, ISTREAM);
      io::Serialize::Read( m_Cutoff, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ROCCurve::Point::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_FalsePositives, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TruePositives, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Cutoff, OSTREAM, INDENT);
      return OSTREAM;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a list of pairs of prediction and test results not yet classified by a threshold
    //! @param TEST_RESULTS list of pairs of prediction ( double) and test outputs( double)
    //! @param THRESHOLD the threshold value to be used for classifying the outputs
    //! @param POSITIVES_ABOVE_THRESHOLD flag whether positives are above threshold
    //! ( classification means that results above the threshold are considered positive = 1)
    ROCCurve::ROCCurve
    (
      storage::List< storage::Pair< double, double> > &TEST_RESULTS,
      const double THRESHOLD,
      bool POSITIVES_ABOVE_THRESHOLD
    ) :
      m_SortedCounts()
    {
      // if positives are above threshold
      if( POSITIVES_ABOVE_THRESHOLD)
      {
        // sort in descending order
        TEST_RESULTS.Sort
        (
          storage::PairBinaryPredicateFirst< double, double>
          (
            util::BinaryFunctionSTLWrapper< std::greater< double> >()
          )
        );
      }
      else
      {
        // sort in ascending order
        TEST_RESULTS.Sort
        (
          storage::PairBinaryPredicateFirst< double, double>
          (
            util::BinaryFunctionSTLWrapper< std::less< double> >()
          )
        );
      }

      // initialize roc curve based on POSITIVES_ABOVE_THRESHOLD
      Initialize
      (
        POSITIVES_ABOVE_THRESHOLD
        ? ClassifyResults( TEST_RESULTS, Comparisons< double>::GetEnums().CreateUnaryPredicate( Comparisons< double>::GetEnums().e_Greater, THRESHOLD))
        : ClassifyResults( TEST_RESULTS, Comparisons< double>::GetEnums().CreateUnaryPredicate( Comparisons< double>::GetEnums().e_Less, THRESHOLD))
      );
    }

    //! @brief constructor from a list of pairs of prediction and test results not yet classified by a threshold
    //! @param TEST_RESULTS list of pairs of prediction ( double) and test outputs( double)
    //! @param THRESHOLD the threshold value to be used for classifying the outputs
    //! ( classification means that results above the threshold are considered positive = 1)
    ROCCurve::ROCCurve
    (
      const storage::List< storage::Pair< double, double> > &TEST_RESULTS,
      const double THRESHOLD
    ) :
      m_SortedCounts()
    {
      Initialize( ClassifyResults( TEST_RESULTS, Comparisons< double>::GetEnums().CreateUnaryPredicate( Comparisons< double>::GetEnums().e_Greater, THRESHOLD)));
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ROCCurve::s_Instance
    (
      GetObjectInstances().AddInstance( new ROCCurve())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the roc curve with a list of pair of predicted values and expected outputs
    void ROCCurve::Initialize( const storage::List< storage::Pair< double, bool> > &RESULTS_CLASSIFIED)
    {
      // reset the counts list
      m_SortedCounts.Reset();

      // if there are no results, return
      if( RESULTS_CLASSIFIED.IsEmpty())
      {
        return;
      }

      // get an iterator range on the results
      storage::List< storage::Pair< double, bool> >::const_iterator result_itr( RESULTS_CLASSIFIED.Begin());

      // add the first result to the counts of true/false
      Point counts;
      counts.Update( result_itr->Second(), result_itr->First());
      ++result_itr;

      // insert the current counts
      m_SortedCounts.PushBack( counts);
      m_SortedCounts.AllocateMemory( RESULTS_CLASSIFIED.GetSize());

      // iterate over sorted results
      for
      (
        storage::List< storage::Pair< double, bool> >::const_iterator
          result_itr_prev( RESULTS_CLASSIFIED.Begin()),
          result_itr_end( RESULTS_CLASSIFIED.End());
        result_itr != result_itr_end;
        ++result_itr, ++result_itr_prev
      )
      {
        // check for new threshold value
        if( result_itr_prev->First() != result_itr->First())
        {
          // insert the current counts
          m_SortedCounts.PushBack( m_SortedCounts.LastElement());
        }

        // initialize the count that corresponds to the classification of the result
        m_SortedCounts.LastElement().Update( result_itr->Second(), result_itr->First());
      }
    }

    //! @brief returns the total number of false positives
    //! @return total number of false positives
    const size_t ROCCurve::GetNumberActualNegatives() const
    {
      BCL_Assert( !m_SortedCounts.IsEmpty(), "There are no results stored in this ROCCurve curve!");
      return m_SortedCounts.LastElement().GetNumberFalsePositives();
    }

    //! @brief returns the total number of true positives
    //! @return total number of true positives
    const size_t ROCCurve::GetNumberActualPositives() const
    {
      BCL_Assert( !m_SortedCounts.IsEmpty(), "There are no results stored in this ROCCurve curve!");
      return m_SortedCounts.LastElement().GetNumberTruePositives();
    }

    //! @brief calculates and returns the area under the curve
    //! @return curve integral
    double ROCCurve::Integral() const
    {
      return Integral( m_SortedCounts.End());
    }

    //! @brief calculates and returns the weighted area under the curve by weighting function
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //! @return entire weighted curve integral
    double ROCCurve::Integral( const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION, const bool ROC_X_AXIS_LOG10_SCALING) const
    {
      return Integral( m_SortedCounts.End(), WEIGHTING_FUNCTION, ROC_X_AXIS_LOG10_SCALING);
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @return curve integral
    double ROCCurve::Integral( const double FRACTION) const
    {
      return Integral( GetEndIteratorForInterval( FRACTION), LinearFunction( 0, 1));
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //!        weighted by a given function
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @return weighted curve integral
    double ROCCurve::Integral
    (
      const double FRACTION,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral( GetEndIteratorForInterval( FRACTION), WEIGHTING_FUNCTION, ROC_X_AXIS_LOG10_SCALING);
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return curve integral
    double ROCCurve::Integral
    (
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        m_SortedCounts.End(),
        LinearFunction( 0, 1),
        ROC_X_AXIS,
        ROC_Y_AXIS,
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return curve integral
    double ROCCurve::Integral
    (
      const double FRACTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        GetEndIteratorForInterval( FRACTION),
        LinearFunction( 0, 1),
        ROC_X_AXIS,
        ROC_Y_AXIS,
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_X_AXIS_FRACTION_CUTOFF x axis value cutoff rather then specifying a number of counts
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return curve integral
    double ROCCurve::Integral
    (
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> &ROC_X_AXIS_FRACTION_CUTOFF,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        m_SortedCounts.End(),
        LinearFunction( 0, 1),
        ROC_X_AXIS,
        ROC_Y_AXIS,
        ROC_X_AXIS_FRACTION_CUTOFF,
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //!        weighted by a given function
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return weighted curve integral
    double ROCCurve::Integral
    (
      const double FRACTION,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        GetEndIteratorForInterval( FRACTION),
        WEIGHTING_FUNCTION,
        ROC_X_AXIS,
        ROC_Y_AXIS,
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
    //!        weighted by a given function
    //!        the contingency matrix functions to plot on x- and y-axis can be specified
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_X_AXIS_FRACTION_CUTOFF x axis cutoff range rather then specifying a number of counts
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return weighted curve integral
    double ROCCurve::Integral
    (
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> &ROC_X_AXIS_FRACTION_CUTOFF,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        m_SortedCounts.End(),
        WEIGHTING_FUNCTION,
        ROC_X_AXIS,
        ROC_Y_AXIS,
        ROC_X_AXIS_FRACTION_CUTOFF,
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates ContingencyMatrix with true/false positive/negatives
    //! @param FRACTION consider the first fraction as being predicted positive, the remainder predicted negative
    ContingencyMatrix ROCCurve::ContingencyMatrixFraction( const double FRACTION) const
    {
      if( m_SortedCounts.IsEmpty())
      {
        return ContingencyMatrix();
      }

      // closed interval [0,1] for valid fraction
      static const Range< double> s_valid_fraction_range( 0.0, 1.0);

      // check that fraction is between 0 and 1
      BCL_Assert
      (
        s_valid_fraction_range.IsWithin( FRACTION),
        "given fraction not in interval " + s_valid_fraction_range.GetString() + ": " + util::Format()( FRACTION)
      );

      // compute the # of predictions desired
      const size_t desired_predictions( std::min( size_t( std::ceil( GetNumberResults() * FRACTION)), GetNumberResults()));

      // number of predictions in first fraction
      storage::Vector< Point>::const_iterator
        itr_fraction_last( m_SortedCounts.Begin()), itr_fraction_end( m_SortedCounts.Last());

      // walk in the list until the # of predicted values is correct
      while
      (
        itr_fraction_last != itr_fraction_end
        && itr_fraction_last->GetNumberPredictedPositives() < desired_predictions
      )
      {
        ++itr_fraction_last;
      }

      // if the total counts is exactly equal to the desired predictions, return the desired contingency matrix
      if( itr_fraction_last->GetNumberPredictedPositives() == desired_predictions)
      {
        return itr_fraction_last->GetContingencyMatrix( m_SortedCounts.LastElement());
      }

      // get the counts at the itr, which is the minimal # of counts greater than or equal to the desired # of predictions
      Point counts_gte_desired( *itr_fraction_last);

      // get the counts just before itr, which is the minimal # of counts less than the desired # of predictions
      Point counts_lt_desired;

      if( itr_fraction_last != m_SortedCounts.Begin())
      {
        // go to the previous value;
        counts_lt_desired = *--itr_fraction_last;
      }

      const size_t total_lt( counts_lt_desired.GetNumberPredictedPositives());
      const size_t total_gte( counts_gte_desired.GetNumberPredictedPositives());

      // compute true and false positives slopes
      const double true_positives_slope
      (
        double( counts_gte_desired.GetNumberTruePositives() - counts_lt_desired.GetNumberTruePositives())
        / double( total_gte - total_lt)
      );
      const double false_positives_slope
      (
        double( counts_gte_desired.GetNumberFalsePositives() - counts_lt_desired.GetNumberFalsePositives())
        / double( total_gte - total_lt)
      );

      // compute true positives and false positives using interpolation between the two values
      const size_t true_positives
      (
        counts_lt_desired.GetNumberTruePositives() + true_positives_slope * ( desired_predictions - total_lt)
      );
      const size_t false_positives
      (
        counts_lt_desired.GetNumberFalsePositives() + false_positives_slope * ( desired_predictions - total_lt)
      );

      // end
      return
        Point( false_positives, true_positives, counts_lt_desired.GetCutoff()).GetContingencyMatrix
        (
          m_SortedCounts.LastElement()
        );
    }

    //! @brief calculates ContingencyMatrix with true/false positive/negatives
    //! @param FRACTION consider the first fraction as being predicted positive, the remainder predicted negative
    double ROCCurve::CutoffFraction( const double FRACTION) const
    {
      if( m_SortedCounts.IsEmpty())
      {
        return 0.0;
      }

      /////////////
      if( FRACTION == 0.0)
      {
        return m_SortedCounts.FirstElement().GetCutoff();
      }

      if( FRACTION == 1.0)
      {
        return m_SortedCounts.LastElement().GetCutoff();
      }
      ////////////////

      // closed interval [0,1] for valid fraction
      static const Range< double> s_valid_fraction_range( 0.0, 1.0);

      // check that fraction is between 0 and 1
      BCL_Assert
      (
        s_valid_fraction_range.IsWithin( FRACTION),
        "given fraction not in interval " + s_valid_fraction_range.GetString() + ": " + util::Format()( FRACTION)
      );

      // compute the # of predictions desired
//      const size_t desired_predictions( std::min( size_t( std::ceil( GetNumberResults() * FRACTION)), GetNumberResults()));
      const float desired_predictions( float( GetNumberResults()) * FRACTION);

      // number of predictions in first fraction
      storage::Vector< Point>::const_iterator
        itr_fraction_last( m_SortedCounts.Begin()), itr_fraction_end( m_SortedCounts.Last());

      // walk in the list until the # of predicted values is correct
      while
      (
        itr_fraction_last != itr_fraction_end
        && float( itr_fraction_last->GetNumberPredictedPositives()) < desired_predictions
      )
      {
        ++itr_fraction_last;
      }

      // if the total counts is exactly equal to the desired predictions, return the desired contingency matrix
      if( itr_fraction_last->GetNumberPredictedPositives() == desired_predictions)
      {
        return itr_fraction_last->GetCutoff();
      }

      // get the counts at the itr, which is the minimal # of counts greater than or equal to the desired # of predictions
      Point counts_gte_desired( *itr_fraction_last);

      // get the counts just before itr, which is the minimal # of counts less than the desired # of predictions
      Point counts_lt_desired;

      if( itr_fraction_last != m_SortedCounts.Begin())
      {
        // go to the previous value;
        counts_lt_desired = *--itr_fraction_last;
      }
      else
      {
          return itr_fraction_last->GetCutoff();
      }
//
//
      // {
        //if( ++itr_fraction_last != m_SortedCounts.End())
       // {
         // counts_lt_desired = *itr_fraction_last;
       // }
      //}

      const float total_lt( float( counts_lt_desired.GetNumberPredictedPositives()));
      const float total_gte( float( counts_gte_desired.GetNumberPredictedPositives()));

      // compute true and false positives slopes
      const double cutoff_slope
      (
        double( counts_gte_desired.GetCutoff() - counts_lt_desired.GetCutoff())
        / double( total_gte - total_lt)
      );
      return cutoff_slope * ( desired_predictions - total_lt) + counts_lt_desired.GetCutoff();
    }

    //! @brief calculates and returns the area under the curve for the counts between provided regions
    //! @param END iterator to the end of the interval for which the integral is going to be calculated
    //! @return curve integral for values between specified regions
    double ROCCurve::Integral
    (
      const storage::Vector< Point>::const_iterator &END
    ) const
    {
      // return integral with no weighting by a function
      return Integral( END, LinearFunction( 0, 1));
    }

    //! @brief calculates and returns the area under the curve plotting False Positive Rate on the x axis
    //!        and True Positive Rate on the y axis for the counts between provided regions
    //!        weighted with a mathematical function
    //! @param END iterator to the end of the interval for which the integral is going to be calculated
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @return weighted curve integral for values between specified regions
    double ROCCurve::Integral
    (
      const storage::Vector< Point>::const_iterator &END,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      return Integral
      (
        END,
        WEIGHTING_FUNCTION,
        &ContingencyMatrix::GetFalsePositiveRate, // X-coordinate
        &ContingencyMatrix::GetTruePositiveRate,  // Y-coordinate
        Range< double>( 0, 1),
        ROC_X_AXIS_LOG10_SCALING
      );
    }

    //! @brief calculates and returns the area under the curve for the counts between provided regions
    //!        weighted with a mathematical function
    //! @param END iterator to the end of the interval for which the integral is going to be calculated
    //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
    //!        accordingly
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_X_AXIS_RANGE range on x axis that specifies a begin and end fraction
    //!        rather then specifying a number of counts through param END
    //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
    //! @return weighted curve integral for values between specified regions,
    //!         the integral is normalized according to the specified range depending on the x axis scaling
    double ROCCurve::Integral
    (
      const storage::Vector< Point>::const_iterator &END,
      const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> ROC_X_AXIS_RANGE,
      const bool ROC_X_AXIS_LOG10_SCALING
    ) const
    {
      BCL_Assert
      (
        ROC_X_AXIS_RANGE.GetMax() > 0 && ROC_X_AXIS_RANGE.GetMin() <= 1,
        "x axis cutoff has to be between 0 < x_axis_cutoff <= 1! current range is " + util::Format()( ROC_X_AXIS_RANGE)
      );

      storage::Vector< Point>::const_iterator counts_itr( m_SortedCounts.Begin());

      // if interval is empty return 0
      if( counts_itr == END)
      {
        BCL_MessageStd
        (
          "This interval asked for ROCCurve curve is empty, therefore returning 0 for integral"
        );
        return double( 0);
      }

      // total count of data points
      const size_t total_count( GetNumberResults() - 1);

      // initialize minimum of roc curve x-axis dependent on enabled log scaling
      const double roc_min
      (
        ROC_X_AXIS_LOG10_SCALING
        ? ROC_X_AXIS_RANGE.GetMin() == 0
          ? std::log10( 1 / double( total_count))
          : std::log10( ROC_X_AXIS_RANGE.GetMin())
        : ROC_X_AXIS_RANGE.GetMin()
      );

      // x axis step size = 1 / (N + P)
      const double step_size
      (
        ROC_X_AXIS_LOG10_SCALING
        ? ( std::log10( ROC_X_AXIS_RANGE.GetMax()) - roc_min) / double( total_count)
        : ( 1.0 - roc_min) / double( total_count)
      );

      // runnning average of integral segments
      RunningAverage< double> integral;

      // map for associating x-axis values to y-axis values
      storage::Map< double, RunningAverage< double> > xy_plot( ComputeXAxisYAxisMapping( ROC_X_AXIS, ROC_Y_AXIS));

      // get the actual end of the x-axis value based on the END iterator provided
      float end_fraction
      (
        END != m_SortedCounts.End()
        ? ROC_X_AXIS( END->GetContingencyMatrix( m_SortedCounts.LastElement()))
        : (
              ROC_X_AXIS_LOG10_SCALING
            ? Pow( 10.0, total_count * step_size + roc_min)
            : total_count * step_size + roc_min
          )
      );
      if( end_fraction > ROC_X_AXIS_RANGE.GetMax())
      {
        end_fraction = ROC_X_AXIS_RANGE.GetMax();
      }

      BCL_Assert( xy_plot.GetSize(), "Roc curve: x axis value map is empty!")

      if( xy_plot.GetSize() == size_t( 1))
      {
        return xy_plot.Begin()->second;
      }

      // iterate over every count and calculate the integral of every integration segment
      for( size_t distance_counter( 0); distance_counter < total_count; ++distance_counter)
      {
        // determine fraction of x axis
        double fraction( distance_counter * step_size);

        // if log scaling is activated
        if( ROC_X_AXIS_LOG10_SCALING)
        {
          fraction = Pow( 10.0, fraction + roc_min);
        }
        else
        {
          fraction += roc_min;
        }

        // stop if x axis range is reached
        if( fraction > end_fraction)
        {
          break;
        }
        else if( fraction < ROC_X_AXIS_RANGE.GetMin())
        {
          continue;
        }

        storage::Map< double, RunningAverage< double> >::const_iterator itr_xy_plot_lower, itr_xy_plot_upper;

        // determine boundaries for every integration segment
        itr_xy_plot_lower = itr_xy_plot_upper = xy_plot.UpperBound( fraction);

        if( itr_xy_plot_upper == xy_plot.End())
        {
          itr_xy_plot_lower = --itr_xy_plot_upper;
        }
        if( itr_xy_plot_lower != xy_plot.Begin())
        {
          --itr_xy_plot_lower;
        }

        if( itr_xy_plot_lower == itr_xy_plot_upper)
        {
          integral += itr_xy_plot_lower->second * WEIGHTING_FUNCTION( fraction);
          continue;
        }

        // final boundary values for iteration segment
        const double upper_x( itr_xy_plot_upper->first);
        const double upper_y( itr_xy_plot_upper->second);
        const double lower_x( itr_xy_plot_lower->first);
        const double lower_y( itr_xy_plot_lower->second);

        double integral_delta( 0);

        // final boundary values when log scaling is enabled
        const double slope( ( upper_y - lower_y) / ( upper_x - lower_x));
        const double y_icept( lower_y - lower_x * slope);

        integral_delta = fraction * slope + y_icept;

        // averaged integral
        integral += integral_delta * WEIGHTING_FUNCTION( fraction);
      }

      // if integral is zero, write out a message
      if( integral == 0)
      {
        BCL_MessageStd( "calculated integral has value 0 !");
        return integral;
      }

      BCL_MessageDbg( "calculated integral: " + util::Format()( integral));
      return integral;
    }

    //! @brief compute the x axis and y axis mapping for computing an interpolated integral under the curve
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
    //! @return x axis and y axis mapping
    storage::Map< double, RunningAverage< double> > ROCCurve::ComputeXAxisYAxisMapping
    (
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS
    ) const
    {
      storage::Map< double, RunningAverage< double> > xy_plot;

      // add data for the 0 point
      {
        ContingencyMatrix con_matrix
        (
          0,                          // true_positives
          0,                          // false_positives
          GetNumberActualPositives(), // false_negatives
          GetNumberActualNegatives()  // true_negatives
        );

        // only add the value if it is defined on both axis!
        const double x_axis_value( ROC_X_AXIS( con_matrix));
        const double y_axis_value( ROC_Y_AXIS( con_matrix));
        if( util::IsDefined( x_axis_value) && util::IsDefined( y_axis_value))
        {
          xy_plot[ x_axis_value] += y_axis_value;
        }
      }

      if( !m_SortedCounts.IsEmpty())
      {
        Point last_point( m_SortedCounts.LastElement());
        // iterate over rates
        for
        (
          storage::Vector< Point>::const_iterator
            counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
          counts_itr != counts_itr_end;
          ++counts_itr
        )
        {
          ContingencyMatrix con_matrix( counts_itr->GetContingencyMatrix( last_point));
          xy_plot[ ROC_X_AXIS( con_matrix)] += ROC_Y_AXIS( con_matrix);
        }
      }

      return xy_plot;
    }

    //! @brief Convert into a map from threshold to contingency matrix
    storage::Map< double, ContingencyMatrix> ROCCurve::ToMap() const
    {
      storage::Map< double, ContingencyMatrix> matrices;

      Point last_point( m_SortedCounts.LastElement());
      // iterate over rates
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        matrices[ counts_itr->GetCutoff()] = counts_itr->GetContingencyMatrix( last_point);
      }
      return matrices;
    }

    //! @brief get the optima for a particular contingency matrix measure
    //! @param MEASURE the measure to optimize
    //! @return an iterator to the maximum point for the particular measure and the maximum value
    std::pair< storage::Vector< ROCCurve::Point>::const_iterator, double> ROCCurve::GetMaxima
    (
      const util::FunctionInterfaceSerializable< ContingencyMatrix, double> &MEASURE
    ) const
    {
      Point last_point( m_SortedCounts.LastElement());
      storage::Vector< ROCCurve::Point>::const_iterator itr_best( m_SortedCounts.Begin());
      double best_val( util::GetUndefined< double>());
      // iterate over points
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        const double local_metric_value( MEASURE( counts_itr->GetContingencyMatrix( last_point)));
        if
        (
          util::IsDefined( local_metric_value)
          && ( !util::IsDefined( best_val) || local_metric_value > best_val)
        )
        {
          best_val = local_metric_value;
          itr_best = counts_itr;
        }
      }
      return std::make_pair( itr_best, best_val);
    }

    //! @brief creates a thinned version of the ca1ounts so that only every Nth( PERIODICITY * n) element is kept
    //! @brief This function is to be used for plotting.
    //! @param PERIODICITY periodicity with which elements should be selected
    //! @return thinned ROC curve counts
    ROCCurve ROCCurve::GetThinnedRocCurvePeriodicity( const size_t PERIODICITY) const
    {
      // check that fraction is between 1 and total size
      BCL_Assert
      (
        PERIODICITY != 0 && PERIODICITY < GetNumberResults(),
        "The provided periodicity: " + util::Format()( PERIODICITY) +
        " cannot be 0 or larger than total number of results: " + util::Format()( GetNumberResults())
      );

      // initialize the new counts list
      storage::Vector< Point> new_counts;

      // initialize iterator to the beginning
      storage::Vector< Point>::const_iterator iterator( m_SortedCounts.Begin());

      // initialize end iterator
      const storage::Vector< Point>::const_iterator end_iterator( m_SortedCounts.End());

      // start the target # of predictions counter
      size_t desired_predictions( PERIODICITY);

      // iterate until you hit the end
      while( iterator != end_iterator)
      {
        // push back this element to new counts
        new_counts.PushBack( *iterator);

        // walk in the list until the # of predicted values is correct
        while( iterator != m_SortedCounts.End() && iterator->GetNumberPredictedPositives() < desired_predictions)
        {
          ++iterator;
        }
        desired_predictions += PERIODICITY;
      }

      // create a new ROCCurve to hold the counts
      ROCCurve roc_with_new_counts;

      // set the counts
      roc_with_new_counts.m_SortedCounts = new_counts;

      // return the roc curve with the new counts
      return roc_with_new_counts;
    }

    //! @brief creates a thinned version of the counts so that only TOTAL_NUMBER_ENTRIES elements are kept
    //! @brief This function is to be used for plotting.
    //! @param TOTAL_NUMBER_ENTRIES periodicity with which elements should be selected
    //! @return thinned ROC curve counts
    ROCCurve ROCCurve::GetThinnedRocCurveTotal( const size_t TOTAL_NUMBER_ENTRIES) const
    {
      return GetThinnedRocCurvePeriodicity( std::max( size_t( 1), m_SortedCounts.GetSize() / TOTAL_NUMBER_ENTRIES));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream from a formatted output
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ROCCurve::ReadPlottingTable( std::istream &ISTREAM)
    {
      // reset the counts
      m_SortedCounts.Reset();

      // read the identifier to make sure the file has correct header
      std::string identifier;
      ISTREAM >> identifier;
      BCL_Assert( identifier == "FALSE_POSITIVE", "The provided identifier is not correct" + identifier);

      // read the remaining of the header
      ISTREAM >> identifier;

      // read the number of examples
      size_t number_examples, example_ctr( 0);
      ISTREAM >> number_examples;

      // while number of examples is not reached
      while( example_ctr < number_examples && !ISTREAM.eof())
      {
        // read a new count pair and insert it into m_SortedCounts
        size_t first, second;
        ISTREAM >> first >> second;
        m_SortedCounts.PushBack( Point( first, second, double( example_ctr)));

        // increase the example ctr
        ++example_ctr;
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream in a formatted output style the result between BEGIN and END
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ROCCurve::WritePlottingTable
    (
      std::ostream &OSTREAM
    ) const
    {
      // print number of elements
      OSTREAM << "FALSE_POSITIVE\tTRUE_POSITIVE\t" << m_SortedCounts.GetSize() << '\n';

      // x axis step size = 1 / (N + P)
      const double step_size( 1 / double( GetNumberResults()));

      // total count of data points
      const double max_fraction( 1);

      ContingencyMatrix con_matrix;

      // iterate over counts
      for( double fraction( step_size); fraction <= max_fraction; fraction += step_size)
      {
        con_matrix = ContingencyMatrixFraction( fraction);

        // output
        OSTREAM << con_matrix.GetNumberFalsePositives() << '\t' << con_matrix.GetNumberTruePositives() << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief write to std::ostream in a formatted output style the result between BEGIN and END
    //! @param OSTREAM output stream to write to
    //! @param FORMAT format the rate
    //! @return output stream which was written to
    std::ostream &ROCCurve::WriteRatePlottingTable
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT
    ) const
    {
      // print number of elements
      OSTREAM << "FALSE_POSITIVE_RATE\tTRUE_POSITIVE_RATE\t"
              << GetNumberActualNegatives() << '\t' << GetNumberActualPositives()
              << '\n';

      // end
      return WriteRatePlottingTableGeneric
      (
        OSTREAM,
        &ContingencyMatrix::GetFalsePositiveRate, // X-coordinate
        &ContingencyMatrix::GetTruePositiveRate,  // Y-coordinate
        Range< double>( 0.0, 1.0),                // X-axis range
        FORMAT
      );
    }

    //! @brief write to std::ostream in a formatted output style the result between BEGIN and END
    //! @param OSTREAM output stream to write to
    //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix, specifies what is plotted
    //!        on x axis of roc curve
    //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix, specifies what is plotted
    //!        on y axis of roc curve
    //! @param ROC_X_AXIS_RANGE x axis value range specifying the x axis value when to begin and stop plotting
    //!        default is set to [ 0.0, 1.0], 0% to 100%
    //! @param FORMAT format the rate
    //! @param PLOTTING_TABLE_TITLE title that describes columns for x axis values and y axis values
    //! @return output stream which was written to
    std::ostream &ROCCurve::WriteRatePlottingTableGeneric
    (
      std::ostream &OSTREAM,
      const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
      const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
      const Range< double> &ROC_X_AXIS_RANGE,
      const util::Format &FORMAT,
      const std::string &PLOTTING_TABLE_TITLE
    ) const
    {
      BCL_Assert
      (
        ROC_X_AXIS_RANGE.GetMin() >= 0 && ROC_X_AXIS_RANGE.GetMax() <= 1,
        "WriteRatePlottingTableGeneric: x axis cutoff range has to be between [ 0.0, 1.0] but is " + util::Format()( ROC_X_AXIS_RANGE)
      );

      if( !PLOTTING_TABLE_TITLE.empty())
      {
        // print title
        OSTREAM << PLOTTING_TABLE_TITLE << '\n';
      }

      // iterate over rates
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        // update running contingency matrix at current rate
        ContingencyMatrix con_matrix( counts_itr->GetContingencyMatrix( m_SortedCounts.LastElement()));

        const double roc_x_axis( ROC_X_AXIS( con_matrix));

        // if x axis cutoff did not enter allowed range yet then continue integration
        if( roc_x_axis < ROC_X_AXIS_RANGE.GetMin())
        {
          continue;
        }

        // if x axis cutoff exceeding allowed maximum range then break integration
        if( roc_x_axis > ROC_X_AXIS_RANGE.GetMax())
        {
          break;
        }

        // output
        OSTREAM << FORMAT( roc_x_axis) << '\t' << FORMAT( ROC_Y_AXIS( con_matrix)) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief write all common contingency matrix measures to the file, return a vector of strings for each measure that was written
    //! @param DATA_FILENAME file name for the table
    //! @param MEASURES a list of measures of interest
    //! @return strings of all measures that were written
    storage::Vector< std::string> ROCCurve::WriteRatePlottingTableComplete
    (
      const std::string &DATA_FILENAME,
      const storage::Vector< std::string> &MEASURES
    ) const
    {
      io::OFStream table;
      io::File::MustOpenOFStream( table, DATA_FILENAME);
      // create all measures
      storage::Map< std::string, ContingencyMatrixMeasures> measure_map;

      // cutoff and local-PPV are special (since they are really unknown by the contingency matrix),
      // so get its alias and check whether it will be used
      size_t cutoff_index( util::GetUndefined< size_t>());
      size_t local_ppv_index( util::GetUndefined< size_t>());
      for( size_t measure( 0); measure < ContingencyMatrixMeasures::s_NumberMeasures; ++measure)
      {
        ContingencyMatrixMeasures this_measure( static_cast< ContingencyMatrixMeasures::Measure>( measure));
        measure_map[ this_measure.GetAlias()] = this_measure;
      }
      storage::Vector< std::string> measure_strings( MEASURES.IsEmpty() ? ContingencyMatrixMeasures::s_NumberMeasures : MEASURES.GetSize());
      storage::Vector< ContingencyMatrixMeasures> all_measures( MEASURES.IsEmpty() ? ContingencyMatrixMeasures::s_NumberMeasures : MEASURES.GetSize());
      if( MEASURES.IsEmpty())
      {
        for( size_t measure( 0); measure < ContingencyMatrixMeasures::s_NumberMeasures; ++measure)
        {
          all_measures( measure) = ContingencyMatrixMeasures( static_cast< ContingencyMatrixMeasures::Measure>( measure));
          measure_strings( measure) = all_measures( measure).GetAlias();
          table << measure_strings( measure) << '\t';
        }
        cutoff_index = size_t( ContingencyMatrixMeasures::e_Cutoff);
        local_ppv_index = size_t( ContingencyMatrixMeasures::e_LocalPPV);
      }
      else
      {
        measure_strings = MEASURES;
        for( size_t measure( 0), n_measures( MEASURES.GetSize()); measure < n_measures; ++measure)
        {
          all_measures( measure) = measure_map[ MEASURES( measure)];
          if( all_measures( measure).GetMeasure() == ContingencyMatrixMeasures::e_Cutoff)
          {
            cutoff_index = measure;
          }
          else if( all_measures( measure).GetMeasure() == ContingencyMatrixMeasures::e_LocalPPV)
          {
            local_ppv_index = measure;
          }
          table << measure_strings( measure) << '\t';
        }
      }
      table << '\n';

      PiecewiseFunction local_ppv;
      if( util::IsDefined( local_ppv_index))
      {
        local_ppv = GetLocalPPVCurve();
      }

      // iterate over rates
      util::Format format;
      for
      (
        storage::Vector< Point>::const_iterator
          counts_itr( m_SortedCounts.Begin()), counts_itr_end( m_SortedCounts.End());
        counts_itr != counts_itr_end;
        ++counts_itr
      )
      {
        // update running contingency matrix at current rate
        ContingencyMatrix con_matrix( counts_itr->GetContingencyMatrix( m_SortedCounts.LastElement()));

        size_t measure_id( 0);
        for
        (
          storage::Vector< ContingencyMatrixMeasures>::const_iterator
            itr_measure( all_measures.Begin()), itr_measure_end( all_measures.End());
          itr_measure != itr_measure_end;
          ++itr_measure, ++measure_id
        )
        {
          if( measure_id != cutoff_index && measure_id != local_ppv_index)
          {
            table << format( ( *itr_measure)( con_matrix)) << '\t';
          }
          else if( measure_id == cutoff_index)
          {
            table << format( counts_itr->GetCutoff()) << '\t';
          }
          else if( measure_id == local_ppv_index)
          {
            table << local_ppv( counts_itr->GetCutoff()) << '\t';
          }
        }
        table << '\n';
      }

      io::File::CloseClearFStream( table);

      // end
      return measure_strings;
    }

    //! @brief get the localized PPV curve
    //! The localized PPV-curve yields a more precise estimate of the local PPVs for a given range
    //! The standard ROC-curve answers the question: what is the PPV of all points predicted at or above a given
    //! cutoff (assuming +parity). The localized ppv-curve, by contrast, provides the PPV of values near a particular
    //! cutoff. This curve is produced by iteratively finding the optimal PPV interest, then removing all
    //! entries (P/N) at or above/below (depending on parity) from the ROC curve. PPV is essentially the only metric
    //! (other than the obscure false discovery rate = 1-PPV) that is independent of negatives and the total number of
    //! positives, and so is the only metric that can be computed locally in this manner.
    //! @return a piecewise function for the PPV curve, which is guaranteed to be monotonic
    PiecewiseFunction ROCCurve::GetLocalPPVCurve() const
    {
      ROCCurve copy( *this);
      const ContingencyMatrixMeasures ppv( ContingencyMatrixMeasures::e_PositivePredictiveValue);

      // true if positive cases should be predicted higher than negative cases
      const bool positives_higher( m_SortedCounts.FirstElement().GetCutoff() > m_SortedCounts.LastElement().GetCutoff());

      // keep track of all x, y points for construction of spline
      storage::Vector< double> xs, ys;

      // keep finding the point with the highest PPV value; then remove all predicted positives contributing to that ppv
      // after recording it in the xs (cutoff) and ys (ppv) vector.
      while( !copy.m_SortedCounts.IsEmpty())
      {
        std::pair< storage::Vector< ROCCurve::Point>::const_iterator, double> optima( copy.GetMaxima( ppv));
        storage::Vector< ROCCurve::Point>::const_iterator itr_optima( optima.first);
        ROCCurve::Point optimal_point( *itr_optima);

        // add this optima to the splined points
        xs.PushBack( copy.CutoffFraction( optimal_point.GetNumberPredictedPositives() / ( 2.0 * copy.GetNumberResults())));

//       if(xs.GetSize() >= size_t(2) && xs.LastElement() < xs(xs.GetSize()-2) && copy.m_SortedCounts.GetSize() <= size_t(5))
//       {
//        BCL_MessageStd(util::Format()(copy.m_SortedCounts));
//       }

        ys.PushBack( optima.second);

        // remove the elements that form this optimal PPV from the curve
        copy.m_SortedCounts.RemoveElements( 0, std::distance( copy.GetSortedCounts().Begin(), itr_optima) + 1);
        // subtract the predicted positives from the copy's points
        for
        (
          storage::Vector< ROCCurve::Point>::iterator
            itr_copy( copy.m_SortedCounts.Begin()), itr_copy_end( copy.m_SortedCounts.End());
          itr_copy != itr_copy_end;
          ++itr_copy
        )
        {
          *itr_copy =
            ROCCurve::Point
            (
              itr_copy->GetNumberFalsePositives() - optimal_point.GetNumberFalsePositives(),
              itr_copy->GetNumberTruePositives() - optimal_point.GetNumberTruePositives(),
              itr_copy->GetCutoff()
            );
        }
      }

      // spline requires sorted points; but if positives are higher, then x is a descending vector
      if( positives_higher)
      {
        std::reverse( xs.Begin(), xs.End());
        std::reverse( ys.Begin(), ys.End());
      }

      // create a damped spline to ensure monotonicity

      // no begin/end slope.
      // Assume that the best ppv on this data is really the best PPV.
      // In general there's no way of knowing what the actual maximum that the model will predict is or what the ppv is
      // at that value, without having access to training dataset predictions
      CubicSplineDamped csvd;
      csvd.Train
      (
        linal::Vector< double>( xs.Begin(), xs.End()),
        linal::Vector< double>( ys.Begin(), ys.End()),
        0.0,
        0.0
      );

      // create a trivial piecewise function. This allows for easy editing later on if we decide to use a different
      // function
      storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > functions;
      functions.PushBack
      (
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        (
          Range< double>
          (
            RangeBorders::e_LeftClosed,
            -std::numeric_limits< double>::max(),
            std::numeric_limits< double>::max(),
            RangeBorders::e_RightClosed
          ),
          util::ShPtr< FunctionInterfaceSerializable< double, double> >( csvd.Clone())
        )
      );

      return PiecewiseFunction( functions);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ROCCurve::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SortedCounts, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ROCCurve::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SortedCounts, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns iterator on the the end of the specified fraction of the list
    //! @param FRACTION top fraction ( between 0 and 1) to be used in consequent operations
    //! @return iterator to ending of the requested interval
    storage::Vector< ROCCurve::Point>::const_iterator
    ROCCurve::GetEndIteratorForInterval( const double FRACTION) const
    {
      // assert that the fraction is between 0 and including 1
      BCL_Assert
      (
        FRACTION > 0.0 && FRACTION <= 1.0,
        "The provided fraction should be between 0 and including 1 and not" + util::Format()( FRACTION)
      );

      // initialize an iterator
      storage::Vector< Point>::const_iterator iterator( m_SortedCounts.Begin());

      // compute the # of predictions desired
      const size_t desired_predictions( size_t( GetNumberResults() * FRACTION));

      // walk in the list until the # of predicted values is correct
      while
      (
        iterator != m_SortedCounts.End() &&
        iterator->GetNumberPredictedPositives() < desired_predictions
      )
      {
        ++iterator;
      }

      // return
      return iterator;
    }

    //! @brief static function to classify provided results and return the classified expected outputs
    //! @param UNCLASSIFIED_RESULTS results that are not yet classified
    //! @param UNARY_CLASSIFIER classifier to be used for classifying expected outputs
    //! @return classified results
    storage::List< storage::Pair< double, bool> > ROCCurve::ClassifyResults
    (
      const storage::List< storage::Pair< double, double> > &UNCLASSIFIED_RESULTS,
      const FunctionInterfaceSerializable< double, bool> &UNARY_CLASSIFIER
    )
    {
      // initialize the list to be returned
      storage::List< storage::Pair< double, bool> > classified_results;

      // iterate over provided results
      for
      (
        storage::List< storage::Pair< double, double> >::const_iterator
          result_itr( UNCLASSIFIED_RESULTS.Begin()), result_itr_end( UNCLASSIFIED_RESULTS.End());
        result_itr != result_itr_end; ++result_itr
      )
      {
        if( util::IsDefined( result_itr->Second()))
        {
          // push back the classified result
          classified_results.PushBack
          (
            storage::Pair< double, bool>( result_itr->First(), UNARY_CLASSIFIER( result_itr->Second()))
          );
        }
      }

      // end
      return classified_results;
    }

  } // namespace math
} // namespace bcl
