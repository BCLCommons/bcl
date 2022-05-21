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

#ifndef BCL_MATH_ROC_CURVE_H_
#define BCL_MATH_ROC_CURVE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "bcl_math_range.h"
#include "function/bcl_function_member_const.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ROCCurve
    //! @brief This class handles storing results in a roc curve
    //! @details This class works on a list of test results where a prediction can be coupled with a real value
    //!
    //! @see @link example_math_roc_curve.cpp @endlink
    //! @author karakam, woetzen, butkiem1
    //! @date Mar 14, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ROCCurve :
      public util::ObjectInterface
    {

    public:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Point
      //! @brief A trivial data class used for storing information needed to construct the roc curve; used in favor of
      //!        storage::Triplet for readability
      //!
      //! @see @link example_math_roc_curve.cpp @endlink
      //! @author mendenjl
      //! @date May 14, 2013
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class Point :
        public util::ObjectInterface
      {
      private:

      //////////
      // data //
      //////////

        size_t m_FalsePositives;  //!< # of false positives
        size_t m_TruePositives;   //!< # of correct predictions on the desired side of the cutoff (most commonly above the cutoff)
        double m_Cutoff;          //!< actual cutoff value

      public:
      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief constructor from members
        //! @param FP # of false positives
        //! @param TP # of true positives
        //! @param CUTOFF actual cutoff value
        Point( const size_t &FP = 0, const size_t &TP = 0, const double &CUTOFF = 0.0);

        //! @brief Clone function
        //! @return pointer to new ROCCurve
        Point *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name of the object behind a pointer or the current object
        //! @return the class name
        const std::string &GetClassIdentifier() const;

        //! @brief get number of true positives
        //! @return number of true positives = TP
        size_t GetNumberTruePositives() const;

        //! @brief get number of false positives
        //! @return number of false positives = FP
        size_t GetNumberFalsePositives() const;

        //! @brief number predicted positives
        //! @return sum of TP and FP = P'
        size_t GetNumberPredictedPositives() const;

        //! @brief get the cutoff
        //! @return the cutoff
        double GetCutoff() const;

      ////////////////
      // operations //
      ////////////////

        //! @brief add true positives to the point
        //! @param RESULT true for true positive, false for false positive
        //! @param CUTOFF the new cutoff
        void Update( const bool &RESULT, const double &CUTOFF);

        //! @brief Construct a contingency matrix, given the final point in the ROC
        //! @param FINAL_POINT the final point in the ROC
        //! @return a contingency matrix for this object
        ContingencyMatrix GetContingencyMatrix( const Point &FINAL_POINT) const;

      protected:

        //! @brief read from std::istream
        //! @param ISTREAM input stream
        //! @return istream which was read from
        std::istream &Read( std::istream &ISTREAM);

        //! @brief write to std::ostream
        //! @param OSTREAM output stream to write to
        //! @return output stream which was written to
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      };

    private:

    //////////
    // data //
    //////////

      //! keeps track of the total count for true positives vs. true negatives
      storage::Vector< Point> m_SortedCounts;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ROCCurve() :
        m_SortedCounts()
      {
      }

      //! @brief constructor from a list of pairs of prediction and boolean classifications of the expected outputs
      //! @param TEST_RESULTS_CLASSIFIED list of test results with their according boolean classification
      ROCCurve( const storage::List< storage::Pair< double, bool> > &TEST_RESULTS_CLASSIFIED) :
        m_SortedCounts()
      {
        Initialize( TEST_RESULTS_CLASSIFIED);
      }

      //! @brief constructor from a list of pairs of prediction and test results not yet classified by a threshold
      //! @param TEST_RESULTS list of pairs of prediction ( double) and test outputs( double)
      //! @param THRESHOLD the threshold value to be used for classifying the outputs
      //! ( classification means that results above the threshold are considered positive = 1)
      ROCCurve
      (
        const storage::List< storage::Pair< double, double> > &TEST_RESULTS,
        const double THRESHOLD
      );

      //! @brief constructor from a list of pairs of prediction and test results not yet classified by a threshold
      //! @param TEST_RESULTS list of pairs of prediction ( double) and test outputs( double)
      //! @param THRESHOLD the threshold value to be used for classifying the outputs
      //! @param POSITIVES_ABOVE_THRESHOLD flag whether positives are above threshold
      //! ( classification means that results above the threshold are considered positive = 1)
      ROCCurve
      (
        storage::List< storage::Pair< double, double> > &TEST_RESULTS,
        const double THRESHOLD,
        bool POSITIVES_ABOVE_THRESHOLD
      );

      //! @brief constructor from a list of pairs of prediction and test results not yet classified by a threshold
      //! @param TEST_RESULTS list of pairs of prediction ( double) and test outputs( double)
      //! @param UNARY_CLASSIFIER classifier to be used for classifying expected outputs
      //! ( classification means that results above the threshold are considered positive = 1)
      ROCCurve
      (
        const storage::List< storage::Pair< double, double> > &TEST_RESULTS,
        const FunctionInterfaceSerializable< double, bool> &UNARY_CLASSIFIER
      ) :
        m_SortedCounts()
      {
        Initialize( ClassifyResults( TEST_RESULTS, UNARY_CLASSIFIER));
      }

      //! @brief Clone function
      //! @return pointer to new ROCCurve
      ROCCurve *Clone() const
      {
        return new ROCCurve( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get sorted counts
      //! @return sorted counts
      const storage::Vector< Point> &GetSortedCounts() const
      {
        return m_SortedCounts;
      }

      //! @brief returns the total number of false positives
      //! @return total number of false positives
      const size_t GetNumberActualNegatives() const;

      //! @brief returns the total number of true positives
      //! @return total number of true positives
      const size_t GetNumberActualPositives() const;

      //! @brief returns the total number of results
      //! return the total number of results
      const size_t GetNumberResults() const
      {
        return GetNumberActualNegatives() + GetNumberActualPositives();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the ROC curve with a list of pair of predicted values and their classifications
      //! @param RESULTS_CLASSIFIED combination of predicted values with classification of true versus false positive
      void Initialize( const storage::List< storage::Pair< double, bool> > &RESULTS_CLASSIFIED);

      //! @brief calculates and returns the area under the curve
      //! @return curve integral
      double Integral() const;

      //! @brief calculates and returns the weighted area under the curve by weighting function
      //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
      //! @return entire weighted curve integral
      double Integral
      (
        const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
        const bool ROC_X_AXIS_LOG10_SCALING = false
      ) const;

      //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
      //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
      //! @return curve integral
      double Integral( const double FRACTION) const;

      //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
      //!        weighted by a given function
      //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
      //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
      //!        accordingly
      //! @param ROC_X_AXIS_LOG10_SCALING whether to use log10 scaling on the x-axis for integral computation
      //! @return weighted curve integral
      double Integral
      (
        const double FRACTION,
        const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
        const bool ROC_X_AXIS_LOG10_SCALING = false
      ) const;

      //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
      //!        the contingency matrix functions to plot on x- and y-axis can be specified
      //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
      //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
      //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
      //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
      //! @return curve integral
      double Integral
      (
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
        const bool ROC_X_AXIS_LOG10_SCALING = bool( false)
      ) const;

      //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
      //!        the contingency matrix functions to plot on x- and y-axis can be specified
      //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
      //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix
      //! @param ROC_Y_AXIS function pointer to a member function in math::ContingencyMatrix
      //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
      //! @return curve integral
      double Integral
      (
        const double FRACTION,
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
        const bool ROC_X_AXIS_LOG10_SCALING = bool( false)
      ) const;

      //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1),
      //!        the contingency matrix functions to plot on x- and y-axis can be specified
      //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
      //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
      //! @param ROC_X_AXIS_FRACTION_CUTOFF x axis value cutoff rather then specifying a number of counts
      //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
      //! @return curve integral
      double Integral
      (
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
        const Range< double> &ROC_X_AXIS_FRACTION_CUTOFF = Range< double>( 0.0, 1.0),
        const bool ROC_X_AXIS_LOG10_SCALING = bool( false)
      ) const;

      //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
      //!        weighted by a given function
      //!        the contingency matrix functions to plot on x- and y-axis can be specified
      //! @param FRACTION top fraction ( between 0 and 1) to be considered when calculating integral
      //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
      //!        accordingly
      //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
      //! @param ROC_X_AXIS function pointer to a member function in math::ContingencyMatrix

      //! @return weighted curve integral
      double Integral
      (
        const double FRACTION,
        const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
        const bool ROC_X_AXIS_LOG10_SCALING = bool( false)
      ) const;

      //! @brief calculates and returns the area under the curve for the given top fraction ( between 0 and 1)
      //!        weighted by a given function
      //!        the contingency matrix functions to plot on x- and y-axis can be specified
      //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
      //!        accordingly
      //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
      //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
      //! @param ROC_X_AXIS_FRACTION_CUTOFF x axis value cutoff rather then specifying a number of counts through param END
      //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
      //! @return weighted curve integral
      double Integral
      (
        const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
        const Range< double> &ROC_X_AXIS_FRACTION_CUTOFF = Range< double>( 0.0, 1.0),
        const bool ROC_X_AXIS_LOG10_SCALING = bool( false)
      ) const;

      //! @brief calculates ContingencyMatrix with true/false positive/negatives
      //! @param FRACTION consider the first fraction as being predicted positive, the remainder predicted negative
      ContingencyMatrix ContingencyMatrixFraction( const double FRACTION) const;

      //! @brief returns the cutoff at which the given FRACTION would be considered positive
      double CutoffFraction( const double FRACTION) const;

      //! @brief Convert into a map from threshold to contingency matrix
      storage::Map< double, ContingencyMatrix> ToMap() const;

      //! @brief get the optima for a particular contingency matrix measure
      //! @param MEASURE the measure to optimize
      //! @return an iterator to the maximum point for the particular measure and the maximum value
      std::pair< storage::Vector< Point>::const_iterator, double> GetMaxima
      (
        const util::FunctionInterfaceSerializable< ContingencyMatrix, double> &MEASURE
      ) const;

    private:

      //! @brief calculates and returns the area under the curve for the counts between provided regions
      //! @param END iterator to the end of the interval for which the integral is going to be calculated
      //! @return curve integral for values between specified regions
      double Integral
      (
        const storage::Vector< Point>::const_iterator &END
      ) const;

      //! @brief calculates and returns the area under the curve for the counts between provided regions
      //!        weighted with a mathematical function
      //! @param END iterator to the end of the interval for which the integral is going to be calculated
      //! @param WEIGHTING_FUNCTION FunctionInterface which weights the ROC curve and allows Integral calculation
      //!        accordingly
      //! @param ROC_X_AXIS_LOG10_SCALING enables log10 scaling of x axis
      //! @return weighted curve integral for values between specified regions
      double Integral
      (
        const storage::Vector< Point>::const_iterator &END,
        const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
        const bool ROC_X_AXIS_LOG10_SCALING = bool( false)
      ) const;

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
      //! @param IS_IDEAL flag indicating whether the returned mapping refers to an ideal (perfect) x-y curve or
      //!        takes the actual values determined by ROC_X_AXIS and ROC_Y_AXIS into account
      //! @return weighted curve integral for values between specified regions
      double Integral
      (
        const storage::Vector< Point>::const_iterator &END,
        const FunctionInterfaceSerializable< double, double> &WEIGHTING_FUNCTION,
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
        const Range< double> ROC_X_AXIS_RANGE = Range< double>( 0.0, 1.0),
        const bool ROC_X_AXIS_LOG10_SCALING = bool( false)
      ) const;

      //! @brief compute the x axis and y axis mapping for computing an interpolated integral under the curve
      //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix
      //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix
      //! @param IS_IDEAL flag indicating whether the returned mapping refers to an ideal (perfect) x-y curve or
      //!        takes the actual values determined by ROC_X_AXIS and ROC_Y_AXIS into account
      //! @return x axis and y axis mapping
      storage::Map< double, RunningAverage< double> > ComputeXAxisYAxisMapping
      (
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS
      ) const;

    public:

      //! @brief get the localized PPV curve
      //! The localized PPV-curve yields a more precise estimate of the local PPVs for a given range
      //! The standard ROC-curve answers the question: what is the PPV of all points predicted at or above a given
      //! cutoff (assuming +parity). The localized ppv-curve, by contrast, provides the PPV of values near a particular
      //! cutoff. This curve is produced by iteratively finding the optimal PPV interest, then removing all
      //! entries (P/N) at or above/below (depending on parity) from the ROC curve. PPV is essentially the only metric
      //! (other than the obscure false discovery rate = 1-PPV) that is independent of negatives and the total number of
      //! positives, and so is the only metric that can be computed locally in this manner.
      //! @return a piecewise function for the PPV curve, which is guaranteed to be monotonic
      PiecewiseFunction GetLocalPPVCurve() const;

      //! @brief creates a thinned version of the counts so that only every Nth( PERIODICITY * n) element is kept
      //! @brief This function is to be used for plotting.
      //! @param PERIODICITY periodicity with which elements should be selected
      //! @return thinned ROC curve counts
      ROCCurve
      GetThinnedRocCurvePeriodicity( const size_t PERIODICITY) const;

      //! @brief creates a thinned version of the counts so that only TOTAL_NUMBER_ENTRIES elements are kept
      //! @brief This function is to be used for plotting.
      //! @param TOTAL_NUMBER_ENTRIES total number of elements to be selected
      //! @return thinned ROC curve counts
      ROCCurve
      GetThinnedRocCurveTotal( const size_t TOTAL_NUMBER_ENTRIES) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream from a formatted output
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadPlottingTable( std::istream &ISTREAM);

      //! @brief write to std::ostream in a formatted output style
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &WritePlottingTable( std::ostream &OSTREAM) const;

      //! @brief write to std::ostream in a formatted output style
      //! @param OSTREAM output stream to write to
      //! @param FORMAT format the rate
      //! @return output stream which was written to
      std::ostream &WriteRatePlottingTable
      (
        std::ostream &OSTREAM,
        const util::Format &FORMAT = util::Format().FFP( 3).W( 5).ForceW()
      ) const;

      //! @brief write to std::ostream in a formatted output style the result between BEGIN and END
      //! @param OSTREAM output stream to write to
      //! @param ROC_X_AXIS function pointer to a member function in ContingencyMatrix, specifies what is plotted
      //!        on x axis of roc curve
      //! @param ROC_Y_AXIS function pointer to a member function in ContingencyMatrix, specifies what is plotted
      //!        on y axis of roc curve
      //! @param ROC_X_AXIS_FRACTION_CUTOFF x axis value cutoff specifying the x axis value when to stop plotting
      //!        default is set to 1.0 or 100%
      //! @param FORMAT format the rate
      //! @param PLOTTING_TABLE_TITLE title that describes columns for x axis values and y axis values
      //! @return output stream which was written to
      std::ostream &WriteRatePlottingTableGeneric
      (
        std::ostream &OSTREAM,
        const function::MemberConst< ContingencyMatrix, double> &ROC_X_AXIS,
        const function::MemberConst< ContingencyMatrix, double> &ROC_Y_AXIS,
        const Range< double> &ROC_X_AXIS_RANGE = Range< double>( 0.0, 1.0),
        const util::Format &FORMAT = util::Format().W( 5).ForceW(),
        const std::string &PLOTTING_TABLE_TITLE = std::string( "")
      ) const;

      //! @brief write all common contingency matrix measures to the file, return a vector of strings for each measure that was written
      //! @param DATA_FILENAME file name for the table
      //! @param MEASURES a list of measures of interest
      //! @return strings of all measures that were written
      storage::Vector< std::string> WriteRatePlottingTableComplete
      (
        const std::string &DATA_FILENAME,
        const storage::Vector< std::string> &MEASURES
      ) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief returns iterator on the the end of the specified fraction of the list
      //! @param FRACTION top fraction ( between 0 and 1) to be used in consequent operations
      //! @return iterator to ending of the requested interval
      storage::Vector< Point>::const_iterator GetEndIteratorForInterval( const double FRACTION) const;

      //! @brief static function to classify provided results and return the classified results
      //! @param UNCLASSIFIED_RESULTS results that are not yet classified
      //! @param UNARY_CLASSIFIER classifier to be used for classifying expected outputs
      //! @return classified results
      static storage::List< storage::Pair< double, bool> > ClassifyResults
      (
        const storage::List< storage::Pair< double, double> > &UNCLASSIFIED_RESULTS,
        const FunctionInterfaceSerializable< double, bool> &UNARY_CLASSIFIER
      );

    }; // class ROCCurve

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_ROC_CURVE_H_
