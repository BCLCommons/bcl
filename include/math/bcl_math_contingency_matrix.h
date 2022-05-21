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

#ifndef BCL_MATH_CONTINGENCY_MATRIX_H_
#define BCL_MATH_CONTINGENCY_MATRIX_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ContingencyMatrix
    //! @brief contingency matrix for analyzing a experiments accuracy, specificity and other measures
    //! @details for a better explanation see
    //! <a href="http://en.wikipedia.org/wiki/Receiver_operating_characteristic" target="blank">
    //! receiver operating characteristics</a>
    //!
    //! @see @link example_math_contingency_matrix.cpp @endlink
    //! @author woetzen, karakam
    //! @date 2008-07-21
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ContingencyMatrix :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //        actual value
      //          p  n   pred
      //         -------
      // p    p' |TP|FP| P'
      // r       |--|--|
      // e    n' |FN|TN| N'
      // d       -------
      // actual   P  N

      size_t m_TruePositives;
      size_t m_FalsePositives;
      size_t m_FalseNegatives;
      size_t m_TrueNegatives;

    public:

      //! @brief default constructor sets everything to 0
      ContingencyMatrix() :
        m_TruePositives( 0),
        m_FalsePositives( 0),
        m_FalseNegatives( 0),
        m_TrueNegatives( 0)
      {
      }

      //! @brief construct from all four values
      //! @param TRUE_POSITIVES number of positive predicted positive instances
      //! @param FALSE_POSITIVES number of positive predicted negative instances
      //! @param FALSE_NEGATIVES number of negative predicted positive instances
      //! @param TRUE_NEGATIVES number of negative predicted negative instances
      ContingencyMatrix
      (
        const size_t TRUE_POSITIVES,
        const size_t FALSE_POSITIVES,
        const size_t FALSE_NEGATIVES,
        const size_t TRUE_NEGATIVES
      ) :
        m_TruePositives( TRUE_POSITIVES),
        m_FalsePositives( FALSE_POSITIVES),
        m_FalseNegatives( FALSE_NEGATIVES),
        m_TrueNegatives( TRUE_NEGATIVES)
      {
      }

      //! @brief Clone function
      //! @return pointer to new ContingencyMatrix
      ContingencyMatrix *Clone() const
      {
        return new ContingencyMatrix( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get number of true positives
      //! @return number of true positives = TP
      size_t GetNumberTruePositives() const
      {
        return m_TruePositives;
      }

      //! @brief get number of false positives
      //! @return number of false positives = FP
      size_t GetNumberFalsePositives() const
      {
        return m_FalsePositives;
      }

      //! @brief get number of false negatives
      //! @return number of false negatives = FN
      size_t GetNumberFalseNegatives() const
      {
        return m_FalseNegatives;
      }

      //! @brief get number of true negatives
      //! @return number of true negatives = TN
      size_t GetNumberTrueNegatives() const
      {
        return m_TrueNegatives;
      }

      //! @brief number actual positives
      //! @return sum of TP and FN = P
      size_t GetNumberActualPositives() const
      {
        return m_TruePositives + m_FalseNegatives;
      }

      //! @brief number actual negatives
      //! @return sum of TN and FP = N
      size_t GetNumberActualNegatives() const
      {
        return m_TrueNegatives + m_FalsePositives;
      }

      //! @brief number predicted positives
      //! @return sum of TP and FP = P'
      size_t GetNumberPredictedPositives() const
      {
        return m_TruePositives + m_FalsePositives;
      }

      //! @brief number predicted negatives
      //! @return sum of TN and TN = N'
      size_t GetNumberPredictedNegatives() const
      {
        return m_FalseNegatives + m_TrueNegatives;
      }

      //! @brief return total
      //! @return TP + FP + FN + TN
      size_t GetTotal() const
      {
        return m_TruePositives + m_FalsePositives + m_FalseNegatives + m_TrueNegatives;
      }

      //! @return true positive rate
      //! @return TPR = TP/P=TP/(TP+FN)
      double GetTruePositiveRate() const
      {
        return double( m_TruePositives) / double( m_TruePositives + m_FalseNegatives);
      }

      //! @return false positive rate
      //! @return FPR = FP/N = FP/(FP+TN)
      double GetFalsePositiveRate() const
      {
        return double( m_FalsePositives) / double( m_FalsePositives + m_TrueNegatives);
      }

      //! @return false negative rate
      //! @return FNR = FN / P = FN/(TP+FN)
      double GetFalseNegativeRate() const
      {
        return double( m_FalseNegatives) / double( m_TruePositives + m_FalseNegatives);
      }

      //! @return false negative rate
      //! @return TNR = TN / N = TN/(FP+TN)
      double GetTrueNegativeRate() const
      {
        return double( m_TrueNegatives) / double( m_FalsePositives + m_TrueNegatives);
      }

      //! @brief get accuracy
      //! @return ACC = ( TP + TN) / ( P + N)
      double GetAccuracy() const
      {
        return double( m_TruePositives + m_TrueNegatives) /
               double( m_TruePositives + m_FalsePositives + m_FalseNegatives + m_TrueNegatives);
      }

      //! @brief get specificity
      //! @return SPC = TN / (FP + TN) = 1 − FPR
      double GetSpecificity() const
      {
        return double( m_TrueNegatives) / double( m_FalsePositives + m_TrueNegatives);
      }

      //! @brief return precision
      //! @return PC = PPV = TP / (TP + FP)
      double GetPrecision() const
      {
        return double( m_TruePositives) / double( m_TruePositives + m_FalsePositives);
      }

      //! @brief return normalized precision taking the recovered percentage rather the actual values into account
      //! @return norm PC = norm PPV = (TP/P) / (TP/P + FP/N)
      double GetPrecisionNormalized() const
      {
        return 2 * ( m_TruePositives / double( GetNumberActualPositives()))
        / double( m_TruePositives / double( GetNumberActualPositives()) + m_FalsePositives / double( GetNumberActualNegatives())) - double( 1);
      }

      //! @brief return recall (TPR - true positive rate)
      //! @return recall = TP / (TP + FN)
      double GetRecall() const
      {
        return GetTruePositiveRate();
      }

      //! @brief get positive predictive value
      //! @return PPV = PC = TP / (TP + FP)
      double GetPositivePredictiveValue() const
      {
        return GetPrecision();
      }

      //! @brief get ideal positive predictive value (e.g. value for a perfect predictor)
      //! @return Ideal PPV
      double GetIdealPositivePredictiveValue() const
      {
        // # predictions = TP + FP
        const double ideal_predicted_positives( m_TruePositives + m_FalsePositives);
        const double actual_positives( m_TruePositives + m_FalseNegatives);
        return ideal_predicted_positives < actual_positives ? 1.0 : actual_positives / ideal_predicted_positives;
      }

      //! @brief get ideal positive predictive value (e.g. value for a perfect predictor)
      //! @return Ideal PPV
      double GetIdealPositivePredictiveValueConsideringFalsePositives() const
      {
        // # predictions = TP + FP
        const double actual_positives( m_TruePositives + m_FalseNegatives);
        return double( actual_positives) / double( actual_positives + m_FalsePositives);
      }

      //! @brief get positive predictive value
      //! @return NPV = TN / (TN + FN)
      double GetNegativePredictiveValue() const
      {
        return double( m_TrueNegatives) / double( std::max( size_t( 1), m_TrueNegatives + m_FalseNegatives));
      }

      //! @brief get false discovery rate
      //! @return FDR = FP / (FP + TP)
      double GetFalseDiscoveryRate() const
      {
        return double( m_FalsePositives) / double( std::max( size_t( 1), m_FalsePositives + m_TruePositives));
      }

      //! @brief get fraction positive predicted
      //! @return FPP = (TP + FP) / (P + N)
      double GetFractionPredictedPositives() const
      {
        return double( GetNumberPredictedPositives()) / double( GetTotal());
      }

      //! @brief get fraction positive predicted normalized
      //! @return norm FPP = (TP/P + FP/N) / (P/P + N/N)
      double GetFractionPredictedPositivesNormalized() const
      {
        return ( m_TruePositives / double( GetNumberActualPositives()) + m_FalsePositives / double( GetNumberActualNegatives())) / double( 2);
      }

      //! @brief get Matthews correlation coefficient
      //! @return MCC = ( TP * TN - FP * FN) / sqrt( P * N * P' * N')
      double GetMatthewsCorrelationCoefficient() const
      {
        // convert all the variables to double.  With the limited precision of size_t's, it is possible for the products
        // in the MCC equation to overflow the range of size_t's
        const double tp( m_TruePositives), fp( m_FalsePositives), tn( m_TrueNegatives), fn( m_FalseNegatives);
        // handle 0-denominator case using max; typically this happens if there are no predicted positives
        // denominator == 0 -> numerator == zero, and taking limits, the correct value for the MCC is 0 in that case
        const double mcc
        (
          ( tp * tn - fp * fn)
          /
          Sqrt( std::max( ( tp + fp) * ( tp + fn) * ( tn + fp) * ( tn + fn), double( 1)))
        );
        // the max() is needed in case any of the terms in the de
        return mcc;
      }

      //! @brief get enrichment
      //! is the ratio between the true predicted positives vs all predicted positive relative
      //! to the original ratio of actual positives vs. the total
      //! @return ER = TP / (TP + FP) * ( P + N) / P
      double GetEnrichment() const
      {
        return double( m_TruePositives) / double( m_TruePositives + m_FalsePositives) *
               double( GetTotal()) / double( GetNumberActualPositives());
      }

      //! @brief get information gain from the contingency matrix
      //! @return the information gain of the contingency matrix
      double GetInformationGainRatio() const
      {
        const double total_values( GetTotal());
        const double true_positive_fraction( double( m_TruePositives) / total_values);
        const double false_positive_fraction( double( m_FalsePositives) / total_values);
        const double false_negative_fraction( double( m_FalseNegatives) / total_values);
        const double true_negative_fraction( double( m_TrueNegatives) / total_values);

        // Compute entropy before split
        const double pre_split
        (
            PLogP( true_positive_fraction + false_negative_fraction)
          + PLogP( false_positive_fraction + true_negative_fraction)
        );

        // Compute entropy of the split
        const double split_entropy
        (
            PLogP( true_positive_fraction + false_positive_fraction)
          + PLogP( false_negative_fraction + true_negative_fraction)
        );

        // split was already perfectly ordered, no information gain
        if( Absolute( split_entropy) <= std::numeric_limits< double>::epsilon())
        {
          return 0.0;
        }

        // compute entropy after the split
        const double matrix_entropy
        (
          PLogP( true_positive_fraction) + PLogP( false_positive_fraction)
          + PLogP( false_negative_fraction) + PLogP( true_negative_fraction)
        );

        // compute information gain ratio
        return ( pre_split + split_entropy - matrix_entropy) / split_entropy;
      }

      //! @brief get enrichment
      //! is the ratio between the true predicted positives vs all predicted positive relative
      //! to the original ratio of actual positives vs. the total
      //! @return ER = TP / (TP + FP) * ( P + N) / P
      double GetOversampledEnrichment( const double &ENRICHMENT_FACTOR) const
      {
        return double
        (
          m_TruePositives / ENRICHMENT_FACTOR) / double( ( m_TruePositives / ENRICHMENT_FACTOR) + m_FalsePositives) *
          double
          (
            ( m_TruePositives / ENRICHMENT_FACTOR) + m_FalsePositives + ( m_FalseNegatives / ENRICHMENT_FACTOR)
            + m_TrueNegatives
          )
          / double( ( m_TruePositives / ENRICHMENT_FACTOR) + ( m_FalseNegatives / ENRICHMENT_FACTOR));
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_TruePositives, ISTREAM);
        io::Serialize::Read( m_FalsePositives, ISTREAM);
        io::Serialize::Read( m_FalseNegatives, ISTREAM);
        io::Serialize::Read( m_TrueNegatives, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_TruePositives, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_FalsePositives, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_FalseNegatives, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_TrueNegatives, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! @brief helper function for information gain ratio
      //! @param VALUE probability of a given event
      //! @return p log p, unless p is 0, then 0
      double PLogP( const double &VALUE) const
      {
        static const double log2_conversion( double( 1.0) / log( double( 2.0)));
        return VALUE > 0.0 ? VALUE * log( VALUE) * log2_conversion : 0.0;
      }

    }; // class ContingencyMatrix

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_CONTINGENCY_MATRIX_H_
