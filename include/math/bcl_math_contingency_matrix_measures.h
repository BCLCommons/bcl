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

#ifndef BCL_MATH_CONTINGENCY_MATRIX_MEASURES_H_
#define BCL_MATH_CONTINGENCY_MATRIX_MEASURES_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_contingency_matrix.h"
#include "function/bcl_function_member_const.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_function_interface_serializable.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ContingencyMatrixMeasures
    //! @brief a class which is allows choice of contingency matrix functions at runtime
    //!
    //! @see @link example_math_contingency_matrix_measures.cpp @endlink
    //! @author mendenjl
    //! @date Apr 08, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ContingencyMatrixMeasures :
      public util::FunctionInterfaceSerializable< ContingencyMatrix, double>
    {
    public:

      //! Enum of different values that can be retrieved by the contingency matrix
      enum Measure
      {
        //! Abbreviations:
        //! TP - True Positives
        //! FP - False Positives
        //! FN - False Negatives
        //! TN - True Negatives
        //! T  - Total = TP + FP + TN + FN
        //! P  - Positives = TP + FN
        //! N  - Negatives = FP + TN
        //! PP - Predicted Positives = TP + FP
        //! PN - Predicted Negatives = TN + FN
        e_TruePositiveRate,                //!< TPR = TP / P
        e_FalsePositiveRate,               //!< FPR = FP / P
        e_FalseNegativeRate,               //!< FNR = 1 - TPR = FN / N
        e_TrueNegativeRate,                //!< TNR = 1 - FPR = TN / N
        e_Accuracy,                        //!< (TP + TN) / T
        e_Specificity,                     //!< Synonym for TNR = TN / N
        e_Precision,                       //!< TP / PP
        e_Recall,                          //!< Synonym for TPR = TP / P
        e_PositivePredictiveValue,         //!< Synonym for precision = TP / PP
        e_IdealPPVRelativeToHitRate,       //!< Max possible precision = Min(1,P/PP)
        e_IdealPPVRelativeToFPR,           //!< P/(P+FP) for plotting against FP or FN sensitive measures
        e_NegativePredictiveValue,         //!< TN / PN
        e_FalseDiscoveryRate,              //!< 1 - Precision = FP / PP
        e_HitRate,                         //!< PP / T
        e_MatthewsCorrelationCoefficient,  //!< (TP * TN - FP * FN) / Sqrt( P * N * PP * PN)
        e_Enrichment,                      //!< Precision / ( P / T)
        e_InformationGain,                 //!< Information gain of the contingency matrix, minus IG in the Positive/Negative split
        e_Cutoff,                          //!< Cutoff used to generate contingency matrix measures
        e_LocalPPV,                        //!< Local-PPV. Best estimate of actual predictive value for predictions in a given range
        s_NumberMeasures
      };

      //! @brief Get the name, description, and function of the given measure
      //! @param MEASURE the measure of interest
      //! @return the short name or abbreviation for the class
      static const storage::Triplet< std::string, std::string, function::MemberConst< ContingencyMatrix, double> > &
        GetMeasureInfo( const Measure &MEASURE);

      Measure m_Measure; //!< Measure that this class calculates on the contingency matrix

      bool m_OptimizationDirection; //!< True if larger is better, false if smaller is better

      //! @brief Measure as string
      //! @param MEASURE the measure
      //! @return the string for the measure
      static const std::string &GetMeasureName( const Measure &MEASURE);

      //! @brief Initialization enum I/O helper
      typedef util::WrapperEnum< Measure, &GetMeasureName, s_NumberMeasures> MeasureEnum;

      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_Instances;

    public:

      //! @brief default constructor
      ContingencyMatrixMeasures() :
        m_Measure( s_NumberMeasures)
      {
      }

      //! @brief default constructor, undefined measure
      ContingencyMatrixMeasures( const Measure &MEASURE);

      //! @brief Clone function
      //! @return pointer to new ContingencyMatrixMeasures
      ContingencyMatrixMeasures *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the measure that this contingency matrix calculator calculates
      const Measure &GetMeasure() const;

      //! @brief determine the optimization direction for this particular measure
      //! @return true if larger values for this metric are better
      const bool &GetOptimizationParity() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief operator() taking a double based on the internal measure
      //! @param MATRIX contingency matrix of interest
      //! @return returns a double based on the selected measure
      double operator()( const ContingencyMatrix &MATRIX) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ContingencyMatrixMeasures

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_CONTINGENCY_MATRIX_MEASURES_H_
