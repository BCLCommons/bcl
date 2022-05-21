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
#include "math/bcl_math_contingency_matrix_measures.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // separate, anonymous namespace to prevent exporting symbols
    namespace
    {
      // add each of the possible instances to the enumerated instances
      util::ObjectInterface *AddInstances()
      {
        // keep a pointer to the last created instance
        util::ObjectInterface *last_instance( NULL);
        for( size_t measure( 0); measure < ContingencyMatrixMeasures::s_NumberMeasures; ++measure)
        {
          last_instance =
            util::Enumerated< util::FunctionInterfaceSerializable< ContingencyMatrix, double> >::AddInstance
            (
              new ContingencyMatrixMeasures( static_cast< ContingencyMatrixMeasures::Measure>( measure))
            );
        }
        return last_instance;
      }
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ContingencyMatrixMeasures::s_Instances( AddInstances());

    //! @brief Get the name, description, and function of the given measure
    //! @param MEASURE the measure of interest
    //! @return the short name or abbreviation for the class
    const storage::Triplet< std::string, std::string, function::MemberConst< ContingencyMatrix, double> > &
      ContingencyMatrixMeasures::GetMeasureInfo( const Measure &MEASURE)
    {
      typedef storage::Triplet< std::string, std::string, function::MemberConst< ContingencyMatrix, double> > t_Info;
      static const t_Info s_info[ s_NumberMeasures + 1] =
      {
        t_Info( "TPR", "True positive rate = TP / P", &math::ContingencyMatrix::GetTruePositiveRate),
        t_Info( "FPR", "False positive rate = FP / P", &math::ContingencyMatrix::GetFalsePositiveRate),
        t_Info( "FNR", "False negative rate = FN / N", &math::ContingencyMatrix::GetFalseNegativeRate),
        t_Info( "TNR", "True negative rate = TN / N", &math::ContingencyMatrix::GetTrueNegativeRate),
        t_Info( "Accuracy", "(TP + TN) / T", &math::ContingencyMatrix::GetAccuracy),
        t_Info( "Specificity", "= TNR = TN / N", &math::ContingencyMatrix::GetTrueNegativeRate),
        t_Info( "Precision", "= TP / PP (predicted positives)", &math::ContingencyMatrix::GetPrecision),
        t_Info( "Recall", "= TPR = TP / P", &math::ContingencyMatrix::GetRecall),
        t_Info( "PPV", "Positive predictive value = Precision = TP / PP (predicted positives)", &math::ContingencyMatrix::GetPositivePredictiveValue),
        t_Info( "Ideal-PPV", "Max possible PPV | precision = Min(1,P/PP).  Not appropriate for plotting against measures that consider FP or FN!", &math::ContingencyMatrix::GetIdealPositivePredictiveValue),
        t_Info( "Ideal-PPV_FPRelative", "P/(P+FP).  Appropriate for plotting against measures that consider FP or FN!", &math::ContingencyMatrix::GetIdealPositivePredictiveValueConsideringFalsePositives),
        t_Info( "NPV", "Negative predictive value = TN / PN (Predicted negatives)", &math::ContingencyMatrix::GetNegativePredictiveValue),
        t_Info( "FDR", "False discovery rate = 1 - precision = FP / PP", &math::ContingencyMatrix::GetFalseDiscoveryRate),
        t_Info( "HitRate", "fraction predicted positive = PP / T", &math::ContingencyMatrix::GetFractionPredictedPositives),
        t_Info( "MCC", "Matthews Correlation Coefficient = (TP * TN - FP * FN) / Sqrt( P * N * PP * PN)", &math::ContingencyMatrix::GetMatthewsCorrelationCoefficient),
        t_Info( "Enrichment", "Ratio of precision to positive rate = Precision / ( P / T)", &math::ContingencyMatrix::GetEnrichment),
        t_Info( "InformationGainRatio", "Information gain ratio, accounts for Positive/Negative split (range 0-1)", &math::ContingencyMatrix::GetInformationGainRatio),
        t_Info( "Cutoff", "The cutoff value at which the contingency matrix was created", NULL),
        t_Info( "LocalPPV", "Local-PPV. Best estimate of actual predictive value for predictions in a given range", NULL),
        t_Info( "", "", &math::ContingencyMatrix::GetTruePositiveRate)
      };
      return s_info[ MEASURE];
    };

    //! @brief Measure as string
    //! @param MEASURE the measure
    //! @return the string for the measure
    const std::string &ContingencyMatrixMeasures::GetMeasureName( const Measure &MEASURE)
    {
      return GetMeasureInfo( MEASURE).First();
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor sets everything to 0
    ContingencyMatrixMeasures::ContingencyMatrixMeasures( const Measure &MEASURE) :
      m_Measure( MEASURE),
      m_OptimizationDirection( true)
    {
      function::MemberConst< ContingencyMatrix, double> function( GetMeasureInfo( m_Measure).Third());
      if( function.IsDefined())
      {
        ContingencyMatrix good( 99, 1, 1, 99), bad( 1, 99, 99, 1);
        m_OptimizationDirection = function( good) > function( bad);
      }
    }

    //! @brief virtual copy constructor
    ContingencyMatrixMeasures *ContingencyMatrixMeasures::Clone() const
    {
      return new ContingencyMatrixMeasures( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ContingencyMatrixMeasures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the measure that this contingency matrix calculator calculates
    const ContingencyMatrixMeasures::Measure &ContingencyMatrixMeasures::GetMeasure() const
    {
      return m_Measure;
    }

    //! @brief determine the optimization direction for this particular measure
    //! @return true if larger values for this metric are better
    const bool &ContingencyMatrixMeasures::GetOptimizationParity() const
    {
      return m_OptimizationDirection;
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &ContingencyMatrixMeasures::GetAlias() const
    {
      return GetMeasureInfo( m_Measure).First();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking a double based on the internal measure
    //! @param MATRIX contingency matrix of interest
    //! @return returns a double based on the selected measure
    double ContingencyMatrixMeasures::operator()( const ContingencyMatrix &MATRIX) const
    {
      return GetMeasureInfo( m_Measure).Third()( MATRIX);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ContingencyMatrixMeasures::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( GetMeasureInfo( m_Measure).Second());
      return serializer;
    }

  } // namespace math
} // namespace bcl
