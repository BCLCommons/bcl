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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "math/bcl_math_contingency_matrix.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_contingency_matrix.cpp
  //!
  //! @author mueller
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathContingencyMatrix :
    public ExampleInterface
  {
  public:

    ExampleMathContingencyMatrix *Clone() const
    {
      return new ExampleMathContingencyMatrix( *this);
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

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // empty contingency matrix
      math::ContingencyMatrix empty_matrix;

      // empty_matrix should be empty
      BCL_ExampleCheck( math::ContingencyMatrix().GetTotal(), 0);

      // constructing contingency matrix from properties
      math::ContingencyMatrix contingency_matrix( 30, 40, 10, 20);

    /////////////////
    // data access //
    /////////////////

      // check number of true positives
      BCL_ExampleCheck( contingency_matrix.GetNumberTruePositives(), 30);

      // check number of false positives
      BCL_ExampleCheck( contingency_matrix.GetNumberFalsePositives(), 40);

      // check number of false negatives
      BCL_ExampleCheck( contingency_matrix.GetNumberFalseNegatives(), 10);

      // check number of true negatives
      BCL_ExampleCheck( contingency_matrix.GetNumberTrueNegatives(), 20);

      // check number of actual positives
      BCL_ExampleCheck( contingency_matrix.GetNumberActualPositives(), 40);

      // check number of actual negatives
      BCL_ExampleCheck( contingency_matrix.GetNumberActualNegatives(), 60);

      // check number of predicted positives
      BCL_ExampleCheck( contingency_matrix.GetNumberPredictedPositives(), 70);

      // check number of predicted negatives
      BCL_ExampleCheck( contingency_matrix.GetNumberPredictedNegatives(), 30);

      // check total number
      BCL_ExampleCheck( contingency_matrix.GetTotal(), 100);

      // check true positive rate
      BCL_ExampleCheck( contingency_matrix.GetTruePositiveRate(), 0.75);

      // check false positive rate
      // note that we need to cast the numerator to double for MinGW, which otherwise considers floating point constants to be of float type
      BCL_ExampleCheck( contingency_matrix.GetFalsePositiveRate(), double( 40.0) / 60.0);

      // check accuracy
      BCL_ExampleCheck( contingency_matrix.GetAccuracy(), 0.5);

      // check specificity
      BCL_ExampleCheck( contingency_matrix.GetSpecificity(), double( 20.0) / 60.0);

      // check GetPrecision
      BCL_ExampleCheck( contingency_matrix.GetPrecision(), double( 3.0) / 7.0);

      // check GetPositivePredictiveValue
      BCL_ExampleCheck( contingency_matrix.GetPositivePredictiveValue(), double( 3.0) / 7.0);

      // check GetNegativePredictiveValue
      BCL_ExampleCheck( contingency_matrix.GetNegativePredictiveValue(), double( 2.0) / 3.0);

      // check GetFalseDiscoveryRate
      BCL_ExampleCheck( contingency_matrix.GetFalseDiscoveryRate(), double( 4.0) / 7.0);

      // check GetFractionPredictedPositives
      BCL_ExampleCheck( contingency_matrix.GetFractionPredictedPositives(), double( 7.0) / 10.0);

      // check GetMatthewsCorrelationCoefficient
      BCL_ExampleCheck( contingency_matrix.GetMatthewsCorrelationCoefficient(), double( 1.0) / math::Sqrt( double( 126)));

      // check GetEnrichment
      BCL_ExampleCheck( contingency_matrix.GetEnrichment(), double( 15.0) / 14.0);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathContingencyMatrix

  const ExampleClass::EnumType ExampleMathContingencyMatrix::s_Instance
  (
    GetExamples().AddEnum( ExampleMathContingencyMatrix())
  );

} // namespace bcl
