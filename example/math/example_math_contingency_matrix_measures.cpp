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
#include "math/bcl_math_contingency_matrix_measures.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_contingency_matrix_measures.cpp
  //!
  //! @author mendenjl
  //! @date Apr 08, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathContingencyMatrixMeasures :
    public ExampleInterface
  {
  public:

    ExampleMathContingencyMatrixMeasures *Clone() const
    {
      return new ExampleMathContingencyMatrixMeasures( *this);
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

      // create
      math::ContingencyMatrixMeasures infogain_measurer( math::ContingencyMatrixMeasures::e_InformationGain);

      // Check the measure type
      BCL_ExampleCheck( infogain_measurer.GetAlias(), "InformationGainRatio");

      // constructing contingency matrix from properties
      math::ContingencyMatrix contingency_matrix( 30, 40, 10, 20);

    /////////////////
    // data access //
    /////////////////

      // check the values returned
      BCL_ExampleCheckWithinTolerance
      (
        infogain_measurer( contingency_matrix),
        contingency_matrix.GetInformationGainRatio(),
        1.0e-5
      );

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

  }; //end ExampleMathContingencyMatrixMeasures

  const ExampleClass::EnumType ExampleMathContingencyMatrixMeasures::s_Instance
  (
    GetExamples().AddEnum( ExampleMathContingencyMatrixMeasures())
  );

} // namespace bcl
