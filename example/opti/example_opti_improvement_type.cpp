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
#include "opti/bcl_opti_improvement_type.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_improvement_type.cpp
  //! @brief The math angle class converts between degrees and radians and shows the unit name
  //!
  //! @author mendenjl
  //! @date Aug 22, 2013
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiImprovementType :
    public ExampleInterface
  {
  public:

    ExampleOptiImprovementType *Clone() const
    {
      return new ExampleOptiImprovementType( *this);
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

      // no constructor

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleCheck( opti::GetImprovementTypeName( opti::e_SmallerIsBetter), "SmallerIsBetter");
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.5, opti::e_SmallerIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.4, opti::e_SmallerIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( 0.4, 0.5, opti::e_SmallerIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.5, opti::e_SmallerEqualIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.4, opti::e_SmallerEqualIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( 0.4, 0.5, opti::e_SmallerEqualIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.5, opti::e_LargerIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.4, opti::e_LargerIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( 0.4, 0.5, opti::e_LargerIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.5, opti::e_LargerEqualIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( 0.5, 0.5, opti::e_SmallerAbsIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( -0.5, 0.5, opti::e_SmallerAbsIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( 1.0, 0.5, opti::e_SmallerAbsIsBetter), false);
      BCL_ExampleCheck( opti::DoesImprove( 0.4, 0.5, opti::e_SmallerAbsIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( 0.4, -0.5, opti::e_SmallerAbsIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( -0.4, 0.5, opti::e_SmallerAbsIsBetter), true);
      BCL_ExampleCheck( opti::DoesImprove( -0.4, -0.5, opti::e_SmallerAbsIsBetter), true);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOptiImprovementType

  const ExampleClass::EnumType ExampleOptiImprovementType::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiImprovementType())
  );

} // namespace bcl
