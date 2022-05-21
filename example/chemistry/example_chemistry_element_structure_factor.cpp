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
#include "chemistry/bcl_chemistry_element_structure_factor.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_sas_data_parameters.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_element_structure_factor.cpp
  //! @details Example class to test structure factor class
  //!
  //! @author putnamdk
  //! @date Jan 29, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryElementStructureFactor :
    public ExampleInterface
  {
  public:

    ExampleChemistryElementStructureFactor *Clone() const
    {
      return new ExampleChemistryElementStructureFactor( *this);
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
      // Crommer Mann Form Factors for the element Carbon
      const double cm_a[] = { 2.310, 1.020, 1.589, 0.865};
      const double cm_b[] = { 20.844, 10.208, 0.569, 51.651};
      const double cm_c( 0.216);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      chemistry::ElementStructureFactor structure_default;

      // construct from Carbon values
      chemistry::ElementStructureFactor structure_param
      (
        cm_a[ 0],
        cm_a[ 1],
        cm_a[ 2],
        cm_a[ 3],
        cm_b[ 0],
        cm_b[ 1],
        cm_b[ 2],
        cm_b[ 3],
        cm_c
      );

      // clone
      const util::ShPtr< math::FunctionInterfaceSerializable< restraint::SasDataParameters, double> > structure_cloned
      (
        structure_param.Clone()
      );
      BCL_ExampleCheck( structure_cloned.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // Example Check for default parameters
      BCL_ExampleCheck( structure_param.GetAValues()( 0), cm_a[0]);
      BCL_ExampleCheck( structure_param.GetAValues()( 1), cm_a[1]);
      BCL_ExampleCheck( structure_param.GetAValues()( 2), cm_a[2]);
      BCL_ExampleCheck( structure_param.GetAValues()( 3), cm_a[3]);
      BCL_ExampleCheck( structure_param.GetBValues()( 0), cm_b[0]);
      BCL_ExampleCheck( structure_param.GetBValues()( 1), cm_b[1]);
      BCL_ExampleCheck( structure_param.GetBValues()( 2), cm_b[2]);
      BCL_ExampleCheck( structure_param.GetBValues()( 3), cm_b[3]);
      BCL_ExampleCheck( structure_param.GetC(), cm_c);

      // Example Check for default constructor without parameters
      BCL_ExampleCheck( util::IsDefined( structure_default.GetAValues()( 0)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetAValues()( 1)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetAValues()( 2)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetAValues()( 3)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetBValues()( 0)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetBValues()( 1)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetBValues()( 2)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetBValues()( 3)), false);
      BCL_ExampleCheck( util::IsDefined( structure_default.GetC()), false);

    ///////////////
    // operators //
    ///////////////

      {
        restraint::SasDataParameters q_value( 0.0);
        // calculate f(x)
        // expected_form_factor is 6.0 with q of zero
        const double calculated_form_factor( structure_param( q_value));
        const double expected_form_factor( 6.0);
        BCL_ExampleIndirectCheckWithinTolerance
        (
          calculated_form_factor, expected_form_factor, 0.001,
          " incorrect result for q = 0"
        );
      }

      {
        // calculate f(x)
        // expected_form_factor is 5.10436 with q of 1.26
        restraint::SasDataParameters q_value2( 1.26);
        const double calculated_form_factor( structure_param( q_value2));
        const double expected_form_factor( 5.10436);
        BCL_ExampleIndirectCheckWithinTolerance
        (
          calculated_form_factor, expected_form_factor, 0.001,
          " incorrect result for q = 1.26"
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      WriteBCLObject( structure_param);
      ReadBCLObject( structure_default);
      BCL_ExampleIndirectCheck( structure_default.GetAValues(), structure_param.GetAValues(), "read and write");
      BCL_ExampleIndirectCheck( structure_default.GetBValues(), structure_param.GetBValues(), "read and write");
      BCL_ExampleIndirectCheck( structure_default.GetC(), structure_param.GetC(), "read and write");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryElementStructureFactor

  const ExampleClass::EnumType ExampleChemistryElementStructureFactor::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryElementStructureFactor())
  );

} // namespace bcl

