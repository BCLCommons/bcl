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
#include "restraint/bcl_restraint_rdc.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"
#include "restraint/bcl_restraint_rdc_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_rdc.cpp
  //!
  //! @author weinerbe
  //! @date Feb 16, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintRDC :
    public ExampleInterface
  {
  public:

    ExampleRestraintRDC *Clone() const
    {
      return new ExampleRestraintRDC( *this);
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
      // read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      restraint::RDC def_construct;

    ////////////////
    // operations //
    ////////////////

      // pushback some restraints
      def_construct.PushBack
      (
        restraint::LocatorCoordinatesHydrogen( 'A', 12, "CA"),
        restraint::LocatorCoordinatesHydrogen( 'A', 12, "HA"),
        1.09,
        -2.5
      );
      def_construct.PushBack
      (
        restraint::LocatorCoordinatesHydrogen( 'A', 12, "N"),
        restraint::LocatorCoordinatesHydrogen( 'A', 12, "H"),
        1.03,
        1.0
      );

      // test NormalizetoNH
      restraint::RDC normalized( def_construct);
      normalized.NormalizetoNH();

      // test AdjustSigns
      restraint::RDC signs( def_construct);
      signs.AdjustSigns();

      // test GenerateAssignment without adjustments
      const restraint::RDCAssignment def_assignment( def_construct.GenerateAssignment( protein_model));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        def_assignment.GetData().FirstElement().Third(), -2.5, 0.001,
        "GenerateAssignment without adjustments"
      );

      // test GenerateAssignment with normalization
      const restraint::RDCAssignment norm_assignment( normalized.GenerateAssignment( protein_model));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        norm_assignment.GetData().FirstElement().Third(), 1.19446, 0.001,
        "GenerateAssignment with normalization"
      );

      // test GenerateAssignment with signs adjusted
      const restraint::RDCAssignment sign_assignment( signs.GenerateAssignment( protein_model));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        sign_assignment.GetData().FirstElement().Third(), 2.5, 0.001,
        "GenerateAssignment with signs adjusted"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( def_construct);
      restraint::RDC read_construct;
      ReadBCLObject( read_construct);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintRDC

  const ExampleClass::EnumType ExampleRestraintRDC::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintRDC())
  );

} // namespace bcl
