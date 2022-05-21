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
#include "chemistry/bcl_chemistry_bond_conformational.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_bond_conformational.cpp
  //!
  //! @author kothiwsk
  //! @date Dec 02, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryBondConformational :
    public ExampleInterface
  {
  public:

    ExampleChemistryBondConformational *Clone() const
    {
      return new ExampleChemistryBondConformational( *this);
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

      // construct from conformation interface
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation
      // create conformation to test bond_conformation
      chemistry::FragmentConformationShared conformation;
      conformation = sdf::FragmentFactory::MakeConformation( input_sdf);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

    /////////////////
    // data access //
    /////////////////

      // test GetBondType
      BCL_ExampleCheck
      (
        conformation.GetAtomsIterator()->GetBonds().FirstElement().GetBondType(),
        chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBondInRing
      );

      // test TargetAtom
      iterate::Generic< const chemistry::AtomConformationalInterface> itr( conformation.GetAtomsIterator());
      BCL_ExampleCheck
      (
        &conformation.GetAtomsIterator()->GetBonds().FirstElement().GetTargetAtom()
        == &*++itr,
        true
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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryBondConformational

  const ExampleClass::EnumType ExampleChemistryBondConformational::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryBondConformational())
  );

} // namespace bcl
