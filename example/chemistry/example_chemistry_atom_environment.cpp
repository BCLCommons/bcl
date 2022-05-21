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
#include "chemistry/bcl_chemistry_atom_environment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_environment.cpp
  //!
  //! @author mueller
  //! @date   08/24/2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomEnvironment :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomEnvironment *Clone() const
    {
      return new ExampleChemistryAtomEnvironment( *this);
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
      const std::string taxol_sdf( AddExampleInputPathToFilename( ExampleInterface::e_Chemistry, "taxol.sdf"));

      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, taxol_sdf);
      chemistry::FragmentConformationShared taxol
        = sdf::FragmentFactory::MakeConformation( sdf::MdlHandler( input));
      io::File::CloseClearFStream( input);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor:
      chemistry::AtomEnvironment default_environment;

      // check that a default-constructed atom-environment is empty
      BCL_ExampleCheck( chemistry::AtomEnvironment().GetEnvironmentAtoms().IsEmpty(), true);

      iterate::Generic< const chemistry::AtomConformationalInterface> itr( taxol.GetAtomsIterator());
      itr.GotoPosition( 3);
      // check that the environment around taxol contains 8 atoms
      BCL_ExampleCheck
      (
        chemistry::AtomEnvironment( *itr, 3).GetEnvironmentAtoms().GetSize(),
        8
      );

    /////////////////
    // data access //
    /////////////////

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

  }; //end ExampleChemistryAtomEnvironment

  const ExampleClass::EnumType ExampleChemistryAtomEnvironment::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomEnvironment())
  );

} // namespace bcl
