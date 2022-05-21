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
#include "chemistry/bcl_chemistry_small_molecule_string_properties_mapped.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_small_molecule_string_properties_mapped.cpp
  //!
  //! @author mendenjl
  //! @date Mar 16, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySmallMoleculeStringPropertiesMapped :
    public ExampleInterface
  {
  public:

    ExampleChemistrySmallMoleculeStringPropertiesMapped *Clone() const
    {
      return new ExampleChemistrySmallMoleculeStringPropertiesMapped( *this);
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
      // construct an implementation that should map values of the property Mol_ID onto TestValues
      const util::Implementation< chemistry::StringPropertyInterface>
      mapper
      (
        "MappedString(file=" + AddExampleInputPathToFilename( e_Chemistry, "input_map.object") + ",property=Cached(Mol_ID))"
      );

      // ensure that each is defined
      BCL_ExampleAssert( mapper.IsDefined(), true);

      // check the descriptor string
      BCL_ExampleCheck( mapper.GetLabel().GetValue(), "MappedString");

      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_map.sdf"));
      // read in ensemble, remove hydrogens
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Remove);
      // close stream
      io::File::CloseClearFStream( input);
      // get each molecule out of the ensemble
      storage::Vector< chemistry::FragmentComplete> fragments( ensemble.Begin(), ensemble.End());

      // get the value that each id should map to in each molecule
      const storage::Vector< std::string> expected_values
      (
        storage::Vector< std::string>::Create( "3.0 14.0 18.0", "5.0 6.0 12.0", "", "")
      );

      for( size_t molecule_number( 0), number_mols( 4); molecule_number < number_mols; ++molecule_number)
      {
        BCL_ExampleCheck( mapper->operator()( fragments( molecule_number)), expected_values( molecule_number));
      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistrySmallMoleculeStringPropertiesMapped

  const ExampleClass::EnumType ExampleChemistrySmallMoleculeStringPropertiesMapped::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySmallMoleculeStringPropertiesMapped())
  );

} // namespace bcl
