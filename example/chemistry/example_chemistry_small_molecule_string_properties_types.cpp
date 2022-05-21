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
#include "chemistry/bcl_chemistry_small_molecule_string_properties_types.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_small_molecule_string_properties_types.cpp
  //!
  //! @author mendenjl
  //! @date Dec 01, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySmallMoleculeStringPropertiesTypes :
    public ExampleInterface
  {
  public:

    ExampleChemistrySmallMoleculeStringPropertiesTypes *Clone() const
    {
      return new ExampleChemistrySmallMoleculeStringPropertiesTypes( *this);
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
      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input);
      const chemistry::ConformationInterface &first_molecule( *ensemble.Begin());
      // close stream
      io::File::CloseClearFStream( input);

      // try out the property that retrieves atom types, bond types, and chirality
      util::Implementation< chemistry::StringPropertyInterface>
        atom_types_retriever( "AtomTypes"), bond_types_retriever( "BondTypes"), chirality_retriever( "Chirality");

      // ensure that each is defined
      BCL_ExampleAssert( atom_types_retriever.IsDefined(), true);
      BCL_ExampleAssert( bond_types_retriever.IsDefined(), true);
      BCL_ExampleAssert( chirality_retriever.IsDefined(), true);

      // check the descriptor string for each
      BCL_ExampleCheck( atom_types_retriever.GetString(), "AtomTypes");
      BCL_ExampleCheck( bond_types_retriever.GetString(), "BondTypes");
      BCL_ExampleCheck( chirality_retriever.GetString(), "Chirality");

      // try each out on a small molecule
      BCL_ExampleCheck
      (
        atom_types_retriever->operator()( first_molecule),
        util::TrimString( first_molecule.GetAtomTypesString())
      );
      BCL_ExampleCheck
      (
        bond_types_retriever->operator()( first_molecule),
        util::TrimString( first_molecule.GetBondTypesString())
      );
      BCL_ExampleCheck
      (
        chirality_retriever->operator()( first_molecule),
        util::TrimString( first_molecule.GetChiralityString())
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistrySmallMoleculeStringPropertiesTypes

  const ExampleClass::EnumType ExampleChemistrySmallMoleculeStringPropertiesTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySmallMoleculeStringPropertiesTypes())
  );

} // namespace bcl
