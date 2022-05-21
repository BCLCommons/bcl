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
#include "descriptor/bcl_descriptor_numeric_string.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_numeric_string.cpp
  //!
  //! @author mendenjl
  //! @date Feb 17, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorNumericString :
    public ExampleInterface
  {
  public:

    ExampleDescriptorNumericString *Clone() const
    {
      return new ExampleDescriptorNumericString( *this);
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
      // read in ensemble, remove hydrogens
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Remove);
      // close stream
      io::File::CloseClearFStream( input);
      // get the first molecule out of the ensemble
      const chemistry::ConformationInterface &first_molecule( *ensemble.Begin());

      // make a property that calculates the molecular mass and returns it as a string
      descriptor::NumericString< chemistry::AtomConformationalInterface> molecular_mass;
      BCL_ExampleCheck( molecular_mass.GetAlias(), "ToString");

      // create the label for the molecular mass property
      util::ObjectDataLabel molecular_mass_label( "ToString(Weight)");
      BCL_ExampleAssert( molecular_mass.TryRead( util::ObjectDataLabel( "Weight"), util::GetLogger()), true);

      descriptor::CheminfoID molecular_mass_impl;
      BCL_ExampleCheck( molecular_mass_impl.TryRead( molecular_mass_label, util::GetLogger()), true);
      // check the resulting label
      BCL_ExampleCheck( molecular_mass.GetLabel(), util::ObjectDataLabel( "Weight"));
      BCL_ExampleCheck( molecular_mass_impl.GetLabel(), molecular_mass_label);

      // create another property that calculates the atomic masses and return them as a string
      descriptor::NumericString< chemistry::AtomConformationalInterface> atomic_masses_retriever;
      util::ObjectDataLabel atomic_mass_label( "ToString(Atom_Mass)");
      BCL_ExampleAssert( atomic_masses_retriever.TryRead( util::ObjectDataLabel( "Atom_Mass"), util::GetLogger()), true);
      // check the resulting label
      BCL_ExampleCheck( util::ObjectDataLabel( "Atom_Mass"), atomic_masses_retriever.GetLabel());

      //TODO
      // check that the molecular mass retriever works
      //      BCL_ExampleCheck( molecular_mass( first_molecule), size_t( 272.2));

      // get the masses of the elements in the first molecule
      // the first molecule has heavy atoms in this order C3O2C17
      float carbon_mass
      (
        chemistry::GetElementTypes().e_Carbon->GetProperty( chemistry::ElementTypeData::e_Mass)
      );
      std::string carbon_mass_str( util::Format().W( 13)( carbon_mass) + " ");
      std::string three_carbon_mass_str( carbon_mass_str + carbon_mass_str + carbon_mass_str);
      std::string seventeen_carbon_mass_str
      (
        three_carbon_mass_str + three_carbon_mass_str + three_carbon_mass_str
        + three_carbon_mass_str + three_carbon_mass_str + carbon_mass_str + carbon_mass_str
      );
      float oxygen_mass
      (
        chemistry::GetElementTypes().e_Oxygen->GetProperty( chemistry::ElementTypeData::e_Mass)
      );
      std::string two_oxygen_mass_str( util::Format().W( 13)( oxygen_mass) + " " + util::Format().W( 13)( oxygen_mass) + " ");
      std::string heavy_atoms_mass_str( three_carbon_mass_str + two_oxygen_mass_str + seventeen_carbon_mass_str);
      std::string atomic_masses( util::TrimString( heavy_atoms_mass_str));

      // check that the retriever gives us the same thing
      BCL_ExampleCheck
      (
        util::TrimString
        (
          std::string( atomic_masses_retriever.CollectValuesOnEachElementOfObject( first_molecule).Begin(), 14 * 22)
        ),
        atomic_masses
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorNumericString

  const ExampleClass::EnumType ExampleDescriptorNumericString::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorNumericString())
  );

} // namespace bcl
