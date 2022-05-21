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
#include "descriptor/bcl_descriptor_numeric.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_numeric.cpp
  //!
  //! @author kothiwsk
  //! @date Dec 04, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorNumeric :
    public ExampleInterface
  {
  public:

    ExampleDescriptorNumeric *Clone() const
    {
      return new ExampleDescriptorNumeric( *this);
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
      descriptor::Numeric< chemistry::AtomConformationalInterface> molecular_mass;
      BCL_ExampleCheck( molecular_mass.GetAlias(), "ToNumeric");

      // create the label for the molecular mass property
      util::ObjectDataLabel molecular_mass_label( "ToNumeric(ToString(Weight),size=1)");
      BCL_ExampleAssert( molecular_mass.TryRead( util::ObjectDataLabel( "(ToString(Weight),size=1)"), util::GetLogger()), true);

      descriptor::CheminfoProperty molecular_mass_impl;
      BCL_ExampleCheck( molecular_mass_impl.TryRead( molecular_mass_label, util::GetLogger()), true);
      // check the resulting label
      BCL_ExampleCheck( molecular_mass.GetLabel(), molecular_mass_label);
      BCL_ExampleCheck( molecular_mass_impl.GetLabel(), molecular_mass_label);

      // create another property that calculates the atomic masses and return them as a string
      descriptor::Numeric< chemistry::AtomConformationalInterface> atomic_masses_retriever;
      util::ObjectDataLabel atomic_mass_label( "ToNumeric(ToString(Atom_Mass),size=1)");
      BCL_ExampleAssert( atomic_masses_retriever.TryRead( util::ObjectDataLabel( "(ToString(Atom_Mass),size=1)"), util::GetLogger()), true);
      // check the resulting label
      BCL_ExampleCheck( atomic_mass_label, atomic_masses_retriever.GetLabel());

      //TODO
      // check that the molecular mass retriever works
      //      BCL_ExampleCheck( molecular_mass( first_molecule), size_t( 272.2));

      // get the masses of the elements in the first molecule
      // the first molecule has heavy atoms in this order C3O2C17
      float carbon_mass
      (
        chemistry::GetElementTypes().e_Carbon->GetProperty( chemistry::ElementTypeData::e_Mass)
      );
      float oxygen_mass
      (
        chemistry::GetElementTypes().e_Oxygen->GetProperty( chemistry::ElementTypeData::e_Mass)
      );
      storage::Vector< float> atomic_masses( 3, carbon_mass);
      atomic_masses.Append( storage::Vector< float>( 2, oxygen_mass));
      atomic_masses.Append( storage::Vector< float>( 17, carbon_mass));

      // check that the retriever gives us the same thing
      BCL_ExampleCheck
      (
        atomic_masses_retriever.CollectValuesOnEachElementOfObject( first_molecule),
        linal::Vector< float>( atomic_masses)
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorNumeric

  const ExampleClass::EnumType ExampleDescriptorNumeric::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorNumeric())
  );

} // namespace bcl
