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
#include "chemistry/bcl_chemistry_pick_atom_random.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_pick_atom_random.cpp
  //!
  //! @author mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryPickAtomRandom :
    public ExampleInterface
  {
  public:

    ExampleChemistryPickAtomRandom *Clone() const
    {
      return new ExampleChemistryPickAtomRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // default constructor
      chemistry::PickAtomRandom pick_atom_random_default;

      // clone
      const util::ShPtr< chemistry::PickAtomRandom> pick_atom_random_clone( pick_atom_random_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_ExampleIndirectCheck( pick_atom_random_clone.IsDefined(), true, "clone");
      BCL_ExampleCheck
      (
        chemistry::PickAtomRandom().GetClassIdentifier(),
        GetStaticClassName< chemistry::PickAtomRandom>()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // setup input stream
      io::IFStream input_sdf;

      // read in molecule
      const std::string filename_in( AddExampleInputPathToFilename( e_Chemistry, "corina_diazepam.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, filename_in);

      chemistry::FragmentComplete diazepam( sdf::FragmentFactory::MakeFragment( input_sdf, sdf::e_Remove));

      iterate::Generic< const chemistry::AtomConformationalInterface> itr_atoms( diazepam.GetAtomsIterator());
      // get the atoms from diazepam into a list
      util::SiPtrList< const chemistry::AtomConformationalInterface> atoms( itr_atoms.Begin(), itr_atoms.End());

      // initialize a storage::Set to hold the bodies that are picked
      storage::Set< util::SiPtr< const chemistry::AtomConformationalInterface> > picked_atoms;

      for( size_t counter( 0); counter < 1000 && picked_atoms.GetSize() < atoms.GetSize(); ++counter)
      {
        picked_atoms.Insert( pick_atom_random_default.Pick( atoms));
      }

      BCL_ExampleIndirectCheck
      (
        picked_atoms.GetSize(), atoms.GetSize(),
        "a truly random picking is expected to pick all heavy atoms in diazepam at least once within 1000 attempts"
      );

      // end
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryPickAtomRandom

  const ExampleClass::EnumType ExampleChemistryPickAtomRandom::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryPickAtomRandom())
  );

} // namespace bcl
