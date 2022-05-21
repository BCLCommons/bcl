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
#include "chemistry/bcl_chemistry_small_molecule_fragment_isomorphism.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_small_molecule_fragment_isomorphism_picked.cpp
  //!
  //! @author kothiwsk
  //! @date Nov 23, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySmallMoleculeFragmentIsomorphism :
    public ExampleInterface
  {
  public:

    ExampleChemistrySmallMoleculeFragmentIsomorphism *Clone() const
    {
      return new ExampleChemistrySmallMoleculeFragmentIsomorphism( *this);
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

      // create input stream for reading molecule of interest which needs to generated
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_molecule_conformations.sdf"));
      // read in molecule
      chemistry::FragmentEnsemble molecule( input, sdf::e_Remove);
      // close stream
      io::File::CloseClearFStream( input);

      // create a tree searcher
      chemistry::SearchFragmentLibraryFromTree search_lib
      (
        ( chemistry::RotamerLibraryFile( chemistry::RotamerLibraryInterface::GetDefault()))
      );

      // get isomorphism between of fragments with molecule , for fragments contained in molecule
      util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> fragment_isomorphisms
      (
        search_lib.FindFragmentsOfMolecule( molecule.GetMolecules().FirstElement())
      );

    /////////////////
    // data access //
    /////////////////

      // test GetFragmentIsomorphism
      BCL_ExampleCheck( fragment_isomorphisms.GetSize(), 16);

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

  }; //end ExampleChemistrySmallMoleculeFragmentIsomorphismPicked

  const ExampleClass::EnumType ExampleChemistrySmallMoleculeFragmentIsomorphism::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySmallMoleculeFragmentIsomorphism())
  );

} // namespace bcl
