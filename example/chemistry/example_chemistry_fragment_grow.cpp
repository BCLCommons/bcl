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
#include "chemistry/bcl_chemistry_fragment_grow.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_pick_atom_random.h"
#include "chemistry/bcl_chemistry_pick_fragment_random.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_grow.cpp
  //!
  //! @author kothiwsk
  //! @date Jun 29, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentGrow :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentGrow *Clone() const
    {
      return new ExampleChemistryFragmentGrow( *this);
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

      // construct from properties
      chemistry::FragmentGrow fragment_grow;

      // clone function
      util::ShPtr< chemistry::FragmentGrow> fragment_grow_clone( fragment_grow.Clone());
      BCL_ExampleCheck( fragment_grow_clone.IsDefined(), true);

      // create input stream for reading a fragment ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
      // close stream
      // creating ShPtr of growfragments
      util::ShPtr< chemistry::FragmentEnsemble> sp_fragment_pool( new chemistry::FragmentEnsemble( input, sdf::e_Remove));
      io::File::CloseClearFStream( input);
      // remove hydrogens to make valency
      sp_fragment_pool->RemoveH();

      chemistry::FragmentGrow fragment_grower
      (
        sp_fragment_pool,
        chemistry::CollectorValence(),
        chemistry::PickAtomRandom(),
        chemistry::PickFragmentRandom()
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

      // create input stream for reading a base molecule
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "hexane.sdf"));
      // read in ensemble
      chemistry::FragmentEnsemble base_mol( input);
      chemistry::FragmentComplete hexane( base_mol.GetMolecules().FirstElement());
      // close stream
      io::File::CloseClearFStream( input);
      // remove hydrogen to make valency
      base_mol.RemoveH();

      // Grow a fragment on hexane
      chemistry::FragmentComplete new_molecule( *fragment_grower( hexane).GetArgument());
      // Write the Generated fragment to a file
      const std::string fragment_grow_out( AddExampleOutputPathToFilename( new_molecule, "fragment_grow.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, fragment_grow_out);
      // write out the constitution as we do not have coordinates
      new_molecule.WriteMDL( output);
      io::File::CloseClearFStream( output);
      std::string correct_filename( fragment_grow_out + ".correct");

      // check if the generated file is right
      // check if the generated file is right
      if
      (
        BCL_ExampleCheck( io::File::FilesMatch( fragment_grow_out, correct_filename), true)
      )
      {
        remove( fragment_grow_out.c_str());
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentGrow

  const ExampleClass::EnumType ExampleChemistryFragmentGrow::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentGrow())
  );

} // namespace bcl
