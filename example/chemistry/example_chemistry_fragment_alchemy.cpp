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
#include "chemistry/bcl_chemistry_fragment_alchemy.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_pick_atom_random.h"
#include "chemistry/bcl_chemistry_pick_fragment_random.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_alchemy.cpp
  //!
  //! @author brownbp1
  //! @date Dec 12, 2019
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentAlchemy :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentAlchemy *Clone() const
    {
      return new ExampleChemistryFragmentAlchemy( *this);
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
      chemistry::FragmentAlchemy fragment_alchemy;

      // clone function
      util::ShPtr< chemistry::FragmentAlchemy> fragment_alchemy_clone( fragment_alchemy.Clone());
      BCL_ExampleCheck( fragment_alchemy_clone.IsDefined(), true);

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
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "d6w.sdf"));

      // read in ensemble
      chemistry::FragmentEnsemble base_mol( input);
      chemistry::FragmentComplete d6w( base_mol.GetMolecules().FirstElement());

      // close stream
      io::File::CloseClearFStream( input);

      // set some options
      fragment_alchemy.SetAllowedElements( storage::Vector< chemistry::ElementType>( 1, chemistry::GetElementTypes().e_Carbon));
      fragment_alchemy.SetMutableAtomIndices( storage::Vector< size_t>( 1, 12));
      fragment_alchemy.SetOVReverseH( true);

      // Transform an atom in d6w
      math::MutateResult< chemistry::FragmentComplete> new_mutate;
      util::ShPtr< chemistry::FragmentComplete> new_molecule( new chemistry::FragmentComplete);
      new_mutate = fragment_alchemy( d6w);
      if( new_mutate.GetArgument().IsDefined())
      {
        new_molecule = new_mutate.GetArgument();
      }

      // Write the Generated fragment to a file
      const std::string fragment_alchemy_out( AddExampleOutputPathToFilename( *new_molecule, "fragment_alchemy.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, fragment_alchemy_out);

      // write out the constitution as we do not have coordinates
      new_molecule->WriteMDL( output);
      io::File::CloseClearFStream( output);
      std::string correct_filename( fragment_alchemy_out + ".correct");

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentAlchemy

  const ExampleClass::EnumType ExampleChemistryFragmentAlchemy::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentAlchemy())
  );

} // namespace bcl
