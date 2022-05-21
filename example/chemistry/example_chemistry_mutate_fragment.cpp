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
#include "chemistry/bcl_chemistry_mutate_fragment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_mutate_fragment.cpp
  //! @author kothiwsk
  //! @date 11/23/2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMutateFragment :
    public ExampleInterface
  {
  public:

    ExampleChemistryMutateFragment *Clone() const
    {
      return new ExampleChemistryMutateFragment( *this);
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

      // create input stream for reading a fragment ensemble
      io::IFStream input;

      // create input stream for reading molecule of interest which needs to generated
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));
      // read in molecule
      chemistry::FragmentEnsemble molecules( input, sdf::e_Remove);
      // close stream
      io::File::CloseClearFStream( input);

      chemistry::FragmentComplete &first_molecule( molecules.GetMolecules().FirstElement());

      // create a tree searcher
      chemistry::SearchFragmentLibraryFromTree search_lib
      (
        ( chemistry::RotamerLibraryFile( chemistry::RotamerLibraryInterface::GetDefault()))
      );

      // get isomorphism between of fragments with molecule , for fragments contained in molecule
      util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> fragment_isomorphisms
      (
        search_lib.FindFragmentsOfMolecule( first_molecule)
      );

      // get the rotamer dihedral bond data for the molecule
      const util::ShPtrVector< chemistry::RotamerDihedralBondData> bond_mapping
      (
        chemistry::SmallMoleculeFragmentMapping().MapFragmentIsomorphisms
        (
          first_molecule, fragment_isomorphisms
        )
      );

      // create a graph where atoms are colored by atom type or chirality and edges by whether edge is a bond ring or not
      chemistry::ConformationGraphConverter graph_maker
      (
        chemistry::ConformationGraphConverter::e_AtomType,
        chemistry::ConfigurationalBondTypeData::e_IsInRing
      );
      const graph::ConstGraph< size_t, size_t> molecule_graph( graph_maker( first_molecule));

      chemistry::MutateChirality mutchi( first_molecule);
      chemistry::MutateFragment first_fragment_mutate
      (
        chemistry::MutateFragment ( bond_mapping.FirstElement(), molecule_graph, first_molecule, mutchi)
      );

      chemistry::FragmentEnsemble new_conformations;
      new_conformations.PushBack( *first_fragment_mutate( first_molecule).GetArgument());

      for( size_t counts( 0); counts < size_t( 5); ++counts)
      {
        new_conformations.PushBack( *first_fragment_mutate( new_conformations.GetMolecules().LastElement()).GetArgument());
      }

      const std::string mutated_conformations( AddExampleOutputPathToFilename( new_conformations, "mutated_conformation.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, mutated_conformations);
      new_conformations.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename( mutated_conformations + ".correct");

      // check if the generated file is right
      if
      (
        BCL_ExampleCheck( io::File::FilesMatchWithinAbsoluteTolerance( mutated_conformations, correct_filename, 0.05), true)
      )
      {
        remove( mutated_conformations.c_str());
      }

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
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryMutateFragment

  const ExampleClass::EnumType ExampleChemistryMutateFragment::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMutateFragment())
  );

} // namespace bcl
