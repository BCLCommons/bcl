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
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_merge_fragment_complete.cpp
  //! @details Tests ChemistryMergeFragmentComplete class
  //! @author kothiwsk
  //! @date Oct 08, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMergeFragmentComplete :
    public ExampleInterface
  {
  public:

    ExampleChemistryMergeFragmentComplete *Clone() const
    {
      return new ExampleChemistryMergeFragmentComplete( *this);
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

      // create input stream for reading a fragment ensemble and read in first molecule
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_E.sdf"));
      // read in ensemble
      util::ShPtr< chemistry::FragmentEnsemble> fragments
      (
        new chemistry::FragmentEnsemble( input, sdf::e_Remove)
      );
      io::File::CloseClearFStream( input);

      // read in the second molecule
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "benzene_prepared.sdf"));
      fragments->Append( chemistry::FragmentEnsemble( input, sdf::e_Remove));
      // close stream
      io::File::CloseClearFStream( input);

      // get reference to first and the second molecules
      chemistry::FragmentComplete &first_mol( fragments->GetMolecules().FirstElement());
      chemistry::FragmentComplete &second_mol( fragments->GetMolecules().LastElement());

      // constructor for MergeFragmentComplete
      chemistry::MergeFragmentComplete merge_molecules;

      // create an ensemble to store the created molecules
      chemistry::FragmentEnsemble merged_molecule_ensemble;

      // create a map and a vector of set to be used for merging fragments
      storage::Vector< storage::Pair< size_t, size_t> > indices_map;
      storage::Map< size_t, size_t> common_indices;

      // create first element to direct merginging of second and first molecule
      indices_map.PushBack( storage::Pair< size_t, size_t>( 1, 0));
      common_indices.Insert( indices_map.LastElement());
      storage::Vector< size_t> appended_atoms;
      chemistry::FragmentComplete merged_molecule_a
      (
        merge_molecules.MergeFragments( second_mol, first_mol, common_indices, appended_atoms).Second()
      );

      // store the merged fragment in ensemble
      merged_molecule_ensemble.PushBack( merged_molecule_a);

      // get another merged molecule by connecting first and second molecule
      chemistry::FragmentComplete merged_molecule_aa
      (
        merge_molecules.MergeFragments
        (
          second_mol, first_mol, chemistry::ConfigurationalBondType( 1), indices_map.FirstElement()
        ).Second()
      );
      merged_molecule_ensemble.PushBack( merged_molecule_aa);

      // create another molecule by merging molecules a and b at multiple overlapping atoms
      indices_map.PushBack( storage::Pair< size_t, size_t>( 0, 1));
      common_indices.Insert( indices_map.LastElement());
      chemistry::FragmentComplete merged_molecule_b
      (
        merge_molecules.MergeFragments( second_mol, first_mol, common_indices, appended_atoms).Second()
      );
      merged_molecule_ensemble.PushBack( merged_molecule_b);

      // get vinyl
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "corina_vinyl.sdf"));
      fragments->Append( chemistry::FragmentEnsemble( input, sdf::e_Remove));
      // close stream
      io::File::CloseClearFStream( input);

      chemistry::FragmentComplete &vinyl( fragments->GetMolecules().LastElement());

      // merge vinyl with first molecule at a single atom
      indices_map.PushBack( storage::Pair< size_t, size_t>( 2, 3));
      common_indices.Insert( indices_map.LastElement());
      chemistry::FragmentComplete merged_molecule_c
      (
        merge_molecules.MergeFragments( first_mol, vinyl, common_indices, appended_atoms).Second()
      );
      merged_molecule_ensemble.PushBack( merged_molecule_c);

      // merge vinyl with first molecule at multiple overlapping atoms
      chemistry::FragmentComplete merged_molecule_ca
      (
        merge_molecules.MergeFragments
        (
          first_mol, vinyl, chemistry::ConfigurationalBondType( 1), indices_map.FirstElement()
        ).Second()
      );

      merged_molecule_ensemble.PushBack( merged_molecule_ca);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      const std::string merge_fragments( AddExampleOutputPathToFilename( merged_molecule_ensemble, "merged_molecules.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, merge_fragments);
      merged_molecule_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename( merge_fragments + ".correct");

      // check if the generated file is right
      if
      (
        BCL_ExampleCheck( io::File::FilesMatch( merge_fragments, correct_filename), true)
      )
      {
        remove( merge_fragments.c_str());
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

  }; //end ExampleChemistryMergeFragmentComplete

  const ExampleClass::EnumType ExampleChemistryMergeFragmentComplete::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMergeFragmentComplete())
  );

} // namespace bcl
