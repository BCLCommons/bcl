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
#include "chemistry/bcl_chemistry_fragment_probability_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_probability_score.cpp
  //!
  //! @author kothiwsk
  //! @date Sep 26, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentProbabilityScore :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentProbabilityScore *Clone() const
    {
      return new ExampleChemistryFragmentProbabilityScore( *this);
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

      // load molecule conformations from the pdb / platinum dataset
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "platinum_dataset_first_ten.sdf"));
      // read in molecule
      chemistry::FragmentEnsemble molecules_pdb( input);
      // close stream
      io::File::CloseClearFStream( input);

      // load molecule conformations from corina
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "platinum_dataset_first_ten.corina.sdf"));
      // read in ensemble
      chemistry::FragmentEnsemble molecules_corina( input);
      // close stream
      io::File::CloseClearFStream( input);

      chemistry::SearchFragmentLibraryFromTree search_lib
      (
        ( chemistry::RotamerLibraryFile( chemistry::RotamerLibraryInterface::GetDefault()))
      );

      // get isomorphism between of fragments with molecule , for fragments contained in molecule

      // expected scores
      storage::Vector< double> expected_scores_pdb
      (
        storage::Vector< double>::Create
        (
          0.553366, 0.16093, 0.469666,
          -0.042187, 0.16591, 0.242851,
          0.458787, 0.450205, 0.389701,
          0.695205
        )
      );

      storage::Vector< double> expected_scores_corina
      (
        storage::Vector< double>::Create
        (
          0.464281, 0.218054, -0.128767,
          -0.0560069, 0.151984, 0.0217566,
          -0.0818386, -0.0212791, 0.0773302,
          0.481244
        )
      );

      // collect scores of all conformations in a container
      storage::Vector< double> conformation_scores_pdb;
      storage::Vector< double> conformation_scores_corina;
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
          itr( molecules_pdb.GetMolecules().Begin()),
          itr_corina( molecules_corina.GetMolecules().Begin()), itr_end( molecules_pdb.GetMolecules().End());
        itr != itr_end;
        ++itr, ++itr_corina
      )
      {
        // find fragment isomorphisms
        util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> fragment_isomorphisms
        (
          search_lib.FindFragmentsOfMolecule( *itr)
        );

        // get the rotamer dihedral bond data for the molecule
        const util::ShPtrVector< chemistry::RotamerDihedralBondData> bond_mapping
        (
          chemistry::SmallMoleculeFragmentMapping().MapFragmentIsomorphisms
          (
            *itr, fragment_isomorphisms
          )
        );

        // object for calculating fragment probability score
        chemistry::FragmentProbabilityScore fragment_probability( bond_mapping, false);

        conformation_scores_pdb.PushBack( fragment_probability( *itr));
        conformation_scores_corina.PushBack( fragment_probability( *itr_corina));
      }

      BCL_ExampleCheckWithinTolerance( conformation_scores_pdb, expected_scores_pdb, 0.01);
      BCL_ExampleCheckWithinTolerance( conformation_scores_corina, expected_scores_corina, 0.01);

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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentProbabilityScore

  const ExampleClass::EnumType ExampleChemistryFragmentProbabilityScore::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentProbabilityScore())
  );

} // namespace bcl
