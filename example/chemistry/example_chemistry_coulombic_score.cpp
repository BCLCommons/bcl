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
#include "chemistry/bcl_chemistry_coulombic_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_coulombic_score.cpp
  //!
  //! @author kothiwsk
  //! @date Jan 05, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryCoulombicScore :
    public ExampleInterface
  {
  public:

    ExampleChemistryCoulombicScore *Clone() const
    {
      return new ExampleChemistryCoulombicScore( *this);
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

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_1_rotamers.sdf"));
      // read in ensemble
      util::ShPtr< chemistry::FragmentEnsemble> fragments
      (
        new chemistry::FragmentEnsemble( input, sdf::e_Remove)
      );
      // close stream
      io::File::CloseClearFStream( input);

      // create input stream for reading molecule of interest which needs to generated
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_1_conformers.sdf"));
      // read in molecule
      chemistry::FragmentEnsemble conformers( input, sdf::e_Remove);
      // close stream
      io::File::CloseClearFStream( input);

      // expected scores
      storage::Vector< double> expected_scores
      (
        storage::Vector< double>::Create( -0.24617, -0.246162, -0.246231, -0.246119, -0.246262)
      );

      descriptor::CheminfoProperty Vcharge( "Atom_Vcharge");
      // object for calculating fragment probability score
      chemistry::CoulombicScore sp_fragment_probability
      (
        Vcharge,
        double( 1.0)
      );

      // collect scores of all conformations in a container
      storage::Vector< double> conformation_scores;
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
          itr( conformers.GetMolecules().Begin()), itr_end( conformers.GetMolecules().End());
        itr != itr_end;
        ++itr
      )
      {
        conformation_scores.PushBack( sp_fragment_probability( *itr));
      }

      BCL_ExampleCheckWithinTolerance( conformation_scores, expected_scores, 0.01);

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

  }; //end ExampleChemistryCoulombicScore

  const ExampleClass::EnumType ExampleChemistryCoulombicScore::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryCoulombicScore())
  );

} // namespace bcl
