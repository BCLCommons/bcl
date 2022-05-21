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
#include "chemistry/bcl_chemistry_fragment_molecule.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_breadth_first_FragmentMolecule.cpp
  //!
  //! @author kothiwsk
  //! @date May 03, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentMolecule :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentMolecule *Clone() const
    {
      return new ExampleChemistryFragmentMolecule( *this);
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

    /////////////////
    // data access //
    /////////////////

      // load a ensembles directly from input streams
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));
      chemistry::FragmentEnsemble ensemble( input, sdf::e_Saturate);
      io::File::CloseClearFStream( input);

      util::ShPtr< chemistry::ConstitutionSet> conformation_set( new chemistry::ConstitutionSet());

      chemistry::FragmentMolecule fragment_molecule( conformation_set, size_t( 100));
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
          itr != itr_end;
        ++itr
      )
      {
        fragment_molecule( *itr);
      }

      const std::string fragments( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "mGluR5_fragments.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, fragments);
      // write out the picked fragmetns
      for
      (
        auto
          itr_ensemble( conformation_set->GetConstitutions().Begin()),
          itr_ensemble_end( conformation_set->GetConstitutions().End());
        itr_ensemble != itr_ensemble_end;
        ++itr_ensemble
      )
      {
        // write the fragment out to a file
        ( *itr_ensemble)->WriteMDL( output);
      }
      io::File::CloseClearFStream( output);
      std::string correct_filename( fragments + ".correct");

      // check if the generated file is right
      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( fragments, correct_filename),
          true,
          "fragmentation algorithm"
        )
      )
      {
        remove( fragments.c_str());
      }
      else
      {
        BCL_MessageStd( "Test file \"" + fragments + "\" did not match expected output in \"" + correct_filename + "\"");
      }

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

  }; //end ExampleChemistryFragmentMolecule

  const ExampleClass::EnumType ExampleChemistryFragmentMolecule::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentMolecule())
  );

} // namespace bcl
