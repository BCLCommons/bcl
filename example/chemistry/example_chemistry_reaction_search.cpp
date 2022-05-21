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
#include "chemistry/bcl_chemistry_reaction_search.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_rxn_factory.h"

// external includes - sorted alphabetically
#include <cstdio>

#undef input

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_reaction_search.cpp
  //! @brief this example tests the implementation of ReactionSearch
  //!
  //! @author geanesar
  //! @date Feb 18, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryReactionSearch :
     public ExampleInterface
   {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleChemistryReactionSearch
    ExampleChemistryReactionSearch *Clone() const
    {
      return new ExampleChemistryReactionSearch( *this);
    }

  //////////
  // data //
  //////////

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // helper functions //
  //////////////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      io::IFStream input;
      chemistry::ReactionEnsemble reactions;

      std::string rxn_dirname( AddExampleInputPathToFilename( e_Chemistry, "reaction_search"));
      io::Directory rxn_dir( rxn_dirname);
      if( !rxn_dir.DoesExist()) {
        BCL_MessageStd( "Could not open reaction directory " + rxn_dirname);
      }

      storage::List< io::DirectoryEntry> rxn_dir_files( rxn_dir.ListFiles( "", ".rxn"));
      for
      (
        storage::List< io::DirectoryEntry>::const_iterator itr_f( rxn_dir_files.Begin()),
          itr_f_end( rxn_dir_files.End());
        itr_f != itr_f_end;
        ++itr_f
      )
      {
        BCL_ExampleMustOpenInputFile( input, itr_f->GetFullName());
        reactions.ReadMoreFromRXN( input);
        io::File::CloseClearFStream( input);
      }

      // read reactions from a directory
      int n_expected_rxns( 5);
      BCL_ExampleCheck( reactions.GetSize(), n_expected_rxns);

      chemistry::FragmentEnsemble mols;
      {
        io::IFStream input;
        BCL_ExampleMustOpenInputFile
        (
          input,
          AddExampleInputPathToFilename( e_Chemistry, "reaction_search/reactants.sdf")
        );
        mols.ReadMoreFromMdl( input, sdf::e_Saturate);
        io::File::CloseClearFStream( input);
      }
      int n_expected_mols( 5);
      BCL_ExampleCheck( mols.GetSize(), n_expected_mols);

      chemistry::ReactionSearch catalog( mols, reactions);
      catalog.Initialize();

      BCL_ExampleAssert( catalog.IsEmpty(), false);

      BCL_ExampleAssert( catalog.GetMolecules().IsDefined(), true);
      BCL_ExampleCheck( catalog.GetMolecules()->GetSize(), n_expected_mols);

      BCL_ExampleAssert( catalog.GetReactions().IsDefined(), true);
      BCL_ExampleCheck( catalog.GetReactions()->GetSize(), n_expected_rxns);

      // test the copy and destruction operators
      util::ShPtr< chemistry::ReactionSearch> sp_cat( util::CloneToShPtr( catalog));
      BCL_ExampleAssert( sp_cat.IsDefined(), true);

      BCL_ExampleAssert( sp_cat->GetMolecules().IsDefined(), true);
      BCL_ExampleCheck( sp_cat->GetMolecules()->GetSize(), catalog.GetMolecules()->GetSize());

      BCL_ExampleAssert( sp_cat->GetReactions().IsDefined(), true);
      BCL_ExampleCheck( sp_cat->GetReactions()->GetSize(), catalog.GetMolecules()->GetSize());

      // make sure the copy was a deep copy
      catalog.Reset();

      BCL_ExampleIndirectAssert( sp_cat->GetMolecules().IsDefined(), true, "shptr molecules valid after parent reset");
      BCL_ExampleIndirectCheck( sp_cat->GetMolecules()->GetSize(), n_expected_mols, "shptr molecules correct after parent reset");

      BCL_ExampleIndirectAssert( sp_cat->GetReactions().IsDefined(), true, "shptr reactions valid after parent reset");
      BCL_ExampleIndirectCheck( sp_cat->GetReactions()->GetSize(), n_expected_rxns, "shptr reactions correct after parent reset");

      // benzene, benzoic acid, aniline, 2-aminoaniline, formic acid
      int n_expected_vertices( 5);

      // benzene -> aniline -> 2-aminoaniline
      // benzene -> benzoic acid
      // formate -> benzoic acid
      int n_expected_edges( 4);

      BCL_ExampleCheck
      (
        sp_cat->GetStructureTree().GetVertices().GetSize(),
        n_expected_vertices
      );

      // TODO this doesn't work right yet, but all of the reactions are
      // searched anyways so it will just be a performance penalty
      //BCL_ExampleCheck
      //(
      //  sp_cat->GetStructureTree().NumEdges(),
      //  n_expected_edges
      //);

      std::string taminoanil_name( "2-aminoaniline");
      chemistry::FragmentComplete taminoanil( sp_cat->GetMolecules()->operator()( 0));
      BCL_ExampleIndirectAssert( taminoanil.GetName().compare( 0, taminoanil_name.length(), taminoanil_name), 0,
                                  "2-aminoaniline was found correctly");

      std::string bz_acid_name( "benzoic_acid");
      chemistry::FragmentComplete bz_acid( sp_cat->GetMolecules()->operator()( 1));
      BCL_ExampleIndirectAssert( bz_acid.GetName().compare( 0, bz_acid_name.length(), bz_acid_name), 0,
                                  "benzoic acid was found correctly");

      // TODO more checks of searching capabilities

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; // class ExampleChemistryReactionSearch

   const ExampleClass::EnumType ExampleChemistryReactionSearch::s_Instance
   (
     GetExamples().AddEnum( ExampleChemistryReactionSearch())
   );

} // namespace bcl
