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

//// (c) Copyright BCL @ Vanderbilt University 2014
//// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
//// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
//// (c) This file is part of the BCL software suite and is made available under license.
//// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
//// (c) For academic and non-profit users:
//// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
//// (c) For commercial users:
//// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
//// (c) For BCL developers at Vanderbilt University:
//// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
//// (c)
//// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the
//// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
//// (c) Vanderbilt University.
//// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
//// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org
//// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)
//
//// include example header
//#include "example.h"
//
//// include the header of the class which this example is for
//#include "chemistry/bcl_chemistry_reaction_worker.h"
//
//// includes from bcl - sorted alphabetically
//#include "chemistry/bcl_chemistry_fragment_ensemble.h"
//#include "io/bcl_io_file.h"
//#include "sdf/bcl_sdf_rxn_factory.h"
//
//// external includes - sorted alphabetically
//#include <cstdio>
//
//namespace bcl
//{
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  //!
//  //! @example example_chemistry_reaction_worker.cpp
//  //! @brief this example tests the implementation of ReactionWorker
//  //!
//  //! @author geanesar, mendenjl
//  //! @date Jan 28, 2015
//  //! @remarks status complete
//  //! @remarks reviewed by nobody on
//  //!
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  class ExampleChemistryReactionWorker :
//     public ExampleInterface
//   {
//
//  //////////////////////////////////
//  // construction and destruction //
//  //////////////////////////////////
//
//  public:
//
//    //! @brief Clone function
//    //! @return pointer to a new ExampleChemistryReactionWorker
//    ExampleChemistryReactionWorker *Clone() const
//    {
//      return new ExampleChemistryReactionWorker( *this);
//    }
//
//  //////////
//  // data //
//  //////////
//
//  /////////////////
//  // data access //
//  /////////////////
//
//    //! @brief returns the class name
//    //! @return the class name as const ref std::string
//    const std::string &GetClassIdentifier() const
//    {
//      return GetStaticClassName( *this);
//    }
//
//  ////////////////
//  // operations //
//  ////////////////
//
//    //! @brief run routine
//    //! this is performing the execution of the example
//    int Run() const
//    {
//
//      io::IFStream input;
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // Construction
//      chemistry::ReactionWorker def;
//      chemistry::ReactionWorker worker;
//
//      // Assignment operator
//      def = worker;
//
//      // Copy constructor
//      util::ShPtr< chemistry::ReactionWorker> sp_rxnworker
//      (
//        util::CloneToShPtr( def)
//      );
//
//      BCL_ExampleCheck( sp_rxnworker.IsDefined(), true);
//
//      // Destructor
//      sp_rxnworker.Reset();
//
//      // Timing test.  Uncomment if you want to test the impact of caching
//      /*      size_t uncached_us( 0);
//      util::Stopwatch query_timer( false);
//      for( size_t i( 0); i < 10000; ++i)
//      {
//        chemistry::ReactionWorker empty_worker;
//        query_timer.Start();
//        size_t op( empty_worker.MatchesReactants( reactant, new_reaction).GetSize());
//        ++op;
//        query_timer.Stop();
//        uncached_us += query_timer.GetTotalTime().GetMicroSeconds();
//        query_timer.Reset();
//      }
//
//      size_t cached_us( 0);
//      worker.MatchesReactants( reactant, new_reaction);
//      for( size_t i( 0); i < 10000; ++i)
//      {
//        query_timer.Start();
//        size_t op( worker.MatchesReactants( reactant, new_reaction).GetSize());
//        ++op;
//        query_timer.Stop();
//        cached_us += query_timer.GetTotalTime().GetMicroSeconds();
//        query_timer.Reset();
//      }
//
//      BCL_MessageStd( "1000x: uncached: " + util::Format()( uncached_us) + ", cached: " + util::Format()( cached_us));
//      */
//
//      //
//      // Reactions.  Must test four categories
//      //   2 reactants, 2 products
//      //   2 reactants, 1 product
//      //   1 reactant, 2 products
//      //   1 reactant, 1 product
//      //
//
//      /*
//      // 2 reactant, 2 product (amide coupling/exchange)
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "amide_coupling.rxn"));
//      sdf::RXNHandler handler_amide( input);
//      io::File::CloseClearFStream( input);
//
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "amide_reactants.sdf"));
//      chemistry::FragmentEnsemble amide_reactants;
//      amide_reactants.ReadMoreFromMdl( input);
//      io::File::CloseClearFStream( input);
//
//      chemistry::ReactionComplete amide_rxn
//      (
//        sdf::RXNFactory::MakeReactionComplete( handler_amide)
//      );
//
//      BCL_MessageStd( "Amide reactants: " + util::Format()( amide_reactants.GetSize()));
//      BCL_MessageStd( "Amide reaction reactants: " + util::Format()( amide_rxn.GetNumberReactants()) + ", "
//          "number products: " + util::Format()( amide_rxn.GetNumberProducts()));
//
//      chemistry::FragmentEnsemble::const_iterator mol_itr( amide_reactants.Begin());
//      BCL_ExampleIndirectCheck( worker.MatchesReactants( *mol_itr, amide_rxn).GetSize(), 1, "Amide reaction, first molecule matches a reactant");
//      ++mol_itr;
//      BCL_ExampleIndirectCheck( worker.MatchesReactants( *mol_itr, amide_rxn).GetSize(), 1, "Amide reaction, second molecule matches a reactant");
//      BCL_ExampleIndirectCheck( worker.MatchesProducts( amide_rxn.GetProduct( 0), amide_rxn).GetSize(), 1, "Amide reaction, first product structure matches itself");
//
//      BCL_MessageStd( "Does reaction work (amide reaction)?");
//      chemistry::FragmentEnsemble amide_res( worker.ExecuteReaction( amide_rxn, amide_reactants));
//      BCL_ExampleIndirectCheck( amide_res.GetSize(), 2, "ExecuteReaction() on amide reaction produced appropriate number of products");
//      */
//
//      // 2 reactant, 2 product (amide coupling/exchange)
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "amide_coupling.rxn"));
//      sdf::RXNHandler handler_amide( input);
//      io::File::CloseClearFStream( input);
//
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "amide_reactants.sdf"));
//      chemistry::FragmentEnsemble amide_reactants;
//      amide_reactants.ReadMoreFromMdl( input);
//      io::File::CloseClearFStream( input);
//
//      chemistry::ReactionComplete amide_rxn
//      (
//        sdf::RXNFactory::MakeReactionComplete( handler_amide)
//      );
//
//      // 2 reactant, 1 product (Click reaction)
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "click.rxn"));
//      sdf::RXNHandler handler_click( input);
//      io::File::CloseClearFStream( input);
//
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "click_reactants.sdf"));
//      chemistry::FragmentEnsemble click_reactants;
//      click_reactants.ReadMoreFromMdl( input);
//      io::File::CloseClearFStream( input);
//
//      chemistry::ReactionComplete click_rxn
//      (
//        sdf::RXNFactory::MakeReactionComplete( handler_click)
//      );
//
//      // 1 reactant, 2 product (retro-diels alder)
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "retro_diels_alder.rxn"));
//      sdf::RXNHandler handler_retroda( input);
//      io::File::CloseClearFStream( input);
//
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "retro_diels_alder_reactants.sdf"));
//      chemistry::FragmentEnsemble retroda_reactants;
//      retroda_reactants.ReadMoreFromMdl( input);
//      io::File::CloseClearFStream( input);
//
//      chemistry::ReactionComplete retroda_rxn
//      (
//        sdf::RXNFactory::MakeReactionComplete( handler_retroda)
//      );
//
//      // 1 reactant, 1 product (sigmatropic rearrangement)
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "sigmatropic_rearrangement.rxn"));
//      sdf::RXNHandler handler_sigmatropic( input);
//      io::File::CloseClearFStream( input);
//
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "sigmatropic_reactants.sdf"));
//      chemistry::FragmentEnsemble sigmatropic_reactants;
//      sigmatropic_reactants.ReadMoreFromMdl( input);
//      io::File::CloseClearFStream( input);
//
//      chemistry::ReactionComplete sigmatropic_rxn
//      (
//        sdf::RXNFactory::MakeReactionComplete( handler_sigmatropic)
//      );
//
//      //
//      // Input and output
//      //
//
//      // Reactant/product matching
//      chemistry::FragmentEnsemble::const_iterator mol_itr( amide_reactants.Begin());
//      BCL_ExampleIndirectCheck( worker.MatchesReactants( *mol_itr, amide_rxn).GetSize(), 1, "Amide reaction, first molecule matches a reactant");
//      ++mol_itr;
//      BCL_ExampleIndirectCheck( worker.MatchesReactants( *mol_itr, amide_rxn).GetSize(), 1, "Amide reaction, second molecule matches a reactant");
//      BCL_ExampleIndirectCheck( worker.MatchesProducts( amide_rxn.GetProduct( 0), amide_rxn).GetSize(), 1, "Amide reaction, first product structure matches itself");
//
//      chemistry::FragmentEnsemble mol_output;
//
//      //
//      // reaction execution
//      //
//
//      BCL_MessageStd( "Amide reaction (2 reactants, 2 products)");
//      chemistry::FragmentEnsemble amide_res( worker.ExecuteReaction( amide_rxn, amide_reactants));
//      BCL_ExampleIndirectCheck( amide_res.GetSize(), 2, "ExecuteReaction() on amide reaction produced appropriate number of products");
//      mol_output.Append( amide_res);
//      amide_res = worker.React( amide_rxn, amide_reactants);
//      BCL_ExampleIndirectCheck( amide_res.GetSize(), 2, "React() gave the same output as ExecuteReaction() for a trivial ordering");
//
//      // test if a 2-reactant reaction can be run as an intramolecular reaction
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "4-amino-butanoic_acid.sdf"));
//      chemistry::FragmentEnsemble mols( input);
//      io::File::CloseClearFStream( input);
//      chemistry::FragmentEnsemble amide_intra_res
//      (
//        worker.ExecuteIntramolecularReaction( amide_rxn, *mols.Begin())
//      );
//      BCL_ExampleIndirectCheck( amide_intra_res.GetSize(), 2, "ExecuteIntramolecularReaction() produced appropriate number of products");
//      mol_output.Append( amide_intra_res);
//
//      // Reorder the reactants to ensure React() can figure the reaction out
//      chemistry::FragmentEnsemble reordered_reactants;
//      chemistry::FragmentEnsemble::const_iterator itr_mol( amide_reactants.End());
//      reordered_reactants.PushBack( *( --itr_mol));
//      reordered_reactants.PushBack( *( --itr_mol));
//
//      //
//      // perform the reactions
//      //
//
//      chemistry::FragmentEnsemble react_res( worker.React( amide_rxn, reordered_reactants));
//      BCL_ExampleIndirectCheck( react_res.GetSize(), 2, "React(), amide reaction, produced the correct number of products after reordering reactants");
//
//      chemistry::FragmentEnsemble click_res( worker.ExecuteReaction( click_rxn, click_reactants));
//      BCL_ExampleIndirectCheck( click_res.GetSize(), 1, "React(), click reaction, produced the correct number of products");
//      mol_output.Append( click_res);
//
//      chemistry::FragmentEnsemble retroda_res( worker.ExecuteReaction( retroda_rxn, retroda_reactants));
//      BCL_ExampleIndirectCheck( retroda_res.GetSize(), 2, "React(), retro diels-alder reaction, produced the correct number of products");
//      mol_output.Append( retroda_res);
//
//      chemistry::FragmentEnsemble sigmatropic_res( worker.ExecuteReaction( sigmatropic_rxn, sigmatropic_reactants));
//      BCL_ExampleIndirectCheck( sigmatropic_res.GetSize(), 1, "React(), sigmatropic reaction, produced the correct number of products");
//      mol_output.Append( sigmatropic_res);
//
//      io::OFStream output;
//      const std::string output_filename( AddExampleOutputPathToFilename( worker, "reaction_worker.sdf"));
//      BCL_ExampleMustOpenOutputFile( output, output_filename);
//
//      mol_output.WriteMDL( output);
//      io::File::CloseClearFStream( output);
//
//      std::string correct_filename( output_filename + ".correct");
//
//      // check if the generated file is right
//      if
//      (
//        BCL_ExampleCheck( io::File::FilesMatch( output_filename, correct_filename), true)
//      )
//      {
//        remove( output_filename.c_str());
//      }
//      else
//      {
//        BCL_MessageStd( output_filename + " did not match " + correct_filename);
//      }
//
//      return 0;
//    }
//
//    static const ExampleClass::EnumType s_Instance;
//
//   }; // class ExampleChemistryReactionWorker
//
//   const ExampleClass::EnumType ExampleChemistryReactionWorker::s_Instance
//   (
//     GetExamples().AddEnum( ExampleChemistryReactionWorker())
//   );
//
//} // namespace bcl
