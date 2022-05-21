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
//#include "chemistry/bcl_chemistry_fragment_ring_swap.h"
//
//// includes from bcl - sorted alphabetically
//#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
//#include "chemistry/bcl_chemistry_collector_valence.h"
//#include "chemistry/bcl_chemistry_pick_atom_random.h"
//#include "chemistry/bcl_chemistry_pick_fragment_random.h"
//#include "io/bcl_io_file.h"
//
//// external includes - sorted alphabetically
//#include <stdio.h>
//
//namespace bcl
//{
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  //!
//  //! @example example_chemistry_fragment_ring_swap.cpp
//  //!
//  //! @author mendenjl
//  //! @date Sep 10, 2019
//  //! @remarks status complete
//  //! @remarks reviewed by nobody on
//  //!
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  class ExampleChemistryFragmentRingSwap :
//    public ExampleInterface
//  {
//  public:
//
//    ExampleChemistryFragmentRingSwap *Clone() const
//    {
//      return new ExampleChemistryFragmentRingSwap( *this);
//    }
//
//  /////////////////
//  // data access //
//  /////////////////
//
//    //! @brief returns class name
//    //! @return the class name as const ref std::string
//    const std::string &GetClassIdentifier() const
//    {
//      return GetStaticClassName( *this);
//    }
//
//    int Run() const
//    {
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // construct from properties
//      chemistry::FragmentRingSwap fragment_grow;
//
//      // clone function
//      util::ShPtr< chemistry::FragmentRingSwap> fragment_grow_clone( fragment_grow.Clone());
//      BCL_ExampleCheck( fragment_grow_clone.IsDefined(), true);
//
//      // create input stream for reading a fragment ensemble
//      io::IFStream input;
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
//      // close stream
//      // creating ShPtr of growfragments
//      util::ShPtr< chemistry::FragmentEnsemble> sp_fragment_pool( new chemistry::FragmentEnsemble( input, sdf::e_Remove));
//      io::File::CloseClearFStream( input);
//      // remove hydrogens to make valency
//      sp_fragment_pool->RemoveH();
//
//      util::ShPtr< chemistry::SearchFragmentLibraryFromTree> tree_search;
//      // create a fragment grower; do not attempt to refine geometry to avoid convoluting changes in sample conformations with
//      // changes in the output of this class/example
//      chemistry::FragmentRingSwap fragment_grower( tree_search, "", chemistry::FragmentComplete(), chemistry::FragmentComplete(), storage::Vector< size_t>(), false);
//      io::OFStream output;
//      std::string fname( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "test_fragment_ring_swap.sdf"));
//      std::string fname_correct
//      (
//        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "test_fragment_ring_swap.correct.sdf")
//      );
//      // output on apple differs due to the std::random_shuffle function being implemented differently on our clang version
//      std::string fname_correct_mac
//      (
//        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "test_fragment_ring_swap.correct.mac.sdf")
//      );
//      BCL_ExampleMustOpenOutputFile
//      (
//        output,
//        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "test_fragment_ring_swap.sdf")
//      );
//      chemistry::FragmentComplete test( *sp_fragment_pool->Begin());
//      for( size_t i( 0), sz( 100); i < sz; ++i)
//      {
//        auto new_test( fragment_grower( test).GetArgument());
//        if( new_test.IsDefined() && new_test->GetNumberAtoms())
//        {
//          test = *new_test;
//          test.WriteMDL( output);
//        }
//      }
//      io::File::CloseClearFStream( output);
//
//      if
//      (
//        io::File::FilesMatchWithinAbsoluteTolerance
//        (
//          fname,
//          fname_correct,
//          0.1
//        )
//      )
//      {
//        BCL_ExampleCheck
//        (
//          io::File::FilesMatchWithinAbsoluteTolerance
//          (
//            fname,
//            fname_correct,
//            0.1
//          ),
//          true
//        );
//        remove( fname.c_str());
//      }
//      else if
//      (
//        BCL_ExampleCheck
//        (
//          io::File::FilesMatchWithinAbsoluteTolerance
//          (
//            fname,
//            fname_correct_mac,
//            0.1
//          ),
//          true
//        )
//      )
//      {
//        remove( fname.c_str());
//      }
//
//    /////////////////
//    // data access //
//    /////////////////
//
//    ///////////////
//    // operators //
//    ///////////////
//
//    ////////////////
//    // operations //
//    ////////////////
//
//    //////////////////////
//    // input and output //
//    //////////////////////
//
//    //////////////////////
//    // helper functions //
//    //////////////////////
//
//      return 0;
//    } // Run
//
//    static const ExampleClass::EnumType s_Instance;
//
//  }; //end ExampleChemistryFragmentRingSwap
//
//  const ExampleClass::EnumType ExampleChemistryFragmentRingSwap::s_Instance
//  (
//    GetExamples().AddEnum( ExampleChemistryFragmentRingSwap())
//  );
//
//} // namespace bcl
