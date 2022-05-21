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
#include "score/bcl_score_protein_model_loop_domain_closure.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_model_loop_domain_closure.cpp
  //!
  //! @author alexanns, fischea
  //! @date January 4, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreProteinModelLoopDomainClosure :
    public ExampleInterface
  {
  public:

    ExampleScoreProteinModelLoopDomainClosure *Clone() const
    {
      return new ExampleScoreProteinModelLoopDomainClosure( *this);
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
      // create protein model
      assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")));

      const assemble::LocatorAA anchor_locator( 'A', 7);
      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel>
      > target_residue_locator( new assemble::LocatorAA( 'A', 8));

      // N terminal SSE locator
      assemble::LocatorSSE n_terminal_sse_locator( 'A', 1, 7);

      // loop segment locators
      storage::List< fold::LocatorLoopSegment> loop_segment_locators
      (
        1, fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 8, 9), false)
      );

      util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > sp_locators
      (
        new util::ShPtrList< fold::LocatorLoopDomain>
        (
          1, util::ShPtr< fold::LocatorLoopDomain>
          (
            new fold::LocatorLoopDomain( loop_segment_locators)
          )
        )
      );

      // create model data
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());
      sp_model_data->Insert( assemble::ProteinModelData::e_LoopDomainLocators, sp_locators);
      model.SetProteinModelData( sp_model_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::ProteinModelLoopDomainClosure def_constr;

      // constructor taking member parameters
      score::ProteinModelLoopDomainClosure param_constr;
      const double correct_score( -3.46931);
      {
        const double score( param_constr( model));
        BCL_MessageDbg( "score is " + util::Format()( score));
        BCL_ExampleCheckWithinTolerance( score, correct_score, 0.00001);
      }

      // clone constructor
      util::ShPtr< score::ProteinModelLoopDomainClosure> clone_constr( param_constr.Clone());
      {
        const double score( clone_constr->operator()( model));
        BCL_MessageDbg( "score is " + util::Format()( score));
        BCL_ExampleCheckWithinTolerance( score, correct_score, 0.00001);
      }

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< score::ProteinModelLoopDomainClosure>(), clone_constr->GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

      // test scoring operator
      {
        const double score( param_constr( model));
        BCL_MessageDbg( "score is " + util::Format()( score));
        BCL_ExampleCheckWithinTolerance( score, correct_score, 0.00001);
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write the score and then read it back in
      WriteBCLObject( param_constr);
      score::ProteinModelLoopDomainClosure score_read;
      ReadBCLObject( score_read);
      // test reading and writing
      {
        const double score( score_read( model));
        BCL_MessageDbg( "score is " + util::Format()( score));
        BCL_ExampleCheckWithinTolerance( score, correct_score, 0.00001);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreProteinModelLoopDomainClosure

  const ExampleClass::EnumType ExampleScoreProteinModelLoopDomainClosure::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelLoopDomainClosure())
  );

} // namespace bcl
