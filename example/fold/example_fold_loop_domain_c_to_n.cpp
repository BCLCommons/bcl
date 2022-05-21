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
#include "fold/bcl_fold_loop_domain_c_to_n.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_loop_domain_c_to_n.cpp
  //!
  //! @author alexanns
  //! @date Jul 6, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldLoopDomainCToN :
    public ExampleInterface
  {
  public:

    ExampleFoldLoopDomainCToN *Clone() const
    {
      return new ExampleFoldLoopDomainCToN( *this);
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
      {
        // create list of loop segments
        storage::List< fold::LocatorLoopSegment> loop_segments;

        // 35-38
        loop_segments.PushBack( fold::LocatorLoopSegment( assemble::LocatorSSE( 'A', 35, 38), false));

        fold::LocatorLoopDomain loop_domain_locator( loop_segments, false);

        const std::string input_pdb( AddExampleInputPathToFilename( e_Biology, "2LZM_rotated_35_38.pdb"));

        util::ShPtr< assemble::ProteinModel> model( Proteins::GetModel( input_pdb).Clone());

        const util::ShPtr< fold::LoopDomain> loop_domain( loop_domain_locator.Locate( *model));

      /////////////////
      // data access //
      /////////////////

      ///////////////
      // operators //
      ///////////////

      ////////////////
      // operations //
      ////////////////

        // test GetResidues
        storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > residues( loop_domain->GetResidues());

        const size_t size( residues.GetSize());
        const size_t expected_size( 4);
        BCL_ExampleCheck( size, expected_size);
        const biol::AABase &last_aa( *residues.LastElement().First());
        BCL_ExampleCheck( last_aa.GetSeqID(), 38);
        BCL_ExampleCheck( last_aa.GetChainID(), 'A');

        // test GetMostProximalLoopSegmentAA
        {
          const util::ShPtr< biol::AABase> &aa( loop_domain->GetMostProximalLoopSegmentAA());
          BCL_ExampleCheck( aa->GetSeqID(), 38);
          BCL_ExampleCheck( aa->GetChainID(), 'A');
        }

        // test GetMostDistalLoopSegment
        {
          const assemble::SSE &sse( loop_domain->GetMostDistalLoopSegment());
          BCL_ExampleCheck( sse.GetFirstAA()->GetSeqID(), 35);
          BCL_ExampleCheck( sse.GetLastAA()->GetSeqID(), 38);
          BCL_ExampleCheck( sse.GetChainID(), 'A');
        }
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

  }; //end ExampleFoldLoopDomainCToN

  const ExampleClass::EnumType ExampleFoldLoopDomainCToN::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldLoopDomainCToN())
  );

} // namespace bcl
