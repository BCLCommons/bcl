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
#include "fold/bcl_fold_mutate_domain_sse_split.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_domain_sse_split.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Sep 3, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateDomainSSESplit :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateDomainSSESplit *Clone() const
    {
      return new ExampleFoldMutateDomainSSESplit( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      assemble::LocatorSSE sse_locator( 'B', 91, 109);

      const util::SiPtr< const assemble::SSE> located_sse( sse_locator.Locate( protein_model));

      // make domain
      assemble::Domain domain
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
        (
          util::ShPtr< assemble::SSE>( located_sse->Clone())
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      fold::MutateDomainSSESplit default_constr;

      // constructor taking member variable parameters
      const std::string scheme( "scheme");
      fold::MutateDomainSSESplit param_constr( scheme, 2.0);
      BCL_ExampleCheck( param_constr.GetScheme(), scheme);

      // clone constructor
      util::ShPtr< fold::MutateDomainSSESplit> clone_constr( param_constr.Clone());

    /////////////////
    // data access //
    /////////////////

      // get scheme
      BCL_ExampleCheck( param_constr.GetScheme(), scheme);

    ///////////////
    // operators //
    ///////////////

      // operator
      {
        BCL_MessageDbg( "test operator");

        // mutate domain
        math::MutateResult< assemble::Domain> mutate_result( param_constr( domain));
        BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);

        // get the sses from this domain
        const util::SiPtrVector< const assemble::SSE> domain_sses( mutate_result.GetArgument()->GetSSEs());
        BCL_ExampleCheck( domain_sses.GetSize(), 2);

        BCL_MessageDbg( "starting sse is " + located_sse->GetIdentification());
        BCL_MessageDbg( "first split sse is " + domain_sses.FirstElement()->GetIdentification());
        BCL_MessageDbg( "second split sse is " + domain_sses.LastElement()->GetIdentification());

        BCL_ExampleCheck( domain_sses.FirstElement()->GetFirstAA()->GetSeqID(), 91);
        BCL_ExampleCheck( domain_sses.FirstElement()->GetLastAA()->GetSeqID(), 102);
        BCL_ExampleCheck( domain_sses.LastElement()->GetFirstAA()->GetSeqID(), 103);
        BCL_ExampleCheck( domain_sses.LastElement()->GetLastAA()->GetSeqID(), 109);
      }
      // operator with 2 residue sse
      {
        BCL_MessageDbg( "test operator with 2 residue sse");
        assemble::SSE short_sse( located_sse->SubSequence( 0, 2), located_sse->GetType());
        assemble::Domain short_domain
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>
          (
            util::ShPtr< assemble::SSE>( short_sse.Clone())
          )
        );
        // mutate domain
         math::MutateResult< assemble::Domain> mutate_result( param_constr( short_domain));
         BCL_ExampleCheck( mutate_result.GetArgument().IsDefined(), true);

         // get the sses from this domain
         const util::SiPtrVector< const assemble::SSE> domain_sses( mutate_result.GetArgument()->GetSSEs());
         BCL_ExampleCheck( domain_sses.GetSize(), 2);

         BCL_MessageDbg( "starting sse is " + located_sse->GetIdentification());
         BCL_MessageDbg( "first split sse is " + domain_sses.FirstElement()->GetIdentification());
         BCL_MessageDbg( "second split sse is " + domain_sses.LastElement()->GetIdentification());

         BCL_ExampleCheck( domain_sses.FirstElement()->GetFirstAA()->GetSeqID(), 91);
         BCL_ExampleCheck( domain_sses.FirstElement()->GetLastAA()->GetSeqID(), 91);
         BCL_ExampleCheck( domain_sses.LastElement()->GetFirstAA()->GetSeqID(), 92);
         BCL_ExampleCheck( domain_sses.LastElement()->GetLastAA()->GetSeqID(), 92);
      }
      // operator with different standard deviation
      {
        BCL_MessageDbg( "test operator with different standard deviation");
        fold::MutateDomainSSESplit mutate( scheme, 200.0);
        // mutate domain
        math::MutateResult< assemble::Domain> mutate_result( mutate( domain));
        BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);

        // get the sses from this domain
        const util::SiPtrVector< const assemble::SSE> domain_sses( mutate_result.GetArgument()->GetSSEs());
        BCL_ExampleCheck( domain_sses.GetSize(), 2);

        BCL_MessageDbg( "starting sse is " + located_sse->GetIdentification());
        BCL_MessageDbg( "first split sse is " + domain_sses.FirstElement()->GetIdentification());
        BCL_MessageDbg( "second split sse is " + domain_sses.LastElement()->GetIdentification());

        BCL_ExampleCheck( domain_sses.FirstElement()->GetFirstAA()->GetSeqID(), 91);
        BCL_ExampleCheck( domain_sses.FirstElement()->GetLastAA()->GetSeqID(), 102);
        BCL_ExampleCheck( domain_sses.LastElement()->GetFirstAA()->GetSeqID(), 103);
        BCL_ExampleCheck( domain_sses.LastElement()->GetLastAA()->GetSeqID(), 109);
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      WriteBCLObject( param_constr);
      fold::MutateDomainSSESplit read_mutate;
      ReadBCLObject( read_mutate);
      {
        BCL_MessageDbg( "test Write and read");
        // mutate domain
        math::MutateResult< assemble::Domain> mutate_result( param_constr( domain));
        BCL_ExampleAssert( mutate_result.GetArgument().IsDefined(), true);

        // get the sses from this domain
        const util::SiPtrVector< const assemble::SSE> domain_sses( mutate_result.GetArgument()->GetSSEs());
        BCL_ExampleCheck( domain_sses.GetSize(), 2);

        BCL_MessageDbg( "starting sse is " + located_sse->GetIdentification());
        BCL_MessageDbg( "first split sse is " + domain_sses.FirstElement()->GetIdentification());
        BCL_MessageDbg( "second split sse is " + domain_sses.LastElement()->GetIdentification());

        BCL_ExampleCheck( domain_sses.FirstElement()->GetFirstAA()->GetSeqID(), 91);
        BCL_ExampleCheck( domain_sses.FirstElement()->GetLastAA()->GetSeqID(), 105);
        BCL_ExampleCheck( domain_sses.LastElement()->GetFirstAA()->GetSeqID(), 106);
        BCL_ExampleCheck( domain_sses.LastElement()->GetLastAA()->GetSeqID(), 109);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateDomainSSESplit

  const ExampleClass::EnumType ExampleFoldMutateDomainSSESplit::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateDomainSSESplit())
  );
  
} // namespace bcl
