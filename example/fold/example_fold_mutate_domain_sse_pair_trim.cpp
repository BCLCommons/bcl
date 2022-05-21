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
#include "fold/bcl_fold_mutate_domain_sse_pair_trim.h"

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
  //! @example example_fold_mutate_domain_sse_pair_trim.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Sep 4, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateDomainSSEPairTrim :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateDomainSSEPairTrim *Clone() const
    {
      return new ExampleFoldMutateDomainSSEPairTrim( *this);
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

      assemble::LocatorSSE sse_locator_b( 'B', 122, 131);

      const util::SiPtr< const assemble::SSE> located_sse_b( sse_locator_b.Locate( protein_model));

      // make domain
      assemble::Domain domain
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::Create
        (
          util::ShPtr< assemble::SSE>( located_sse->Clone()), util::ShPtr< assemble::SSE>( located_sse_b->Clone())
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      fold::MutateDomainSSEPairTrim default_constr;

      // constructor taking member variable parameters
      storage::Map< biol::SSType, size_t> min_sse_size;
      min_sse_size =
        storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().HELIX, 5));
      fold::MutateDomainSSEPairTrim param_constr( 4, min_sse_size, "scheme");

      // clone constructor
      util::ShPtr< fold::MutateDomainSSEPairTrim> clone_constr( param_constr.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // test operator
      {
        BCL_MessageDbg( "test operator");
        const math::MutateResult< assemble::Domain> mutate( param_constr( domain));
        BCL_ExampleAssert( mutate.GetArgument().IsDefined(), true);
        const assemble::Domain &mutated_domain( *mutate.GetArgument());
        const util::SiPtrVector< const assemble::SSE> sses( mutated_domain.GetSSEs());
        BCL_ExampleAssert( sses.GetSize(), 2);

        BCL_MessageDbg( "starting sse a is " + located_sse->GetIdentification());
        BCL_MessageDbg( "starting sse b is " + located_sse_b->GetIdentification());
        BCL_MessageDbg( "ending sse a is " + sses.FirstElement()->GetIdentification());
        BCL_MessageDbg( "ending sse b is " + sses.LastElement()->GetIdentification());
        BCL_ExampleCheck( sses.FirstElement()->GetIdentification(), "HELIX B   91 ASP <==>  105 ASN");
        BCL_ExampleCheck( sses.LastElement()->GetIdentification(), "HELIX B  122 ASN <==>  131 ASP");
      }

      // test operator don't use all possible residues to trim
      {
        BCL_MessageDbg( "test operator don't use all possible residues to trim");
        assemble::LocatorSSE sse_locator( 'B', 289, 300);

        const util::SiPtr< const assemble::SSE> located_sse( sse_locator.Locate( protein_model));

        assemble::LocatorSSE sse_locator_b( 'B', 304, 313);

        const util::SiPtr< const assemble::SSE> located_sse_b( sse_locator_b.Locate( protein_model));

        // make domain
        assemble::Domain domain
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::Create
          (
            util::ShPtr< assemble::SSE>( located_sse->Clone()), util::ShPtr< assemble::SSE>( located_sse_b->Clone())
          )
        );
        fold::MutateDomainSSEPairTrim param_constr( 5, min_sse_size, "scheme");
        const math::MutateResult< assemble::Domain> mutate( param_constr( domain));
        BCL_ExampleAssert( mutate.GetArgument().IsDefined(), true);
        const assemble::Domain &mutated_domain( *mutate.GetArgument());
        const util::SiPtrVector< const assemble::SSE> sses( mutated_domain.GetSSEs());
        BCL_ExampleAssert( sses.GetSize(), 2);

        BCL_MessageDbg( "starting sse a is " + located_sse->GetIdentification());
        BCL_MessageDbg( "starting sse b is " + located_sse_b->GetIdentification());
        BCL_MessageDbg( "ending sse a is " + sses.FirstElement()->GetIdentification());
        BCL_MessageDbg( "ending sse b is " + sses.LastElement()->GetIdentification());
        BCL_ExampleCheck( sses.FirstElement()->GetIdentification(), "HELIX B  289 ALA <==>  297 GLY");
        BCL_ExampleCheck( sses.LastElement()->GetIdentification(), "HELIX B  304 SER <==>  313 VAL");
      }

      // test operator with possibility to remove too many residues from sse
      {
        BCL_MessageDbg( "test operator with possibility to remove too many residues from sse");
        assemble::LocatorSSE sse_locator( 'B', 289, 300);

        const util::SiPtr< const assemble::SSE> located_sse( sse_locator.Locate( protein_model));

        assemble::LocatorSSE sse_locator_b( 'B', 304, 313);

        const util::SiPtr< const assemble::SSE> located_sse_b( sse_locator_b.Locate( protein_model));

        // make domain
        assemble::Domain domain
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::Create
          (
            util::ShPtr< assemble::SSE>( located_sse->Clone()), util::ShPtr< assemble::SSE>( located_sse_b->Clone())
          )
        );
        fold::MutateDomainSSEPairTrim param_constr( 10, min_sse_size, "scheme");
        const math::MutateResult< assemble::Domain> mutate( param_constr( domain));
        BCL_ExampleAssert( mutate.GetArgument().IsDefined(), true);
        const assemble::Domain &mutated_domain( *mutate.GetArgument());
        const util::SiPtrVector< const assemble::SSE> sses( mutated_domain.GetSSEs());
        BCL_ExampleAssert( sses.GetSize(), 2);

        BCL_MessageDbg( "starting sse a is " + located_sse->GetIdentification());
        BCL_MessageDbg( "starting sse b is " + located_sse_b->GetIdentification());
        BCL_MessageDbg( "ending sse a is " + sses.FirstElement()->GetIdentification());
        BCL_MessageDbg( "ending sse b is " + sses.LastElement()->GetIdentification());
        BCL_ExampleCheck( sses.FirstElement()->GetIdentification(), "HELIX B  289 ALA <==>  294 HIS");
        BCL_ExampleCheck( sses.LastElement()->GetIdentification(), "HELIX B  308 VAL <==>  313 VAL");
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      WriteBCLObject( *clone_constr);
      fold::MutateDomainSSEPairTrim read_mutate;
      ReadBCLObject( read_mutate);

      {
        BCL_MessageDbg( "test operator");
        const math::MutateResult< assemble::Domain> mutate( read_mutate( domain));
        BCL_ExampleAssert( mutate.GetArgument().IsDefined(), true);
        const assemble::Domain &mutated_domain( *mutate.GetArgument());
        const util::SiPtrVector< const assemble::SSE> sses( mutated_domain.GetSSEs());
        BCL_ExampleAssert( sses.GetSize(), 2);

        BCL_MessageDbg( "starting sse a is " + located_sse->GetIdentification());
        BCL_MessageDbg( "starting sse b is " + located_sse_b->GetIdentification());
        BCL_MessageDbg( "ending sse a is " + sses.FirstElement()->GetIdentification());
        BCL_MessageDbg( "ending sse b is " + sses.LastElement()->GetIdentification());
        BCL_ExampleCheck( sses.FirstElement()->GetIdentification(), "HELIX B   91 ASP <==>  105 ASN");
        BCL_ExampleCheck( sses.LastElement()->GetIdentification(), "HELIX B  122 ASN <==>  131 ASP");
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateDomainSSEPairTrim

  const ExampleClass::EnumType ExampleFoldMutateDomainSSEPairTrim::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateDomainSSEPairTrim())
  );
  
} // namespace bcl
