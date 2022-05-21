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
#include "fold/bcl_fold_mutate_sse_type.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_sse_type.cpp
  //! @brief this examples demonstrates how sse types (and consequently the conformation) can be mutated
  //!
  //! @author woetzen
  //! @date Jun 16, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateSSEType :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateSSEType *Clone() const
    {
      return new ExampleFoldMutateSSEType( *this);
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
      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // get the SSEs
      util::ShPtr< assemble::SSE>
        sp_helix( Proteins::GetSSE( pdb_filename, 'A', 23, 34, biol::GetAAClasses().e_AABackBone)),
        sp_strand( Proteins::GetSSE( pdb_filename, 'A', 1, 7, biol::GetAAClasses().e_AABackBone));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from ss types
      fold::MutateSSEType mutate_sse_type
      (
        storage::Set< biol::SSType>::Create( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND, biol::GetSSTypes().COIL)
      );

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( mutate_sse_type.GetClassIdentifier(), GetStaticClassName< fold::MutateSSEType>());

      // scheme
      BCL_ExampleCheck( mutate_sse_type.GetScheme(), fold::MutateSSEType::GetDefaultScheme());

    ///////////////
    // operators //
    ///////////////

      storage::Map< biol::SSType, size_t> sstype_count;

      // mutate an sse
      size_t count( 0);
      while( count < 40)
      {
        sp_helix = mutate_sse_type( *sp_helix).GetArgument();
        ++sstype_count[ sp_helix->GetType()];
        ++count;
      }

      BCL_ExampleIndirectCheck
      (
        sstype_count[ biol::GetSSTypes().HELIX] > 10 && sstype_count[ biol::GetSSTypes().STRAND] && sstype_count[ biol::GetSSTypes().COIL],
        true,
        "not every ss type appeared at least 10 times: " + util::Format()( sstype_count)
      );

      util::ShPtr< assemble::SSE> sp_new_ss( mutate_sse_type( *sp_strand).GetArgument());
      BCL_ExampleCheck( sp_new_ss->GetType() != sp_strand->GetType(), true);

      Proteins::WriteSSEToPDB
      (
        sp_new_ss, AddExampleOutputPathToFilename( mutate_sse_type, "1ubi_1-7E_to_HC.pdb")
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // check read write
      WriteBCLObject( mutate_sse_type);
      fold::MutateSSEType mutate_read( storage::Set< biol::SSType>( biol::GetSSTypes().COIL));
      ReadBCLObject( mutate_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateSSEType

  const ExampleClass::EnumType ExampleFoldMutateSSEType::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateSSEType())
  );

} // namespace bcl
