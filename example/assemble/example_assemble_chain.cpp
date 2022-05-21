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
#include "assemble/bcl_assemble_chain.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_chain.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleChain :
    public ExampleInterface
  {
  public:

    ExampleAssembleChain *Clone() const
    {
      return new ExampleAssembleChain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::Chain def_construct;
      BCL_ExampleIndirectCheck
      (
        def_construct.GetData().IsEmpty() && !def_construct.GetSequence().IsDefined(),
        true, "default constructor"
      );

      // test constructor from a sequence and domain
      BCL_MessageStd( "testing constructor from a sequence and domain");
      const util::ShPtr< biol::AASequence> sp_sequence( protein_model.GetChain( 'A')->GetSequence());
      const assemble::Domain source_domain( protein_model.GetSSEsAsDomain());
      const assemble::Chain domain_construct( sp_sequence, source_domain);

      // test constructor from ShPtr of AAsequence
      assemble::Chain sequence_construct( sp_sequence);
      BCL_ExampleIndirectCheck
      (
        sequence_construct.GetData().IsEmpty() && sequence_construct.GetSequence().IsDefined(),
        true, "constructor from AAsequence"
      );

      // test copy constructor
      assemble::Chain copy_construct( domain_construct);
      BCL_ExampleIndirectCheck
      (
        domain_construct.GetSequence() == copy_construct.GetSequence() &&
        domain_construct.GetData().InternalData() == copy_construct.GetData().InternalData(),
        true, "copy constructor"
      );

      // test Clone function
      util::ShPtr< assemble::Chain> clone_construct( domain_construct.Clone());
      BCL_ExampleIndirectCheck
      (
        domain_construct.GetSequence() == clone_construct->GetSequence() &&
        domain_construct.GetData().InternalData() == clone_construct->GetData().InternalData(),
        true, "clone constructor"
      );

      // test HardCopy function
      util::ShPtr< assemble::Chain> hard_copy_construct( domain_construct.HardCopy());
      BCL_ExampleIndirectCheck
      (
        domain_construct.GetSequence()->Sequence() == hard_copy_construct->GetSequence()->Sequence() &&
        domain_construct.GetData().GetSize() == hard_copy_construct->GetData().GetSize() &&
        ( **domain_construct.GetData().Begin()).Sequence() == ( **hard_copy_construct->GetData().Begin()).Sequence(),
        true, "hard copy constructor"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetChainID
      BCL_ExampleCheck( domain_construct.GetChainID(), 'A');

      //check SetChainID
      copy_construct.SetChainID( 'F');
      BCL_ExampleIndirectCheck( copy_construct.GetChainID(), 'F', "SetChainID");

      //check GetSequence
      BCL_ExampleCheck( domain_construct.GetSequence(), sp_sequence);

    ////////////////
    // operations //
    ////////////////

      assemble::Chain insert_single_sse( sequence_construct);
      assemble::Chain insert_multiple_sses( sequence_construct);
      assemble::Chain insert_new_domain( sequence_construct);

      // check Insert with NewSSE
      BCL_ExampleCheck
      (
        insert_single_sse.Insert( util::ShPtr< assemble::SSE>( protein_model.GetSSEs().FirstElement()->Clone())), true
      );

      // check Insert with NewSSEVector
      const util::SiPtrVector< const assemble::SSE> sses( protein_model.GetSSEs());
      util::ShPtrVector< assemble::SSE> sse_vector;
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        sse_vector.PushBack( util::ShPtr< assemble::SSE>( ( *sse_itr)->Clone()));
      }
      BCL_ExampleCheck( insert_multiple_sses.Insert( sse_vector), true);

      // check Insert with NewDomain
      BCL_ExampleCheck( insert_new_domain.Insert( source_domain), true);

      // check AddLoops
      copy_construct.AddLoops( true, false);
      BCL_ExampleIndirectCheck
      (
        copy_construct.GetSSEs( biol::GetSSTypes().COIL).GetSize(), 5, "AddLoops"
      );

    ///////////////
    // operators //
    ///////////////

      // check = operator
      def_construct = domain_construct;
      BCL_ExampleIndirectCheck
      (
        domain_construct.GetSequence() == def_construct.GetSequence() &&
        domain_construct.GetData().InternalData() == def_construct.GetData().InternalData(),
        true, "= operator"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( domain_construct);

      // read the object back in
      assemble::Chain chain_read;
      ReadBCLObject( chain_read);

      BCL_ExampleIndirectCheck( domain_construct.GetNumberSSEs(), chain_read.GetNumberSSEs(), "read and write");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleChain

  const ExampleClass::EnumType ExampleAssembleChain::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleChain())
  );

} // namespace bcl

