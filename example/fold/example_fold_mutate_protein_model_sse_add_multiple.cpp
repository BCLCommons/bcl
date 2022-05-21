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
#include "fold/bcl_fold_mutate_protein_model_sse_add_multiple.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_placement_domain_using_fold_template.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_add_multiple.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSEAddMultiple :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSEAddMultiple *Clone() const
    {
      return new ExampleFoldMutateProteinModelSSEAddMultiple( *this);
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
      BCL_MessageStd( "testing default constructor");
      fold::MutateProteinModelSSEAddMultiple def_construct;

      // create the SSE pool
      assemble::SSEPool sse_pool( protein_model.GetSSEs());

      // create the placement
      const fold::PlacementDomainUsingFoldTemplate domain_placement;

      // test constructor from a SSEPool and PlacementDomainInterface
      BCL_MessageStd( "testing constructor from sse pool and placement");
      const fold::MutateProteinModelSSEAddMultiple domain_construct
      (
        util::ShPtr< assemble::SSEPool>( sse_pool.Clone()), domain_placement
      );

      // test constructor from a SSEPool and a ShPtr to a PlacementDomainInterface
      BCL_MessageStd( "testing constructor from sse pool and a ShPtr to a placement");
      const fold::MutateProteinModelSSEAddMultiple domain_interface_construct
      (
        util::ShPtr< assemble::SSEPool>( sse_pool.Clone()),
        util::ShPtr< fold::PlacementDomainInterface>( domain_placement.Clone())
      );

      // test copy constructor
      BCL_MessageStd( "testing copy constructor");
      const fold::MutateProteinModelSSEAddMultiple copy_construct( domain_construct);

      // test clone function
      BCL_MessageStd( "testing clone function");
      const util::ShPtr< fold::MutateProteinModelSSEAddMultiple> clone_construct( domain_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::fold::MutateProteinModelSSEAddMultiple");
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelSSEAddMultiple>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< fold::MutateProteinModelSSEAddMultiple>() +
        " but should give " + correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< fold::MutateProteinModelSSEAddMultiple>() == clone_construct->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_construct->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ///////////////
    // operators //
    ///////////////

      // check operator ()
      BCL_MessageStd( "testing () operator");

      assemble::ProteinModel empty_model( protein_model.GetEmptyChains());
      const assemble::ProteinModel mutated_model( *domain_construct( empty_model).GetArgument());
      const size_t known_number_of_sses( 5);
      BCL_Example_Check
      (
        mutated_model.GetSSEs().GetSize() == known_number_of_sses,
        "The mutated model has " + util::Format()( mutated_model.GetSSEs().GetSize()) + " SSEs but should have " +
        util::Format()( known_number_of_sses)
      );

      // check operator () with small template
      BCL_MessageStd( "testing () operator with small template");
      const fold::PlacementDomainUsingFoldTemplate small_template( -1);
      const fold::MutateProteinModelSSEAddMultiple add_using_small_template
      (
        util::ShPtr< assemble::SSEPool>( sse_pool.Clone()), small_template
      );
      const assemble::ProteinModel small_template_model( *add_using_small_template( empty_model).GetArgument());
      BCL_Example_Check
      (
        small_template_model.GetSSEs().GetSize() == known_number_of_sses - 1,
        "The mutated model has " + util::Format()( small_template_model.GetSSEs().GetSize()) +
        " SSEs but should have " + util::Format()( known_number_of_sses - 1)
      );

      // check operator () with large template
      BCL_MessageStd( "testing () operator with large template");
      const fold::PlacementDomainUsingFoldTemplate large_template( 1);
      const fold::MutateProteinModelSSEAddMultiple add_using_large_template
      (
        util::ShPtr< assemble::SSEPool>( sse_pool.Clone()), large_template
      );
      const assemble::ProteinModel large_template_model( *add_using_large_template( empty_model).GetArgument());
      BCL_Example_Check
      (
        large_template_model.GetSSEs().GetSize() == known_number_of_sses,
        "The mutated model has " + util::Format()( large_template_model.GetSSEs().GetSize()) +
        " SSEs but should have " + util::Format()( known_number_of_sses)
      );

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write mutated_sse_add_multiple.pdb");
      Proteins::WriteModelToPDB
      (
        mutated_model, AddExampleOutputPathToFilename( domain_construct, "mutated_sse_add_multiple.pdb")
      );

      // test = operator
      def_construct = domain_construct;

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( domain_construct);

      // read the object back in
      fold::MutateProteinModelSSEAddMultiple mutate_read;
      ReadBCLObject( mutate_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSEAddMultiple

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSEAddMultiple::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSEAddMultiple())
  );

} // namespace bcl
