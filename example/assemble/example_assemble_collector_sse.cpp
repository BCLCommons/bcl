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
#include "assemble/bcl_assemble_collector_sse.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_sse.cpp
  //!
  //! @author loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorSSE :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorSSE *Clone() const
    {
      return new ExampleAssembleCollectorSSE( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_unpaired_strand.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

      // create SiPtrVector "sse_vector" and initialize with all strands of "protein_model"
      util::SiPtrVector< const assemble::SSE> sse_vector( protein_model.GetSSEs());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::CollectorSSE collector;

    ////////////////
    // operations //
    ////////////////

      // test Collect function
      BCL_MessageStd( "test Collect function");
      util::SiPtrList< const assemble::SSE> sse_collection( collector.Collect( protein_model));

      // make sure that the correct number of SSEs are collected
      BCL_Example_Check
      (
        sse_collection.GetSize() == sse_vector.GetSize(),
        "number of SSEs collected is " + util::Format()( sse_collection.GetSize()) +
        " but should be" + util::Format()( sse_vector.GetSize())
      );

      // create ss_types vector
      storage::Set< biol::SSType> ss_types;
      ss_types.InsertElement( biol::GetSSTypes().HELIX);

      // collect helices
      util::SiPtrVector< const assemble::SSE> sse_vector_helix( protein_model.GetSSEs( ss_types));

      // test constructor from a set of SSTypes
      BCL_MessageStd
      (
        "test constructor from a set of SSTypes"
      );
      assemble::CollectorSSE collector_b( ss_types);

      // check that the constructor worked properly
      BCL_Example_Check
      (
        collector_b.GetSSTypes().GetSize() == 1 && *collector_b.GetSSTypes().Begin() == biol::GetSSTypes().HELIX,
        " the constructor that takes a Set of SSTypes did not work correctly"
      );

      BCL_MessageStd( "test Collect function");
      util::SiPtrList< const assemble::SSE> sse_collection_b( collector_b.Collect( protein_model));
      // make sure that the correct number of SSEs are collected
      BCL_Example_Check
      (
        sse_collection_b.GetSize() == sse_vector_helix.GetSize(),
        "number of SSEs collected is " + util::Format()( sse_collection_b.GetSize())
        + " but should be " + util::Format()( sse_vector_helix.GetSize())
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorSSE

  const ExampleClass::EnumType ExampleAssembleCollectorSSE::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorSSE())
  );

} // namespace bcl
