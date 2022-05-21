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
#include "descriptor/bcl_descriptor_aa_sse_info.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_sse_info.cpp
  //!
  //! @author teixeipl, mendenjl
  //! @date Mar 09, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAASSEInfo :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAASSEInfo *Clone() const
    {
      return new ExampleDescriptorAASSEInfo( *this);
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
      // construct an implementation that will read proteins from the biol directory
      const std::string pdb_fasta_dir( AddExampleInputPathToFilename( e_Biology, ""));

      // create an object to read proteins, as either pdbs or fastas, from the directory
      util::Implementation
      <
        assemble::RetrieveProteinModelWithCache
      > retriever_all( "SequenceDirectory( " + pdb_fasta_dir + ")");

      // retrieve 1b43A
      assemble::ProteinModelWithCache protein( *retriever_all->Retrieve( "1b43A"));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////
      // create a descriptor that will retrieve the SSE size
      descriptor::AASSEInfo sse_size_retriever( descriptor::AASSEInfo::e_AASSESize);
      // call set object with protein
      sse_size_retriever.SetObject( protein);

      // create a descriptor that will retrieve the SSE ID
      descriptor::AASSEInfo sse_id_retriever( descriptor::AASSEInfo::e_AASSEID);
      // call set object with protein
      sse_id_retriever.SetObject( protein);

      // create a descriptor that will retrieve the AA position
      descriptor::AASSEInfo sse_pos_retriever( descriptor::AASSEInfo::e_AAPositionInSSE);
      // call set object with protein
      sse_pos_retriever.SetObject( protein);

      // create a descriptor that will retrieve the distance of an AA from the center of its SSE
      descriptor::AASSEInfo sse_dist_from_center_retriever( descriptor::AASSEInfo::e_AADistanceFromSSECenter);
      // call set object with protein
      sse_dist_from_center_retriever.SetObject( protein);

      // clone
      util::ShPtr< descriptor::AASSEInfo> sp_aa_info( sse_size_retriever.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( sse_size_retriever.GetSizeOfFeatures(), 1);

    ///////////////
    // operators //
    ///////////////

      descriptor::Iterator< biol::AABase> itr( sse_size_retriever.GetType(), protein);
      // operator()
      BCL_ExampleCheck
      (
        sse_size_retriever( itr),
        linal::Vector< float>( size_t( 1), 63)
      );
      itr.GotoPosition( 63);
      BCL_ExampleCheck
      (
        sse_size_retriever( itr),
        linal::Vector< float>( size_t( 1), 21)
      );

      descriptor::Iterator< biol::AABase> itr_id( sse_id_retriever.GetType(), protein);
      // operator()
      BCL_ExampleCheck
      (
        sse_id_retriever( itr_id),
        linal::Vector< float>( size_t( 1), 0.0)
      );
      itr_id.GotoPosition( 63);
      BCL_ExampleCheck
      (
        sse_id_retriever( itr_id),
        linal::Vector< float>( size_t( 1), 1.0)
      );

      descriptor::Iterator< biol::AABase> itr_pos( sse_pos_retriever.GetType(), protein);
      // operator()
      BCL_ExampleCheck
      (
        sse_pos_retriever( itr_pos),
        linal::Vector< float>( size_t( 1), 0.0)
      );
      itr_pos.GotoPosition( 64);
      BCL_ExampleCheck
      (
        sse_pos_retriever( itr_pos),
        linal::Vector< float>( size_t( 1), 1 / float( 20.0))
      );

      descriptor::Iterator< biol::AABase> itr_dist( sse_dist_from_center_retriever.GetType(), protein);
      // operator()
      BCL_ExampleCheck
      (
        sse_dist_from_center_retriever( itr_dist),
        linal::Vector< float>( size_t( 1), 1.0)
      );
      itr_dist.GotoPosition( 64);
      BCL_ExampleCheck
      (
        sse_dist_from_center_retriever( itr_dist),
        linal::Vector< float>( size_t( 1), -( ( 1 / 20.0) - 0.5) * 2)
      );
    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAASSEInfo

  const ExampleClass::EnumType ExampleDescriptorAASSEInfo::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAASSEInfo())
  );

} // namespace bcl
