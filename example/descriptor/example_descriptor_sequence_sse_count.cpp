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
#include "descriptor/bcl_descriptor_sequence_sse_count.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_sequence_sse_count.cpp
  //!
  //! @author mendenjl
  //! @date Feb 13, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorSequenceSSECount :
    public ExampleInterface
  {
  public:

    ExampleDescriptorSequenceSSECount *Clone() const
    {
      return new ExampleDescriptorSequenceSSECount( *this);
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

      // create a descriptor that will retrieve the SSE counts both inside and outside the membrane
      descriptor::SequenceSSECount sse_counts_retriever( false, false);
      BCL_ExampleCheck
      (
        sse_counts_retriever.TryRead( util::ObjectDataLabel( "(method=OCTOPUS,min helix size=10)"), util::GetLogger()),
        true
      );
      // call set object with protein
      sse_counts_retriever.SetObject( protein);
      // create a descriptor that will retrieve the SSE counts only inside the membrane
      descriptor::SequenceSSECount sse_membrane_counts_retriever( true, false);
      BCL_ExampleCheck
      (
        sse_membrane_counts_retriever.TryRead
        (
          util::ObjectDataLabel( "(method=OCTOPUS,membrane method=OCTOPUS,min helix size=10)"),
          util::GetLogger()
        ),
        true
      );
      // call set object with protein
      sse_counts_retriever.SetObject( protein);
      sse_membrane_counts_retriever.SetObject( protein);

    /////////////////
    // data access //
    /////////////////

      // test the GetLength function
      BCL_ExampleCheck( sse_counts_retriever.GetSizeOfFeatures(), 3);
      BCL_ExampleCheck( sse_membrane_counts_retriever.GetSizeOfFeatures(), 3);

    ///////////////
    // operators //
    ///////////////

      descriptor::Iterator< biol::AABase> itr( sse_counts_retriever.GetType(), protein);
      // operator()
      BCL_ExampleCheck
      (
        sse_counts_retriever( itr),
        linal::MakeVector< float>( 1.0, 0.0, 2.0)
      );
      BCL_ExampleCheck
      (
        sse_membrane_counts_retriever( itr),
        linal::MakeVector< float>( 1.0, 0.0, 0.0)
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

  }; //end ExampleDescriptorSequenceSSECount

  const ExampleClass::EnumType ExampleDescriptorSequenceSSECount::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorSequenceSSECount())
  );

} // namespace bcl
