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
#include "descriptor/bcl_descriptor_aa_distance_from_aa_type.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_distance_from_aa_type.cpp
  //!
  //! @author mendenjl
  //! @date May 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAADistanceFromAAType :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAADistanceFromAAType *Clone() const
    {
      return new ExampleDescriptorAADistanceFromAAType( *this);
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
      descriptor::AADistanceFromAAType proline_retriever( false);
      // find the two nearest prolines, up to 40 AAs away
      BCL_ExampleCheck
      (
        proline_retriever.TryRead
        (
          util::ObjectDataLabel( "(type=PROLINE,max=40,direction=Center,size=2)"),
          util::GetLogger()
        ),
        true
      );
      // call set object with protein
      proline_retriever.SetObject( protein);

    /////////////////
    // data access //
    /////////////////

      // test  the GetLength function
      BCL_ExampleCheck( proline_retriever.GetSizeOfFeatures(), 2);

    ///////////////
    // operators //
    ///////////////

      descriptor::Iterator< biol::AABase> itr( proline_retriever.GetType(), protein);
      // operator()
      BCL_ExampleCheck
      (
        proline_retriever( itr),
        linal::MakeVector< float>( 3.0, 9.0)
      );
      itr.GotoPosition( 63);
      BCL_ExampleCheck
      (
        proline_retriever( itr),
        linal::MakeVector< float>( 11.0, 18.0)
      );
      itr.GotoPosition( 121);
      BCL_ExampleCheck
      (
        proline_retriever( itr),
        linal::MakeVector< float>( 23.0, 28.0)
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

  }; //end ExampleDescriptorAADistanceFromAAType

  const ExampleClass::EnumType ExampleDescriptorAADistanceFromAAType::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAADistanceFromAAType())
  );

} // namespace bcl
