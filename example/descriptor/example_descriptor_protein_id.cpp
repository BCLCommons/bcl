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
#include "descriptor/bcl_descriptor_protein_id.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_protein_id.cpp
  //!
  //! @author mendenjl
  //! @date Feb 15, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorProteinId :
    public ExampleInterface
  {
  public:

    ExampleDescriptorProteinId *Clone() const
    {
      return new ExampleDescriptorProteinId( *this);
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
      // setup

      // construct an implementation that will read proteins from the biol directory
      const std::string pdb_fasta_dir( AddExampleInputPathToFilename( e_Biology, ""));

      // create an object to read proteins, as either pdbs or fastas, from the directory
      util::Implementation
      <
        assemble::RetrieveProteinModelWithCache
      > retriever_all( "ProteinDirectory( " + pdb_fasta_dir + ")");

      // retrieve 1ubi
      assemble::ProteinModelWithCache ubiquiton( *retriever_all->Retrieve( "1ubi"));

      // create an iterator for ubiquiton
      descriptor::Iterator< biol::AABase> itr_ubiquiton
      (
        descriptor::Type( 0, false, descriptor::Type::e_Symmetric),
        ubiquiton
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a descriptor that will retrieve the protein id
      descriptor::ProteinId id_retriever;

      // call set object with ubiquiton
      id_retriever.SetObject( ubiquiton);

    /////////////////
    // data access //
    /////////////////

      // test the GetSizeOfFeatures function
      BCL_ExampleCheck( id_retriever.GetSizeOfFeatures(), 5);
      BCL_ExampleCheck( id_retriever.GetType().GetDimension(), 0);

    ///////////////
    // operators //
    ///////////////

      // operator()
      BCL_ExampleCheck
      (
        id_retriever( itr_ubiquiton),
        linal::Vector< char>( id_retriever.GetNormalSizeOfFeatures(), "1ubi ")
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAASSTMPrediction

  const ExampleClass::EnumType ExampleDescriptorProteinId::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorProteinId())
  );

} // namespace bcl
