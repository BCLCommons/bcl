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
#include "descriptor/bcl_descriptor_aa_pair_probability.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_pair_probability.cpp
  //!
  //! @author mendenjl
  //! @date Feb 15, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAAPairProbability :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAAPairProbability *Clone() const
    {
      return new ExampleDescriptorAAPairProbability( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////
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
        descriptor::Type( 2, false, descriptor::Type::e_Symmetric),
        ubiquiton
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // retrieve an implementation using the pair probability
      util::Implementation< descriptor::Base< biol::AABase, float> > pam100( "AAPairProbability(method=pam100)");
      pam100->SetObject( ubiquiton);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( pam100->GetSizeOfFeatures(), 1);
      BCL_ExampleCheck( pam100->GetType().GetDimension(), 2);
      BCL_ExampleCheck( pam100->GetType().ConsiderRepeatedObjects(), false);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check the descriptions produced for amino acids 1-2, 1-3, 1-4
      BCL_ExampleCheckWithinAbsTolerance( pam100->operator ()( itr_ubiquiton)( 0), -0.2, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( pam100->operator ()( ++itr_ubiquiton)( 0), 0.1, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( pam100->operator ()( ++itr_ubiquiton)( 0), -0.1, 0.001);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

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

  }; //end ExampleDescriptorAAPairProbability

  const ExampleClass::EnumType ExampleDescriptorAAPairProbability::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAAPairProbability())
  );

} // namespace bcl
