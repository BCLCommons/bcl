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
#include "descriptor/bcl_descriptor_aa_ss_tm_prediction.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_retrieve_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_ss_tm_prediction.cpp
  //!
  //! @author mendenjl
  //! @date Feb 15, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAASSTMPrediction :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAASSTMPrediction *Clone() const
    {
      return new ExampleDescriptorAASSTMPrediction( *this);
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
        descriptor::Type( 1, false, descriptor::Type::e_Symmetric),
        ubiquiton
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a descriptor that will retrieve PSIPRED_SS2 information (secondary structure only)
      descriptor::AASSTMPrediction psipred( "", false, true);
      BCL_ExampleAssert
      (
        psipred.TryRead
        (
          util::ObjectDataLabel
          (
            "",
            psipred.GetAlias(),
            storage::Vector< util::ObjectDataLabel>::Create( util::ObjectDataLabel( "method", sspred::GetMethods().e_PSIPRED.GetName()))
          ),
          util::GetLogger()
        ),
        true
      );

      // create a descriptor that will retrieve jufo SS information
      descriptor::AASSTMPrediction jufo( "", false, true);
      jufo.SetMethod( sspred::GetMethods().e_JUFO);

      // create a descriptor that will retrieve jufo SS & TM information
      descriptor::AASSTMPrediction jufo9d( sspred::GetMethods().e_JUFO9D, true, true);
      jufo9d.SetMethod( sspred::GetMethods().e_JUFO);

      // create a descriptor that will retrieve jufo SS & TM information
      descriptor::AASSTMPrediction undefined_method( "", true, true);

      // call set object with ubiquiton
      psipred.SetObject( ubiquiton);
      jufo.SetObject( ubiquiton);
      jufo9d.SetObject( ubiquiton);
      undefined_method.SetObject( ubiquiton);

    /////////////////
    // data access //
    /////////////////

      // test the GetSizeOfFeatures function
      BCL_ExampleCheck( psipred.GetSizeOfFeatures(), 3);
      BCL_ExampleCheck( jufo9d.GetSizeOfFeatures(), 9);

      // get type
      BCL_ExampleCheck( jufo9d.GetType().GetDimension(), 1);

    ///////////////
    // operators //
    ///////////////

      // operator()
      BCL_ExampleCheckWithinTolerance
      (
        jufo( itr_ubiquiton),
        linal::MakeVector< float>( 0.006, 0.391, 0.603),
        1.0e-6
      );

      // check that the undefined method returns a vector of NANs
      BCL_ExampleCheck( undefined_method( itr_ubiquiton).IsDefined(), false);
      BCL_ExampleCheck( undefined_method( itr_ubiquiton).GetSize(), 9);

      // try it for psipred too
      BCL_ExampleCheckWithinTolerance
      (
        psipred( itr_ubiquiton),
        linal::MakeVector< float>( 0.001, 0.014, 0.990).SetToSum( 1.0),
        1.0e-6
      );

      // increment, make sure results are still valid
      ++itr_ubiquiton;
      BCL_ExampleIndirectCheckWithinTolerance
      (
        psipred( itr_ubiquiton),
        linal::MakeVector< float>( 0.000, 0.845, 0.196).SetToSum( 1.0),
        1.0e-6,
        "Incrementing iterator changes the output"
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

  const ExampleClass::EnumType ExampleDescriptorAASSTMPrediction::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAASSTMPrediction())
  );

} // namespace bcl
