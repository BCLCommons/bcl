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
#include "sspred/bcl_sspred_method_handler.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "sspred/bcl_sspred_jufo.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_method_handler.cpp
  //!
  //! @author weinerbe, teixeipl
  //! @date Oct 15, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredMethodHandler :
    public ExampleInterface
  {
  public:

    ExampleSspredMethodHandler *Clone() const
    {
      return new ExampleSspredMethodHandler( *this);
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
      // get a protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test ReadPredictionsForAA
      std::stringstream jufo_aa_pred( "   1 M U   0.603  0.006  0.391\n");
      const linal::Vector3D correct_pred( 0.006, 0.391, 0.603);
      util::ShPtr< biol::AABase> sp_first_aa( protein_model.GetSequences().FirstElement()->GetFirstAA());
      storage::Set< sspred::Method> initial_set( sp_first_aa->GetSSPredictions().GetKeys());
      sspred::MethodHandler::ReadPredictionsForAA
      (
        jufo_aa_pred,
        *sp_first_aa,
        sspred::GetMethods().e_JUFO
      );
      BCL_ExampleCheck
      (
        sp_first_aa->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction(),
        correct_pred
      );

      // open Jufo prediction file
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubiA.jufo"));

      // get aasequence
      util::ShPtr< biol::AASequence> sp_aa_sequence( protein_model.GetChains().FirstElement()->GetSequence());

      // test ReadPredictionsForAASequence
      sspred::MethodHandler::ReadPredictionsForAASequence
      (
        read,
        *sp_aa_sequence,
        sspred::GetMethods().e_JUFO
      );
      BCL_ExampleCheck
      (
        sp_aa_sequence->GetLastAA()->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction(),
        linal::Vector3D( 0.053, 0.068, 0.879)
      );

      // close read stream
      io::File::CloseClearFStream( read);

      // test PossibleFileEntries
      storage::Vector< io::DirectoryEntry> result_vector
      (
        sspred::MethodHandler::PossibleFileEntries
        (
          sspred::GetMethods().e_JUFO,
          'A',
          "1ubi",
          AddExampleInputPathToFilename( e_Biology, "")
        )
      );

      BCL_ExampleCheck( result_vector.GetSize(), 1);

      // test AvailablePredictionFiles
      storage::Set< sspred::Method> temp_set
      (
        sspred::MethodHandler::AvailablePredictionFiles
        (
          storage::Set< sspred::Method>( sspred::GetMethods().Begin(), sspred::GetMethods().End()),
          'A',
          "1ubi",
          AddExampleInputPathToFilename( e_Biology, "")
        )
      );
      BCL_ExampleCheck( temp_set.IsEmpty(), false);
      storage::Set< sspred::Method> full_set( temp_set);
      full_set.InsertElements( initial_set);

      // test ReadPredictionsForAASequence
      sspred::MethodHandler::ReadPredictionsForAASequence
      (
        temp_set,
        *sp_aa_sequence,
        "1ubi",
        AddExampleInputPathToFilename( e_Biology, "")
      );
      BCL_ExampleCheck
      (
        sp_aa_sequence->GetLastAA()->GetSSPredictions().GetKeys(),
        full_set
      );

      // get a new protein model
      assemble::ProteinModel new_protein_model
      (
        pdb::Factory().ProteinModelFromPDBFilename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

      // test ReadPredictionsForProteinModel
      sspred::MethodHandler::ReadPredictionsForProteinModel
      (
        temp_set,
        new_protein_model,
        "1ubi",
        AddExampleInputPathToFilename( e_Biology, "")
      );
      BCL_ExampleCheck
      (
        new_protein_model.GetSequences().FirstElement()->GetLastAA()->GetSSPredictions().GetSize(),
        temp_set.GetSize()
      );

      // Open write stream
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( sspred::MethodHandler(), "1ubi_test_jufo.txt"));

      // test WritePredictionsForAASequence
      sspred::MethodHandler::WritePredictionsForAASequence
      (
        write,
        *sp_aa_sequence,
        sspred::GetMethods().e_JUFO
      );

      // close write stream
      io::File::CloseClearFStream( write);

    //////////////////////
    // helper functions //
    //////////////////////

      // make a jufo prediction
      sspred::JUFO jufo_prediction( linal::Vector3D( 1, 0, 0));

      // test InitializePredictionsForAASequence
      sspred::MethodHandler::InitializePredictionsForAASequence
      (
        sspred::GetMethods().e_JUFO,
        *sp_aa_sequence,
        jufo_prediction
      );
      BCL_ExampleCheck
      (
        sp_aa_sequence->GetLastAA()->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction(),
        jufo_prediction.GetThreeStatePrediction()
      );

      // make a jufo prediction
      sspred::JUFO new_jufo_prediction( linal::Vector3D( 0, 1, 0));

      // test SetPredictionsForSubSequence
      sspred::MethodHandler::SetPredictionsForSubSequence
      (
        sspred::GetMethods().e_JUFO,
        *sp_aa_sequence,
        new_jufo_prediction,
        5,
        5
      );
      BCL_ExampleCheck
      (
        sp_aa_sequence->GetAA( 5)->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction(),
        new_jufo_prediction.GetThreeStatePrediction()
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredMethodHandler

  const ExampleClass::EnumType ExampleSspredMethodHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredMethodHandler())
  );

} // namespace bcl
