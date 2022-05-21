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
#include "sspred/bcl_sspred_jufo.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_jufo.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredJUFO :
    public ExampleInterface
  {
  public:

    ExampleSspredJUFO *Clone() const
    {
      return new ExampleSspredJUFO( *this);
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

      // attach read to a fasta file
      io::IFStream read;
      BCL_MessageStd( "read fasta: " + AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));

      // create a sequence and read in the fasta sequence
      biol::AASequence seq_a( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // open the blast profile file for reading
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.ascii"));

      // read the blast profile and close the stream
      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq_a);
      io::File::CloseClearFStream( read);

      // make a hard copy of the sequence
      util::ShPtr< biol::AASequence> sp_seq_b( seq_a.HardCopy());
      biol::AASequence seq_b( *sp_seq_b);

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      sspred::JUFO jufo_default;
      BCL_ExampleCheck( sspred::JUFO().GetThreeStatePrediction(), sspred::MethodInterface::GetDefaultPredictionVector());

      // test default constructor
      linal::Vector3D jufo_vector( 0.9, 0.1, 0.0);
      sspred::JUFO jufo_from_vector( jufo_vector);
      BCL_ExampleCheck( sspred::JUFO( jufo_vector).GetThreeStatePrediction(), jufo_vector);

    /////////////////
    // data access //
    /////////////////

      // check file extension
      BCL_ExampleCheck( sspred::JUFO().GetFileExtension(), ".jufo");

    ////////////////
    // operations //
    ////////////////

      // test three state prediction
      BCL_MessageStd( "Testing GetNineStatePrediction()");
      linal::Matrix< double> jufo_matrix
      (
        sspred::MethodInterface::ConvertThreeStateToNineState( jufo_vector, biol::GetEnvironmentTypes().e_Solution)
      );
      BCL_ExampleCheck( jufo_from_vector.GetNineStatePrediction(), jufo_matrix);

      // test ReadPredictionsForAASequence function
      BCL_MessageStd( "Testing ReadPredictionsForAASequence function");

      // attach read to jufo file
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.jufo"));

      // read the jufo for that sequence
      sspred::JUFO().ReadPredictionsForAASequence( read, seq_b);

      // run jufo for this sequence
      BCL_MessageStd( "Testing CalculateJUFO function");
      sspred::JUFO::Calculate( seq_a);
      BCL_MessageStd( "CalculateJUFO completed");

      // create iterators to the beginning and end of sequences
      biol::AASequence::const_iterator seq_a_itr( seq_a.Begin());
      biol::AASequence::const_iterator seq_b_itr( seq_b.Begin());
      const biol::AASequence::const_iterator seq_a_itr_end( seq_a.End());
      const biol::AASequence::const_iterator seq_b_itr_end( seq_b.End());

      // iterate over both sequences
      for
      (
        ;
        seq_a_itr != seq_a_itr_end && seq_b_itr != seq_b_itr_end;
        ++seq_a_itr, ++seq_b_itr
      )
      {

        // TODO: Once the new version of JUFO is integrated into BCL, then we would need to update all the jufo
        // files not only for 1fms, but also other proteins, and also update the ascii files with the new db and
        // then we'll be able to compare jufo values from file and calculated on the go. If there are questions, please
        // contact Mert
//        BCL_Example_Check
//        (
//          ExampleClass::ExampleResult::e_Trivial,
//          math::EqualWithinTolerance
//          (
//            ( *seq_a_itr)->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction(),
//            ( *seq_b_itr)->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()
//          ),
//          "The vectors differ significantly from each other\n" +
//            util::Format()( ( *seq_a_itr)->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction()) +
//            "\nvs\n" +
//            util::Format()( ( *seq_b_itr)->GetSSPrediction( sspred::GetMethods().e_JUFO)->GetThreeStatePrediction())
//        );
      }

      // output the calculated JUFO for this sequence to a file
      BCL_MessageStd( "Writing out calculated JUFO to 1fms_created.jufo file");
      // create the output file for writing calculate jufo values for this sequence
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleInputPathToFilename( e_Biology, "1fms_created.jufo"));

      // use the sspred to output the jufo
      sspred::MethodHandler::WritePredictionsForAASequence( write, seq_a, sspred::GetMethods().e_JUFO);
      io::File::CloseClearFStream( write);

      // alert user
      BCL_MessageStd( "jufo file has been created 1fms_created.jufo");

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      WriteBCLObject( jufo_from_vector);
      // read the object back in
      sspred::JUFO jufo_read;
      ReadBCLObject( jufo_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_ExampleIndirectCheck
      (
        jufo_from_vector.GetThreeStatePrediction(),
        jufo_read.GetThreeStatePrediction(),
        "I/O"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredJUFO

  const ExampleClass::EnumType ExampleSspredJUFO::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredJUFO())
  );

} // namespace bcl
