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
#include "contact/bcl_contact_prediction_map.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_prediction_map.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactPredictionMap :
    public ExampleInterface
  {
  public:

    ExampleContactPredictionMap *Clone() const
    {
      return new ExampleContactPredictionMap( *this);
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
      // create AAsequence sequence
      io::IFStream read;
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
      BCL_MessageStd( "read pdb: " + pdb_filename);
      util::ShPtr< biol::AASequence> sp_sequence_full
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AACaCb).GetChain( 'A')->GetSequence()
      );
      // read jufo prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1IE9A.jufo"));
      sspred::MethodHandler::ReadPredictionsForAASequence( read, *sp_sequence_full, sspred::GetMethods().e_JUFO);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1IE9A.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *sp_sequence_full);
      io::File::CloseClearFStream( read);

      // create a subsequence
      util::ShPtr< biol::AASequence> sp_sequence( new biol::AASequence( sp_sequence_full->SubSequence( 0, 40)));

      // create a chain with the sequence
      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( sp_sequence));

      // create a protein model
      assemble::ProteinModel protein_model( sp_chain);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the PredictionMap
      BCL_MessageStd( "Calling default constructor")
      contact::PredictionMap prediction_map_empty;

      BCL_MessageStd( "Creating the prediction map for a chain of 1IE9");
      contact::PredictionMap prediction_map_a( *sp_chain);

      BCL_MessageStd( "Creating the prediction map for a model of 1IE9");
      contact::PredictionMap prediction_map_model( protein_model);

      BCL_MessageStd( "Calling clone")
      util::ShPtr< contact::PredictionMap> sp_prediction_map( prediction_map_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_MessageStd
      (
        "This class has the following identifier" + prediction_map_a.GetClassIdentifier()
      );

      // example check
      BCL_Example_Check
      (
        prediction_map_a.GetClassIdentifier() == "bcl::contact::PredictionMap",
        "The class identifier is wrong : " + prediction_map_a.GetClassIdentifier()
      );

      // get the sequences
      BCL_MessageStd
      (
        "Number of sequences stored: " + util::Format()( prediction_map_a.GetSequences().GetSize())
      );

      // example check
      BCL_ExampleCheck( prediction_map_a.GetSequences().GetSize(), 1);

      // iterate over every sequence
      for
      (
        util::ShPtrVector< biol::AASequence>::const_iterator seq_itr( prediction_map_a.GetSequences().Begin()),
          seq_itr_end( prediction_map_a.GetSequences().End());
        seq_itr != seq_itr_end;
        ++seq_itr
      )
      {
        // output the sequence and the chain id
        BCL_MessageStd
        (
          " Sequence with chain ID " + util::Format()( ( *seq_itr)->GetChainID()) + "\n" + ( *seq_itr)->Sequence()
        );
      }

      // check the chain id
      BCL_ExampleAssert( prediction_map_a.GetSequence( 'A').IsDefined(), true);

      // check the chain id
      BCL_ExampleCheck( prediction_map_a.GetSequence( 'A')->GetChainID(), 'A');

      // get a specific sequences
      BCL_MessageStd
      (
        "Sequence of the chain A\n" + prediction_map_a.GetSequence( 'A')->Sequence()
      );

    ////////////////
    // operations //
    ////////////////

      // test the GetVector function for various residue couples
      BCL_MessageStd( "Printing out predictions for different residue couples :");
      storage::Vector< storage::VectorND< 2, size_t> > indices;
      indices.PushBack( storage::VectorND< 2, size_t>( 10, 20));
      indices.PushBack( storage::VectorND< 2, size_t>( 20, 10));
      indices.PushBack( storage::VectorND< 2, size_t>( 20, 30));
      indices.PushBack( storage::VectorND< 2, size_t>( 30, 20));
      indices.PushBack( storage::VectorND< 2, size_t>( 10, 30));
      indices.PushBack( storage::VectorND< 2, size_t>( 30, 10));

      // iterate over index pairs
      for
      (
        storage::Vector< storage::VectorND< 2, size_t> >::const_iterator index_itr( indices.Begin()),
          index_itr_end( indices.End());
        index_itr != index_itr_end; ++index_itr
      )
      {
        // initialize ptr to amino acids
        const util::SiPtr< const biol::AABase> aa_ptr_a( sp_sequence->GetAA( index_itr->First()));
        const util::SiPtr< const biol::AABase> aa_ptr_b( sp_sequence->GetAA( index_itr->Second()));

        // form aa pair
        const storage::VectorND< 2, util::SiPtr< const biol::AAData> > aa_pair
        (
          aa_ptr_a->GetData(),
          aa_ptr_b->GetData()
        );

        // get the predictions
        const storage::Pair< linal::Vector< double>, double> predictions( prediction_map_a.GetPredictions( aa_pair));

        // output the predictions
        BCL_MessageStd
        (
          util::Format().W( 5)( aa_ptr_a->GetSeqID()) + " " +
          aa_ptr_a->GetType()->GetOneLetterCode() + " " +
          util::Format().W( 5)( aa_ptr_b->GetSeqID()) + " " +
          aa_ptr_b->GetType()->GetOneLetterCode() + " " +
          util::Format().W( 5).FFP( 3)( predictions.First()( contact::GetTypes().HELIX_HELIX)) + " " +
          util::Format().W( 5).FFP( 3)( predictions.First()( contact::GetTypes().HELIX_SHEET)) + " " +
          util::Format().W( 5).FFP( 3)( predictions.First()( contact::GetTypes().SHEET_HELIX)) + " " +
          util::Format().W( 5).FFP( 3)( predictions.First()( contact::GetTypes().STRAND_STRAND)) + " " +
          util::Format().W( 5).FFP( 3)( predictions.First()( contact::GetTypes().SHEET_SHEET)) + " " +
          util::Format().W( 5).FFP( 3)( predictions.Second())
        );
      } // end index itr

    //////////////////////
    // input and output //
    //////////////////////

      // output the prediction map to a file
      const std::string out_filename( AddExampleOutputPathToFilename( prediction_map_a, "1IE9A.contact_prediction_map"));
      BCL_MessageStd( "Writing the the map to file " + out_filename);
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, out_filename);

      prediction_map_a.WritePredictionMap( write);
      io::File::CloseClearFStream( write);

      // read the prediction map back from the file
      BCL_MessageStd( "Reading a new contact map from file");
      BCL_ExampleMustOpenInputFile( read, out_filename);

      contact::PredictionMap prediction_map_b( sp_sequence);
      prediction_map_b.ReadPredictionMap( read);

      BCL_MessageStd( "Comparing values between two maps");

      // iterate over every residue
      for
      (
        biol::AASequence::const_iterator aa_itr_a( sp_sequence->Begin()),
          aa_itr_end( sp_sequence->End());
        aa_itr_a != aa_itr_end;
        ++aa_itr_a
      )
      {
        // versus every residue
        for
        (
          biol::AASequence::const_iterator aa_itr_b( sp_sequence->Begin());
          aa_itr_b != aa_itr_end;
          ++aa_itr_b
        )
        {

          // create pointers to amino acids
          const storage::VectorND< 2, util::SiPtr< const biol::AAData> > aa_pair
          (
            ( *aa_itr_a)->GetData(),
            ( *aa_itr_b)->GetData()
          );

          // get the predictions
          const storage::Pair< linal::Vector< double>, double> predictions_a( prediction_map_a.GetPredictions( aa_pair));
          const storage::Pair< linal::Vector< double>, double> predictions_b( prediction_map_b.GetPredictions( aa_pair));

          // if both predictions are defined
          if
          (
            predictions_a.Second() != util::GetUndefined< double>() &&
            predictions_b.Second() == util::GetUndefined< double>()
          )
          {
            // assert that either both are undefined or both are equal
            BCL_Example_Check
            (
              math::EqualWithinTolerance
              (
                predictions_a.Second(),
                predictions_b.Second(),
                0.01
              ),
              "Numbers are not same, check the reading protocol for the error!!!" +
              util::Format()( ( *aa_itr_a)->GetSeqID()) + " vs " + util::Format()( ( *aa_itr_b)->GetSeqID()) +
              " prediction created \n" +
              util::Format()( predictions_a) +
              " prediction read    \n" +
              util::Format()( predictions_b)
            );
          };
        } // end aa_itr_b
      } // end aa_itr_a

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleContactPredictionMap

  const ExampleClass::EnumType ExampleContactPredictionMap::s_Instance
  (
    GetExamples().AddEnum( ExampleContactPredictionMap())
  );

} // namespace bcl

