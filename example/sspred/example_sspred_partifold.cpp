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
#include "sspred/bcl_sspred_partifold.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_partifold.cpp
  //!
  //! @author koehlej, weinerbe
  //! @date Sep 13, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredPARTIFOLD :
    public ExampleInterface
  {
  public:

    ExampleSspredPARTIFOLD *Clone() const
    {
      return new ExampleSspredPARTIFOLD( *this);
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
      // read in fasta sequence
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1b43A.fasta"));
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      sspred::PARTIFOLD def_construct;
      const linal::Vector3D def_prediction( 0.0, 0.0, 1.0);
      const linal::Vector3D strand( 0.0, 1.0, 0.0);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), def_prediction, 0.001);

    /////////////////
    // data access //
    /////////////////

      // test GetFileExtension
      const std::string extension( ".mfe");
      BCL_ExampleCheck( def_construct.GetFileExtension(), extension);

    ////////////////
    // operations //
    ////////////////

      // read from file
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1b43A.mfe"));
      def_construct.ReadPredictionsForAASequence( read, seq);
      const linal::Vector3D correct_pred( 0.0, 1.0, 0.0);
      const sspred::PARTIFOLD pred_construct( strand);

      // test GetThreeStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetThreeStatePrediction(), strand, 0.001);

      // test GetNineStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetNineStatePrediction().GetRow( 0), strand, 0.001);

      // last residue before TM-strand, should be coil
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 105)->GetSSPrediction( sspred::GetMethods().e_PARTIFOLD)->GetThreeStatePrediction(),
        def_prediction,
        0.001
      );

      // first residue in TM-strand
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 106)->GetSSPrediction( sspred::GetMethods().e_PARTIFOLD)->GetThreeStatePrediction(),
        correct_pred,
        0.001
      );

      // last residue in TM-strand
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 116)->GetSSPrediction( sspred::GetMethods().e_PARTIFOLD)->GetThreeStatePrediction(),
        correct_pred,
        0.001
      );

      // first residue after TM-strand, should be coil
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 117)->GetSSPrediction( sspred::GetMethods().e_PARTIFOLD)->GetThreeStatePrediction(),
        def_prediction,
        0.001
      );

      // iterate over seq
      for
      (
        util::ShPtrVector< biol::AABase>::const_iterator aa_itr( seq.Begin()), aa_itr_end( seq.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        BCL_MessageStd
        (
          ( *aa_itr)->GetIdentification() + ": " +
          ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_PARTIFOLD)->GetOneStateSSPrediction().GetName()
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // input and output
      WriteBCLObject( pred_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), strand, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredPARTIFOLD

  const ExampleClass::EnumType ExampleSspredPARTIFOLD::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredPARTIFOLD())
  );

} // namespace bcl
