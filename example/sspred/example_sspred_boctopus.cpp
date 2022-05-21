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
#include "sspred/bcl_sspred_boctopus.h"

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
  //! @example example_sspred_boctopus.cpp
  //!
  //! @author mendenjl
  //! @date Aug 25, 2014
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredBOCTOPUS :
    public ExampleInterface
  {
  public:

    ExampleSspredBOCTOPUS *Clone() const
    {
      return new ExampleSspredBOCTOPUS( *this);
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
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1a0tP.fasta"));
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      sspred::BOCTOPUS def_construct;
      const linal::Vector3D def_prediction( 0, 0, 1);
      const linal::Vector3D strand( 0.0, 1.0, 0.0);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), def_prediction, 0.001);

    /////////////////
    // data access //
    /////////////////

      // test GetFileExtension
      const std::string extension( ".bocto_topo");
      BCL_ExampleCheck( def_construct.GetFileExtension(), extension);

    ////////////////
    // operations //
    ////////////////

      // read from file
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1a0tP.bocto_topo"));
      def_construct.ReadPredictionsForAASequence( read, seq);
      const linal::Vector3D correct_pred( 0.0, 1.0, 0.0);
      const sspred::BOCTOPUS pred_construct( 'M');

      // test GetThreeStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetThreeStatePrediction(), strand, 0.001);

      // test GetNineStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetNineStatePrediction().GetRow( 0), strand, 0.001);

      // last residue before predicted TM-strand, should be coil
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 3)->GetSSPrediction( sspred::GetMethods().e_BOCTOPUS)->GetThreeStatePrediction(),
        def_prediction,
        0.001
      );

      // first residue in predicted TM-strand
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 4)->GetSSPrediction( sspred::GetMethods().e_BOCTOPUS)->GetThreeStatePrediction(),
        correct_pred,
        0.001
      );

      // last residue in predicted TM-strand
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 14)->GetSSPrediction( sspred::GetMethods().e_BOCTOPUS)->GetThreeStatePrediction(),
        correct_pred,
        0.001
      );

      // first residue after predicted TM-strand, should be coil
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 15)->GetSSPrediction( sspred::GetMethods().e_BOCTOPUS)->GetThreeStatePrediction(),
        def_prediction,
        0.001
      );

    //////////////////////
    // input and output //
    //////////////////////

      // input and output
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( pred_construct, sspred::BOCTOPUS()), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredBOCTOPUS

  const ExampleClass::EnumType ExampleSspredBOCTOPUS::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredBOCTOPUS())
  );

} // namespace bcl
