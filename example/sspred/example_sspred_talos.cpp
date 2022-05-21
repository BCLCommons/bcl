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
#include "sspred/bcl_sspred_talos.h"

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
  //! @example example_sspred_talos.cpp
  //!
  //! @author weinerbe
  //! @date Jun 22, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredTALOS :
    public ExampleInterface
  {
  public:

    ExampleSspredTALOS *Clone() const
    {
      return new ExampleSspredTALOS( *this);
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
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubiA.fasta"));
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      sspred::TALOS def_construct;
      const linal::Vector3D def_prediction( 0, 0, 1);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), def_prediction, 0.001);

      // test constructor from prediction
      const linal::Vector3D helix( 0.1, 0.8, 0.1);
      const linal::Vector< double> helix_math( 3, helix.Begin());
      const sspred::TALOS pred_construct( helix);

    /////////////////
    // data access //
    /////////////////

      // test GetFileExtension
      const std::string extenstion( "SS.tab");
      BCL_ExampleCheck( def_construct.GetFileExtension(), extenstion);

    ////////////////
    // operations //
    ////////////////

      // test GetThreeStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetThreeStatePrediction(), helix, 0.001);

      // test GetNineStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetNineStatePrediction().GetRow( 2), helix_math, 0.001);

      // read from file
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubiASS.tab"));
      def_construct.ReadPredictionsForAASequence( read, seq);
      const linal::Vector3D correct_pred( 0.0, 0.798, 0.202);
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 1)->GetSSPrediction( sspred::GetMethods().e_TALOS)->GetThreeStatePrediction(),
        correct_pred,
        0.001
      );

    //////////////////////
    // input and output //
    //////////////////////

      // input and output
      WriteBCLObject( pred_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), helix, 0.001);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredTALOS

  const ExampleClass::EnumType ExampleSspredTALOS::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredTALOS())
  );

} // namespace bcl
