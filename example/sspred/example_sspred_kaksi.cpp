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
#include "sspred/bcl_sspred_kaksi.h"

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
  //! @example example_sspred_kaksi.cpp
  //!
  //! @author mendenjl
  //! @date Jul 15, 2014
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredKaksi :
    public ExampleInterface
  {
  public:

    ExampleSspredKaksi *Clone() const
    {
      return new ExampleSspredKaksi( *this);
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
      sspred::Kaksi def_construct;
      const linal::Vector3D def_prediction( 0, 0, 1);
      const linal::Vector3D helix( 1.0, 0.0, 0.0);
      const linal::Vector3D strand( 0.0, 1.0, 0.0);
      const linal::Vector3D coil( 0.0, 0.0, 1.0);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), def_prediction, 0.001);

    /////////////////
    // data access //
    /////////////////

      // test GetFileExtension
      const std::string extension( ".kaksi");
      BCL_ExampleCheck( def_construct.GetFileExtension(), extension);

    ////////////////
    // operations //
    ////////////////

      // read from file
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1b43.kaksi"));
      def_construct.ReadPredictionsForAASequence( read, seq);
      const linal::Vector3D correct_pred( 1.0, 0.0, 0.0);
      const sspred::Kaksi pred_construct( helix);

      // test GetThreeStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetThreeStatePrediction(), helix, 0.001);

      // test GetNineStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetNineStatePrediction().GetRow( 2), helix, 0.001);

      // A residue in the long helix
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 62)->GetSSPrediction( sspred::GetMethods().e_Kaksi)->GetThreeStatePrediction(),
        helix,
        0.001
      );

      // Another residue in the long helix
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 63)->GetSSPrediction( sspred::GetMethods().e_Kaksi)->GetThreeStatePrediction(),
        helix,
        0.001
      );

      // Last residue in strand
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 79)->GetSSPrediction( sspred::GetMethods().e_Kaksi)->GetThreeStatePrediction(),
        strand,
        0.001
      );

      // Post strand, in coil
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 84)->GetSSPrediction( sspred::GetMethods().e_Kaksi)->GetThreeStatePrediction(),
        coil,
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
          ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_Kaksi)->GetOneStateSSPrediction().GetName()
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // input and output
      WriteBCLObject( pred_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), helix, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredKaksi

  const ExampleClass::EnumType ExampleSspredKaksi::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredKaksi())
  );

} // namespace bcl
