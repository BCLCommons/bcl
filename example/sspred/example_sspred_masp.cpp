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
#include "sspred/bcl_sspred_masp.h"

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
  //! @example example_sspred_masp.cpp
  //!
  //! @author mendenjl
  //! @date May 15, 2014
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredMASP :
    public ExampleInterface
  {
  public:

    ExampleSspredMASP *Clone() const
    {
      return new ExampleSspredMASP( *this);
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
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1lghA.fasta"));
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      sspred::MASP def_construct;
      const linal::Vector3D def_prediction( 0, 0, 1);
      const linal::Vector3D helix( 1.0, 0.0, 0.0);
      BCL_ExampleCheckWithinTolerance( def_construct.GetThreeStatePrediction(), def_prediction, 0.001);

    /////////////////
    // data access //
    /////////////////

      // test GetFileExtension
      const std::string extension( ".masp");
      BCL_ExampleCheck( def_construct.GetFileExtension(), extension);

    ////////////////
    // operations //
    ////////////////

      // read from file
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1lghA.masp"));
      def_construct.ReadPredictionsForAASequence( read, seq);
      const linal::Vector3D correct_pred( 1.0, 0.0, 0.0);
      const sspred::MASP pred_construct( linal::Vector3D( 1.0, 0.0, 0.0), linal::Vector3D( 0.0));

      // test GetThreeStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetThreeStatePrediction(), helix, 0.001);

      // test GetNineStatePrediction
      BCL_ExampleCheckWithinTolerance( pred_construct.GetNineStatePrediction().GetRow( 0), helix, 0.001);

      // last residue before TM-helix, should be coil
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 1)->GetSSPrediction( sspred::GetMethods().e_MASP)->GetThreeStatePrediction(),
        linal::Vector3D( 0.164, 0.13, 0.706),
        0.001
      );

      // first residue in TM-helix
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 15)->GetSSPrediction( sspred::GetMethods().e_MASP)->GetThreeStatePrediction(),
        correct_pred,
        0.001
      );

      // last residue in TM-helix
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 39)->GetSSPrediction( sspred::GetMethods().e_MASP)->GetThreeStatePrediction(),
        correct_pred,
        0.001
      );

      // first residue after TM-helix, should be coil
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 40)->GetSSPrediction( sspred::GetMethods().e_MASP)->GetThreeStatePrediction(),
        linal::Vector3D( 0.575, 0.093, 0.332),
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
          ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_MASP)->GetOneStateSSPrediction().GetName()
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

  }; //end ExampleSspredMASP

  const ExampleClass::EnumType ExampleSspredMASP::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredMASP())
  );

} // namespace bcl
