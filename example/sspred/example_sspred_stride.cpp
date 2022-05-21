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
#include "sspred/bcl_sspred_stride.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_stride.cpp
  //!
  //! @author mendenjl
  //! @date Jun 18, 2013
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredStride :
    public ExampleInterface
  {
  public:

    ExampleSspredStride *Clone() const
    {
      return new ExampleSspredStride( *this);
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
      sspred::Stride strand_construct( biol::GetSSTypes().STRAND);
      sspred::Stride helix_construct( biol::GetSSTypes().e_HelixLeftAlpha);
      const linal::Vector3D def_prediction( 0, 0, 1);
      const linal::Vector3D helix( 1.0, 0.0, 0.0);
      const linal::Vector3D strand( 0.0, 1.0, 0.0);
      BCL_ExampleCheckWithinTolerance( strand_construct.GetThreeStatePrediction(), strand, 0.001);
      BCL_ExampleCheckWithinTolerance( helix_construct.GetThreeStatePrediction(), helix, 0.001);

    /////////////////
    // data access //
    /////////////////

      // test GetFileExtension
      const std::string extension( ".stride");
      BCL_ExampleCheck( strand_construct.GetFileExtension(), extension);

    ////////////////
    // operations //
    ////////////////

      // read from file
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubiA.stride"));
      strand_construct.ReadPredictionsForAASequence( read, seq);

      // middle residue in strand
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 12)->GetSSPrediction( sspred::GetMethods().e_Stride)->GetThreeStatePrediction(),
        strand,
        0.001
      );

      // first residue in helix
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 22)->GetSSPrediction( sspred::GetMethods().e_Stride)->GetThreeStatePrediction(),
        helix,
        0.001
      );

      // last residue in TM-helix
      BCL_ExampleCheckWithinTolerance
      (
        seq.GetAA( 33)->GetSSPrediction( sspred::GetMethods().e_Stride)->GetThreeStatePrediction(),
        helix,
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
          ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_Stride)->GetOneStateSSPrediction().GetName()
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // input and output
      WriteBCLObject( helix_construct);
      sspred::Stride def_construct;
      ReadBCLObject( def_construct);
      BCL_ExampleCheckWithinTolerance
      (
        def_construct.GetThreeStatePrediction(),
        helix_construct.GetThreeStatePrediction(),
        0.001
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredStride

  const ExampleClass::EnumType ExampleSspredStride::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredStride())
  );

} // namespace bcl
