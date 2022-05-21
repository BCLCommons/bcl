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
#include "biol/bcl_biol_exposure_prediction.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_exposure_prediction.cpp
  //!
  //! @author weinerbe, lib14
  //! @date May 12, 2014
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolExposurePrediction :
    public ExampleInterface
  {
  public:

    ExampleBiolExposurePrediction *Clone() const
    {
      return new ExampleBiolExposurePrediction( *this);
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
      BCL_MessageStd( "read fasta: " + AddExampleInputPathToFilename( e_Biology, "2K73A.fasta"));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "2K73A.fasta"));

      // create a sequence and read in the fasta sequence
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // open the blast profile file for reading
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "2K73A.uniref50.ascii5"));

      // read the blast profile and close the stream
      biol::BlastProfileHandler::ReadProfileForAASequence( read, seq);
      io::File::CloseClearFStream( read);
      util::ShPtr< biol::AASequence> seq_hardcopy( seq.HardCopy());

    ////////////////
    // operations //
    ////////////////

      // calculate exposure for sequence
      biol::ExposurePrediction::Calculate( seq);

      // check protomeric contact number prediction
      biol::ExposurePrediction::Calculate( seq, biol::ExposurePrediction::e_ProtomericContactNumber);
      BCL_ExampleCheckWithinAbsTolerance( seq.GetAA( 45)->GetExposurePrediction(), 1.349, 0.001);

      // check another exposure type
      biol::ExposurePrediction::Calculate( seq, biol::ExposurePrediction::e_OligomericContactNumber);
      BCL_ExampleCheckWithinAbsTolerance( seq.GetAA( 45)->GetExposurePrediction(), 1.752, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // write out the predictions
      io::OFStream write;
      const std::string filename( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "2K73A.exposure"));
      BCL_ExampleMustOpenOutputFile( write, filename);
      biol::ExposurePrediction::WritePredictions( write, seq);
      io::File::CloseClearFStream( write);

      // read in the predictions
      BCL_ExampleMustOpenInputFile( read, filename);
      biol::ExposurePrediction::ReadPredictions( read, *seq_hardcopy);
      io::File::CloseClearFStream( read);
      BCL_ExampleCheckWithinAbsTolerance( seq_hardcopy->GetAA( 45)->GetExposurePrediction(), 1.752, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolExposurePrediction

  const ExampleClass::EnumType ExampleBiolExposurePrediction::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolExposurePrediction())
  );

} // namespace bcl
