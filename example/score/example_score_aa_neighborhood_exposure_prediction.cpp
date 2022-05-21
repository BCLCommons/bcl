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
#include "score/bcl_score_aa_neighborhood_exposure_prediction.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "biol/bcl_biol_exposure_prediction.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_neighborhood_exposure_prediction.cpp
  //!
  //! @author weinerbe, lib14
  //! @date May 12, 2014
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAANeighborhoodExposurePrediction :
    public ExampleInterface
  {
  public:

    ExampleScoreAANeighborhoodExposurePrediction *Clone() const
    {
      return new ExampleScoreAANeighborhoodExposurePrediction( *this);
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
      // read in protein model
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"))
      );
      util::ShPtr< biol::AASequence> sp_sequence( protein_model.GetChain( 'A')->GetSequence());

      // read in blast profile
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "2K73A.uniref50.ascii5"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *sp_sequence);
      io::File::CloseClearFStream( read);

      biol::ExposurePrediction exposure_prediction;
      exposure_prediction.Calculate( *sp_sequence);

      // get buried residue
      assemble::LocatorAA locator( 'A', 46);
      util::SiPtr< const biol::AABase> sp_aa( locator.Locate( protein_model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct membrane
      const util::SiPtr< const biol::Membrane> sp_membrane;

      // construct neighbor count
      assemble::AANeighborCount neighbor_count;

      // construct exposure score
      const score::AANeighborhoodExposurePrediction exposure_score;

    ///////////////
    // operators //
    ///////////////

      // test () operator
      const double expected_score( 4.392);
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        exposure_score
        (
          assemble::AANeighborList
          (
            *sp_aa,
            protein_model.GetAminoAcids(),
            neighbor_count.GetDistanceCutoff(),
            neighbor_count.GetMinimalSequenceSeparation(),
            false
          ),
          sp_membrane
        ),
        expected_score,
        0.01,
        "() operator"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAANeighborhoodExposurePrediction

  const ExampleClass::EnumType ExampleScoreAANeighborhoodExposurePrediction::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAANeighborhoodExposurePrediction())
  );

} // namespace bcl
