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
#include "score/bcl_score_protein_model_inverted.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "score/bcl_score_aa_neighborhood_exposure.h"
#include "score/bcl_score_protein_model_aa_neighborhood.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_model_inverted.cpp
  //!
  //! @author karakam
  //! @date Oct 10, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreProteinModelInverted :
    public ExampleInterface
  {
  public:

    ExampleScoreProteinModelInverted *Clone() const
    {
      return new ExampleScoreProteinModelInverted( *this);
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
      // construct min sse sizes map
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 4;
      min_sse_sizes[ biol::GetSSTypes().COIL] = 999;

      // get 1UBI model
      const std::string filename_a( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel native_model
      (
        Proteins::GetModel( filename_a, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // get 1UBI model
      const std::string filename_b( AddExampleInputPathToFilename( e_Biology, "1ubi_inverted.pdb"));
      assemble::ProteinModel native_inverted_model
      (
        Proteins::GetModel( filename_b, biol::GetAAClasses().e_AABackBone)
      );

      // construct score
      util::ShPtr< score::ProteinModel> sp_score
       (
         new score::ProteinModelAANeighborhood
         (
           util::CloneToShPtr( score::AANeighborhoodExposure( assemble::AANeighborCount()))
         )
       );

      // initialize inverter
      util::ShPtr< assemble::ProteinModelInverter> sp_inverter( new assemble::ProteinModelInverter());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      score::ProteinModelInverted score_default;

      //! construct from an inverter and scoring function
      const std::string scheme( "neighbor_inverted");
      score::ProteinModelInverted score_inverted( sp_score, sp_inverter, scheme);
      BCL_ExampleCheck( score_inverted.GetScheme(), scheme);

    ///////////////
    // operators //
    ///////////////

      // test operator
      BCL_MessageStd( "Testing operator()");

      // store the expected scores
      const double expected_value( sp_score->operator ()( native_inverted_model));
      BCL_MessageStd( "score of the native inverted model: " + util::Format()( expected_value));
      const double calculated_value( score_inverted( native_model));
      BCL_MessageStd( "calculated score: " + util::Format()( calculated_value));
      BCL_ExampleCheck( expected_value, calculated_value);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreProteinModelInverted

  const ExampleClass::EnumType ExampleScoreProteinModelInverted::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelInverted())
  );

} // namespace bcl
