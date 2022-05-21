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
#include "score/bcl_score_aa_assignment_property.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_assignment_property.cpp
  //!
  //! @author heinzes1, woetzen
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAAssignmentProperty :
    public ExampleInterface
  {
  public:

    ExampleScoreAAAssignmentProperty *Clone() const
    { return new ExampleScoreAAAssignmentProperty( *this);}

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

    /////////////////
    // preparation //
    /////////////////

      io::IFStream read;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      biol::AASequence seq2( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      BCL_MessageStd
      (
        "aatype of first aa of 1fms_ and first aa of 1f5mA " +
        seq1.GetFirstAA()->GetType()->GetThreeLetterCode() + " " + seq2.GetFirstAA()->GetType()->GetThreeLetterCode()
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::AAAssignmentProperty default_constructor;

      storage::Map< biol::AATypeData::PropertyTypeEnum, util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> > > assignment_score_functions;
      for( size_t type( biol::AATypeData::e_NaturalPrevalence); type < biol::AATypeData::s_NumberPropertyTypes; ++type)
      {
        assignment_score_functions[ biol::AATypeData::PropertyType( type)] = util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::PropertyType( type)));
      }

    ///////////////
    // operators //
    ///////////////

      // calculated scores
      storage::Map< biol::AATypeData::PropertyTypeEnum, double> calculated_scores;

      // iterate through constructed scores
      for
      (
        storage::Map< biol::AATypeData::PropertyTypeEnum, util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> > >::const_iterator
          score_itr( assignment_score_functions.Begin()), score_itr_end( assignment_score_functions.End());
        score_itr != score_itr_end;
        ++score_itr
      )
      {
        calculated_scores[ score_itr->first] = score_itr->second->operator ()( *seq1.GetFirstAA(), *seq2.GetFirstAA());
      }

      storage::Map< biol::AATypeData::PropertyTypeEnum, double> expected_scores;

      expected_scores[ biol::AATypeData::e_NaturalPrevalence                     ] = -0.015;
      expected_scores[ biol::AATypeData::e_StericalParameter                     ] = -1.31;
      expected_scores[ biol::AATypeData::e_Polarizability                        ] = -0.06;
      expected_scores[ biol::AATypeData::e_Volume                                ] = -1.6;
      expected_scores[ biol::AATypeData::e_Hydrophobicity                        ] = -0.04;
      expected_scores[ biol::AATypeData::e_IsoelectricPoint                      ] = -0.37;
      expected_scores[ biol::AATypeData::e_HelixProbability                      ] = -0.07;
      expected_scores[ biol::AATypeData::e_StrandProbability                     ] = -0.13;
      expected_scores[ biol::AATypeData::e_TransferFreeEnergyWhimleyWhite        ] = -0.69;
      expected_scores[ biol::AATypeData::e_TransferFreeEnergyEngelmanSeitzGoldman] = -0.4;
      expected_scores[ biol::AATypeData::e_TransferFreeEnergyKyteDoolittle       ] = -0.4;
      expected_scores[ biol::AATypeData::e_TransferFreeEnergyEisenberg           ] = -0.66;
      expected_scores[ biol::AATypeData::e_FreeEnergyCore                        ] = -0.135;
      expected_scores[ biol::AATypeData::e_FreeEnergyTransition                  ] = -0.004;
      expected_scores[ biol::AATypeData::e_FreeEnergySolution                    ] = -0.146;
      expected_scores[ biol::AATypeData::e_SASA                                  ] = -37.884;

      // iterate through expected scores
      for
      (
        storage::Map< biol::AATypeData::PropertyTypeEnum, double>::const_iterator
          itr( expected_scores.Begin()), itr_end( expected_scores.End());
        itr != itr_end;
        ++itr
      )
      {
        const double calculated_score( calculated_scores.Find( itr->first)->second);

        BCL_MessageStd
        (
          "assignment score " + util::Format()( itr->first) + " is " + util::Format()( calculated_score)
        );

        BCL_ExampleIndirectCheck
        (
          math::EqualWithinTolerance( calculated_score, itr->second), true,
          " Assignment score " + util::Format()( itr->first) + " calculated: " + util::Format()( calculated_score) +
          " expected: " + util::Format()( itr->second)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( default_constructor);
      // create instance of class "AAAssignmentProperty" and read from file
      score::AAAssignmentProperty read_score;
      ReadBCLObject( read_score);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAAssignmentProperty

  const ExampleClass::EnumType ExampleScoreAAAssignmentProperty::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAAssignmentProperty())
  );

} // namespace bcl
