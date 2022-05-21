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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "score/bcl_score_aa_neighborhood_exposure_prediction.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list.h"
#include "biol/bcl_biol_aa_base.h"
#include "command/bcl_command_flag_static.h"
#include "fold/bcl_fold_default_flags.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AANeighborhoodExposurePrediction::s_Instance
    (
      util::Enumerated< AANeighborhoodInterface>::AddInstance( new AANeighborhoodExposurePrediction())
    );

    //! @brief return command line flag for using exposure score
    //! @return command line flag for using exposure score
    util::ShPtr< command::FlagInterface> &AANeighborhoodExposurePrediction::GetFlagScoreExposure()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "score_exposure",
          "\tflag to enable use of exposure scores, requires a prediction file given by " +
            fold::DefaultFlags::GetFlagReadSequenceDataPath()->GetName() + " flag"
        )
      );
      // end
      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AANeighborhoodExposurePrediction::AANeighborhoodExposurePrediction()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AANeighborhoodExposurePrediction
    AANeighborhoodExposurePrediction *AANeighborhoodExposurePrediction::Clone() const
    {
      return new AANeighborhoodExposurePrediction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborhoodExposurePrediction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the minimal sequence separation
    //! @return minimal sequence separation for neighbors of that exposure score between amino acids in the same chain
    size_t AANeighborhoodExposurePrediction::GetMinimalSequenceSeparation() const
    {
      return GetAANeighborCount().GetMinimalSequenceSeparation();
    }

    //! @brief access to the distance cutoff
    //! @return distance cutoff above which the neighbor does not have influence on the score anymore
    double AANeighborhoodExposurePrediction::GetDistanceCutoff() const
    {
      return GetAANeighborCount().GetDistanceCutoff();
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AANeighborhoodExposurePrediction::GetScheme() const
    {
      static const std::string s_scheme( "exposure_prediction");
      return s_scheme;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &AANeighborhoodExposurePrediction::GetAlias() const
    {
      static const std::string s_name( "ExposurePrediction");
      return s_name;
    }
  ///////////////
  // operators //
  ///////////////

    //! @brief calculates the score for the passed neighborhood
    //! @param NEIGHBOR_LIST AA neighbor list
    //! @param MEMBRANE membrane object
    //! @return exposure score based on the prediction
    double AANeighborhoodExposurePrediction::operator()
    (
      const assemble::AANeighborList &NEIGHBOR_LIST,
      const util::SiPtr< const biol::Membrane> &MEMBRANE
    ) const
    {
      // get predicted neighbor count
      const double predicted_exposure( NEIGHBOR_LIST.GetCenterAminoAcid()->GetExposurePrediction());

      // if undefined, return score of 0
      if( !util::IsDefined( predicted_exposure))
      {
        return 0.0;
      }

      // return the absolute neighbor count difference
      return math::Absolute( GetAANeighborCount()( NEIGHBOR_LIST) - predicted_exposure);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AANeighborhoodExposurePrediction::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AANeighborhoodExposurePrediction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

    //! @brief write detailed scheme and values to OSTREAM
    //! @param NEIGHBOR_LIST AA neighbor list
    //! @param MEMBRANE membrane object
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &AANeighborhoodExposurePrediction::WriteDetailedSchemeAndValues
    (
      const assemble::AANeighborList &NEIGHBOR_LIST,
      const util::SiPtr< const biol::Membrane> &MEMBRANE,
      std::ostream &OSTREAM
    ) const
    {
      // get the predicted neighbor count
      const double predicted_exposure( NEIGHBOR_LIST.GetCenterAminoAcid()->GetExposurePrediction());

      // get the calculated neighbor count
      const double model_exposure( GetAANeighborCount()( NEIGHBOR_LIST));

      // compare calculated neighbor count to predicted neighbor count
      const double this_score( math::Absolute( model_exposure - predicted_exposure));

      // write information bcl
      OSTREAM << NEIGHBOR_LIST.GetCenterAminoAcid()->GetSeqID() << '\t' << predicted_exposure << '\t'
              << model_exposure << '\t' << this_score;

      // end
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AANeighborhoodExposurePrediction::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.SetClassDescription( "Scores the deviation of exposure in the model from given exposure.");

      return serializer;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gets the static AANeighborCount object
    //! @return the static AANeighborCount object
    const assemble::AANeighborCount &AANeighborhoodExposurePrediction::GetAANeighborCount()
    {
      // create neighbor count object
      static const assemble::AANeighborCount s_neighbor_count;

      return s_neighbor_count;
    }

  } // namespace score
} // namespace bcl
