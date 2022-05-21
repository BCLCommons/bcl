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

#ifndef BCL_SCORE_AA_NEIGHBORHOOD_EXPOSURE_PREDICTION_H_
#define BCL_SCORE_AA_NEIGHBORHOOD_EXPOSURE_PREDICTION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "command/bcl_command.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_neighborhood_interface.h"
#include "command/bcl_command_flag_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborhoodExposurePrediction
    //! @brief score the neighborhood based on the exposure prediction of the center AA
    //! @details Compares the predicted neighbor counts with the calculated neighbor count from BCL folded model
    //!
    //! @see @link example_score_aa_neighborhood_exposure_prediction.cpp @endlink
    //! @author weinerbe, lib14
    //! @date March 20, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborhoodExposurePrediction :
      public AANeighborhoodInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief return command line flag for using exposure score
      //! @return command line flag for using exposure score
      static util::ShPtr< command::FlagInterface> &GetFlagScoreExposure();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AANeighborhoodExposurePrediction();

      //! @brief Clone function
      //! @return pointer to new AANeighborhoodExposurePrediction
      AANeighborhoodExposurePrediction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation for neighbors of that exposure score between amino acids in the same chain
      size_t GetMinimalSequenceSeparation() const;

      //! @brief access to the distance cutoff
      //! @return distance cutoff above which the neighbor does not have influence on the score anymore
      double GetDistanceCutoff() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculates the score for the passed neighborhood
      //! @param NEIGHBOR_LIST AA neighbor list
      //! @param MEMBRANE membrane object
      //! @return exposure score based on the prediction
      double operator()
      (
        const assemble::AANeighborList &NEIGHBOR_LIST,
        const util::SiPtr< const biol::Membrane> &MEMBRANE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param NEIGHBOR_LIST AA neighbor list
      //! @param MEMBRANE membrane object
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::AANeighborList &NEIGHBOR_LIST,
        const util::SiPtr< const biol::Membrane> &MEMBRANE,
        std::ostream &OSTREAM
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

//      //! @brief gets the linear regression function used to calculate rSASA from NC
//      //! @return the linear regression function used to calculate rSASA from NC
//      static const math::LinearFunction &GetNCTorSASAFunction();

      //! @brief gets the static AANeighborCount object
      //! @return the static AANeighborCount object
      static const assemble::AANeighborCount &GetAANeighborCount();

    }; // class AANeighborhoodExposurePrediction

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_AA_NEIGHBORHOOD_EXPOSURE_PREDICTION_H_
