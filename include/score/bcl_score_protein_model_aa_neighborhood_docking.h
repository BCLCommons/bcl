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

#ifndef BCL_SCORE_PROTEIN_MODEL_AA_NEIGHBORHOOD_DOCKING_H_
#define BCL_SCORE_PROTEIN_MODEL_AA_NEIGHBORHOOD_DOCKING_H_

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_neighborhood_interface.h"
#include "bcl_score_protein_model.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace score
  {
    /////////////////////////////////////////////
    //! @class ProteinModelAANeighborhoodDocking
    //! @brief Scores neighborhood environment for docking models
    //!
    //! @see @link example_score_protein_model_aa_neighborhood_docking.cpp @endlink
    //! @author lib14
    //! @modified April 20, 2018
    //!
    /////////////////////////////////////////////

    class BCL_API ProteinModelAANeighborhoodDocking :
        public ProteinModel
    {

    private:

      //
      util::Implementation< AANeighborhoodInterface> m_ScoreAANeighborhood;

      //
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      ProteinModelAANeighborhoodDocking();

      //! @param SCORE_AA_NEIGHBORHOOD
      //! @param SCHEME
      ProteinModelAANeighborhoodDocking
      (
        const AANeighborhoodInterface &SCORE_AA_NEIGHBORHOOD,
        const std::string &SCHEME
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new ProteinModelAANeighborhood copied from this one
      ProteinModelAANeighborhoodDocking *Clone() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the sum of exposures of all amino acids for the given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return the sum of exposures of all amino acids for the given ProteinModel
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief reads the exposure predictions from a file
      //! @param ISTREAM stream to read from
      //! @param PROTEIN_MODEL protein model to contain predictions
      static void ReadPredictions( std::istream &ISTREAM, assemble::ProteinModel &PROTEIN_MODEL);

      //! @brief reads the exposure predictions from a file
      //! @param ISTREAM stream to read from
      //! @param SEQUENCE sequence to contain predictions
      static void ReadPredictions
      (
        const storage::Vector< storage::Pair< int, std::string> > &CHAIN_EXPOSURE,
        biol::AASequence &SEQUENCE
      );

    }; // end of class ProteinModelAANeighborhoodDocking
  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_AA_NEIGHBORHOOD_DOCKING_H_
