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

#ifndef BCL_FOLD_DEFAULT_SCORES_H_
#define BCL_FOLD_DEFAULT_SCORES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "sspred/bcl_sspred.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_scores.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DefaultScores
    //! @brief scores used in default fold protocol
    //! @details class for initializing scores used in default fold protocol as well as forming the ScoreWeightSet with
    //! the weights
    //!
    //! @see @link example_fold_default_scores.cpp @endlink
    //! @author karakam
    //! @date Nov 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DefaultScores :
      public util::ObjectInterface
    {

    public:

    //////////
    // data //
    //////////

      Score e_ScoreAAPairClash;                //!< amino acid pair clash score based on first side chain atom distance
      Score e_ScoreAAPairDistance;             //!< amino acid pair distance potential based on first side chain atom distance
      Score e_ScoreAAPairHiResClash;           //!< amino acid hi-resolution clashing term, considers Ca-Cb vector and angles
      Score e_ScoreAAPairSCInteraction;        //!< amino acid side-chain interaction score
      Score e_ScoreAANeighborCount;            //!< amino acid neighbor count score based on counting neighbors within a radius
      Score e_ScoreAANeighborCountEntropy;     //!< amino acid neighbor count score for amino acids that are not in the model
      Score e_ScoreLoop;                       //!< loop score, evaluating likelihood of euclidean vs. sequence distance
      Score e_ScoreLoopAngle;                  //!< loop angle score, evaluating likelihood of angle between SSEs
      Score e_ScoreLoopClosure;                //!< loop closure score penalizing, if euclidean distance cannot be bridged by residues
      Score e_ScoreLoopClosureGradient;        //!< loop closure score penalizing, if euclidean distance cannot be bridged by residues, uses a wider transition region to max penalty
      Score e_ScorePhiPsi;                     //!< phi psi score, depending on the type of the SSE
      Score e_ScoreRadiusOfGyration;           //!< radius of gyration evaluating model compactness
      Score e_ScoreContactOrder;               //!< contact order evaluating the complexity of model
      Score e_ScoreSSEFragmentPairPacking;     //!< SSE fragment pairing evaluating the relative geometry of two SSE fragments
      Score e_ScoreSSELoopClash;               //!< Naive linear loop clashes with other SSEs; set to 0 for now
      Score e_ScoreStrandFragmentPairing;      //!< SSE strand fragment pairing evaluating proximity of strand fragments facing each other with the hydrogen bond donors and acceptors
      Score e_ScoreSymmetry;                   //!< score symmetry
      Score e_ScoreSSECompleteness;            //!< # of sses still missing from the model
      Score e_ScoreSSETripletChirality;        //!< SSE triplets in sequence space propensity to be in contact
      Score e_ScoreSSEDerivedTripletChirality; //!< SSE triplets in sequence space propensity to be in contact
      Score e_ScoreSSEContactType;             //!< SSE pair contact type
      Score e_ScoreSSEAdjacentContact;         //!< Scores propensity of adjacent SSEs to be in contact
      Score e_ScoreSSEOrientation;             //!< Score orientations of SSEs that are in contact, by # of SSEs in between
      Score e_ScoreSSEInteractionWeight;       //!< Scores SSE interaction weight

      storage::Map< sspred::Method, storage::VectorND< 2, Score> > m_SSPredScores; //!< score secondary structure model assignment relative to given predictions

      util::ShPtr< assemble::ProteinModelInverter> m_ProteinInverter; //!< protein model inverter

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      DefaultScores();

      //! @brief Clone function
      //! @return pointer to new DefaultScores
      DefaultScores *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief access to static instance of this class
      static DefaultScores &GetInstance();

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the protein inverter
      //! @return the protein inverter
      const util::ShPtr< assemble::ProteinModelInverter> &GetProteinInverter() const;

      //! @brief get the sspred scores
      //! the sspred scores
      const storage::Map< sspred::Method, storage::VectorND< 2, Score> > &GetSSPredScores() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initialize the scores and add them to Scores enumerator
      void InitializeScores();

      //! @brief modify the score weight set
      //! @param SCORE_WEIGHT_SET Score weight set
      void ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const;

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

    }; // class DefaultScores

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_DEFAULT_SCORES_H_
