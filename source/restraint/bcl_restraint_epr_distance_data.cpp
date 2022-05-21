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
#include "restraint/bcl_restraint_epr_distance_data.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "score/bcl_score_restraint_atom_attraction.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "score/bcl_score_restraint_distance_spin_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score EPRDistanceData::e_ScoreEPRDistanceRestraint( fold::GetScores().e_Undefined);
    fold::Score EPRDistanceData::e_ScoreEPRDistanceUpperPenalty( fold::GetScores().e_Undefined);
    fold::Score EPRDistanceData::e_ScoreEPRDistanceLowerPenalty( fold::GetScores().e_Undefined);

    const util::SiPtr< const score::RestraintAtomDistance> EPRDistanceData::s_SpinLabelScore
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintDistanceSpinLabel(),
          1.0,
          score::RestraintDistanceSpinLabel::GetDefaultScheme()
        )
      )
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EPRDistanceData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new EPRDistanceData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    EPRDistanceData::EPRDistanceData() :
      m_Handler( GetDefaultHandler()),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new EPRDistanceData
    EPRDistanceData *EPRDistanceData::Clone() const
    {
      return new EPRDistanceData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &EPRDistanceData::GetAlias() const
    {
      static const std::string s_name( "DistanceEPR");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EPRDistanceData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default handler
    //! @return default handler for this class
    const HandlerAtomDistanceAssigned &EPRDistanceData::GetDefaultHandler()
    {
      static const HandlerAtomDistanceAssigned s_handler( EPRDistanceData::GetDefaultExtension());
      return s_handler;
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &EPRDistanceData::GetDefaultExtension()
    {
      static const std::string s_extension( ".epr_cst_bcl");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &EPRDistanceData::GetAtomDistanceRestraints() const
    {
      return *m_Restraints;
    }

    //! @brief gives the scoring object that is used to score this type of restraint
    //! @return the scoring object that is used to score this type of restraint
    const util::ShPtrVector< score::ProteinModel> &EPRDistanceData::GetScores() const
    {
      // initialize scores
      static util::ShPtrVector< score::ProteinModel> scores;

      // if the vector has not been filled
      if( scores.IsEmpty())
      {
        // initialize scores then add them to vector
        EPRDistanceData epr_data;
        epr_data.InitializeScores();
        scores.PushBack( *e_ScoreEPRDistanceRestraint);
        scores.PushBack( *e_ScoreEPRDistanceLowerPenalty);
        scores.PushBack( *e_ScoreEPRDistanceUpperPenalty);
      }

      // end
      return scores;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void EPRDistanceData::InitializeScores()
    {
      if( !e_ScoreEPRDistanceRestraint.IsDefined())
      {
        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler.ReadRestraintsFromFile());
        }
        util::ShPtr< score::RestraintAtomDistance> score( s_SpinLabelScore->Clone());
        score->SetRestraints( m_Restraints);
        e_ScoreEPRDistanceRestraint = fold::GetScores().AddScore( score);

        // create Histogram "histogram" which will be used to hold the histogram of SL-CB distances
        math::Histogram histogram;

        // create IFStream "read"
        io::IFStream read;

        // open "read" and bind it to the histogram file containing SL-CB distances
        io::File::MustOpenIFStream
        (
          read,
          score::Score::AddHistogramPath( score::RestraintDistanceSpinLabel::GetDefaultHistogramFilename())
        );

        // read in from "read" into "histogram"
        read >> histogram;

        // close and clear read stream
        io::File::CloseClearFStream( read);

        // score for the transition on the left (negative x-axis) side of the kb potential
        util::ShPtr< score::RestraintAtomDistance> lower_penalty
        (
          new score::RestraintAtomDistance
          (
            score::RestraintAtomAttraction
            (
              score::RestraintAtomAttraction::GetDefaultDepthRange(),
              score::RestraintAtomAttraction::GetDefaultLeftEndWell( histogram),
              score::RestraintAtomAttraction::GetDefaultTransitionWidth(),
              true
            ),
            1.0,
            "epr_lower_penalty",
            m_Restraints
          )
        );
        util::Enumerated< score::RestraintAtomDistance>::AddInstance( lower_penalty->Clone());
        e_ScoreEPRDistanceLowerPenalty = fold::GetScores().AddScore( lower_penalty);

        // score for the transition on the right (right x-axis) side of the kb potential
        util::ShPtr< score::RestraintAtomDistance> upper_penalty
        (
          new score::RestraintAtomDistance
          (
            score::RestraintAtomAttraction
            (
              score::RestraintAtomAttraction::GetDefaultDepthRange(),
              score::RestraintAtomAttraction::GetDefaultRightEndWell( histogram),
              score::RestraintAtomAttraction::GetDefaultTransitionWidth(),
              false
            ),
            1.0,
            "epr_upper_penalty",
            m_Restraints
          )
        );
        util::Enumerated< score::RestraintAtomDistance>::AddInstance( upper_penalty->Clone());
        e_ScoreEPRDistanceUpperPenalty = fold::GetScores().AddScore( upper_penalty);
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void EPRDistanceData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreEPRDistanceRestraint, 5);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreEPRDistanceUpperPenalty, 5);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreEPRDistanceLowerPenalty, 5);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void EPRDistanceData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRDistanceData::GetSerializer() const
    {
      io::Serializer serial( m_Handler.GetSerializer());
      serial.SetClassDescription( "Allows use of electron paramagnetic resonance-based accessibility data measurements");
      return serial;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EPRDistanceData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Handler, ISTREAM);
      io::Serialize::Read( m_Restraints, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &EPRDistanceData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
