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
#include "restraint/bcl_restraint_epr_accessibility_data.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_scores.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_accessibility.h"
#include "score/bcl_score_accessibility_hydrophobic_moment_magnitude.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentNiEDDA( fold::GetScores().e_Undefined);
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentOxygen( fold::GetScores().e_Undefined);
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentMagnitudeNiEDDA( fold::GetScores().e_Undefined);
    fold::Score EPRAccessibilityData::e_ScoreExposureMomentMagnitudeOxygen( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EPRAccessibilityData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new EPRAccessibilityData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    EPRAccessibilityData::EPRAccessibilityData
    (
      const HandlerAccessibilityAA &HANDLER
    ) :
      m_Handler( HANDLER),
      m_HydrophobicMomentWindowSizes( size_t( 5), size_t( 5), size_t( 5)),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new EPRAccessibilityData
    EPRAccessibilityData *EPRAccessibilityData::Clone() const
    {
      return new EPRAccessibilityData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &EPRAccessibilityData::GetAlias() const
    {
      static const std::string s_name( "AccessibilityEPR");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EPRAccessibilityData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &EPRAccessibilityData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".access_bcl");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const HandlerAccessibilityAA &EPRAccessibilityData::GetDefaultHandler()
    {
      static const HandlerAccessibilityAA s_handler
      (
        util::ShPtr< assemble::AAExposureInterface>( new assemble::AANeighborVector())
      );
      return s_handler;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void EPRAccessibilityData::InitializeScores()
    {
      // read restraints from the restraint file if necessary
      if( !m_Restraints.IsDefined() && m_Handler.Exists())
      {
        m_Restraints = util::CloneToShPtr( m_Handler.ReadRestraintsFromFile());
      }
      const size_t &helix_window( m_HydrophobicMomentWindowSizes( biol::GetSSTypes().HELIX));
      const size_t &strand_window( m_HydrophobicMomentWindowSizes( biol::GetSSTypes().STRAND));
      const size_t &coil_window( m_HydrophobicMomentWindowSizes( biol::GetSSTypes().COIL));

      if( !e_ScoreExposureMomentNiEDDA.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_NiEDDA);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMoment
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentNiEDDA = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score,
              m_Restraints,
              score::Accessibility::GetDefaultScheme() + environment.GetString()
            )
          )
        );
      }
      if( !e_ScoreExposureMomentOxygen.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_Oxygen);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMoment
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentOxygen = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score, m_Restraints, score::Accessibility::GetDefaultScheme() + environment.GetString()
            )
          )
        );
      }
      if( !e_ScoreExposureMomentMagnitudeNiEDDA.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_NiEDDA);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMomentMagnitude
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentMagnitudeNiEDDA = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score,
              m_Restraints,
              score::Accessibility::GetDefaultScheme() + "_magn_" + environment.GetString()
            )
          )
        );
      }
      if( !e_ScoreExposureMomentMagnitudeOxygen.IsDefined())
      {
        const AccessibilityAA::EnvironmentEnum environment( AccessibilityAA::e_Oxygen);
        const util::ShPtr< score::AccessibilityHydrophobicMoment> score
        (
          new score::AccessibilityHydrophobicMomentMagnitude
          (
            environment,
            pdb::Factory::GetSSETypeMinSizes( helix_window, strand_window, coil_window)
          )
        );
        e_ScoreExposureMomentMagnitudeOxygen = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::Accessibility
            (
              score, m_Restraints, score::Accessibility::GetDefaultScheme() + "_magn_" + environment.GetString()
            )
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void EPRAccessibilityData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void EPRAccessibilityData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EPRAccessibilityData::GetSerializer() const
    {
      io::Serializer serial( m_Handler.GetSerializer());
      serial.SetClassDescription( "Allows use of electron paramagnetic resonance-based accessibility data measurements");

      serial.AddInitializer
      (
        "helix window size",
        "Number of residues considered when computing hydrophobic moment in helices",
        io::Serialization::GetAgent( &m_HydrophobicMomentWindowSizes( 0)),
        "5"
      );

      serial.AddInitializer
      (
        "strand window size",
        "Number of residues considered when computing hydrophobic moment in strands",
        io::Serialization::GetAgent( &m_HydrophobicMomentWindowSizes( 1)),
        "5"
      );
      serial.AddInitializer
      (
        "loop window size",
        "Number of residues considered when computing hydrophobic moment in loops",
        io::Serialization::GetAgent( &m_HydrophobicMomentWindowSizes( 2)),
        "5"
      );
      return serial;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EPRAccessibilityData::Read( std::istream &ISTREAM)
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
    std::ostream &EPRAccessibilityData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
