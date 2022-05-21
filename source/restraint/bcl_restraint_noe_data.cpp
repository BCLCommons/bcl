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
#include "restraint/bcl_restraint_noe_data.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "fold/bcl_fold_default_mutates.h"
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_placement_sse_distance_restraint.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "score/bcl_score_restraint_noe_attraction.h"
#include "score/bcl_score_restraint_noe_knowledge_based.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score NOEData::e_ScoreNOERestraint( fold::GetScores().e_Undefined);
    fold::Score NOEData::e_ScoreNOEPenalty( fold::GetScores().e_Undefined);

    //! the scores used internally
    const util::SiPtr< const score::RestraintAtomDistance> NOEData::s_ScoreNOERestraint
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintNoeKnowledgeBased(),
          1.0,
          "noe_restraint"
        )
      )
    );
    const util::SiPtr< const score::RestraintAtomDistance> NOEData::s_ScoreNOEPenalty
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintNoeAttraction(),
          1.0,
          "noe_penalty"
        )
      )
    );

    // initialize mutate
    fold::Mutate NOEData::e_MutateAddSSENOE( fold::GetMutates().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> NOEData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new NOEData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NOEData::NOEData() :
      m_Handler( GetDefaultHandler()),
      m_Restraints()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    NOEData::NOEData( const HandlerBase< util::ShPtrVector< AtomDistance> > &HANDLER) :
      m_Handler( HANDLER),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new NOEData
    NOEData *NOEData::Clone() const
    {
      return new NOEData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NOEData::GetAlias() const
    {
      static const std::string s_name( "NOE");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NOEData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &NOEData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".noe_star");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const nmr::StarNOEHandler &NOEData::GetDefaultHandler()
    {
      static const nmr::StarNOEHandler s_handler( ".noe_star");
      return s_handler;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &NOEData::GetAtomDistanceRestraints() const
    {
      return *m_Restraints;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void NOEData::InitializeScores()
    {
      if( !e_ScoreNOERestraint.IsDefined())
      {
        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }
        util::ShPtr< score::RestraintAtomDistance> noe_restraint( s_ScoreNOERestraint->Clone());
        noe_restraint->SetRestraints( m_Restraints);
        util::ShPtr< score::RestraintAtomDistance> noe_penalty( s_ScoreNOEPenalty->Clone());
        noe_penalty->SetRestraints( m_Restraints);
        e_ScoreNOERestraint = fold::GetScores().AddScore( noe_restraint);
        e_ScoreNOEPenalty = fold::GetScores().AddScore( noe_penalty);
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void NOEData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreNOERestraint, 5);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreNOEPenalty, 5);
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void NOEData::InitializeMutates()
    {
      if( !e_MutateAddSSENOE.IsDefined())
      {
        // random picker
        const assemble::PickSSERandom picker_sse_random;

        // pool picker
        const find::PickCriteriaWrapper
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        > picker_sse_pool( picker_sse_random);

        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }

        // placement
        const fold::PlacementSSEDistanceRestraint placement( m_Restraints);

        // add SSE next to SSE weighted by # of restraints
        e_MutateAddSSENOE = fold::GetMutates().AddMutate
        (
          fold::MutateProteinModelSSEAdd
          (
            picker_sse_pool,
            placement,
            fold::MutateTree::GetMutateTypeName( fold::MutateTree::e_Add) + "_sse_next_to_sse_noe"
          )
        );
      }
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void NOEData::ModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      // get the probability of normal SSE add
      const double original_prob
      (
        MUTATE_TREE.GetMutateProbability( fold::DefaultMutates::GetInstance().e_AddSSENextToSSE)
      );

      // set the original probability to zero
      MUTATE_TREE.SetMutateProbability( fold::DefaultMutates::GetInstance().e_AddSSENextToSSE, 0.0);

      // set the NOE probability
      MUTATE_TREE.SetMutateProbability( e_MutateAddSSENOE, original_prob);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void NOEData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NOEData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription
      (
        "Nuclear overhauser effect (NOE) restraints; allows for placement of SSEs to satisfy NOE data."
      );
      serial.AddInitializer
      (
        "",
        "Handler for reading NOEs",
        io::Serialization::GetAgent( &m_Handler),
        util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > >( GetDefaultHandler()).GetString()
      );
      serial.AddInitializer
      (
        "spin label length",
        "# of bonds the spin label is away from the backbone",
        io::Serialization::GetAgent( &score::RestraintNMRDistanceInterface::GetSpinLabelLength()),
        "6"
      );
      return serial;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &NOEData::ReadRestraints( std::istream &ISTREAM)
    {
      // read the restraints
      m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraints( ISTREAM));

      // endStage
      return ISTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &NOEData::Read( std::istream &ISTREAM)
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
    std::ostream &NOEData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
