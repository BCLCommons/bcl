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
#include "restraint/bcl_restraint_pre_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "nmr/bcl_nmr_star_noe_handler.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "score/bcl_score_restraint_energy_well.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize scores
    fold::Score PREData::e_ScorePRERestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PREData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new PREData())
    );

    const util::SiPtr< const score::RestraintAtomDistance> PREData::s_Score
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          score::RestraintEnergyWell( score::RestraintEnergyWell::e_PRE),
          1.0,
          "pre_restraint"
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    PREData::PREData() :
      m_Handler( GetDefaultHandler()),
      m_Restraints()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    PREData::PREData( const HandlerBase< util::ShPtrVector< AtomDistance> > &HANDLER) :
      m_Handler( HANDLER),
      m_Restraints()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PREData
    PREData *PREData::Clone() const
    {
      return new PREData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &PREData::GetAlias() const
    {
      static const std::string s_name( "PRE");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PREData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &PREData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".pre_star");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const nmr::StarNOEHandler &PREData::GetDefaultHandler()
    {
      static const nmr::StarNOEHandler s_handler( ".pre_star");
      return s_handler;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &PREData::GetAtomDistanceRestraints() const
    {
      return *m_Restraints;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void PREData::InitializeScores()
    {
      if( !e_ScorePRERestraint.IsDefined())
      {
        // read in restraints
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }

        util::ShPtr< score::RestraintAtomDistance> score( s_Score->Clone());
        score->SetRestraints( m_Restraints);
        e_ScorePRERestraint = fold::GetScores().AddScore( score);
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void PREData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScorePRERestraint, 5);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void PREData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PREData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription( "Paramagenetic relaxation enhancement based restraints");
      serial.AddInitializer
      (
        "",
        "Handler for reading PREs",
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

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PREData::Read( std::istream &ISTREAM)
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
    std::ostream &PREData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Handler, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
