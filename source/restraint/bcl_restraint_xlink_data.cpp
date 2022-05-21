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
#include "restraint/bcl_restraint_xlink_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "score/bcl_score_restraint_xlink.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> XlinkData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new XlinkData())
    );

    //! score for evaluating the agreement of a protein model with cross-linking restraints
    fold::Score XlinkData::e_ScoreXlinkRestraint( fold::GetScores().e_Undefined);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    XlinkData::XlinkData() :
      m_Handler( GetDefaultExtension())
    {
    }

    //! @brief returns a pointer to a new XlinkData
    //! @return pointer to a new XlinkData
    XlinkData *XlinkData::Clone() const
    {
      return new XlinkData( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &XlinkData::GetAlias() const
    {
      static const std::string s_name( "Xlink");
      return s_name;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &XlinkData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the default file extension of files containing cross-link restraints
    //! @return the default file extension of files containing cross-link restraints
    const std::string &XlinkData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".xlink_bcl");
      return s_extension;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void XlinkData::InitializeScores()
    {
      // read in restraints
      if( !m_Restraints.IsDefined() && m_Handler.Exists())
      {
        m_Restraints = util::CloneToShPtr( m_Handler.ReadRestraintsFromFile());
      }
      if( !e_ScoreXlinkRestraint.IsDefined())
      {
        e_ScoreXlinkRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>( new score::RestraintXlink( m_Restraints))
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void XlinkData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreXlinkRestraint, 10);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void XlinkData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer XlinkData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription( "Cross linking restraints");
      serial.AddInitializer
      (
        "",
        "Handler for reading cross-link data",
        io::Serialization::GetAgent( &m_Handler),
        HandlerAtomDistanceAssigned( GetDefaultExtension()).GetString()
      );
      return serial;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from an input stream
    //! @param ISTREAM input stream to read members from
    //! @return returns the input stream
    std::istream &XlinkData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Restraints, ISTREAM);

      return ISTREAM;
    }

    //! @brief writes members into an output stream
    //! @param OSTREAM output stream to write members into
    //! @INDENT number of indentations to use
    //! @return returns the output stream
    std::ostream &XlinkData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
