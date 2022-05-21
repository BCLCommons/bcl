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
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "restraint/bcl_restraint_contact_data.h"
#include "score/bcl_score_restraint_atom_attraction.h"
#include "score/bcl_score_restraint_atom_distance.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize distance cutoff for contact
    const double ContactData::s_ContactDistance( 8.0);

    // initialize score
    fold::Score ContactData::e_ScoreContactRestraint( fold::GetScores().e_Undefined);

    //! contact restraint atom distance score
    const util::SiPtr< const score::RestraintAtomDistanceAssignment> ContactData::s_ScoreContactRestraintAssignment
    (
      util::Enumerated< score::RestraintAtomDistanceAssignment>::AddInstance
      (
        new score::RestraintAtomAttraction
        (
          score::RestraintAtomAttraction::GetDefaultDepthRange(),
          0,
          s_ContactDistance,
          true,
          score::RestraintAtomAttraction::GetDefaultScheme() + "_" + ContactData().GetAlias()
        )
      )
    );

    const util::SiPtr< const score::RestraintAtomDistance> ContactData::s_ScoreContactDistance
    (
      util::Enumerated< score::RestraintAtomDistance>::AddInstance
      (
        new score::RestraintAtomDistance
        (
          *s_ScoreContactRestraintAssignment,
          1.0,
          "contact_restraint"
        )
      )
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ContactData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new ContactData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ContactData::ContactData() :
      m_Handler( ".RR", "-SASACCN", 0.0, 10.0, 8.0),
      m_Restraints(),
      m_Fraction( 1.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ContactData
    ContactData *ContactData::Clone() const
    {
      return new ContactData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ContactData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ContactData::GetAlias() const
    {
      static const std::string s_name( "Contact");
      return s_name;
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &ContactData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".RR");
      return s_extension;
    }

    //! @brief gives reference to the specific type of data this restraint uses
    //! @return gives reference to the specific type of data this restraint uses
    const util::ShPtrVector< AtomDistance> &ContactData::GetAtomDistanceRestraints() const
    {
      return m_Restraints;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void ContactData::InitializeScores()
    {
      if( !e_ScoreContactRestraint.IsDefined())
      {
        // read in restraints
        if( m_Restraints.IsEmpty())
        {
          m_Restraints = m_Handler.ReadRestraintsFromFile();
        }

        e_ScoreContactRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintAtomDistance
            (
              *s_ScoreContactRestraintAssignment,
              m_Fraction,
              "contact_restraint",
              util::CloneToShPtr( m_Restraints)
            )
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void ContactData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreContactRestraint, 500);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ContactData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads restraints formatted for this restraint type from an istream
    //! @return istream restraints formatted for this restraint type were read from
    std::istream &ContactData::ReadRestraints( std::istream &ISTREAM)
    {
      m_Restraints = m_Handler.ReadRestraints( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ContactData::GetSerializer() const
    {
      io::Serializer serial( m_Handler.GetSerializer());
      serial.SetClassDescription( "Predicted or known AA or atom contacts (within 8A by default)");
      serial.AddInitializer
      (
        "fraction",
        "fraction of contacts that the protein model is expected to satisfy. "
        "The best scoring contacts (eq. to this fraction of the contacts in the file) will be used when computing the "
        "score. Useful when contacts are less than certain and the folding process is expected to satisfy only a "
        "fraction of them",
        io::Serialization::GetAgent( &m_Fraction),
        "1.0"
      );
      return serial;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ContactData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Restraints, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ContactData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
