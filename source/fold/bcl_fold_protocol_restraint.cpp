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
#include "fold/bcl_fold_protocol_restraint.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_mutates.h"
#include "restraint/bcl_restraint_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolRestraint::ProtocolRestraint()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolRestraint
    ProtocolRestraint *ProtocolRestraint::Clone() const
    {
      return new ProtocolRestraint( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolRestraint &ProtocolRestraint::GetInstance()
    {
      static ProtocolRestraint s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolRestraint::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolRestraint::GetAlias() const
    {
      static const std::string s_name( "ProtocolRestraint");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolRestraint::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol for using restraints in folding run");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolRestraint::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector
      (
        util::ShPtrVector< command::FlagInterface>::Create
        (
          restraint::GetFlagRestraintsTypes(),
          restraint::GetFlagRestraintsFilePrefix()
        )
      );

      // end
      return s_all_flags_vector;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolRestraint::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolRestraint::InitializeScores()
    {
      // get the restraint types
      const storage::Set< restraint::Type> restraint_types
      (
        restraint::GetFlagRestraintsTypes()->GetObjectSet< restraint::Type>()
      );

      // iterate through the restraint types and add their scores to scores
      for
      (
        storage::Set< restraint::Type>::const_iterator
          restraint_itr( restraint_types.Begin()), restraint_itr_end( restraint_types.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        // copy the restraint data so that the scores can be set
        restraint::Type sp_restraint_data( *restraint_itr);
        sp_restraint_data->InitializeScores();
      }
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolRestraint::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      // get the restraint types
      const storage::Set< restraint::Type> restraint_types
      (
        restraint::GetFlagRestraintsTypes()->GetObjectSet< restraint::Type>()
      );

      // iterate through the restraint types to modify their score weights
      for
      (
        storage::Set< restraint::Type>::const_iterator
          restraint_itr( restraint_types.Begin()), restraint_itr_end( restraint_types.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        ( *restraint_itr)->ModifyScoreWeightSet( SCORE_WEIGHT_SET);
      }
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolRestraint::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolRestraint::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolRestraint::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolRestraint::InitializeMutates()
    {
      // get the restraint types
      const storage::Set< restraint::Type> restraint_types
      (
        restraint::GetFlagRestraintsTypes()->GetObjectSet< restraint::Type>()
      );

      // iterate through the restraint types and add their mutates to mutates
      for
      (
        storage::Set< restraint::Type>::const_iterator
          restraint_itr( restraint_types.Begin()), restraint_itr_end( restraint_types.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        // copy the restraint data so that the mutates can be set
        restraint::Type sp_restraint_data( *restraint_itr);
        sp_restraint_data->InitializeMutates();
      }
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolRestraint::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      // get the restraint types
      const storage::Set< restraint::Type> restraint_types
      (
        restraint::GetFlagRestraintsTypes()->GetObjectSet< restraint::Type>()
      );

      // let each restraint modify the mutate tree
      for
      (
        storage::Set< restraint::Type>::const_iterator
          restraint_itr( restraint_types.Begin()), restraint_itr_end( restraint_types.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        ( *restraint_itr)->ModifyMutateTree( MUTATE_TREE);
      }
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolRestraint::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolRestraint::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolRestraint::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "protocol for using restraints such as NOE, distance restraints"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolRestraint::GetReadMe() const
    {
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "This protocol adapts the BCL::Fold method to function with experimental restraints. To see what is available "
        "use the -help flag and view the description for restraint_types.\n"
        "When using restraint folding in a publication, and for more detailed information regarding the method, please "
        "cite the publication describing the application's development which is currently in preparation. Refer "
        "to www.meilerlab.org for future details.\n\n"
        "Specific flags:\n"
        "-restraint_types #Types of restraints available to use\n"
        "-restraint_prefix #Prefix for the restraint file(s). BCL will add the \".\" and the appropriate file "
        "extension defined by the type of restraint used (e.g. NMR NOE restraints have the extension, \".noe_star\""
      );

      // end
      return s_readme;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolRestraint::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolRestraint::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
