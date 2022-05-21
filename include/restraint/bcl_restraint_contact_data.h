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

#ifndef BCL_RESTRAINT_CONTACT_DATA_H_
#define BCL_RESTRAINT_CONTACT_DATA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_atom_distance.h"
#include "bcl_restraint_handler_atom_distance_assigned.h"
#include "bcl_restraint_interface.h"
#include "fold/bcl_fold_scores.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ContactData
    //! @brief Reads, stores, and writes contact data using CASP ROLL format
    //! @details This class allows for the usage of contact data as restraints for folding/scoring.  The format is used
    //!          as part of the CASP ROLL experiment. http://www.predictioncenter.org/casprol/index.cgi
    //!
    //! @see @link example_restraint_contact_data.cpp @endlink
    //! @author weinerbe
    //! @date Mar 19, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ContactData :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! the handler that will be used to read and write the restraints
      HandlerAtomDistanceAssigned m_Handler;

      //! the atom distance restraint data
      util::ShPtrVector< AtomDistance> m_Restraints;

      //! expected fraction of contacts that the protein model should satisfy
      double m_Fraction;

      //! distance cutoff for contact
      static const double s_ContactDistance;

    public:

      //! contact restraint atom distance score
      static const util::SiPtr< const score::RestraintAtomDistanceAssignment> s_ScoreContactRestraintAssignment;

      static const util::SiPtr< const score::RestraintAtomDistance>           s_ScoreContactDistance;

      //! contact restraint
      static fold::Score e_ScoreContactRestraint;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ContactData();

      //! @brief Clone function
      //! @return pointer to new ContactData
      ContactData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get the default file extension
      //! @return the default file extension
      const std::string &GetDefaultExtension() const;

      //! @brief gives reference to the specific type of data this restraint uses
      //! @return gives reference to the specific type of data this restraint uses
      const util::ShPtrVector< AtomDistance> &GetAtomDistanceRestraints() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initialize the scores and add them to Scores enumerator
      void InitializeScores();

      //! @brief sets the weights of scores in a weight set
      //! @param SCORE_WEIGHT_SET the score weight set that will be modified
      void ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const;

      //! @brief get the mutate tree associated with this protocol
      //! @return the mutate tree associated with this protocol
      util::ShPtr< fold::MutateTree> GetMutateTree() const
      {
        util::ShPtr< fold::MutateTree> sp_mutate_tree( new fold::MutateTree());
        ModifyMutateTree( *sp_mutate_tree);
        return sp_mutate_tree;
      }

      //! @brief merges this protocol's mutate tree into given mutate tree
      //! @param MUTATE_TREE tree into which to merge this protocol's tree
      void MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief reads restraints formatted for this restraint type from an istream
      //! @return istream restraints formatted for this restraint type were read from
      std::istream &ReadRestraints( std::istream &ISTREAM);

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ContactData

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_CONTACT_DATA_H_
