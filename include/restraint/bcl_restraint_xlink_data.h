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

#ifndef BCL_RESTRAINT_XLINK_DATA_H_
#define BCL_RESTRAINT_XLINK_DATA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_atom_distance.h"
#include "bcl_restraint_handler_atom_distance_assigned.h"
#include "bcl_restraint_interface.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "command/bcl_command_flag_interface.h"
#include "score/bcl_score_protein_model.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class XlinkData
    //! @brief Restraint protocol that incorporates distance restraint obtained from cross-linking experiments into the
    //! de novo folding algorithm.
    //!
    //! @see @link example_restraint_xlink_data.cpp @endlink
    //! @author fischea
    //! @date June 17, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API XlinkData :
      public Interface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! score for evaluating the agreement of a protein model with cross-linking restraints
      static fold::Score e_ScoreXlinkRestraint;

    private:

      //! the atom distance restraint data
      util::ShPtr< util::ShPtrVector< AtomDistance> > m_Restraints;

      //! the handler that will be used to read and write the restraints
      HandlerAtomDistanceAssigned m_Handler;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default consturctor
      XlinkData();

      //! @brief returns a pointer to a new XlinkData
      //! @return pointer to a new XlinkData
      XlinkData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the default file extension of files containing cross-link restraints
      //! @return the default file extension of files containing cross-link restraints
      const std::string &GetDefaultExtension() const;

      //! @brief gives reference to the specific type of data this restraint uses
      //! @return gives reference to the specific type of data this restraint uses
      const util::ShPtrVector< AtomDistance> &GetAtomDistanceRestraints() const
      {
        return *m_Restraints;
      }

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

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief reads members from an input stream
      //! @param ISTREAM input stream to read members from
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes members into an output stream
      //! @param OSTREAM output stream to write members into
      //! @INDENT number of indentations to use
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class XlinkData

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_XLINK_DATA_H_
