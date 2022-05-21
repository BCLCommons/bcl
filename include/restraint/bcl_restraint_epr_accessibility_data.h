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

#ifndef BCL_RESTRAINT_EPR_ACCESSIBILITY_DATA_H_
#define BCL_RESTRAINT_EPR_ACCESSIBILITY_DATA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "fold/bcl_fold.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_accessibility_profile.h"
#include "bcl_restraint_handler_accessibility_aa.h"
#include "bcl_restraint_handler_base.h"
#include "bcl_restraint_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EPRAccessibilityData
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_restraint_EPRAccessibilityData.cpp @endlink
    //! @author alexanns
    //! @date Jan 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EPRAccessibilityData :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! the handler that will be used to read and write the restraints
      HandlerAccessibilityAA m_Handler;

      //! hydrophobic moment window sizes for sse types helix, strand, coil
      storage::VectorND< 3, size_t> m_HydrophobicMomentWindowSizes;

      //! the atom distance restraint data
      util::ShPtr< AccessibilityProfile> m_Restraints;

    public:

      //! EPR accessibility hydrophobic moment for NiEDDA
      static fold::Score e_ScoreExposureMomentNiEDDA;

      //! EPR accessibility hydrophobic moment for Oxygen
      static fold::Score e_ScoreExposureMomentOxygen;

      //! EPR accessibility hydrophobic moment for NiEDDA
      static fold::Score e_ScoreExposureMomentMagnitudeNiEDDA;

      //! EPR accessibility hydrophobic moment for Oxygen
      static fold::Score e_ScoreExposureMomentMagnitudeOxygen;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking member variable parameters
      //! @param HANDLER the handler that will be used to read and write the restraints
      EPRAccessibilityData
      (
        const HandlerAccessibilityAA &HANDLER = GetDefaultHandler()
      );

      //! @brief Clone function
      //! @return pointer to new EPRAccessibilityData
      EPRAccessibilityData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the default file extension
      //! @return the default file extension
      const std::string &GetDefaultExtension() const;

      //! @brief get the default restraint handler
      //! @return the default restraint handler
      static const HandlerAccessibilityAA &GetDefaultHandler();

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

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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

    }; // class EPRAccessibilityData

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_EPR_ACCESSIBILITY_DATA_H_
