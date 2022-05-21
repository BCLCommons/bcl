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

#ifndef BCL_FOLD_PROTOCOL_INTERFACE_H_
#define BCL_FOLD_PROTOCOL_INTERFACE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "command/bcl_command.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "mc/bcl_mc.fwd.hh"
#include "pdb/bcl_pdb.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_mutate_tree.h"
#include "bcl_fold_mutates.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProtocolInterface
    //! @brief interface for defining a fold protocol
    //! @details interface class that defines all necessary functions that has to be overwritten for initializing a new
    //! folding protocol. The functions that have to be overwritten include
    //!  - GetAllFlags - returns all unique flags for this protocol
    //!  - ModifyStartModel - how to change the starting model in each iteration of folding
    //!  - InitializeScores - add unique scores for this protocol to Scores enum
    //!  - ModifyScoreWeightSet - modify the score weight set accordingly
    //!  - InitiailizeMutates - add unique mutates for this protocol to Mutates enum
    //!  - ModifyMutateTree - modify the mutate tree accordingly
    //!  - ModifyTerminate - modify the terminate object accordingly
    //!  - GetDescription - return a one sentence description of this protocol
    //!  - GetReadMe - return a detail readme for this protocol
    //!
    //! @remarks example unnecessary
    //! @author karakam, fischea
    //! @date Nov 11, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProtocolInterface :
      public virtual util::ObjectInterface
    {

    public:

    ///////////
    // flags //
    ///////////

      //! @brief returns all flags that are specialized for this protocol
      //! @return all flags that are specialized for this protocol
      virtual const util::ShPtrVector< command::FlagInterface> &GetAllFlags() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief modifies the start model
      //! @param START_MODEL Start model to be modified
      virtual void ModifyStartModel( assemble::ProteinModel &START_MODEL) const = 0;

      //! @brief initialize the scores and add them to Scores enumerator
      virtual void InitializeScores() = 0;

      //! @brief modify the score weight set
      //! @param SCORE_WEIGHT_SET Score weight set
      virtual void ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const = 0;

      //! @brief modify the terminate object
      //! @param CRITERION which will be modified by protocols
      //! @param STAGE Stage in which this terminate will be used
      virtual void ModifyCriterion
      (
        opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
        const mc::Stage &STAGE
      ) const = 0;

      //! @brief modify the printer object
      //! @param PRINTER which will be modified by protocols
      //! @param STAGE Stage in which this terminate will be used
      virtual void ModifyPrinter
      (
        mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
        const mc::Stage &STAGE
      ) const = 0;

      //! @brief modify the pdb factory object
      //! @param FACTORY pdb factory to be modified
      //! @param STAGE Stage in which this terminate will be used
      virtual void ModifyFactory
      (
        util::ShPtr< pdb::Factory> &FACTORY,
        const mc::Stage &STAGE
      ) const = 0;

      //! @brief initialize the mutates and add them to Mutates enumerator
      virtual void InitializeMutates() = 0;

      //! @brief modify the mutate tree used
      //! @param MUTATE_TREE MutateTree to be modified
      virtual void ModifyMutateTree( MutateTree &MUTATE_TREE) const = 0;

      //! @brief reset the protocol (i.e. empty cache of scoring function)
      virtual void Reset()
      {
      }

      //! @brief get the mutate tree associated with this protocol
      //! @return the mutate tree associated with this protocol
      virtual util::ShPtr< MutateTree> GetMutateTree() const = 0;

      //! @brief merges this protocol's mutate tree into given mutate tree
      //! @param MUTATE_TREE tree into which to merge this protocol's tree
      virtual void MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const = 0;

    ////////////
    // readme //
    ////////////

      //! @brief returns string containing short description of the protocol
      //! @return string containing short description of the protocol
      virtual const std::string &GetDescription() const = 0;

      //! @brief returns readme information
      //! @return string containing information about application
      virtual const std::string &GetReadMe() const = 0;

    }; // class ProtocolInterface

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PROTOCOL_INTERFACE_H_
