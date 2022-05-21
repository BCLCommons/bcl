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

#ifndef BCL_RESTRAINT_INTERFACE_H_
#define BCL_RESTRAINT_INTERFACE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "score/bcl_score.fwd.hh"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_mutates.h"
#include "fold/bcl_fold_protocol_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Interface
    //! @brief interface class for restraint objects
    //! @details none
    //!
    //!
    //! @remarks example unnecessary
    //! @author alexanns
    //! @date Jul 29, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Interface :
      public virtual util::SerializableInterface
    {
    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief Clone function
      virtual Interface *Clone() const = 0;

    ///////////////////////////////////////////////////////////////////
    // functions from ProtocolInterface not needed by all restraints //
    ///////////////////////////////////////////////////////////////////

      //! @brief initialize the mutates and add them to Mutates enumerator
      virtual void InitializeScores()
      {
      }

      //! @brief sets the weights of scores in a weight set
      //! @param SCORE_WEIGHT_SET the score weight set that will be modified
      virtual void ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
      {
      }

      //! @brief initialize the mutates and add them to Mutates enumerator
      virtual void InitializeMutates()
      {
      }

      //! @brief modify the mutate tree used
      //! @param MUTATE_TREE MutateTree to be modified
      virtual void ModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
      {
      }

      //! @brief get the mutate tree associated with this protocol
      //! @return the mutate tree associated with this protocol
      virtual util::ShPtr< fold::MutateTree> GetMutateTree() const = 0;

      //! @brief merges this protocol's mutate tree into given mutate tree
      //! @param MUTATE_TREE tree into which to merge this protocol's tree
      virtual void MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const = 0;

    }; // class Interface

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_INTERFACE_H_
