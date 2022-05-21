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

#ifndef BCL_FOLD_PLACEMENT_SSE_NEXT_TO_SSE_H_
#define BCL_FOLD_PLACEMENT_SSE_NEXT_TO_SSE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "find/bcl_find.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementSSENextToSSE
    //! @brief class places an SSE next to an SSE already in the protein model
    //! @details This class first picks a neighbor SSE in the given protein model, next to which the given SSE will
    //! be placed, and then based on the SSTypes of these two SSEs it picks a random axis and random distance close
    //! to a preferred distance for that contact type and places the SSE.
    //!
    //! @see @link example_fold_placement_sse_next_to_sse.cpp @endlink
    //! @author karakam
    //! @date Apr 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementSSENextToSSE :
      public PlacementInterface< assemble::SSE, assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

    private:

      //! Locator to locate the neighbor sse
      util::Implementation
      <
        find::LocatorCriteriaInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>
      > m_NeighborSSELocator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PlacementSSENextToSSE();

      //! @brief constructor from a LocatorCriteriaInterface and a deviation value
      //! @param LOCATOR_CRITERIA LocatorCriteriaInterface to be used
      PlacementSSENextToSSE
      (
        const find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE
        > &LOCATOR_CRITERIA
      );

      //! @brief constructor from a ShPtr to LocatorCriteriaInterface
      //! @param SP_LOCATOR_CRITERIA ShPtr to LocatorCriteriaInterface to be used
      PlacementSSENextToSSE
      (
        const util::ShPtr
        <
          find::LocatorCriteriaInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>
        > &SP_LOCATOR_CRITERIA
      );

      //! @brief Clone function
      //! @return pointer to new PlacementSSENextToSSE
      PlacementSSENextToSSE *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief generate placement for the given SSE at a random orientation wrt to a located SSE from the given model
      //! @param SELECTED_SSE SiPtr to SSE to be placed
      //! @param PROTEIN_MODEL to which the SSE is going to be added
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SELECTED_SSE,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief generate placement for the given SSE at a random orientation wrt to a located SSE from the given model
      //! @param SELECTED_SSE SiPtr to SSE to be placed
      //! @param NEIGHBOR_SSE reference SSE
      static storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SELECTED_SSE,
        const assemble::SSE &NEIGHBOR_SSE
      );

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PlacementSSENextToSSE

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_PLACEMENT_SSE_NEXT_TO_SSE_H_
