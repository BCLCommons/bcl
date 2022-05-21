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

#ifndef BCL_FOLD_PLACEMENT_DOMAIN_INTERFACE_H_
#define BCL_FOLD_PLACEMENT_DOMAIN_INTERFACE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementDomainInterface
    //! @brief interface class for placing a Domain into a protein model
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Feb 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementDomainInterface :
      public util::ObjectInterface
    {

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief determines transformation matrices for placing SSEs from a Domain into a protein model
      //! @param DOMAIN_TO_PLACE domain containing SSEs to be placed
      //! @param PROTEIN_MODEL model that domain will be placed into
      //! @return transformation matrices for placing SSEs from a Domain into a protein model
      virtual storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool> Place
      (
        const assemble::DomainInterface &DOMAIN_TO_PLACE, const assemble::ProteinModel &PROTEIN_MODEL
      ) const = 0;

      //! @brief determines transformation matrices for placing SSEs from a Domain into a protein model
      //! @param DOMAIN_TO_PLACE domain containing SSEs to be placed
      //! @return transformation matrices for placing SSEs from a Domain into a protein model
      virtual storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool> Place
      (
        const assemble::DomainInterface &DOMAIN_TO_PLACE
      ) const = 0;

    }; // class PlacementDomainInterface

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PLACEMENT_DOMAIN_INTERFACE_H_
