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

#ifndef BCL_MODEL_COLLECT_FEATURES_INTERFACE_H_
#define BCL_MODEL_COLLECT_FEATURES_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectFeaturesInterface
    //! @brief interface for methods of collecting features by score
    //!
    //! @tparam t_DataType can be float, double, int
    //!
    //! @author mendenjl
    //! @date Oct 31, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectFeaturesInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone the CollectFeaturesInterface
      //! @return pointer to new CollectFeaturesInterface
      virtual CollectFeaturesInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! Collect the feature indices
      //! @param ARGUMENT scores for the features
      //! @return returns the selected feature indices
      virtual storage::Vector< size_t> Collect( const linal::Vector< float> &ARGUMENT) const = 0;

    ///////////////
    // operators //
    ///////////////

    }; // class CollectFeaturesInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_COLLECT_FEATURES_INTERFACE_H_
