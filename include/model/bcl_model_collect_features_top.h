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

#ifndef BCL_MODEL_COLLECT_FEATURES_TOP_H_
#define BCL_MODEL_COLLECT_FEATURES_TOP_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_collect_features_interface.h"
#include "find/bcl_find_collector_interface.h"
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectFeaturesTop
    //! @brief reduces a data set to a series of cluster centers that cover the same space
    //! There may be a smaller data set that covers the same space, which may be obtained using k-means
    //! This method should be faster, and can be used on its own or to provide the initial cluster points
    //! for the k-means algorithm
    //!
    //! @author mendenjl
    //! @see @link example_model_collect_features_top.cpp @endlink
    //! @date Apr 03, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectFeaturesTop :
      public CollectFeaturesInterface
    {
    private:

    //////////
    // data //
    //////////

      size_t m_NumberToKeep; //!< number of features to keep

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from parameter
      //! @param NUMBER_TO_KEEP number of features to keep
      explicit CollectFeaturesTop( const size_t &NUMBER_TO_KEEP);

      //! @brief Clone function
      //! @return pointer to new CollectFeaturesTop
      CollectFeaturesTop *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief collect the top m_NumberToKeep value indices by score
      //! @param SCORES any values
      //! @return top m_NumberToKeep value indices
      storage::Vector< size_t> Collect( const linal::Vector< float> &SCORES) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class CollectFeaturesTop

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_COLLECT_FEATURES_TOP_H_

