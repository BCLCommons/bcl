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

#ifndef BCL_MODEL_FEATURE_INTERFACE_H_
#define BCL_MODEL_FEATURE_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureInterface
    //! @brief this interface is for features
    //!
    //! @remarks example unnecessary
    //! @author loweew, woetzen, butkiem1, mendenjl
    //! @date 09/14/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class FeatureInterface :
      public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief get size of DataSet
      //! @return the number of data items in the DataSet
      virtual size_t GetSize() const = 0;

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Matrix
      virtual const t_DataType *Begin() const = 0;

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Matrix
      virtual const t_DataType *End() const = 0;

      //! @brief bool indicating whether or not feature has ownership
      //! @return bool
      virtual bool HasOwnership() const = 0;

    }; // template FeatureInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_INTERFACE_H_
