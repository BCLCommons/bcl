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

#ifndef BCL_DESCRIPTOR_ATOM_AROMATICITY_AXES_H_
#define BCL_DESCRIPTOR_ATOM_AROMATICITY_AXES_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomAromaticityAxes
    //! @brief is a class which calculates the
    //!
    //! @see @link example_descriptor_atom_aromaticity_axes.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Sep 06, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomAromaticityAxes :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    public:

      //! the instance of this class is created in chemistry::AtomProperties, so no s_Instance is necessary

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      AtomAromaticityAxes *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 3;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_NeverCache;
      }
    }; // class AtomAromaticityAxes

  } // namespace descriptor
} // namespace bcl

#endif //BCL_DESCRIPTOR_ATOM_AROMATICITY_AXES_H_
