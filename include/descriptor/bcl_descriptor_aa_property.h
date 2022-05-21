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

#ifndef BCL_DESCRIPTOR_AA_PROPERTY_H_
#define BCL_DESCRIPTOR_AA_PROPERTY_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "biol/bcl_biol_aa_base.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAProperty
    //! @brief Returns the values for selected PropertyType for a given amino acid
    //!
    //! @see @link example_descriptor_aa_property.cpp @endlink
    //! @author teixeipl
    //! @date Feb 6, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAProperty :
      public BaseElement< biol::AABase, float>
    {

    private:

    //////////
    // data //
    //////////

      //! set of PropertyType to be used
      biol::AATypeData::PropertyTypeEnum m_Property;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a property to use
      //! @param PROPERTY property to be used
      AAProperty( const biol::AATypeData::PropertyType &PROPERTY);

      //! @brief Clone function
      //! @return pointer to new AABlastProfile
      AAProperty *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class AAProperty

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_PROPERTY_H_
