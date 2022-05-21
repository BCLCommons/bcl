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

#ifndef BCL_DESCRIPTOR_ELEMENT_TYPE_PROPERTY_RETRIEVER_H_
#define BCL_DESCRIPTOR_ELEMENT_TYPE_PROPERTY_RETRIEVER_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_element_type_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ElementTypePropertyRetriever
    //! @brief Get the value of an ElementTypeData::Properties for all atoms in a small molecule
    //!
    //! @see @link example_descriptor_element_type_property_retriever.cpp @endlink
    //! @author mendenjl, kothiwsk
    //! @date Feb 11, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ElementTypePropertyRetriever :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    private:

    //////////
    // data //
    //////////

      std::string                 m_Alias;       //!< alias for this property, depends on the property chosen
      chemistry::ElementTypeData::PropertyEnum  m_Property;    //!< property to retrieve from each atom in the small molecule

    public:

      //! the instance of this class is created in chemistry::AtomProperties, so no s_Instance is necessary

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ElementTypePropertyRetriever();

      //! @brief constructor from a property
      //! @param PROPERTY the property to retrieve
      explicit ElementTypePropertyRetriever( const chemistry::ElementTypeData::Properties &PROPERTY);

      //! @brief virtual copy constructor
      ElementTypePropertyRetriever *Clone() const;

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
        return 1;
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

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    }; // class ElementTypePropertyRetriever

  } // namespace descriptor
} // namespace bcl

#endif //BCL_DESCRIPTOR_ELEMENT_TYPE_PROPERTY_RETRIEVER_H_
