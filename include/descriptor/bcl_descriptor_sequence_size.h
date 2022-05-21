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

#ifndef BCL_DESCRIPTOR_SEQUENCE_SIZE_H_
#define BCL_DESCRIPTOR_SEQUENCE_SIZE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SequenceSize
    //! @brief Returns the sequence size as either a floating point number or a string
    //!
    //! @tparam t_DataType type within sequence being iterated over
    //! @tparam t_ReturnType elemental type returned (char or float)
    //!
    //! @see @link example_descriptor_sequence_size.cpp @endlink
    //! @author mendenjl
    //! @date Feb 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class SequenceSize :
      public BaseSequence< t_DataType, t_ReturnType>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SequenceSize
      SequenceSize< t_DataType, t_ReturnType> *Clone() const;

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
      //! @note only one value is needed if this is a numeric descriptor, for char descriptors, assume that 99999 is the
      //! @note max, so 5 characters is sufficient
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< t_ReturnType> &STORAGE);

    }; // class SequenceSize

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< chemistry::AtomConformationalInterface, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSize< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEQUENCE_SIZE_H_
