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

#ifndef BCL_DESCRIPTOR_PARTIAL_H_
#define BCL_DESCRIPTOR_PARTIAL_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_type.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Partial
    //! @brief Select specific parts of other descriptors
    //!
    //! @see @link example_descriptor_partial.cpp @endlink
    //! @author mendenjl
    //! @date Feb 13, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Partial :
      public Base< t_DataType, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      //! descriptor to consider only a part of
      util::Implementation< Base< t_DataType, t_ReturnType> > m_Property;
      storage::Vector< size_t> m_Columns; //!< Columns to select of that property
      size_t m_MaxColumn; //!< Max(m_Columns)

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Partial() :
        m_MaxColumn( 0)
      {
      }

      //! @brief virtual copy constructor
      Partial *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      Type::Symmetry GetSymmetry() const
      {
        return m_Property->GetSymmetry();
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const
      {
        return m_Property->ConsiderRepeatedElements();
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
      void SetDimensionHook();

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, t_ReturnType> > GetInternalDescriptors();

    private:

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< t_ReturnType> &STORAGE
      );

    }; // class Partial

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< chemistry::AtomConformationalInterface, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Partial< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_PARTIAL_H_
