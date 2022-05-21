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

#ifndef BCL_DESCRIPTOR_COMBINE_H_
#define BCL_DESCRIPTOR_COMBINE_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

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
    //! @class Combine
    //! @brief combines the output of descriptors into one vector
    //!
    //! @see @link example_descriptor_combine.cpp @endlink
    //! @author mendenjl
    //! @date Jan 15, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_ReturnType>
    class Combine :
      public Base< t_DataType, t_ReturnType>
    {
    private:

    //////////
    // data //
    //////////

      //! descriptors interfaces which results to combine
      storage::Vector< util::Implementation< Base< t_DataType, t_ReturnType> > > m_Properties;

      //! generic type of the descriptor
      Type m_Type;

      //! a running sum of the total features per object
      size_t m_FeatureColumnsPerObject;

      CachePreference m_CachePreference;

    public:

      typedef typename
        storage::Vector< util::Implementation< Base< t_DataType, t_ReturnType> > >::const_iterator const_iterator;
      typedef typename storage::Vector< util::Implementation< Base< t_DataType, t_ReturnType> > >::iterator iterator;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Combine( const CachePreference &CACHE_PREF = e_IgnoreCache);

      //! @brief comparison operator
      bool operator ==( const Combine &COMBINE) const;

      //! @brief virtual copy constructor
      Combine *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the type of this descriptor
      //! @return the type of this descriptor (should ignore dimension setting)
      Type GetType() const
      {
        return m_Type;
      }

      //! @brief const iterator begin
      const_iterator Begin() const
      {
        return m_Properties.Begin();
      }

      //! @brief const iterator end
      const_iterator End() const
      {
        return m_Properties.End();
      }

      //! @brief iterator begin
      iterator Begin()
      {
        return m_Properties.Begin();
      }

      //! @brief iterator end
      iterator End()
      {
        return m_Properties.End();
      }

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief push back a SmallMoleculeProperty
      //! @param PROPERTY a SmallMoleculeProperty
      void PushBack( const util::Implementation< Base< t_DataType, t_ReturnType> > &PROPERTY);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief Get the code / label set for the combine with sizes of each property
      //! @return the code / label for the combine with sizes of each property
      //! the feature code set
      model::FeatureLabelSet GetLabelsWithSizes() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    protected:

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return m_CachePreference;
      }

      //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
      void SetDimensionHook();

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_FeatureColumnsPerObject;
      }

      //! @brief get the normal dimension for this descriptor
      //! @return the normal dimension for this descriptor
      size_t GetNormalDimension() const
      {
        return m_Type.GetDimension();
      }

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

    }; // class Combine

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< chemistry::AtomConformationalInterface, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< chemistry::AtomConformationalInterface, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< biol::AABase, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< biol::AABase, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< biol::Mutation, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< biol::Mutation, float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< char, char>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Combine< char, float>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_COMBINE_H_
