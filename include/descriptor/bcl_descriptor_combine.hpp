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

// include header of this class
#include "bcl_descriptor_combine.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"
#include "model/bcl_model_feature_label_set.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType, typename t_ReturnType>
    Combine< t_DataType, t_ReturnType>::Combine( const CachePreference &CACHE_PREF) :
      m_FeatureColumnsPerObject( 0),
      m_CachePreference( CACHE_PREF)
    {
    }

    //! @brief comparison operator
    template< typename t_DataType, typename t_ReturnType>
    bool Combine< t_DataType, t_ReturnType>::operator ==( const Combine< t_DataType, t_ReturnType> &COMBINE) const
    {
      // check whether both objects have the same members
      return m_Properties == COMBINE.m_Properties;
    }

    //! @brief virtual copy constructor
    template< typename t_DataType, typename t_ReturnType>
    Combine< t_DataType, t_ReturnType> *Combine< t_DataType, t_ReturnType>::Clone() const
    {
      return new Combine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Combine< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Combine< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_name( "Combine");
      return s_name;
    }

    //! @brief push back a SmallMoleculeProperty
    //! @param PROPERTY a SmallMoleculeProperty
    template< typename t_DataType, typename t_ReturnType>
    void Combine< t_DataType, t_ReturnType>::PushBack
    (
      const util::Implementation< Base< t_DataType, t_ReturnType> > &PROPERTY
    )
    {
      // add the property to the vector
      m_Properties.PushBack( PROPERTY);

      // copy the type to see if it changes when generalized
      Type old_type( m_Type);
      m_Type.GeneralizeToHandle( PROPERTY->GetType());

      // check whether the dimension changed
      if( old_type.GetDimension() != m_Type.GetDimension())
      {
        // recalculate m_FeatureColumnsPerObject and reset the dimension for all sub-objects
        m_FeatureColumnsPerObject = 0;
        for( iterator itr( m_Properties.Begin()), itr_end( m_Properties.End()); itr != itr_end; ++itr)
        {
          ( *itr)->SetDimension( m_Type.GetDimension());
          m_FeatureColumnsPerObject += ( *itr)->GetSizeOfFeatures();
        }
        Base< t_DataType, t_ReturnType>::SetDimension( m_Type.GetDimension());
      }
      else
      {
        // add the number of features per small molecule to the running total
        m_Properties.LastElement()->SetDimension( m_Type.GetDimension());
        m_FeatureColumnsPerObject += m_Properties.LastElement()->GetSizeOfFeatures();
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetDimension is called
    template< typename t_DataType, typename t_ReturnType>
    void Combine< t_DataType, t_ReturnType>::SetDimensionHook()
    {
      m_FeatureColumnsPerObject = 0;
      for( iterator itr( m_Properties.Begin()), itr_end( m_Properties.End()); itr != itr_end; ++itr)
      {
        if( itr->IsDefined())
        {
          m_FeatureColumnsPerObject += ( *itr)->GetSizeOfFeatures();
        }
      }
      m_Type = Type( this->GetDimensionSetting(), m_Type.ConsiderRepeatedObjects(), m_Type.GetSymmetry());
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType, typename t_ReturnType>
    bool Combine< t_DataType, t_ReturnType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_FeatureColumnsPerObject = 0;
      m_Type = Type();
      m_CachePreference = e_PreferCache;

      // determine the final dimension
      for( const_iterator itr( m_Properties.Begin()), itr_end( m_Properties.End()); itr != itr_end; ++itr)
      {
        m_Type.GeneralizeToHandle( ( *itr)->GetType());
      }
      Base< t_DataType, t_ReturnType>::SetDimension( m_Type.GetDimension());
      m_FeatureColumnsPerObject = 0;

      for( iterator itr( m_Properties.Begin()), itr_end( m_Properties.End()); itr != itr_end; ++itr)
      {
        if( itr->IsDefined())
        {
          m_FeatureColumnsPerObject += ( *itr)->GetSizeOfFeatures();
          m_CachePreference = m_CachePreference == e_NeverCache || ( *itr)->GetCachePreference() == e_NeverCache
                              ? e_NeverCache
                              : m_CachePreference;
        }
      }
      return true;
    }

    //! @brief Get the code / label set for the combine with sizes of each property
    //! @return the code / label for the combine with sizes of each property
    //! the feature code set
    template< typename t_DataType, typename t_ReturnType>
    model::FeatureLabelSet Combine< t_DataType, t_ReturnType>::GetLabelsWithSizes() const
    {
      storage::Vector< size_t> sizes;

      // allocate enough extra memory to hold all the arguments
      sizes.AllocateMemory( m_Properties.GetSize());

      // Iterate through all the internal properties.
      for( const_iterator itr( m_Properties.Begin()), itr_end( m_Properties.End()); itr != itr_end; ++itr)
      {
        // append size for that PropertiesInterface
        sizes.PushBack( ( *itr)->GetSizeOfFeatures());
      }

      return model::FeatureLabelSet
             (
               GetAlias(),
               Base< t_DataType, t_ReturnType>::GetLabel().GetArguments(),
               sizes,
               util::Implementation< util::ImplementationInterface>
               (
                 util::Implementation< Base< t_DataType, t_ReturnType> >()
               )
             );
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer Combine< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Array of descriptors");
      parameters.AddInitializer( "", "descriptors to concatenate", io::Serialization::GetAgent( &m_Properties));

      return parameters;
    } // GetParameters

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType, typename t_ReturnType>
    iterate::Generic< Base< t_DataType, t_ReturnType> > Combine< t_DataType, t_ReturnType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, t_ReturnType> >( m_Properties.Begin(), m_Properties.End());
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType, typename t_ReturnType>
    void Combine< t_DataType, t_ReturnType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< t_ReturnType> &STORAGE
    )
    {
      size_t offset( 0);
      for( iterator itr( m_Properties.Begin()), itr_end( m_Properties.End()); itr != itr_end; ++itr)
      {
        const size_t n_features( ( *itr)->GetSizeOfFeatures());
        const linal::VectorConstReference< t_ReturnType> result( ( **itr)( ITR));
        STORAGE.CreateSubVectorReference( n_features, offset).CopyValues( result);
        if( !result.IsDefined())
        {
          BCL_MessageVrb
          (
            "Undefined descriptor: " + ( *itr)->GetLabel().ToString()
            + " at position: " + util::Format()( ITR.GetPosition())
          );
        }
        offset += n_features;
      }
    }

  } // namespace descriptor
} // namespace bcl
