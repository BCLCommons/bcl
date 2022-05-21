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

#ifndef BCL_DESCRIPTOR_BASE_HPP_
#define BCL_DESCRIPTOR_BASE_HPP_

// include the header of this class
#include "bcl_descriptor_base.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_cache_map.h"
#include "bcl_descriptor_iterator.h"
#include "bcl_descriptor_sequence_interface.h"
#include "io/bcl_io_fixed_line_width_writer.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "type/bcl_type_compare.h"
#include "type/bcl_type_enable_if.h"
#include "type/bcl_type_is_a.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

// Uncomment the next line to verify that the cache is being used properly; it validates every retrieval from the
// cache by recomputing the descriptor, and so can be used to locate improper cache sharing
//#define VALIDATE_DESCRIPTOR_CACHING

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType, typename t_ReturnType>
    Base< t_DataType, t_ReturnType>::Base() :
      m_DimensionSetting( util::GetUndefined< size_t>())
    {
    }

    //! @brief constructor from dimension
    template< typename t_DataType, typename t_ReturnType>
    Base< t_DataType, t_ReturnType>::Base( const size_t &DIMENSION) :
      m_DimensionSetting( util::GetUndefined< size_t>())
    {
      SetDimension( DIMENSION);
    }

    //! @brief copy constructor
    template< typename t_DataType, typename t_ReturnType>
    Base< t_DataType, t_ReturnType>::Base( const Base< t_DataType, t_ReturnType> &BASE) :
      m_DimensionSetting( BASE.m_DimensionSetting),
      m_Storage( BASE.m_Storage),
      m_SequencePtr( BASE.m_SequencePtr),
      m_CachedDescriptor( BASE.m_CachedDescriptor),
      m_CacheLabel( BASE.m_CacheLabel)
    {
      if( m_Storage.GetSize() && BASE.m_CachedDescriptor.Begin() == BASE.m_Storage.Begin())
      {
        // in some cases, the cached descriptor may refer to the storage
        m_CachedDescriptor =
          linal::MatrixConstReference< t_ReturnType>
          (
            m_CachedDescriptor.GetNumberRows(),
            m_CachedDescriptor.GetNumberCols(),
            m_Storage.Begin()
          );
      }
    }

    //! @brief get the type of this descriptor
    //! @return the type of this descriptor (should ignore dimension setting)
    template< typename t_DataType, typename t_ReturnType>
    Type Base< t_DataType, t_ReturnType>::GetType() const
    {
      return Type( GetNormalDimension(), ConsiderRepeatedElements(), GetSymmetry());
    }

    //! @brief get the type of this descriptor
    //! @return the type of this descriptor (should ignore dimension setting)
    template< typename t_DataType, typename t_ReturnType>
    Type Base< t_DataType, t_ReturnType>::GetEffectiveType() const
    {
      if( !util::IsDefined( m_DimensionSetting))
      {
        return GetType();
      }
      return Type( m_DimensionSetting, ConsiderRepeatedElements(), GetSymmetry());
    }

    //! @brief get the cache preference under the current dimension setting (e.g. m_DimensionSetting)
    //! @return the cache preference for the descriptor
    template< typename t_DataType, typename t_ReturnType>
    CachePreference Base< t_DataType, t_ReturnType>::GetCachePreference() const
    {
      // if no dimension setting was given, or the native dimension or higher dimension was given,
      // just return the normal cache preference
      return m_DimensionSetting >= GetNormalDimension() ? GetNormalCachePreference() : e_PreferCache;
    }

  ////////////////
  // operators //
  ////////////////

    class BCL_API DescriptorHelper
    {
    public:

      //! @brief helper function to cache something if the CacheMap contains linal::Vector< t_Returntype>
      template< typename t_DataType, typename t_ReturnType>
      static typename type::EnableIf
      <
        type::Compare< linal::Vector< t_ReturnType>, CacheMap::value_type>::e_Same,
        linal::MatrixConstReference< t_ReturnType>
      >::Type FromCache
      (
        Base< t_DataType, t_ReturnType> &BASE,
        const Type &TYPE
      )
      {
        const size_t native_dimension( BASE.GetNormalDimension());
        const size_t feature_size( BASE.GetSizeOfFeatures());
        const SequenceInterface< t_DataType> &sequence( *BASE.m_SequencePtr);
        const size_t dimension_setting( BASE.m_DimensionSetting);

        // look for this descriptor within the cache
        const util::ObjectDataLabel &cache_label( BASE.GetCacheLabel());
        BCL_MessageVrb
        (
          "FromCache: NativeDimension " + util::Format()( native_dimension)
          + " Dim setting: " + util::Format()( dimension_setting)
          + " cache label: " + cache_label.ToString()
          + " feat size: " + util::Format()( feature_size)
          + " normal feat size: " + util::Format()( BASE.GetNormalSizeOfFeatures())
        );
        util::SiPtr< const linal::Vector< t_ReturnType> > cache_ptr( sequence.FindInCache( cache_label));
        if( !cache_ptr.IsDefined())
        {
          if( !dimension_setting && native_dimension)
          {
            if( cache_label.GetNumberArguments() == size_t( 1))
            {
              // it is possible that the values are cached for higher dimensions, in which case a summation can be performed
              cache_ptr = sequence.FindInCache( *cache_label.Begin());
            }
            if( cache_ptr.IsDefined())
            {
              const size_t n_features( TYPE.GetNumberFeatures( sequence.GetSize()));
              if( cache_ptr->GetSize() != feature_size * n_features)
              {
                BCL_MessageVrb
                (
                  "Found: " + BASE.m_CacheLabel.ToString()
                  + " in cache with wrong # features, uncaching"
                );
                // bad feature size, remove it from the cache
                sequence.Uncache( *cache_label.Begin());
              }
              else
              {
                // higher-dimensional result found, sum it up and add it to the cache
                linal::Vector< t_ReturnType> sum_result( feature_size, cache_ptr->Begin());
                const t_ReturnType *cached_val( cache_ptr->Begin() + feature_size);
                t_ReturnType *storage_row( sum_result.Begin());
                for( size_t row_num( 1); row_num < n_features; ++row_num)
                {
                  for( size_t i( 0); i < feature_size; ++i, ++cached_val)
                  {
                    storage_row[ i] += *cached_val;
                  }
                }
                return
                  linal::MatrixConstReference< t_ReturnType>
                  (
                    size_t( 1),
                    feature_size,
                    sequence.Cache( cache_label, sum_result).Begin()
                  );
              }
            }
          }
        }
        else // if( cache_ptr.IsDefined())
        {
          // descriptor was found
          if( !dimension_setting || !native_dimension)
          {
            if( cache_ptr->GetSize() != feature_size)
            {
              BCL_MessageVrb
              (
                "Found: " + BASE.m_CacheLabel.ToString()
                + " in cache with wrong # features, uncaching"
              );
              // bad feature size, remove it from the cache
              sequence.Uncache( cache_label);
            }
            else
            {
              return linal::MatrixConstReference< t_ReturnType>( size_t( 1), feature_size, cache_ptr->Begin());
            }
          }
          else if( dimension_setting == native_dimension)
          {
            const size_t n_features( TYPE.GetNumberFeatures( sequence.GetSize()));
            if( cache_ptr->GetSize() != n_features * feature_size)
            {
              BCL_MessageVrb
              (
                "Found: " + BASE.m_CacheLabel.ToString()
                + " in cache with wrong # features, uncaching"
              );
              // bad feature size, remove it from the cache
              sequence.Uncache( cache_label);
            }
            else
            {
              return linal::MatrixConstReference< t_ReturnType>( n_features, feature_size, cache_ptr->Begin());
            }
          }
          else if( dimension_setting == size_t( 1))
          {
            const size_t n_features( TYPE.GetNumberFeatures( sequence.GetSize()));
            const size_t native_feature_size( feature_size / dimension_setting);
            if( cache_ptr->GetSize() != n_features * native_feature_size)
            {
              BCL_MessageVrb
              (
                "Found: " + BASE.m_CacheLabel.ToString()
                + " in cache with wrong # features, uncaching"
              );
              // bad feature size, remove it from the cache
              sequence.Uncache( cache_label);
            }
            else
            {
              return linal::MatrixConstReference< t_ReturnType>( n_features, native_feature_size, cache_ptr->Begin());
            }
          }
        }

        // descriptor was not found, calculate it and cache it
        const size_t native_feature_size( BASE.GetNormalSizeOfFeatures());

        BCL_MessageDbg
        (
          BASE.m_CacheLabel.ToString()
          + " was not in cache, recalculating; native feature size: " + util::Format()( native_feature_size)
        );

        // create a native iterator
        Iterator< t_DataType> itr( TYPE, sequence);
        if( !native_dimension)
        {
          linal::Vector< t_ReturnType> storage( native_feature_size, t_ReturnType( 0));
          linal::VectorReference< t_ReturnType> storage_ref( storage);
          // simple sequence descriptor
          BASE.RecalculateImpl( itr, storage_ref);

          BCL_MessageDbg
          (
            "Caching scalar descriptor: " + BASE.m_CacheLabel.ToString()
            + " with # values = " + util::Format()( storage.GetSize())
          );
          return
            linal::MatrixConstReference< t_ReturnType>
            (
              size_t( 1),
              feature_size,
              sequence.Cache( cache_label, storage).Begin()
            );
        }
        else if( !dimension_setting)
        {
          linal::Vector< t_ReturnType> storage( native_feature_size, t_ReturnType( 0));
          linal::Vector< t_ReturnType> local_storage( storage);
          // sum the results, but do not store the lower dimensional results
          for( linal::VectorReference< t_ReturnType> local_storage_ref( local_storage); itr.NotAtEnd(); ++itr)
          {
            local_storage = t_ReturnType( 0);
            BASE.RecalculateImpl( itr, local_storage_ref);
            storage += local_storage;
          }
          BCL_MessageDbg
          (
            "Caching summed scalar descriptor: " + BASE.m_CacheLabel.ToString()
            + " with # values = " + util::Format()( storage.GetSize())
          );
          return
            linal::MatrixConstReference< t_ReturnType>
            (
              size_t( 1),
              feature_size,
              sequence.Cache( cache_label, storage).Begin()
            );
        }
        else if( native_dimension == dimension_setting || native_dimension == size_t( 1))
        {
          // cache all the results at the native dimension
          const size_t n_features( itr.GetSize());
          const size_t native_feature_size( BASE.GetNormalSizeOfFeatures());
          BCL_MessageDbg
          (
            "Caching descriptor: " + BASE.m_CacheLabel.ToString()
            + " with native_feature_size = " + util::Format()( native_feature_size)
            + " with dimension setting: " + util::Format()( dimension_setting)
            + " with total # features: " + util::Format()( n_features)
          );
          // calculate all the values for this descriptor and throw them in the cache
          linal::Vector< t_ReturnType> full_storage( n_features * native_feature_size, t_ReturnType( 0));
          t_ReturnType *itr_storage( full_storage.Begin());
          for( ; itr.NotAtEnd(); ++itr, itr_storage += native_feature_size)
          {
            // calculate the descriptor for this iterator
            linal::VectorReference< t_ReturnType> this_row( native_feature_size, itr_storage);

            // calculate the descriptor
            BASE.RecalculateImpl( itr, this_row);
          }
          return
            linal::MatrixConstReference< t_ReturnType>
            (
              n_features,
              native_feature_size,
              sequence.Cache( cache_label, full_storage).Begin()
            );
        }

        return linal::MatrixConstReference< t_ReturnType>();
      }

      template< typename t_DataType, typename t_ReturnType>
      static typename type::EnableIf
      <
        type::Compare< linal::Vector< t_ReturnType>, CacheMap::value_type>::e_Different,
        linal::MatrixConstReference< t_ReturnType>
      >::Type FromCache( Base< t_DataType, t_ReturnType> &BASE, const Type &TYPE)
      {
        return linal::MatrixConstReference< t_ReturnType>();
      }
    };

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest
    //! @return generated description for given argument
    template< typename t_DataType, typename t_ReturnType>
    linal::VectorConstReference< t_ReturnType> Base< t_DataType, t_ReturnType>::Recalculate
    (
      const Iterator< t_DataType> &ITR
    )
    {
      // check dimensions
      if( ITR.GetType().GetDimension() != m_DimensionSetting && GetType().GetDimension() != size_t( 0))
      {
        SetDimension( ITR.GetType().GetDimension());
      }

      const size_t native_dimension( GetNormalDimension());

      m_Storage = t_ReturnType( 0);

      linal::VectorReference< t_ReturnType> storage_ref( m_Storage);
      // create a native iterator
      if( native_dimension == m_DimensionSetting)
      {
        // simple sequence descriptor
        RecalculateImpl( ITR, storage_ref);
      }
      else if( !native_dimension)
      {
        const Iterator< t_DataType> itr_base( GetType(), *m_SequencePtr);
        RecalculateImpl( itr_base, storage_ref);
      }
      else if( !m_DimensionSetting)
      {
        // sum the results, but do not store the lower dimensional results
        linal::Vector< t_ReturnType> local_storage( m_Storage);
        linal::VectorReference< t_ReturnType> local_storage_ref( local_storage);
        for( Iterator< t_DataType> itr_sum( GetType(), *m_SequencePtr); itr_sum.NotAtEnd(); ++itr_sum)
        {
          local_storage = t_ReturnType( 0);
          RecalculateImpl( itr_sum, local_storage_ref);
          m_Storage += local_storage;
        }
      }
      else if( native_dimension == size_t( 1))
      {
        // cache all the results at the native dimension
        const size_t native_feature_size( GetNormalSizeOfFeatures());

        // calculate all the values for this descriptor and throw them in the cache
        Iterator< t_DataType> itr_element( GetType(), *m_SequencePtr);
        linal::VectorReference< t_ReturnType> storage_ref( m_Storage);
        size_t sub_size( 0);
        for( size_t dimension( 0); dimension < m_DimensionSetting; ++dimension, sub_size += native_feature_size)
        {
          itr_element.GotoPosition( ITR( dimension).GetPosition());
          // calculate the descriptor for this iterator
          linal::VectorReference< t_ReturnType> this_row
          (
            storage_ref.CreateSubVectorReference( native_feature_size, sub_size)
          );

          // test that the descriptor was calculated successfully
          RecalculateImpl( itr_element, this_row);
        }
      }
      else
      {
        BCL_Exit
        (
          "Changing from dimension " + util::Format()( native_dimension)
          + " to " + util::Format()( m_DimensionSetting) + " not yet implemented",
          -1
        );
      }

      return linal::VectorConstReference< t_ReturnType>( m_Storage);
    }

    //! @brief operator to calculate the descriptors for a given sequence iterator position, considering cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest
    //! @return generated description for given argument
    template< typename t_DataType, typename t_ReturnType>
    linal::VectorConstReference< t_ReturnType> Base< t_DataType, t_ReturnType>::operator()( const Iterator< t_DataType> &ITR)
    {
      BCL_Assert( m_SequencePtr.IsDefined(), GetClassIdentifier() + " operator() called before SetObject");

      // check dimensions
      if( ITR.GetType().GetDimension() != m_DimensionSetting && GetType().GetDimension() != size_t( 0))
      {
        SetDimension( ITR.GetType().GetDimension());
      }
      if( m_CachedDescriptor.GetNumberRows() == size_t( 0) && GetCachePreference() == e_PreferCache)
      {
        Type this_type
        (
          GetNormalDimension(),
          ITR.GetType().ConsiderRepeatedObjects(),
          ITR.GetType().GetSymmetry()
        );
        m_CachedDescriptor = DescriptorHelper::FromCache( *this, this_type);
      }

      // try to retrieve it from the cache, and if this descriptor should be cached, calculate and cache it
      const size_t n_cached_rows( m_CachedDescriptor.GetNumberRows());
      if( n_cached_rows)
      {
        const size_t n_cached_cols( m_CachedDescriptor.GetNumberCols());

        // this should never happen, since we check column sizes earlier, but better to check
        BCL_Assert
        (
          n_cached_cols == GetNormalSizeOfFeatures(),
          "Wrong number columns in cached output; should have been "
          + util::Format()( GetNormalSizeOfFeatures())
          + " but was " + util::Format()( n_cached_cols) + " for " + this->GetString()
        );

        if( GetType().GetDimension() == size_t( 0) || m_DimensionSetting == size_t( 0))
        {
          BCL_Assert
          (
            n_cached_rows == 1,
            "Wrong number rows (" + util::Format()( n_cached_rows) + ", should be 1) for "
            + this->GetString() + " " + m_CacheLabel.ToString()
          );
#ifdef VALIDATE_DESCRIPTOR_CACHING
          if
          (
            !math::EqualWithinTolerance
            (
              linal::VectorConstReference< t_ReturnType>( n_cached_cols, m_CachedDescriptor[ 0]),
              Recalculate( ITR),
              1.0e-4,
              1.0e-4
            )
          )
          {
            BCL_MessageCrt
            (
              "Cached result was wrong for " + m_CacheLabel.ToString()
              + " should have been: " + util::Format()( Recalculate( ITR))
              + " but was " + util::Format()( linal::VectorConstReference< t_ReturnType>( n_cached_cols, m_CachedDescriptor[ 0]))
            );
          }
#endif
          return linal::VectorConstReference< t_ReturnType>( n_cached_cols, m_CachedDescriptor[ 0]);
        }
        else if( GetSizeOfFeatures() == GetNormalSizeOfFeatures())
        {
          BCL_Assert
          (
            ITR.GetSize() == n_cached_rows,
            "Wrong number rows, had " + util::Format()( n_cached_rows) + " should be " + util::Format()( ITR.GetSize())
            + " for " + this->GetString()
          );
#ifdef VALIDATE_DESCRIPTOR_CACHING
          if
          (
            !math::EqualWithinTolerance
            (
              linal::VectorConstReference< t_ReturnType>( n_cached_cols, m_CachedDescriptor[ ITR.GetPosition()]),
              Recalculate( ITR),
              1.0e-4,
              1.0e-4
            )
          )
          {
            BCL_MessageCrt
            (
              "Cached result was wrong for " + m_CacheLabel.ToString()
              + " should have been: " + util::Format()( Recalculate( ITR))
              + " but was " + util::Format()( linal::VectorConstReference< t_ReturnType>( n_cached_cols, m_CachedDescriptor[ ITR.GetPosition()]))
            );
          }
#endif
          return linal::VectorConstReference< t_ReturnType>( n_cached_cols, m_CachedDescriptor[ ITR.GetPosition()]);
        }
        else
        {
          BCL_Assert( GetNormalDimension() == size_t( 1), "Should have been an element-wise descriptor");
          t_ReturnType *storage_row( m_Storage.Begin());
          for( size_t i( 0); i < m_DimensionSetting; ++i)
          {
            const t_ReturnType *cached_row( m_CachedDescriptor[ ITR( i).GetPosition()]);
            storage_row = std::copy( cached_row, cached_row + n_cached_cols, storage_row);
          }
          #ifdef VALIDATE_DESCRIPTOR_CACHING
          if
          (
            !math::EqualWithinTolerance
            (
              linal::VectorConstReference< t_ReturnType>( n_cached_cols * m_DimensionSetting, m_Storage.Begin()),
              Recalculate( ITR),
              1.0e-4,
              1.0e-4
            )
          )
          {
            BCL_MessageCrt
            (
              "Cached result was wrong for " + m_CacheLabel.ToString()
              + " should have been: " + util::Format()( Recalculate( ITR))
              + " but was " + util::Format()( linal::VectorConstReference< t_ReturnType>( n_cached_cols * m_DimensionSetting, m_Storage.Begin()))
            );
          }
#endif
          return linal::VectorConstReference< t_ReturnType>( n_cached_cols * m_DimensionSetting, m_Storage.Begin());
        }
      }

      // descriptor was not available, just calculate it
      return Recalculate( ITR);
    }

    //! @brief get the logical name for the complete sequence; evaluates to e.g. sequence, molecule, or string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Base< t_DataType, t_ReturnType>::GetObjectName()
    {
      static const std::string s_name
      (
        type::IsA< t_DataType, chemistry::AtomConformationalInterface>::value
        ? "Molecule"
        : type::IsA< t_DataType, char>::value
          ? "String"
          : type::IsA< t_DataType, biol::AABase>::value
            ? "Sequence"
            : "MutationSet"
      );
      return s_name;
    }

    //! @brief get the logical name for each element of the sequence; evaluates to e.g. amino acid, atom, or character
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Base< t_DataType, t_ReturnType>::GetElementName()
    {
      static const std::string s_name
      (
        type::IsA< t_DataType, chemistry::AtomConformationalInterface>::value
        ? "Atom"
        : type::IsA< t_DataType, char>::value
          ? "Character"
          : type::IsA< t_DataType, biol::AABase>::value
            ? "Amino Acid"
            : "Mutation"
      );
      return s_name;
    }

    //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be appended
    //! @param SEQ sequence of interest
    //! @return generated description for given argument
    //! returns GetNormalSizeOfFeatures() * SEQ.GetSize() values (GetNormalSizeOfFeatures for each element)
    template< typename t_DataType, typename t_ReturnType>
    linal::Vector< t_ReturnType> Base< t_DataType, t_ReturnType>::CollectValuesOnEachElementOfObject
    (
      const SequenceInterface< t_DataType> &SEQ
    )
    {
      SetObject( SEQ);
      Type type( 1, ConsiderRepeatedElements(), GetSymmetry());
      this->SetDimension( 1);
      Iterator< t_DataType> itr( type, SEQ);
      linal::Vector< t_ReturnType> return_vector( GetNormalSizeOfFeatures() * itr.GetSize());
      t_ReturnType *itr_return( return_vector.Begin());

      for( ; itr.NotAtEnd(); ++itr)
      {
        const linal::VectorConstReference< t_ReturnType> &result_this_element( this->operator()( itr));
        itr_return = std::copy( result_this_element.Begin(), result_this_element.End(), itr_return);
      }
      return return_vector;
    }

    //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be appended
    //! @param SEQ sequence of interest
    //! @return generated description for given argument
    //! returns GetNormalSizeOfFeatures() * SEQ.GetSize() values (GetNormalSizeOfFeatures for each element)
    template< typename t_DataType, typename t_ReturnType>
    linal::Vector< t_ReturnType> Base< t_DataType, t_ReturnType>::CollectValuesOnEachElementOfObject
    (
      const SequenceInterface< t_DataType> &SEQ
    ) const
    {
      util::ShPtr< Base< t_DataType, t_ReturnType> > clone( this->Clone());
      return clone->CollectValuesOnEachElementOfObject( SEQ);
    }

    //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be summed
    //! @param SEQ sequence of interest
    //! @return generated description for given argument
    //! returns GetNormalSizeOfFeatures(); if the descriptor is multi-dimensional; results are summed
    template< typename t_DataType, typename t_ReturnType>
    linal::Vector< t_ReturnType>
      Base< t_DataType, t_ReturnType>::SumOverObject( const SequenceInterface< t_DataType> &SEQ)
    {
      this->SetDimension( 0);
      SetObject( SEQ);
      Iterator< t_DataType> itr( Type(), SEQ);
      return linal::Vector< t_ReturnType>( this->operator()( itr));
    }

    //! @brief operator to calculate the descriptors for a given sequence; all inner descriptors will be summed
    //! @param SEQ sequence of interest
    //! @return generated description for given argument
    //! returns GetNormalSizeOfFeatures(); if the descriptor is multi-dimensional; results are summed
    template< typename t_DataType, typename t_ReturnType>
    linal::Vector< t_ReturnType>
      Base< t_DataType, t_ReturnType>::SumOverObject( const SequenceInterface< t_DataType> &SEQ) const
    {
      util::ShPtr< Base< t_DataType, t_ReturnType> > clone( this->Clone());
      clone->SetDimension( 0);
      clone->SetObject( SEQ);
      Iterator< t_DataType> itr( Type(), SEQ);
      return linal::Vector< t_ReturnType>( clone->operator()( itr));
    }

    //! @brief helper function to write help for all instances of this class
    //! @param STREAM stream to write the help to
    //! @param INSTANCES all instances of the class
    template< typename t_DataType, typename t_ReturnType>
    io::FixedLineWidthWriter &Base< t_DataType, t_ReturnType>::WriteInstancesHelp
    (
      io::FixedLineWidthWriter &STREAM,
      const storage::Map< std::string, util::OwnPtr< Base> > &INSTANCES,
      bool FULL_HELP
    )
    {
      // write a short description
      STREAM << "choose any " << ' ' << GetObjectName() << " / " << GetElementName()
             << ' ' << ( type::Compare< t_ReturnType, float>::e_Same ? "Numeric" : "String") << " descriptor :";
      STREAM.AddIndent( 2);

      storage::Vector< storage::Vector< storage::Vector< util::SiPtr< const Base> > > >
        interfaces_by_dimension_type
        (
          size_t( 2),
          storage::Vector< storage::Vector< util::SiPtr< const Base> > >( size_t( util::FunctionalType::s_NumberTypes))
        );

      // write the default data label of each instance out to STREAM
      for
      (
        typename storage::Map< std::string, util::OwnPtr< Base> >::const_iterator
          itr( INSTANCES.Begin()),
          itr_end( INSTANCES.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->first == "")
        {
          continue;
        }
        size_t type
        (
          itr->second->InferFunctionalType( GetStaticClassName< util::Implementation< Base> >())
        );
        const bool defined_dimension( itr->second->DimensionIsWellDefined());
        if( !defined_dimension)
        {
          if( itr->second->DimensionIsContextual())
          {
            interfaces_by_dimension_type( interfaces_by_dimension_type.GetSize() - 2)( type).PushBack( *itr->second);
          }
          else
          {
            interfaces_by_dimension_type.LastElement()( type).PushBack( *itr->second);
          }
          continue;
        }
        const size_t normal_dimension( itr->second->GetNormalDimension());
        if( interfaces_by_dimension_type.GetSize() < normal_dimension + 3)
        {
          interfaces_by_dimension_type.InsertElements
          (
            interfaces_by_dimension_type.GetSize() - 2,
            storage::Vector< storage::Vector< util::SiPtr< const Base> > >( size_t( util::FunctionalType::s_NumberTypes)),
            normal_dimension + 3 - interfaces_by_dimension_type.GetSize()
          );
        }
        interfaces_by_dimension_type( normal_dimension)( type).PushBack( *itr->second);
      }
      static const size_t n_explicit_tuple_names( 4);
      static const std::string s_tuple_names[ n_explicit_tuple_names] =
      {
        GetObjectName(),
        GetElementName(),
        GetElementName() + " pair",
        GetElementName() + " triplet"
      };
      for( size_t dimension( 0), max_dimension( interfaces_by_dimension_type.GetSize()); dimension < max_dimension; ++dimension)
      {
        // check whether there are any descriptors with this dimension
        bool found_one( false);
        const storage::Vector< storage::Vector< util::SiPtr< const Base> > >
          &interfaces_of_dimension( interfaces_by_dimension_type( dimension));
        for( size_t type_num( 0); type_num < util::FunctionalType::s_NumberTypes; ++type_num)
        {
          if( interfaces_of_dimension( type_num).GetSize())
          {
            found_one = true;
            break;
          }
        }
        if( !found_one)
        {
          continue;
        }

        STREAM.NewLineIndent();
        STREAM.NewLineIndent();
        if( dimension + 2 < max_dimension)
        {
          STREAM.WriteHeading
          (
            "Descriptors of "
            + (
                dimension < n_explicit_tuple_names
                ? s_tuple_names[ dimension] + "s"
                : GetElementName() + " " + util::Format()( dimension) + "-tuples"
              ),
            '*',
            true
          );
          STREAM.AddIndent( 2);
          if( dimension)
          {
            STREAM.NewLineIndent();
            STREAM.SetAutoIndentExtra( 0);
            STREAM << "These can be converted into a " << s_tuple_names[ 0] << "-level descriptor using "
                   << s_tuple_names[ 0]
                   << "Sum(X), where X is any descriptor listed below";
            STREAM.SetAutoIndentExtra( 2);
          }
        }
        else if( dimension + 2 == max_dimension)
        {
          STREAM.WriteHeading
          (
            "Descriptors that can be natively computed for both "
            + s_tuple_names[ 0] + " and " + s_tuple_names[ 1],
            '*',
            true
          );
          STREAM.AddIndent( 2);
          STREAM.NewLineIndent();
          STREAM.SetAutoIndentExtra( 0);
          STREAM << "By default, these are computed for each " << s_tuple_names[ 1]
                 << ", except when using GenerateDataset with a " << s_tuple_names[ 0] << "-level result descriptor. "
                 << s_tuple_names[ 0] << "-level descriptors can be obtained instead by using " << s_tuple_names[ 0]
                 << "Sum(X), where X is any descriptor listed below";
          STREAM.SetAutoIndentExtra( 2);
        }
        else if( dimension + 1 == max_dimension)
        {
          STREAM.WriteHeading( "General-purpose descriptor operations", '*', true);
          STREAM.AddIndent( 2);
          STREAM.NewLineIndent();
          STREAM.SetAutoIndentExtra( 0);
          STREAM << " These can be used for any type of descriptor ("
                 << s_tuple_names[ 0] << "-level, " << s_tuple_names[ 1]
                 << "-level, etc, or other general purpose descriptors)";
          STREAM.SetAutoIndentExtra( 2);
        }
        STREAM.PopIndent();
        STREAM.NewLineIndent();
        for( size_t type_num( 0); type_num < util::FunctionalType::s_NumberTypes; ++type_num)
        {
          const storage::Vector< util::SiPtr< const Base> > &impls_of_type( interfaces_of_dimension( type_num));
          if( impls_of_type.IsEmpty())
          {
            continue;
          }
          STREAM.NewLineIndent();
          STREAM << util::FunctionalType::GetTypeName( util::FunctionalType::Type( type_num));
          STREAM.AddIndent( 2);
          STREAM.NewLineIndent();
          for
          (
            typename storage::Vector< util::SiPtr< const Base> >::const_iterator
              itr( impls_of_type.Begin()), itr_end( impls_of_type.End());
            itr != itr_end;
            ++itr
          )
          {
            if( !FULL_HELP)
            {
              util::Implementation< Base>( ( *itr)->Clone()).GetCompleteSerializer().WriteBriefHelp( STREAM);
              STREAM.NewLineIndent();
            }
            else
            {
              STREAM.WriteOnOneLine( "* " + ( *itr)->GetAlias() + " : ");
              ( *itr)->WriteHelp( STREAM);
            }
          }
          STREAM.PopIndent();
        }
      }
      STREAM.PopIndent();
      if( util::Enumerated< Base>::GetDefaultImplementation() != INSTANCES.End())
      {
        STREAM << "\nOther strings will be interpreted as follows:\n";
        util::Enumerated< Base>::GetDefaultImplementation()->second->GetCompleteSerializer().WriteBriefHelp( STREAM);
      }
      return STREAM;
    }

    //! @brief set the sequence object and reset the storage
    //! @param SEQUENCE the sequence object of interest
    template< typename t_DataType, typename t_ReturnType>
    void Base< t_DataType, t_ReturnType>::SetObject( const SequenceInterface< t_DataType> &SEQUENCE)
    {
      BCL_MessageDbg( "Set object called on " + GetString());
      if( !util::IsDefined( m_DimensionSetting))
      {
        SetDimension( GetType().GetDimension());
      }
      m_SequencePtr = util::ToSiPtr( SEQUENCE);
      m_CachedDescriptor = linal::MatrixConstReference< t_ReturnType>();
      // call set object on all internally-held descriptors
      for( iterate::Generic< Base< t_DataType, t_ReturnType> > itr( GetInternalDescriptors()); itr.NotAtEnd(); ++itr)
      {
        itr->SetObject( SEQUENCE);
      }

      SetObjectHook();
    }

    //! @brief override the dimension associated with this descriptor
    //! @param NEW_DIMENSION the new dimension to use
    template< typename t_DataType, typename t_ReturnType>
    void Base< t_DataType, t_ReturnType>::SetDimension( const size_t &NEW_DIMENSION)
    {
      // do not try to set an undefined dimension
      BCL_Assert( util::IsDefined( NEW_DIMENSION), "Tried to set undefined dimension!");
      const size_t native_dimension( GetType().GetDimension());
      {
        // for wrapper type descriptors (for which InjectDimensions is true), it is unnecessary to check for
        // whether dimensional changes would be valid
        if( !InjectDimensions())
        {
          BCL_Assert
          (
            !NEW_DIMENSION || native_dimension == NEW_DIMENSION || native_dimension < size_t( 2),
            "Cannot use a descriptor that takes " + util::Format()( native_dimension)
            + this->GetElementName() + "s with " + util::Format()( NEW_DIMENSION) + " " + this->GetElementName() + "s"
          );
        }
        if( native_dimension == size_t( 0))
        {
          // from a sequence-descriptors perspective, the dimension is always zero
          m_DimensionSetting = 0;
        }
        else
        {
          m_DimensionSetting = NEW_DIMENSION;
        }
      }

      // call set object on all internally-held descriptors
      if( InjectDimensions())
      {
        for( iterate::Generic< Base< t_DataType, t_ReturnType> > itr( GetInternalDescriptors()); itr.NotAtEnd(); ++itr)
        {
          itr->SetDimension( NEW_DIMENSION);
        }
      }

      SetDimensionHook();

      size_t desired_feature_size( GetNormalSizeOfFeatures());

      // elementwise descriptors are repeated when used in higher-dimensional contexts
      if( !InjectDimensions() && m_DimensionSetting > size_t( 1) && GetNormalDimension() == size_t( 1))
      {
        desired_feature_size *= m_DimensionSetting;
      }
      BCL_MessageDbg
      (
        this->GetString() + " desired feature size: " + util::Format()( desired_feature_size) +
        " SetDimension called with: " + util::Format()( NEW_DIMENSION)
      );
      if( m_Storage.GetSize() != desired_feature_size)
      {
        m_Storage = linal::Vector< t_ReturnType>( desired_feature_size, t_ReturnType( 0));
      }
      else
      {
        m_Storage = t_ReturnType( 0);
      }

      // update the cache label
      if( GetCachePreference() == e_PreferCache)
      {
        // the cache-label is non-trivial to produce because it requires that the actual object be recognized
        // as a dynamic object (otherwise we do not get this type's alias involved)
        io::Serializer serializer( this->GetSerializer());
        serializer.SetCommandLineIdentifier( this->GetAlias());
        serializer.SetTypeAndFinalize( util::DataType::e_DynamicObject);

        if( !m_DimensionSetting && native_dimension)
        {
          m_CacheLabel = util::ObjectDataLabel( "", this->GetObjectName() + "Sum", serializer.GetLabel());
        }
        else
        {
          m_CacheLabel = serializer.GetLabel();
        }
      }
      else
      {
        m_CacheLabel = util::ObjectDataLabel();
      }

      // a change in dimension would require m_DescriptorCache to be updated, so reset the cache
      m_CachedDescriptor = linal::MatrixConstReference< t_ReturnType>();
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    template< typename t_DataType, typename t_ReturnType>
    CachePreference Base< t_DataType, t_ReturnType>::GetNormalCachePreference() const
    {
      // while we avoid const casts at virtually any cost, the problem is here that GetInternalDescriptors is in a lot
      // of classes, and so it seems better to have this here rather than requiring all other classes to overload a
      // const version of it as well
      auto itr( const_cast< Base< t_DataType, t_ReturnType> *>( this)->GetInternalDescriptors());
      for( ; itr.NotAtEnd(); ++itr)
      {
        if( itr->GetCachePreference() == e_NeverCache)
        {
          return e_NeverCache;
        }
      }
      return itr.GetSize() == size_t( 0) ? e_IgnoreCache : e_PreferCache;
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType, typename t_ReturnType>
    iterate::Generic< Base< t_DataType, t_ReturnType> > Base< t_DataType, t_ReturnType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, t_ReturnType> >();
    }

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_BASE_H_
