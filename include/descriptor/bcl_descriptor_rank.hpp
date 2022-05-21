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
#include "bcl_descriptor_rank.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_combine.h"
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, accepts bool of whether to auto-compute mean
    template< typename t_DataType>
    Rank< t_DataType>::Rank( const bool &ASCENDING) :
      m_Ascending( ASCENDING)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    Rank< t_DataType> *Rank< t_DataType>::Clone() const
    {
      return new Rank( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &Rank< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &Rank< t_DataType>::GetAlias() const
    {
      static const std::string s_NameAsc( "RankAsc"), s_NameDesc( "RankDesc");

      return m_Ascending ? s_NameAsc : s_NameDesc;
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    template< typename t_DataType>
    size_t Rank< t_DataType>::GetNormalDimension() const
    {
      return m_Property->GetType().GetDimension();
    }

    //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t Rank< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_Property->GetNormalSizeOfFeatures();
    }

    //! @brief set properties to be logarithmized
    //! @param PROPERTY SmallMoleculeProperty
    template< typename t_DataType>
    void Rank< t_DataType>::SetProperty( const util::Implementation< Base< t_DataType, float> > &PROPERTY)
    {
      m_Property = PROPERTY;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer Rank< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "rank property values in " +
        std::string
        (
          m_Ascending
          ? "ascending order. e.g. 0.2 0.21 0.5 0.1 -> 1 2 3 0"
          : "descending order. e.g. 0.2 0.21 0.5 0.1 -> 2 1 0 3"
        )
        + "Duplicate values receive the same rank (averaged) assigned, so 1 1 4 -> "
        + std::string( m_Ascending ? "0.5 0.5 2" : "1.5 1.5 0")
      );

      parameters.AddInitializer
      (
        "",
        "descriptor to rank",
        io::Serialization::GetAgent( &m_Property)
      );
      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > Rank< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Property, &m_Property + 1);
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    template< typename t_DataType>
    Type::Symmetry Rank< t_DataType>::GetSymmetry() const
    {
      return m_Property->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    template< typename t_DataType>
    bool Rank< t_DataType>::ConsiderRepeatedElements() const
    {
      return m_Property->GetType().ConsiderRepeatedObjects();
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void Rank< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // instantiate the property as a vector with indices that correspond to atoms
      linal::VectorConstReference< float> ref( m_Property->operator ()( ITR));

      // maps containing the values of containers as key and occurences as value
      // the automatic sorting of a map takes care of the ranking order of elements
      storage::Map< float, float> ranks_map;

      // add values of container and update their respective count
      for( const float *itr_ref( ref.Begin()), *itr_ref_end( ref.End()); itr_ref != itr_ref_end; ++itr_ref)
      {
        ranks_map[ *itr_ref] += 1.0;
      }

      // counter for ranks of container
      float rank( 0);

      // iterate over all values in rank map and determine the respective rank
      // elements with same rank share the average rank
      if( m_Ascending)
      {
        for
        (
          storage::Map< float, float>::iterator
            itr_map( ranks_map.Begin()), itr_map_end( ranks_map.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          // count for key in rank map
          float count( itr_map->second);
          // adjust rank in map and compute average rank
          itr_map->second = rank + ( count + 1.0) * 0.5;
          // increment rank
          rank += count;
        }
      }
      else
      {
        for
        (
          storage::Map< float, float>::reverse_iterator
            itr_map( ranks_map.ReverseBegin()), itr_map_end( ranks_map.ReverseEnd());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          // count for key in rank map
          float count( itr_map->second);
          // adjust rank in map and compute average rank
          itr_map->second = rank + ( count + 1.0) * 0.5;
          // increment rank
          rank += count;
        }
      }

      // fill final ranking containers
      float *itr_storage( STORAGE.Begin());
      for
      (
        const float *itr_ref( ref.Begin()), *itr_ref_end( ref.End());
        itr_ref != itr_ref_end;
        ++itr_ref, ++itr_storage
      )
      {
        *itr_storage = ranks_map[ *itr_ref];
      }
    }

  } // namespace descriptor
} // namespace bcl
