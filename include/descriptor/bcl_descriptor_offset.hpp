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
#include "bcl_descriptor_offset.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_power_spectrum.h"
#include "bcl_descriptor_window.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType>
    Offset< t_DataType>::Offset( const bool &REFLECTING, const int &OFFSET) :
      m_Descriptor(),
      m_Offset( OFFSET),
      m_Reflecting( REFLECTING)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Offset
    template< typename t_DataType>
    Offset< t_DataType> *Offset< t_DataType>::Clone() const
    {
      return new Offset( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Offset< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &Offset< t_DataType>::GetAlias() const
    {
      static const std::string s_reflecting_name( "ReflectingOffset"), s_name( "Offset");
      return m_Reflecting ? s_reflecting_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t Offset< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetNormalSizeOfFeatures();
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer Offset< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes the descriptor at a specified offset. If the end of the sequence is reached before the offset, the "
        + std::string( m_Reflecting ? " the iterator reflects and continues " : "last value ") + " is returned"
      );
      parameters.AddInitializer
      (
        "",
        "descriptor of interest",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      parameters.AddInitializer
      (
        "offset",
        "Desired # of " + GetStaticClassName< t_DataType>() + "away to retrieve the descriptor for; sign implies directionality",
        io::Serialization::GetAgent( &m_Offset)
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void Offset< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // copy the iterator
      iterate::Generic< const t_DataType> itr_seq( ELEMENT);

      // iterate until the offset is reached
      int local_offset( m_Offset);
      while( local_offset)
      {
        while( local_offset > 0 && itr_seq.NotAtEnd())
        {
          ++itr_seq;
          --local_offset;
        }
        if( !itr_seq.NotAtEnd())
        {
          ++local_offset;
          --itr_seq;
          if( !m_Reflecting)
          {
            // stop for non-reflecting offsets
            break;
          }
          local_offset = -local_offset;
        }
        while( local_offset < 0 && itr_seq.GetPosition())
        {
          ++local_offset;
          --itr_seq;
        }
        if( !itr_seq.GetPosition() && local_offset)
        {
          local_offset = -local_offset;
          if( !m_Reflecting)
          {
            // stop for non-reflecting offsets
            break;
          }
        }
      }

      // create a new iterator at this position
      Iterator< t_DataType> itr( itr_seq);

      STORAGE.CopyValues( m_Descriptor->operator()( itr));
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > Offset< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool Offset< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( !m_Descriptor.IsDefined() || m_Descriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Offset only works with element-wise descriptors that return a single value";
        return false;
      }

      return true;
    }

  } // namespace descriptor
} // namespace bcl
