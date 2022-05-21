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
#include "bcl_descriptor_window_conditional_average.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor from whether to be springy at the end and whether to do a cumulative average
    template< typename t_DataType>
    WindowConditionalAverage< t_DataType>::WindowConditionalAverage( const bool &AVERAGE) :
      m_MaxSize( 1),
      m_InternalDescriptorSize( 1),
      m_Alignment( e_Center),
      m_Stride( 1),
      m_Average( AVERAGE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new WindowConditionalAverage
    template< typename t_DataType>
    WindowConditionalAverage< t_DataType> *WindowConditionalAverage< t_DataType>::Clone() const
    {
      return new WindowConditionalAverage( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &WindowConditionalAverage< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &WindowConditionalAverage< t_DataType>::GetAlias() const
    {
      static const std::string s_window_name( "WindowConditionalAverage"), s_sum_name( "WindowConditionalSum");
      return m_Average ? s_window_name : s_sum_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t WindowConditionalAverage< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_InternalDescriptorSize;
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
    io::Serializer WindowConditionalAverage< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the average of one descriptor so long as the conditional descriptor returns the same value"
      );
      parameters.AddInitializer
      (
        "size",
        "absolute maximum size of the window",
        io::Serialization::GetAgent( &m_MaxSize)
      );
      parameters.AddInitializer
      (
        "condition",
        "Descriptor that should return the same value so long as the averaging should continue; typically a comparison",
        io::Serialization::GetAgent( &m_ConditionalDescriptor)
      );
      parameters.AddInitializer
      (
        "",
        "actual descriptor to average",
        io::Serialization::GetAgent( &m_AveragedDescriptor)
      );
      parameters.AddInitializer
      (
        "alignment",
        "Direction(s) to allow travel when computing the conditional average",
        io::Serialization::GetAgent( &m_Alignment),
        "Center"
      );
      parameters.AddOptionalInitializer
      (
        "stride",
        "distance between adjacent t_DataTypes that are averaged",
        io::Serialization::GetAgentWithMin( &m_Stride, size_t( 1))
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void WindowConditionalAverage< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // copy the iterator
      Iterator< t_DataType> itr_seq( ELEMENT);

      // determine the conditioned value
      linal::Vector< float> condition_desired_value( m_ConditionalDescriptor->operator()( itr_seq));

      math::RunningAverage< linal::Vector< float> > averager;

      const size_t initial_position( itr_seq.GetPosition());

      // collect values to the left, while the conditional descriptor is unchanged
      if( m_Alignment != e_Right)
      {
        for
        (
          size_t number_averaged_values( 0), position( initial_position);
          number_averaged_values < m_MaxSize;
          ++number_averaged_values, --itr_seq, --position
        )
        {
          const linal::VectorConstReference< float> conditional_reference( m_ConditionalDescriptor->operator()( itr_seq));
          if( conditional_reference != condition_desired_value)
          {
            break;
          }
          if( !( number_averaged_values % m_Stride))
          {
            averager += m_AveragedDescriptor->operator()( itr_seq);
          }
          if( !position)
          {
            break;
          }
        }
      }

      if( m_Alignment != e_Left)
      {
        // collect values to the right, while the conditional descriptor is unchanged
        itr_seq = Iterator< t_DataType>( ELEMENT);
        size_t number_averaged_values( 0);
        if( m_Alignment != e_Right)
        {
          ++itr_seq;
          ++number_averaged_values;
        }
        for
        (
          ;
          number_averaged_values < m_MaxSize && itr_seq.NotAtEnd();
          ++number_averaged_values, ++itr_seq
        )
        {
          const linal::VectorConstReference< float> conditional_reference( m_ConditionalDescriptor->operator()( itr_seq));
          if( conditional_reference != condition_desired_value)
          {
            break;
          }
          if( !( number_averaged_values % m_Stride))
          {
            averager += m_AveragedDescriptor->operator()( itr_seq);
          }
        }
      }
      if( averager.GetWeight())
      {
        STORAGE.CopyValues( averager.GetAverage());
        if( !m_Average)
        {
          STORAGE *= averager.GetWeight();
        }
      }
      else
      {
        STORAGE = float( 0.0);
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > WindowConditionalAverage< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_ConditionalDescriptor, &m_AveragedDescriptor + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool WindowConditionalAverage< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_InternalDescriptorSize = m_AveragedDescriptor->GetSizeOfFeatures();
      if( m_AveragedDescriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Expected an element-wise descriptor for averaged descriptor";
        return false;
      }
      if( m_ConditionalDescriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Expected an element-wise descriptor for conditioned descriptor";
        return false;
      }

      return true;
    }

  } // namespace descriptor
} // namespace bcl
