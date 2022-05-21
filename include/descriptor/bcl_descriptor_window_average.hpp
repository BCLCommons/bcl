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
#include "bcl_descriptor_window_average.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor from whether to be springy at the end and whether to do a cumulative average
    template< typename t_DataType>
    WindowAverage< t_DataType>::WindowAverage( const bool &REFLECTIVE, const bool &CUMULATIVE) :
      m_Window( REFLECTIVE),
      m_Cumulative( CUMULATIVE),
      m_InternalDescriptorSize( 1)
    {
    }

    //! @brief Clone function
    //! @return pointer to new WindowAverage
    template< typename t_DataType>
    WindowAverage< t_DataType> *WindowAverage< t_DataType>::Clone() const
    {
      return new WindowAverage( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &WindowAverage< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &WindowAverage< t_DataType>::GetAlias() const
    {
      static const std::string s_window_name( "WindowAverage"), s_reflective_name( "ReflectingWindowAverage");
      static const std::string s_cma_window_name( "CumulativeWindowAverage"),
                               s_cma_reflective_name( "CumulativeReflectingWindowAverage");
      return
        m_Cumulative
        ? ( m_Window.IsReflective() ? s_cma_reflective_name : s_cma_window_name)
        : ( m_Window.IsReflective() ? s_reflective_name : s_window_name);
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t WindowAverage< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_Cumulative ? m_WindowWeights.GetSize() * m_InternalDescriptorSize : m_InternalDescriptorSize;
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
    io::Serializer WindowAverage< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters( m_Window.GetSerializer());
      parameters.SetClassDescription
      (
        "computes descriptor weighted window averages for " + GetStaticClassName< t_DataType>() +
        std::string
        (
          m_Window.IsReflective()
          ? " reflecting at boundaries"
          : " respecting the boundaries"
        )
      );
      parameters.AddInitializer
      (
        "weighting",
        "method of generating weights for the window",
        io::Serialization::GetAgent( &m_WindowWeightsCreator)
      );

      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void WindowAverage< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      const linal::VectorConstReference< float> &window( m_Window( Iterator< t_DataType>( ELEMENT)));
      float *itr_results( STORAGE.Begin());
      for( size_t feature( 0); feature < m_InternalDescriptorSize; ++feature)
      {
        if( !m_Cumulative)
        {
          math::RunningAverage< float> average;
          for
          (
            const float *itr_win( window.Begin() + feature), *itr_win_end( window.End()), *itr_weight( m_WindowWeights.Begin());
            itr_win < itr_win_end;
            itr_win += m_InternalDescriptorSize, ++itr_weight
          )
          {
            if( util::IsDefined( *itr_win))
            {
              average.AddWeightedObservation( *itr_win, *itr_weight);
            }
          }
          *itr_results = average.GetAverage();
          ++itr_results;
        }
        else
        {
          // reset all bin averagers
          for
          (
            storage::Vector< math::RunningAverage< float> >::iterator
              itr( m_BinAverages.Begin()), itr_end( m_BinAverages.End());
            itr != itr_end;
            ++itr
          )
          {
            itr->Reset();
          }
          // for each position in the window
          storage::Vector< size_t>::const_iterator itr_bin_position( m_Bins.Begin());
          for
          (
            const float
              *itr_win( window.Begin() + feature), *itr_win_end( window.End()), *itr_weight( m_WindowWeights.Begin());
            itr_win < itr_win_end;
            itr_win += m_InternalDescriptorSize, ++itr_weight, ++itr_bin_position
          )
          {
            if( util::IsDefined( *itr_win))
            {
              m_BinAverages( *itr_bin_position).AddWeightedObservation( *itr_win, *itr_weight);
            }
          }
          // accumulate results
          *itr_results = m_BinAverages.Begin()->GetAverage();
          ++itr_results;
          for
          (
            storage::Vector< math::RunningAverage< float> >::iterator
              itr( m_BinAverages.Begin() + 1), itr_prev( m_BinAverages.Begin()), itr_end( m_BinAverages.End());
            itr < itr_end;
            ++itr, ++itr_prev, ++itr_results
          )
          {
            itr->AddWeightedObservation( itr_prev->GetAverage(), itr_prev->GetWeight());
            *itr_results = itr->GetAverage();
          }
        }
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > WindowAverage< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Window, &m_Window + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool WindowAverage< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( !m_Window.ReadInitializerSuccessHook( LABEL, ERR_STREAM))
      {
        return false;
      }
      // get the window's relative indices
      storage::Vector< int> window_relative_ids( m_Window.GetRelativeIndices());
      math::Absolute( window_relative_ids);

      m_InternalDescriptorSize = m_Window.GetSizeOfFeatures() / window_relative_ids.GetSize();

      // get the coefficients for the half window
      linal::Vector< float> half_window_coefficients
      (
        m_WindowWeightsCreator->operator()( window_relative_ids.LastElement() + 1)
      );

      m_WindowWeights = linal::Vector< float>( window_relative_ids.GetSize());
      for( size_t window_position( 0), size( window_relative_ids.GetSize()); window_position < size; ++window_position)
      {
        m_WindowWeights( window_position) = half_window_coefficients( window_relative_ids( window_position));
      }
      if( m_Cumulative)
      {
        m_Bins = storage::Vector< size_t>( window_relative_ids.Begin(), window_relative_ids.End());
        m_BinAverages.Resize( window_relative_ids.LastElement() + 1);
      }

      return true;
    }

  } // namespace descriptor
} // namespace bcl
