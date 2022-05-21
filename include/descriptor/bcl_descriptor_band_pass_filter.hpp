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
#include "bcl_descriptor_band_pass_filter.h"
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
    BandPassFilter< t_DataType>::BandPassFilter( const bool &REFLECTING) :
      m_Descriptor(),
      m_Size( 0),
      m_HighPassPeriod( 1),
      m_LowPassPeriod( 1000),
      m_Alignment( e_Center),
      m_Reflect( REFLECTING)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BandPassFilter
    template< typename t_DataType>
    BandPassFilter< t_DataType> *BandPassFilter< t_DataType>::Clone() const
    {
      return new BandPassFilter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &BandPassFilter< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &BandPassFilter< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "BandPassFilter"), s_reflecting_name( "ReflectingBandPassFilter");
      return m_Reflect ? s_reflecting_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t BandPassFilter< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor.IsDefined() ? m_Descriptor->GetSizeOfFeatures() : 1;
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
    io::Serializer BandPassFilter< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Bandpass a single element of a windowed signal. "
        "If only maxp is given, functions as a high pass filter  (only high frequency signals preserved). "
        "If only minp is given, functions as a low pass filter (only low frequency signals preserved).  "
      );
      parameters.AddInitializer
      (
        "",
        "descriptor to use for computing the amplitude",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      parameters.AddInitializer
      (
        "size",
        "desired size of the window (0 = just central value, 1 = central value and neighboring values, ...) used to "
        "compute the periodogram. Larger values decrease noise but include more non-local information and are slower. "
        "Typically a value somewhat larger than the maximum period of interest is used",
        io::Serialization::GetAgentWithMin( &m_Size, size_t( 1))
      );
      parameters.AddInitializer
      (
        "alignment",
        "Alignment of the window.  Use center to consider the window around this element; left to consider the window"
        " up until this element, and right to consider the window following this element. ",
        io::Serialization::GetAgent( &m_Alignment),
        "Center"
      );
      parameters.AddInitializer
      (
        "minp",
        "Minimal period desired; <= 2  for a high-pass filter",
        io::Serialization::GetAgent( &m_HighPassPeriod),
        "1"
      );
      parameters.AddInitializer
      (
        "maxp",
        "Maximum period desired; >= window size  for a low-pass filter",
        io::Serialization::GetAgent( &m_LowPassPeriod),
        "10000"
      );
      parameters.AddInitializer
      (
        "weighting",
        "method of generating weights for the windowing function; rectangular emphasizes longer-range frequencies, while "
        "triangular emphasizes the local oscillation",
        io::Serialization::GetAgent( &m_WindowWeightsCreator),
        "Rectangular"
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void BandPassFilter< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // If the descriptor only returns a single value, this function reduces to:
      // STORAGE = linal::ScalarProduct( m_Filter, m_UndefinedReplacer( Iterator< t_DataType>( ELEMENT)));
      // copy the iterator and call the internal (window) descriptor
      const linal::VectorConstInterface< float> &vals_window( m_UndefinedReplacer( Iterator< t_DataType>( ELEMENT)));
      const float *itr_descriptor_window( vals_window.Begin());
      for
      (
        const float *itr_filter( m_Filter.Begin()), *itr_filter_end( m_Filter.End());
        itr_filter != itr_filter_end;
        ++itr_filter
      )
      {
        for
        (
          float *itr_storage( STORAGE.Begin()), *itr_storage_end( STORAGE.End());
          itr_storage != itr_storage_end;
          ++itr_storage, ++itr_descriptor_window
        )
        {
          *itr_storage += *itr_filter * *itr_descriptor_window;
        }
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > BandPassFilter< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType>
    void BandPassFilter< t_DataType>::SetObjectHook()
    {
      m_UndefinedReplacer.SetObject( *this->GetCurrentObject());
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool BandPassFilter< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_Descriptor.IsDefined() && m_Descriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "BandPassFilter only works with element-wise descriptors";
        return false;
      }

      storage::Vector< util::ObjectDataLabel> window_args;
      window_args.PushBack( util::ObjectDataLabel( "size", util::Format()( m_Size)));
      window_args.PushBack( util::ObjectDataLabel( "", m_Descriptor.GetLabel()));
      window_args.PushBack( util::ObjectDataLabel( "alignment", m_Alignment.GetLabel()));
      Window< t_DataType, float> window_descriptor( m_Reflect);
      if( !window_descriptor.TryRead( util::ObjectDataLabel( "", "", window_args), ERR_STREAM))
      {
        return false;
      }

      // Internal descriptor size
      const size_t internal_descriptor_size( m_Descriptor->GetSizeOfFeatures());

      // setup windowing function
      const size_t total_window_size( window_descriptor.GetSizeOfFeatures() / internal_descriptor_size);

      // get the relative indices for the window descriptor
      storage::Vector< int> window_indices( window_descriptor.GetRelativeIndices());

      m_UndefinedReplacer =
        ReplaceUndefinedValues< t_DataType>
        (
          util::Implementation< Base< t_DataType, float> >( window_descriptor),
          util::Implementation< Base< t_DataType, float> >()
        );

      m_Filter = linal::Vector< float>( total_window_size, float( 1));
      const float two_pi( 2.0f * math::g_Pi);

      const bool is_highpass( m_HighPassPeriod <= 2.0);
      const bool is_lowpass( m_LowPassPeriod > total_window_size);
      if( is_highpass && is_lowpass)
      {
        // trivial case, no frequencies to be filtered, pass signal through unchanged
        m_Filter = float( 0);
        m_Filter( 0) = float( 1);
        return true;
      }
      else if( is_highpass)
      {
        // accepting all frequencies above the low pass
        const float frequency( two_pi / m_LowPassPeriod);
        for( size_t i( 0); i < total_window_size; ++i)
        {
          const int w( window_indices( i));
          if( !w)
          {
            m_Filter( i) = 1.0 - 2.0 / m_LowPassPeriod;
          }
          else
          {
            m_Filter( i) = -sin( frequency * w) / ( math::g_Pi * w);
          }
        }
      }
      else if( is_lowpass)
      {
        // accepting all frequencies below the high pass
        const float frequency( two_pi / m_HighPassPeriod);
        for( size_t i( 0); i < total_window_size; ++i)
        {
          const int w( window_indices( i));
          if( !w)
          {
            m_Filter( i) = 2.0 / m_HighPassPeriod;
          }
          else
          {
            m_Filter( i) = sin( frequency * w) / ( math::g_Pi * w);
          }
        }
      }
      else
      {
        // true band pass filter
        const float high_frequency( two_pi / m_HighPassPeriod);
        const float low_frequency( two_pi / m_LowPassPeriod);
        for( size_t i( 0); i < total_window_size; ++i)
        {
          const int w( window_indices( i));
          if( !w)
          {
            m_Filter( i) = ( high_frequency - low_frequency) / math::g_Pi;
          }
          else
          {
            m_Filter( i) = ( sin( high_frequency * w) - sin( low_frequency * w)) / ( math::g_Pi * w);
          }
        }
      }

      // multiply by window weights
      linal::Vector< float> half_window_coefficients( m_WindowWeightsCreator->operator()( m_Size + 1));
      for( size_t i( 0); i < total_window_size; ++i)
      {
        m_Filter( i) *= half_window_coefficients( math::Absolute( window_indices( i)));
      }
      //m_Filter /= linal::AbsSum( m_Filter);
      const float expected_sum( is_lowpass ? 1.0 : 0.0);
      m_Filter -= ( m_Filter.Sum() - expected_sum) / float( total_window_size);
      //m_Filter.Normalize();

      return true;
    }

  } // namespace descriptor
} // namespace bcl
