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
#include "bcl_descriptor_power_spectrum_sequence_width.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_fourier_analysis_window_creator.h"
#include "bcl_descriptor_power_spectrum.h"
#include "bcl_descriptor_window.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType>
    PowerSpectrumSequenceWidth< t_DataType>::PowerSpectrumSequenceWidth() :
      m_Descriptor(),
      m_NumBins( 1),
      m_FreqsPerBin( 1),
      m_MinPeriod( 2),
      m_Alignment( e_Center),
      m_BandPass( false),
      m_BandPassOffset( false),
      m_MaxPeriods( 48),
      m_Window( true),
      m_CurrentWindowIsUpToDate( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PowerSpectrumSequenceWidth
    template< typename t_DataType>
    PowerSpectrumSequenceWidth< t_DataType> *PowerSpectrumSequenceWidth< t_DataType>::Clone() const
    {
      return new PowerSpectrumSequenceWidth( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &PowerSpectrumSequenceWidth< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &PowerSpectrumSequenceWidth< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "PowerSpectrumSequenceWidth");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t PowerSpectrumSequenceWidth< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return 2 + int( m_BandPass) + int( m_BandPassOffset);
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
    io::Serializer PowerSpectrumSequenceWidth< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "computes the spectral power and phase of oscillations in a given variable");
      parameters.AddInitializer
      (
        "",
        "descriptor to use for computing the amplitude",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      parameters.AddInitializer
      (
        "number bins",
        "Number of different frequency bands to assign each sequence",
        io::Serialization::GetAgentWithMin( &m_NumBins, size_t( 1))
      );
      parameters.AddInitializer
      (
        "frequencies per bin",
        "Number of frequency bands per bin",
        io::Serialization::GetAgentWithMin( &m_FreqsPerBin, size_t( 1)),
        "1"
      );
      parameters.AddInitializer
      (
        "min period",
        "Minimum period to allow.  If the sequence is too short, the result will be 0 padded",
        io::Serialization::GetAgentWithMin( &m_MinPeriod, size_t( 2))
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
        "bandpass",
        "True to include the band-pass filtered values in the output (not normalized)",
        io::Serialization::GetAgent( &m_BandPass),
        "False"
      );
      parameters.AddInitializer
      (
        "bandpass offset",
        "True to include the band-pass filtered values in the output, offset by the mean value",
        io::Serialization::GetAgent( &m_BandPassOffset),
        "False"
      );
      parameters.AddInitializer
      (
        "weighting",
        "method of generating weights for the windowing function; rectangular emphasizes longer-range frequencies, while "
        "triangular emphasizes the local oscillation",
        io::Serialization::GetAgent( &m_WindowWeightsCreator),
        "Hamming"
      );
      parameters.AddInitializer
      (
        "max periods",
        "for each frequency; maximum number of periods to consider before zero-ing out the contribution; helps improve"
        "localization of high-frequencies using large windows when set to 2-3",
        io::Serialization::GetAgent( &m_MaxPeriods),
        "2"
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void PowerSpectrumSequenceWidth< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( !m_CurrentWindowIsUpToDate)
      {
        UpdateCurrentWindow( ELEMENT);
      }

      // get the current position
      const size_t current_pos( ELEMENT.GetPosition());

      // get the sequence size
      const size_t seqsize( ELEMENT.GetSize());

      // get the position in the window array
      const size_t window_array_position( m_Alignment != e_Left ? current_pos : seqsize - current_pos - 1);

      const linal::Matrix< float> &wavelets( m_WaveletsAndPeriods->First());

      // compute a reference to the internal array
      linal::VectorConstReference< float> reference
      (
        wavelets.GetNumberCols(),
        m_CurrentWindow[ window_array_position]
      );

      // handle the constant amplitude
      STORAGE( 0) = m_WaveletsAndPeriods->First().GetRow( 0) * reference;

      const float cos_component( wavelets.GetRow( m_CurrentMaximumBin * 2) * reference);
      const float sin_component( wavelets.GetRow( m_CurrentMaximumBin * 2 + 1) * reference);

      // handle periods that have both sin and cos components
      STORAGE( 1) = math::Sqrt( math::Sqr( cos_component) + math::Sqr( sin_component));

      // add phases, if desired
      if( m_BandPass || m_BandPassOffset)
      {
        // add the band-passed values, which are just the cosines
        if( m_BandPass)
        {
          STORAGE( 2) = cos_component;
        }
        if( m_BandPassOffset)
        {
          STORAGE( 2 + m_BandPass) = cos_component + STORAGE( 0);
        }
      }
    }

    //! @brief Function to update the current window; called if working on a new object and this descriptor was not
    //!        already cached
    //! @param ELEMENT the element of the sequence of interest
    template< typename t_DataType>
    void PowerSpectrumSequenceWidth< t_DataType>::UpdateCurrentWindow( const iterate::Generic< const t_DataType> &ELEMENT)
    {
      const size_t sequence_size( ELEMENT.GetSize());

      // copy the iterator
      Iterator< t_DataType> itr_seq( ELEMENT.Begin());

      // get the values for the entire sequence
      m_Window.SetRadius( sequence_size - 1);
      m_Window.SetAlignment( e_Right);
      const bool is_monodirectional( m_Alignment == e_Left || m_Alignment == e_Right);

      // determine the target periods
      storage::Vector< float> target_periods;
      target_periods.AllocateMemory( m_FreqsPerBin * m_NumBins + 1);
      target_periods.PushBack( 1.0);
      const float freq_increment( 1.0 / float( m_FreqsPerBin));
      for( size_t integer( 1); integer <= m_NumBins; ++integer)
      {
        const size_t start_period( sequence_size / integer);
        if( start_period < m_MinPeriod)
        {
          break;
        }
        for( size_t sub_freq( 0); sub_freq < m_FreqsPerBin; ++sub_freq)
        {
          const float candidate_period( float( sequence_size) / ( float( integer) + sub_freq * freq_increment));
          if( candidate_period < float( m_MinPeriod))
          {
            break;
          }
          target_periods.PushBack( candidate_period);
        }
      }

      {
        linal::VectorConstReference< float> reference_whole_sequence( m_Window( itr_seq));

        // create a copy of the reference vector, reversed;
        linal::Vector< float> reference_reversed( reference_whole_sequence);
        std::reverse( reference_reversed.Begin(), reference_reversed.End());

        // determine the size of the final vector and setup m_CurrentWindow to reflect that
        m_CurrentWindow
          = linal::Vector< float>( size_t( is_monodirectional ? 2 * sequence_size - 1 : 3 * sequence_size - 2));

        if( m_Alignment == e_Right)
        {
          // copy the vectors into the current window in strictly ascending order, accounting for reflection
          std::copy
          (
            reference_reversed.Begin() + 1,
            reference_reversed.End(),
            std::copy
            (
              reference_whole_sequence.Begin(),
              reference_whole_sequence.End(),
              m_CurrentWindow.Begin()
            )
          );
        }
        else
        {
          // copy left facing window, then right, then left again for center window reflection
          float *itr_copy
          (
            std::copy
            (
              reference_whole_sequence.Begin(),
              reference_whole_sequence.End(),
              std::copy
              (
                reference_reversed.Begin(),
                reference_reversed.End() - 1,
                m_CurrentWindow.Begin()
              )
            )
          );

          if( m_Alignment != e_Left)
          {
            std::copy
            (
              reference_reversed.Begin() + 1,
              reference_reversed.End(),
              itr_copy
            );
          }
        }
      }

      // now get window indices
      m_Window.SetAlignment( is_monodirectional ? m_Alignment.GetEnum() : e_JufoCenter);

      if( !m_WaveletsAndPeriods.IsDefined() || m_WaveletsAndPeriods->First().GetNumberCols() != sequence_size)
      {
        // update the wavelets for this object
        m_WaveletsAndPeriods =
          FourierAnalysisWindowCreator::CreateWindowMatrix
          (
            is_monodirectional ? m_Alignment.GetEnum() : e_JufoCenter,
            sequence_size - 1,
            target_periods,
            m_MaxPeriods,
            m_WindowWeightsCreator,
            1.0,
            sequence_size < 300 // only cache for reasonably small sequence sizes
          );
        m_TempCoefficients = linal::Vector< float>( m_WaveletsAndPeriods->First().GetNumberRows(), float( 0));
      }
      m_CurrentWindowIsUpToDate = true;

      const size_t n_freq_calc( m_TempCoefficients.GetSize() / 2);

      linal::Vector< float> amplitudes( n_freq_calc);

      // now walk through the sequence, one vector at a time, accumulate coefficient magnitudes
      for( size_t i( 0); i < sequence_size; ++i)
      {
        // compute a reference to the internal array
        linal::VectorConstReference< float> reference( m_WaveletsAndPeriods->First().GetNumberCols(), m_CurrentWindow[ i]);
        m_TempCoefficients = float( 0.0);
        linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( m_TempCoefficients, m_WaveletsAndPeriods->First(), reference);

        for( size_t j( 0), component( 0); j < n_freq_calc; ++j, component += 2)
        {
          amplitudes( j) +=
            math::Sqrt( math::Sqr( m_TempCoefficients( component)) + math::Sqr( m_TempCoefficients( component + 1)));
        }
      }

      // determine the maximum amplitude of any oscillation
      float maximum_bin_amp( 0);
      size_t max_bin( 0);
      for( size_t j( 1); j < n_freq_calc; ++j)
      {
        if( amplitudes( j) > maximum_bin_amp)
        {
          max_bin = j;
          maximum_bin_amp = amplitudes( j);
        }
      }

      m_CurrentMaximumBin = max_bin;

      // determine the oscillation of the given frequency
//      float max_frequency( target_periods( max_bin));
//      BCL_MessageStd
//      (
//        "Max frequency: " + util::Format()( max_frequency)
//        + " for seq of size " + util::Format()( sequence_size)
//        + " had amplitude: " + util::Format()( maximum_bin_amp)
//      );
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > PowerSpectrumSequenceWidth< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType>
    void PowerSpectrumSequenceWidth< t_DataType>::SetObjectHook()
    {
      m_CurrentWindowIsUpToDate = false;
      m_Window.SetObject( *this->GetCurrentObject());
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool PowerSpectrumSequenceWidth< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if
      (
        m_Descriptor.IsDefined()
        &&
        ( m_Descriptor->GetType().GetDimension() != size_t( 1) || m_Descriptor->GetSizeOfFeatures() != size_t( 1))
      )
      {
        ERR_STREAM << "PowerSpectrumSequenceWidth only works with element-wise descriptors that return a single value";
        return false;
      }

      m_Window = Window< t_DataType, float>( true, m_Descriptor);

      return true;
    }

  } // namespace descriptor
} // namespace bcl
