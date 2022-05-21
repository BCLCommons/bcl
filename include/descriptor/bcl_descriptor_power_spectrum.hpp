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
#include "bcl_descriptor_power_spectrum.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_fourier_analysis_window_creator.h"
#include "bcl_descriptor_window.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType>
    PowerSpectrum< t_DataType>::PowerSpectrum( const bool &REFLECT) :
      m_Descriptor(),
      m_Size( 0),
      m_Alignment( e_Center),
      m_BandPass( false),
      m_BandPassOffset( false),
      m_MaxPeriods( 2),
      m_Reflect( REFLECT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PowerSpectrum
    template< typename t_DataType>
    PowerSpectrum< t_DataType> *PowerSpectrum< t_DataType>::Clone() const
    {
      return new PowerSpectrum( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &PowerSpectrum< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &PowerSpectrum< t_DataType>::GetAlias() const
    {
      static const std::string s_reflecting_name( "ReflectingPowerSpectrum"), s_name( "PowerSpectrum");
      return m_Reflect ? s_reflecting_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t PowerSpectrum< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return ( 1 + int( m_BandPass) + int( m_BandPassOffset)) * m_PeriodicBins.GetSize();
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
    io::Serializer PowerSpectrum< t_DataType>::GetSerializer() const
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
        "size",
        "desired size of the window (0 = just central value, 1 = central value and neighboring values, ...) used to "
        "compute the spectrogram. Larger values decrease noise but include more non-local information. "
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
      parameters.AddOptionalInitializer
      (
        "bins",
        "Contains the ranges of periods of interest to bin together; by default, all calculated amplitudes are returned "
        " as separate values / phases; but this option allows selection of the ranges of periods to bin amplitudes & phases across",
        io::Serialization::GetAgent( &m_PeriodicBins)
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
    void PowerSpectrum< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // copy the iterator
      Iterator< t_DataType> itr_seq( ELEMENT);

      // compute the window
      linal::Vector< float> window_result( m_UndefinedReplacer( itr_seq));

      // linal::ElementwiseMultiply( window_result, m_WindowWeights);

      // compute components
      window_result = m_WaveletsAndPeriods->First() * window_result;

      //BCL_Debug( window_result);

      STORAGE = float( 0.0);

      // handle periods that have both sin and cos components
      for( size_t i( 0), component( 0), n_freq( m_PeriodicBinMapping.GetSize()); i < n_freq; ++i, component += 2)
      {
        if( util::IsDefined( m_PeriodicBinMapping( i)))
        {
          STORAGE( m_PeriodicBinMapping( i)) +=
            math::Sqrt( math::Sqr( window_result( component)) + math::Sqr( window_result( component + 1)));
        }
      }

      // add phases, if desired
      if( m_BandPass || m_BandPassOffset)
      {
        const size_t base_offset( m_PeriodicBins.GetSize());
        size_t offset( base_offset);
        for( size_t i( 0), component( 0), n_freq( m_PeriodicBinMapping.GetSize()); i < n_freq; ++i, component += 2)
        {
          if( util::IsDefined( m_PeriodicBinMapping( i)))
          {
            STORAGE( base_offset + m_PeriodicBinMapping( i)) += window_result( component);
          }
        }
        if( m_BandPass)
        {
          for( size_t i( 0), n_freq( m_PeriodicBins.GetSize()); i < n_freq; ++i)
          {
            STORAGE( offset + i) = STORAGE( offset + i) / std::max( STORAGE( i), float( 1.0e-5));
          }
          offset += m_PeriodicBins.GetSize();
        }
        if( m_BandPassOffset)
        {
          const float mean( window_result( 0));
          for( size_t i( 0), component( 0), n_freq( m_PeriodicBinMapping.GetSize()); i < n_freq; ++i, component += 2)
          {
            if( util::IsDefined( m_PeriodicBinMapping( i)))
            {
              STORAGE( offset + m_PeriodicBinMapping( i)) += window_result( component);
            }
          }
          for( size_t i( 0), n_freq( m_PeriodicBins.GetSize()); i < n_freq; ++i)
          {
            STORAGE( offset + i) += mean;
          }
        }
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > PowerSpectrum< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType>
    void PowerSpectrum< t_DataType>::SetObjectHook()
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
    bool PowerSpectrum< t_DataType>::ReadInitializerSuccessHook
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
        ERR_STREAM << "PowerSpectrum only works with element-wise descriptors that return a single value";
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

      // setup windowing function
      const size_t window_size( window_descriptor.GetSizeOfFeatures());

      m_UndefinedReplacer =
        ReplaceUndefinedValues< t_DataType>
        (
          util::Implementation< Base< t_DataType, float> >( window_descriptor),
          util::Implementation< Base< t_DataType, float> >()
        );

      // create the periods
      storage::Vector< float> target_periods;
      size_t n_periods( window_size * 2 + 1);
      target_periods.PushBack( 1.0);
      target_periods.PushBack( 2.0);
      const float freq_inc( float( window_size - 2) / ( 2.0 * window_size * ( n_periods - 2)));
      const float freq_offset( 1.0 / float( window_size) - freq_inc);
      for( size_t period_number( 2); period_number < n_periods; ++period_number)
      {
        target_periods.PushBack( 1.0 / float( ( n_periods - period_number) * freq_inc + freq_offset));
      }
      target_periods.AllocateMemory( n_periods);

      m_WaveletsAndPeriods =
        FourierAnalysisWindowCreator::CreateWindowMatrix
        (
          m_Alignment,
          m_Size,
          target_periods,
          m_MaxPeriods,
          m_WindowWeightsCreator,
          0.75,
          true
        );
      const storage::Vector< float> &targetted_periods( m_WaveletsAndPeriods->Second());
      n_periods = targetted_periods.GetSize();

      if( m_PeriodicBins.IsEmpty())
      {
        // binning was not given, assign default bins
        // unless the user specifically requests it, the average is never given
        for( size_t period_number( 0); period_number < n_periods; ++period_number)
        {
          m_PeriodicBins.PushBack
          (
            math::Range< float>( targetted_periods( period_number), targetted_periods( period_number))
          );
          m_PeriodicBinMapping.PushBack( period_number);
        }
      }
      else
      {
        m_PeriodicBinMapping.Resize( n_periods);
        m_PeriodicBinMapping.SetAllElements( util::GetUndefined< size_t>());
        // set up the periodic bin mapping
        for( size_t period_number( 0); period_number < n_periods; ++period_number)
        {
          const float period_of_interest( targetted_periods( period_number));
          for( size_t bin_number( 0), n_bins( m_PeriodicBins.GetSize()); bin_number < n_bins; ++bin_number)
          {
            if( m_PeriodicBins( bin_number).IsWithin( period_of_interest))
            {
              m_PeriodicBinMapping( period_number) = bin_number;
              break;
            }
          }
        }
      }

      return true;
    }

  } // namespace descriptor
} // namespace bcl
