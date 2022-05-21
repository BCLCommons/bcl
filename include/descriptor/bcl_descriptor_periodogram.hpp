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
#include "bcl_descriptor_periodogram.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_fourier_analysis_window_creator.h"
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
    Periodogram< t_DataType>::Periodogram( const bool &REFLECTING) :
      m_Descriptor(),
      m_Size( 0),
      m_Alignment( e_Center),
      m_BandPass( false),
      m_Reflect( REFLECTING),
      m_ReturnDominantPeriod( false),
      m_NumberNonTrivialPeriodRanges( 0),
      m_MaxPeriods( 2),
      m_WindowWeightsCreator()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Periodogram
    template< typename t_DataType>
    Periodogram< t_DataType> *Periodogram< t_DataType>::Clone() const
    {
      return new Periodogram( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Periodogram< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &Periodogram< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "Periodogram"), s_reflecting_name( "ReflectingPeriodogram");
      return m_Reflect ? s_reflecting_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t Periodogram< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return ( 1 + int( m_BandPass)) * m_PeriodicBins.GetSize()
             + m_ReturnDominantPeriod * m_NumberNonTrivialPeriodRanges;
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
    io::Serializer Periodogram< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "computes amplitudes of oscillations in a given property at desired frequencies");
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
        "compute the periodogram. Larger values decrease noise but include more non-local information. "
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
        "maxperiod",
        "True to include the period that had the maximum amplitude for each range composed of multiple frequencies",
        io::Serialization::GetAgent( &m_ReturnDominantPeriod),
        "False"
      );
      parameters.AddInitializer
      (
        "bins",
        "Contains the ranges of periods of interest to bin together",
        io::Serialization::GetAgent( &m_PeriodicBins)
      );
      parameters.AddInitializer
      (
        "bin resolution",
        "Contains the resolution of bins, should be one value per bin",
        io::Serialization::GetAgent( &m_Resolution)
      );
      parameters.AddInitializer
      (
        "weighting",
        "method of generating weights for the windowing function; rectangular emphasizes longer-range frequencies, while "
        "triangular emphasizes the local oscillation",
        io::Serialization::GetAgent( &m_WindowWeightsCreator),
        "Rectangular"
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

    //! @brief helper function; returns all local maxima positions in a vector
    //! @param VECTOR the vector of values to consider
    //! @return all local maxima positions in a vector
    storage::Vector< size_t> GetMaximaPositions( const linal::VectorConstInterface< float> &VECTOR)
    {
      if( VECTOR.IsEmpty())
      {
        return storage::Vector< size_t>();
      }
      storage::Vector< size_t> positions;
      positions.PushBack( 0);
      //BCL_Debug( VECTOR);
      for( size_t current_pos( 1), sz( VECTOR.GetSize()); current_pos < sz; ++current_pos)
      {
        if( VECTOR( current_pos) >= VECTOR( current_pos - 1))
        {
          if( current_pos - 1 == positions.LastElement())
          {
            positions.LastElement() = current_pos;
          }
          else if( VECTOR( current_pos) > VECTOR( current_pos - 1))
          {
            positions.PushBack( current_pos);
          }
        }
      }
      return positions;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void Periodogram< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // copy the iterator
      Iterator< t_DataType> itr_seq( ELEMENT);

      // compute the window
      linal::Vector< float> window_result( m_UndefinedReplacer( itr_seq));

      // BCL_Debug( min_amp);

      // basic algorithm is: choose the dominant frequency for each bin at each step
      // Use least squares spectral analysis if the inverse is reasonably well conditioned, otherwise
      // remove the most significant frequencies and move on
      m_Components = 0.0;
      m_LocalComponents = 0.0;
      const storage::Vector< float> &target_periods( m_WaveletsAndPeriods->Second());
      const size_t n_periods( m_WaveletsAndPeriods->Second().GetSize());
      for
      (
        size_t dominant_component_id( 0), max_components( 2 * m_PeriodicBins.GetSize());
        dominant_component_id < max_components;
        ++dominant_component_id
      )
      {
        // use the storage vector to hold all amplitudes
        m_LocalComponents = m_WaveletsAndPeriods->First() * window_result;
        for( size_t period_number( 0); period_number < n_periods; ++period_number)
        {
          // compute the amplitude and store
          m_LocalAmplitudes( period_number)
            = math::Sqrt( math::Sqr( m_LocalComponents( 2 * period_number)) + math::Sqr( m_LocalComponents( 2 * period_number + 1)));
        }

        // BCL_Debug( max_primary_amp);
        // BCL_Debug( min_amp);
        // get all local maxima out of the spectra
        storage::Vector< size_t> direct_local_maxima( GetMaximaPositions( m_LocalAmplitudes));

        //BCL_Debug( direct_local_maxima);
        // remove the last component, if any
        if( direct_local_maxima.IsEmpty())
        {
          break;
        }
//        if( direct_local_maxima.LastElement() == target_periods.GetSize() - 1 && direct_local_maxima.GetSize() > size_t( 1))
//        {
//          // aliasing bias; often undue amounts of weight gets transferred to the last bin due to the welch window
//          // so never use the last frequency if there are any others available
//          direct_local_maxima.PopBack();
//        }
        {
          float max_amp( 0), max_sine_amp( 0), max_cos_amp( 0);
          size_t max_period_number( util::GetUndefined< size_t>());
          // for each period
          for( size_t maxim_number( 0), n_maxima( direct_local_maxima.GetSize()); maxim_number < n_maxima; ++maxim_number)
          {
            // compute the amplitude and store
            float local_amp( m_LocalAmplitudes( direct_local_maxima( maxim_number)));
            if( local_amp > max_amp)
            {
              max_amp = local_amp;
              max_period_number = direct_local_maxima( maxim_number);
              max_sine_amp = m_LocalComponents( 2 * max_period_number + 1);
              max_cos_amp = m_LocalComponents( 2 * max_period_number);
            }
          }
          if( util::IsDefined( max_period_number))
          {
            m_Components( 2 * max_period_number) += max_cos_amp;
            m_Components( 2 * max_period_number + 1) += max_sine_amp;
            const float *cos_window( m_CosWindowsUnweighted( max_period_number).Begin());
            const float *sin_window( m_SinWindowsUnweighted( max_period_number).Begin());
//            BCL_MessageDbg
//            (
//              "Removing component: " + util::Format()( target_periods( max_period_number))
//              + " amp: " + util::Format()( max_amp)
//            );
            // BCL_Debug( m_Components);
            BCL_Assert
            (
              m_CosWindowsUnweighted( max_period_number).GetSize() == window_result.GetSize(),
              "Too short cos window, was " + util::Format()( m_CosWindowsUnweighted( max_period_number).GetSize())
              + " should be " + util::Format()( window_result.GetSize())
            );
            BCL_Assert
            (
              m_SinWindowsUnweighted( max_period_number).GetSize() == window_result.GetSize(),
              "Too short sin window, was " + util::Format()( m_SinWindowsUnweighted( max_period_number).GetSize())
              + " should be " + util::Format()( window_result.GetSize())
            );

            for( size_t i( 0), signal_size( window_result.GetSize()); i < signal_size; ++i)
            {
              window_result( i) -= max_sine_amp * sin_window[ i] + max_cos_amp * cos_window[ i];
            }
            // BCL_Debug( window_result);
          }
          else
          {
            break;
          }
        }
      }

      float mean( 0);
      if( target_periods( 0) < float( 1.5))
      {
        mean = m_Components( 0);
      }
      else
      {
        mean = window_result.Sum() / float( window_result.GetSize());
      }
      {
        //BCL_Debug( m_Components);
        storage::Vector< float>::const_iterator
          itr_periods( target_periods.Begin()), itr_periods_end( target_periods.End());
        size_t period_number( 0), bin_number( 0);
        size_t dominant_period_index( m_PeriodicBins.GetSize() * ( 1 + int( m_BandPass)) - 1);
        for
        (
          storage::Vector< math::Range< float> >::const_iterator
            itr_bin( m_PeriodicBins.Begin()), itr_bin_end( m_PeriodicBins.End());
          itr_bin != itr_bin_end;
          ++itr_bin, ++bin_number
        )
        {
          const size_t start_period_number( period_number);
          for( ; itr_periods != itr_periods_end && itr_bin->IsWithin( *itr_periods); ++itr_periods, ++period_number)
          {
            STORAGE( bin_number) +=
              math::Sqrt( math::Sqr( m_Components( 2 * period_number)) + math::Sqr( m_Components( 2 * period_number + 1)));
            if( m_BandPass)
            {
              STORAGE( bin_number + m_PeriodicBins.GetSize()) += m_Components( 2 * period_number);
            }
          }
          if( m_BandPass)
          {
            STORAGE( bin_number + m_PeriodicBins.GetSize()) += mean;
          }
          if( m_ReturnDominantPeriod && period_number - start_period_number > 1)
          {
            float actual_dominant_period( target_periods( start_period_number));
            float actual_dominant_amplitude( 0);
            // determine the dominant period on this interval
            for( size_t i( start_period_number); i < period_number; ++i)
            {
              const float local_amp( math::Sqrt( math::Sqr( m_Components( 2 * i)) + math::Sqr( m_Components( 2 * i + 1))));
              if( local_amp > actual_dominant_amplitude)
              {
                actual_dominant_period = target_periods( i);
                actual_dominant_amplitude = local_amp;
              }
            }
            STORAGE( ++dominant_period_index) = actual_dominant_period;
          }
        }
      }
      //BCL_MessageDbg( "Remainder: " + util::Format()( window_result));
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > Periodogram< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType>
    void Periodogram< t_DataType>::SetObjectHook()
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
    bool Periodogram< t_DataType>::ReadInitializerSuccessHook
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
        ERR_STREAM << "Periodogram only works with element-wise descriptors that return a single value";
        return false;
      }

      // determine all target periods
      if( m_PeriodicBins.GetSize() != m_Resolution.GetSize())
      {
        ERR_STREAM << "One resolution level is needed per bin";
        return false;
      }

      // ensure that the ranges are sorted
      {
        storage::Vector< math::Range< float> > bins_copy( m_PeriodicBins);
        bins_copy.Sort( std::less< math::Range< float> >());
        if( !( bins_copy == m_PeriodicBins))
        {
          ERR_STREAM << "Periodic bins must be sorted";
          return false;
        }
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

      // store the target periods, temporarily
      storage::Vector< float> target_periods;

      // determine target periods based on periodic bins and resolution
      {
        storage::Vector< float>::const_iterator itr_res( m_Resolution.Begin());
        m_NumberNonTrivialPeriodRanges = 0;
        // carry-over frequency; used to ensure distinct frequencies are sampled
        float carry_over_freq( 0);
        for
        (
          storage::Vector< math::Range< float> >::const_iterator
            itr_bin( m_PeriodicBins.Begin()), itr_bin_end( m_PeriodicBins.End());
          itr_bin != itr_bin_end;
          ++itr_bin, ++itr_res
        )
        {
          const float resolution( *itr_res);
          if( !resolution)
          {
            carry_over_freq = 0.0;
            target_periods.PushBack( itr_bin->GetMiddle());
            continue;
          }

          // iterate from bin min to bin max, in steps of size resolution
          math::Range< float> closed_range( itr_bin->CloseBorders());
          size_t initial_target_periods_sz( target_periods.GetSize());

          if( carry_over_freq)
          {
            closed_range = math::Range< float>( closed_range.GetMin() + carry_over_freq, closed_range.GetMax());
            carry_over_freq = 0.0;
          }

          storage::Vector< math::Range< float> >::const_iterator itr_bin_next( itr_bin + 1);
          if( itr_bin_next != itr_bin_end)
          {
            // adjust the min/max frequency
            math::Range< float> next_closed_range( itr_bin_next->CloseBorders());
            const float res_next( *( itr_res + 1));
            const float this_width( closed_range.GetWidth());
            const float next_width( next_closed_range.GetWidth());
            const float this_effective_resolution( std::min( resolution, this_width));
            const float next_effective_resolution( std::min( next_width, res_next));
            const float ave_resolution( 0.5 * ( this_effective_resolution + next_effective_resolution));
            const float overlap( closed_range.GetMax() - next_closed_range.GetMin() + ave_resolution);
            if( overlap > 0.0)
            {
              const float width_ratio_this( this_width / ( this_width + next_width));
              const float this_resolution_distribution( overlap * width_ratio_this);
              carry_over_freq = overlap * ( 1.0 - width_ratio_this);
              closed_range =
                math::Range< float>
                (
                  closed_range.GetMin(),
                  closed_range.GetMax() - this_resolution_distribution
                );
            }
          }

          for( float min( closed_range.GetMin()), max( closed_range.GetMax()); min <= max; min += resolution)
          {
            target_periods.PushBack( min);
          }
          if( target_periods.GetSize() > initial_target_periods_sz + 1)
          {
            ++m_NumberNonTrivialPeriodRanges;
          }
        }
      }

      m_UndefinedReplacer =
        ReplaceUndefinedValues< t_DataType>
        (
          util::Implementation< Base< t_DataType, float> >( window_descriptor),
          util::Implementation< Base< t_DataType, float> >()
        );

      m_WaveletsAndPeriods
        = FourierAnalysisWindowCreator::CreateWindowMatrix
          (
            m_Alignment,
            m_Size,
            target_periods,
            m_MaxPeriods,
            m_WindowWeightsCreator,
            1.0,
            true
          );

      // precalculate cosine and sine windows at all desired frequences
      m_CosWindowsUnweighted.Resize( target_periods.GetSize());
      m_SinWindowsUnweighted.Resize( target_periods.GetSize());

      size_t period_number( 0);
      for
      (
        storage::Vector< float>::const_iterator
          itr_period( target_periods.Begin()), itr_period_end( target_periods.End());
        itr_period != itr_period_end;
        ++itr_period, ++period_number
      )
      {
        storage::Pair< linal::VectorConstReference< float>, linal::VectorConstReference< float> > cos_sin_window
        (
          FourierAnalysisWindowCreator::GetUnweightedCosSinCoefficients( m_Alignment, m_Size, *itr_period)
        );
        m_CosWindowsUnweighted( period_number) = cos_sin_window.First();
        m_SinWindowsUnweighted( period_number) = cos_sin_window.Second();
      }

      // compute cross correlations between all windows
      size_t total_windows( target_periods.GetSize() * 2);
      m_Components = m_LocalComponents = linal::Vector< float>( total_windows, float( 0.0));
      m_LocalAmplitudes = linal::Vector< float>( total_windows / 2, float( 0.0));

      return true;
    }

  } // namespace descriptor
} // namespace bcl
