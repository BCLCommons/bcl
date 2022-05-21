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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_fourier_analysis_window_creator.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_constants.h"
#include "linal/bcl_linal_vector_operations.h"
#include "sched/bcl_sched_mutex.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief set up the window matrix
    //! @param ALIGNMENT alignment to use for the matrix
    //! @param SIZE size of the window to use
    //! @param TARGET_PERIODS periods of interest of interest
    //! @param MAX_PERIODS maximum # of periods to insert into the signal; also controls windowing
    //! @param WEIGHTING object that creates the weights vector
    //! @param MAX_OVERLAP the maximum overlap to allow for any two adjacent signals
    //!        Values >= 1 ensure that all frequencies are included, even if they are not orthogonal
    //! @param CACHEABLE true to allow caching the matrix and retrieval of the matrix from the cache, if it is available
    //! @return wavelets matrix and actual periods represented by each pair of rows in the matrix
    FourierAnalysisWindowCreator::WaveletsFrequenciesPtr FourierAnalysisWindowCreator::CreateWindowMatrix
    (
      const WindowAlignmentEnum &ALIGNMENT,
      const size_t &SIZE,
      const storage::Vector< float> &TARGET_PERIODS,
      const size_t &MAX_PERIODS,
      const util::Implementation< WindowWeightingInterface> &WEIGHTING,
      const float &MAX_OVERLAP,
      const bool &CACHEABLE
    )
    {
      //! Mutex to protect access to s_Inverses
      static sched::Mutex s_wavelets_mutex;

      //! Map from initializer string/type to the inverse spectral matrix; saves recalculation and storing of inverses
      //! and target periods
      typedef storage::Map
      <
        std::string, //! WindowSize;Alignment;Weighting;MaxPeriods;MaxOverlap;TargetPeriods
        util::OwnPtr
        <
          storage::Pair
          <
            linal::Matrix< float>,    //! Wavelets; arranged in ascending order of frequency, interweaved cos-sin
            storage::Vector< float> //! Frequencies represented by each row-pair in the matrix
          >
        >
      > MapType;

      static MapType s_wavelets;

      // if the window matrix is to be cacheable, first look for it in the map
      std::string query_string;
      if( CACHEABLE)
      {
        s_wavelets_mutex.Lock();
        query_string =
          ALIGNMENT.GetString() + ';' + util::Format()( SIZE) + ';'
          + util::Format()( MAX_PERIODS) + ';' + util::Format()( std::min( MAX_OVERLAP, float( 1.0))) + ';'
          + WEIGHTING.GetString() + ';' + util::Format()( TARGET_PERIODS);
        MapType::iterator itr_map( s_wavelets.Find( query_string));
        if( itr_map != s_wavelets.End())
        {
          s_wavelets_mutex.Unlock();
          // return a non-owning pointer to the mapped-type
          return WaveletsFrequenciesPtr( &*itr_map->second, false);
        }
      }

      // create the window ids mapping
      Window< char, float> window( false, Constants< char, float>( 1.0));
      window.SetRadius( SIZE);
      window.SetAlignment( ALIGNMENT);
      storage::Vector< int> window_ids( window.GetRelativeIndices());

      // get the maximum absolute window position, this controls the window coefficients
      const size_t last_window_place( math::Absolute( window_ids.LastElement()));
      // get the actual window size
      const size_t window_size( window_ids.GetSize());

      // detect central alignment
      const bool is_center( ALIGNMENT == e_JufoCenter || ALIGNMENT == e_Center);

      // get the coefficients for the half window
      linal::Vector< float> window_weights( window_size);

      size_t n_periods( TARGET_PERIODS.GetSize());

      linal::Vector< float> cos_window( window_size), sin_window( window_size);
      linal::Vector< float> unweighted_cos_window( window_size), unweighted_sin_window( window_size);

      storage::Vector< float> self_scalar_product;
      self_scalar_product.AllocateMemory( n_periods);

      linal::Matrix< float> window_matrix( 2 * n_periods, window_size, float( 0.0));

      storage::Vector< float> new_target_periods;
      new_target_periods.AllocateMemory( n_periods);

      size_t row_number( 0);
      for( size_t period_number( 0); period_number < n_periods; ++period_number)
      {
        const float &period( TARGET_PERIODS( period_number));
        const float frequency( 2.0f * math::g_Pi / period);
        size_t max_position( MAX_PERIODS * period);
        if( period > 1.1)
        {
          max_position = std::min( max_position, last_window_place);
          linal::Vector< float> half_window_coefficients
          (
            WEIGHTING->operator()( max_position + 1)
          );
          for( size_t window_position( 0); window_position < window_size; ++window_position)
          {
            if( size_t( math::Absolute( window_ids( window_position))) <= max_position)
            {
              window_weights( window_position) = half_window_coefficients( math::Absolute( window_ids( window_position)));
            }
            else
            {
              window_weights( window_position) = 0;
            }
          }
        }
        else
        {
          // period = 1; use normal weight coefficients
          max_position = last_window_place;
          linal::Vector< float> half_window_coefficients
          (
            WEIGHTING->operator()( last_window_place + 1)
          );
          //BCL_Debug( half_window_coefficients);
          for( size_t window_position( 0); window_position < window_size; ++window_position)
          {
            window_weights( window_position) = half_window_coefficients( math::Absolute( window_ids( window_position)));
          }
        }
        for( size_t i( 0); i < window_size; ++i)
        {
          const float position( frequency * window_ids( i));
          unweighted_cos_window( i) = cos( position);
          cos_window( i) = unweighted_cos_window( i) * window_weights( i);
          unweighted_sin_window( i) = sin( position);
          sin_window( i) = unweighted_sin_window( i) * window_weights( i);
        }

        cos_window /= linal::ScalarProduct( unweighted_cos_window, cos_window);
        if( period <= 2.0001)
        {
          sin_window = float( 0.0);
        }
        else
        {
          const float effective_window_size( ( is_center ? 2 : 1) * max_position + 1);
          const float cos_offset( cos_window.Sum() / effective_window_size);
          sin_window /= linal::ScalarProduct( unweighted_sin_window, sin_window);
          const float sin_offset( sin_window.Sum() / effective_window_size);
          for( size_t i( 0); i < window_size; ++i)
          {
            if( size_t( math::Absolute( window_ids( i))) <= max_position)
            {
              cos_window( i) -= cos_offset;
              sin_window( i) -= sin_offset;
            }
          }
        }
        const float cos_sqnorm( 0.5 * linal::ScalarProduct( cos_window, cos_window));
        const float sin_sqnorm( 0.5 * linal::ScalarProduct( sin_window, sin_window));

        if( MAX_OVERLAP < 1.0)
        {
          float mx( 0);
          // compute the approximate covariance between the cos window and previous windows, ditto for the sine window
          // only add this frequency if product / ave sum is less than 0.8.  This ensures relatively little issues with
          // orthogonality
          for( size_t i( 0); i < row_number; ++i)
          {
            mx = std::max( double( mx), math::Absolute( linal::ScalarProduct( cos_window, window_matrix.GetRow( i)) / ( cos_sqnorm + self_scalar_product( i))));
            if( period < 2.0001)
            {
              mx = std::max( double( mx), math::Absolute( linal::ScalarProduct( sin_window, window_matrix.GetRow( i)) / ( sin_sqnorm + self_scalar_product( i))));
            }
          }

          if( mx > MAX_OVERLAP)
          {
            continue;
          }
        }

        // copy cosine and sine window
        std::copy( cos_window.Begin(), cos_window.End(), window_matrix[ row_number++]);
        std::copy( sin_window.Begin(), sin_window.End(), window_matrix[ row_number++]);
        new_target_periods.PushBack( period);
        self_scalar_product.PushBack( cos_sqnorm);
        self_scalar_product.PushBack( sin_sqnorm);
      }
      if( row_number < window_matrix.GetNumberRows())
      {
        window_matrix = linal::Matrix< float>( row_number, window_size, window_matrix.Begin());
      }

      // create the return matrix
      WaveletsFrequenciesPtr return_value
      (
        new storage::Pair
        <
          linal::Matrix< float>,
          storage::Vector< float>
        >
        (
          window_matrix,
          new_target_periods
        )
      );

      // the following code allows for more exact determination of frequencies, however, it comes at the cost of
      // substantially higher computation and a substantial risk of ill-conditioned matrices.  Likewise, it should not
      // normally be used
      // window_matrix.MoorePenroseInverse();
      // window_matrix.Transpose();

      if( !CACHEABLE)
      {
        // not cacheable, return the own pointer directly
        return return_value;
      }

      // insert the pointer into the map; return a non-owning copy
      MapType::iterator itr_return( s_wavelets.Insert( std::make_pair( query_string, return_value)).first);
      s_wavelets_mutex.Unlock();
      // return a non-owning pointer to the matrix and periods
      return MapType::mapped_type( &*itr_return->second, false);
    }

    //! @brief get unweighted sin/cosine vectors
    //! @param ALIGNMENT alignment to use; determines meaning of each index in the vector
    //! @param SIZE radius of the window to use
    //! @param PERIOD period desired
    //! @return VectorND2 to linal::VectorConstRef to Cosine (first) and Sine (second) unweighted components
    storage::Pair< linal::VectorConstReference< float>, linal::VectorConstReference< float> >
      FourierAnalysisWindowCreator::GetUnweightedCosSinCoefficients
      (
        const WindowAlignmentEnum &ALIGNMENT,
        const size_t &SIZE,
        const float &PERIOD
      )
    {
      //! Mutex to protect access to s_windows
      static sched::Mutex s_windows_mutex;

      //! Map from initializer string/type to the unweighted cosine/sine windows
      //! and target periods
      typedef storage::Map
      <
        std::string, //! WindowSize;Alignment;Weighting;MaxPeriods;MaxOverlap;TargetPeriods
        storage::Pair< linal::Vector< float>, linal::Vector< float> >
      > MapType;

      static MapType s_windows;

      // if the window matrix is to be cacheable, first look for it in the map
      const std::string query_string( ALIGNMENT.GetString() + ';' + util::Format()( SIZE) + ';' + util::Format()( PERIOD));
      s_windows_mutex.Lock();
      MapType::iterator itr_map( s_windows.Find( query_string));
      if( itr_map != s_windows.End())
      {
        s_windows_mutex.Unlock();
        // return a non-owning pointer to the mapped-type
        return
          storage::Pair< linal::VectorConstReference< float>, linal::VectorConstReference< float> >
          (
            linal::VectorConstReference< float>( itr_map->second.First()),
            linal::VectorConstReference< float>( itr_map->second.Second())
          );
      }

      // create the window ids mapping
      Window< char, float> window( false, Constants< char, float>( 1.0));
      window.SetRadius( SIZE);
      window.SetAlignment( ALIGNMENT);
      storage::Vector< int> window_ids( window.GetRelativeIndices());
      const size_t total_window_size( window_ids.GetSize());
      storage::Pair< linal::Vector< float>, linal::Vector< float> >
        cos_sine_window
        (
          linal::Vector< float>( total_window_size, 0.0f),
          linal::Vector< float>( total_window_size, 0.0f)
        );
      linal::Vector< float> &cos_window( cos_sine_window.First());
      linal::Vector< float> &sin_window( cos_sine_window.Second());
      const float frequency( 2.0f * math::g_Pi / PERIOD);
      if( PERIOD > 2.0001)
      {
        for( size_t i( 0); i < total_window_size; ++i)
        {
          const float position( frequency * window_ids( i));
          cos_window( i) = cos( position);
          sin_window( i) = sin( position);
        }
      }
      else
      {
        for( size_t i( 0); i < total_window_size; ++i)
        {
          cos_window( i) = cos( frequency * window_ids( i));
        }
      }

      itr_map = s_windows.Insert( std::make_pair( query_string, cos_sine_window)).first;
      s_windows_mutex.Unlock();
      // return a non-owning pointer to the mapped-type
      return
        storage::Pair< linal::VectorConstReference< float>, linal::VectorConstReference< float> >
        (
          linal::VectorConstReference< float>( itr_map->second.First()),
          linal::VectorConstReference< float>( itr_map->second.Second())
        );
    }

  } // namespace descriptor
} // namespace bcl
