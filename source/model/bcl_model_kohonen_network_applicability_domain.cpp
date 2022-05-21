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
#include "model/bcl_model_kohonen_network_applicability_domain.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_tertiary_function_job_with_data.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> KohonenNetworkApplicabilityDomain::s_Instance
    (
      GetObjectInstances().AddInstance( new KohonenNetworkApplicabilityDomain())
    );

    //! @brief constructor from parameters
    //! @param NETWORK the network to copy
    //! @param SHARE_DISTANCE_METRIC Whether to share the distance metric across all nodes
    //! @param MODEL_RETRIEVER retriever for the models
    KohonenNetworkApplicabilityDomain::KohonenNetworkApplicabilityDomain
    (
      const KohonenNetworkAverage &NETWORK,
      const bool &SHARE_DISTANCE_METRIC,
      const util::Implementation< RetrieveInterface> &MODEL_RETRIEVER
    ) :
      KohonenNetworkAverage( NETWORK),
      m_NodesSharedDistanceMetric( SHARE_DISTANCE_METRIC),
      m_ModelRetriever( MODEL_RETRIEVER)
    {
      ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief setup all internal splines
    //! @param NORMALIZED_DATA all the training data, already normalized
    //! @param PREVIOUS_WINNERS previous rounds node winners, will be updated
    void KohonenNetworkApplicabilityDomain::SetupSplines
    (
      const FeatureDataSetInterface< float> &NORMALIZED_DATA,
      const storage::Vector< size_t> &PREVIOUS_WINNERS,
      const bool &SHARE_DISTANCE_METRIC
    )
    {
      m_NodesSharedDistanceMetric = SHARE_DISTANCE_METRIC;
      linal::MatrixConstReference< float> matrix( NORMALIZED_DATA.GetMatrix());

      if( !m_NodesSharedDistanceMetric)
      {
        // different distance metric for each node of the map
        const size_t n_nodes( KohonenNetworkAverage::m_Codebook.GetSize());

        // create a vector of references, one per node in the kohonen map, with all the closest points to that node
        storage::Vector< storage::Vector< linal::VectorConstReference< float> > > closest_neighbor_map( n_nodes);

        for( size_t i( 0), n_features( NORMALIZED_DATA.GetNumberFeatures()); i < n_features; ++i)
        {
          closest_neighbor_map( PREVIOUS_WINNERS( i)).PushBack( matrix.GetRow( i));
        }
        m_Splines.Resize( n_nodes);
        for( size_t node_id( 0); node_id < n_nodes; ++node_id)
        {
          m_Splines( node_id).Resize( 1);
          m_Splines( node_id)( 0) =
           MapData
           (
             closest_neighbor_map( node_id),
             KohonenNetworkAverage::m_Codebook( node_id).GetFeatureVector()
           );
        }
      }
      else
      {
        m_Splines.Resize( 1);
        m_Splines( 0).Resize( 1);

        const size_t number_points( NORMALIZED_DATA.GetNumberFeatures());
        const storage::Vector< KohonenNode> &nodes( KohonenNetworkAverage::m_Codebook);

        // compute all distances
        linal::Vector< double> distances( number_points);
        for( size_t point_id( 0); point_id < number_points; ++point_id)
        {
          distances( point_id)
            = linal::Distance( matrix.GetRow( point_id), nodes( PREVIOUS_WINNERS( point_id)).GetFeatureVector());
        }

        // sort the vector
        std::sort( distances.Begin(), distances.End());

        // extract some representative points from the series 0, 1/100, 1/8, 1/4, 3/8, 1/2, 3/4, 7/8, 15/16, 31/32, 63/64, 127/128, 1
        // interpolate between closest points. These points are chosen with a heavy tail because distances from kohonen
        // nodes have a high slope towards the end of this distribution
        int n_interpolation_points( 0);
        for( size_t x( number_points); x; x >>= 1)
        {
          ++n_interpolation_points;
        }
        n_interpolation_points = std::max( n_interpolation_points, 2) + 5;
        linal::Vector< double> cumul_dist_function( n_interpolation_points, double( 0.0));
        linal::Vector< double> rep_distances( n_interpolation_points, double( 0.0));
        for( int i( 7), denominator( 4); i < n_interpolation_points; ++i, denominator <<= 1)
        {
          cumul_dist_function( i) = double( denominator - 1) / double( denominator);
        }
        cumul_dist_function( 0) = 0.0;
        cumul_dist_function( 1) = 0.01;
        cumul_dist_function( 2) = 0.125;
        cumul_dist_function( 3) = 0.25;
        cumul_dist_function( 4) = 0.375;
        cumul_dist_function( 5) = 0.5;
        cumul_dist_function( 6) = 0.625;

        cumul_dist_function( n_interpolation_points - 1) = 1.0;
        if( !distances.IsEmpty())
        {
          rep_distances( n_interpolation_points - 1) = distances.Last();
        }
        else
        {
          rep_distances( n_interpolation_points - 1) = 0.01;
        }
        for( int i( 1); i < n_interpolation_points - 1; ++i)
        {
          const size_t lower_bound( cumul_dist_function( i) * number_points);
          const double divided_lo( double( lower_bound) / double( number_points));
          rep_distances( i) = distances( lower_bound);
          if( divided_lo == lower_bound || lower_bound + 1 == number_points)
          {
            continue;
          }
          const size_t upper_bound( lower_bound + 1);
          const double dy( 1.0 / double( number_points));
          const double dx( distances( upper_bound) - distances( lower_bound));
          // interpolation based on nearest points
          rep_distances( i) += ( cumul_dist_function( i) - divided_lo) * dx / dy;
        }

        // train the spline
        m_Splines( 0)( 0).Train( rep_distances, cumul_dist_function, 0, 0);
      }
    } // Calibrate

    //! @brief setup all internal splines
    //! @param SPLINES splines, one for each node, or just one (for whole network)
    void KohonenNetworkApplicabilityDomain::SetupSplines
    (
      const storage::Vector< storage::Vector< math::CubicSplineDamped> > &SPLINES
    )
    {
      m_Splines = SPLINES;
      m_NodesSharedDistanceMetric = m_Splines.GetSize() == size_t( 1);
    }

    //! @brief set the model that this object should call to obtain the un-transformed prediction
    //! @param MODELS the models to call to compute the underlying prediction
    void KohonenNetworkApplicabilityDomain::SetupModels
    (
      const util::Implementation< RetrieveInterface> &MODEL_RETRIEVER
    )
    {
      m_ModelRetriever = MODEL_RETRIEVER;
    }

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURES rescaled features
    //! @return predicted result vector using a model
    FeatureDataSet< float> KohonenNetworkApplicabilityDomain::PredictWithoutRescaling
    (
      const FeatureDataSetInterface< float> &FEATURES
    ) const
    {
      if( m_ModelRetriever.IsDefined() && m_Models.IsEmpty())
      {
        m_Models = m_ModelRetriever->RetrieveEnsemble();
      }
      storage::Vector< size_t> winning_nodes;
      storage::Vector< sched::Mutex> muteces;
      GetWinningNodeIndices
      (
        FEATURES,
        math::Range< size_t>( 0, FEATURES.GetNumberFeatures()),
        storage::Vector< size_t>(),
        winning_nodes,
        muteces
      );

      const linal::MatrixConstReference< float> features_mat( FEATURES.GetMatrix());

      // create a results matrix
      FeatureDataSet< float> results_fds
      (
        FEATURES.GetNumberFeatures(),
        GetNumberOutputs(),
        float( 0.0)
      );
      linal::Matrix< float> &results( results_fds.GetRawMatrix());

      if( !m_ModelRetriever.IsDefined())
      {
        const bool distance_based_single_result
        (
          GetNumberOutputs() > size_t( 1) && m_Splines( 0).GetSize() == size_t( 1)
        );
        for
        (
          size_t result_n( 0), n_results( distance_based_single_result ? size_t( 1) : GetNumberOutputs());
          result_n < n_results;
          ++result_n
        )
        {
          for
          (
            size_t feature_id( 0), number_features( FEATURES.GetNumberFeatures());
            feature_id < number_features;
            ++feature_id
          )
          {
            // get the feature from the dataset
            const linal::VectorConstReference< float> feature( features_mat.GetRow( feature_id));

            if( !feature.IsDefined())
            {
              if( distance_based_single_result)
              {
                // all results will be equal
                results.GetRow( feature_id) = 0.0;
              }
              else
              {
                results( feature_id, result_n) = 0.0;
              }
              continue;
            }

            // get the corresponding node in the kohonen map
            const linal::VectorConstInterface< float> &node( m_Codebook( winning_nodes( feature_id)).GetFeatureVector());

            // get the corresponding spline
            const math::CubicSplineDamped &spline
            (
              m_Splines( m_NodesSharedDistanceMetric ? 0 : winning_nodes( feature_id))( result_n)
            );

            // distance between kohonen map node and feature
            const float distance( linal::Distance( node, feature));

            // get spline value. Handle fencepost cases too
            const float spline_fract
            (
              distance > spline.GetXValues().Last()
              ? distance / spline.GetXValues().Last()
              : float( std::max( spline( distance), double( 0.0)))
            );

            if( distance_based_single_result)
            {
              // all results will be equal
              results.GetRow( feature_id) = spline_fract;
            }
            else
            {
              results( feature_id, result_n) = spline_fract;
            }
          }
        }
      }
      else
      {
        storage::Vector< linal::Matrix< float> > all_results( m_Models.GetSize());
        // iterate through models, average results
        // create a results matrix
        math::RunningAverageSD< linal::Matrix< float> > ave_results, ave_local_ppv;
        FeatureDataSet< float> fds_descaled( FEATURES);
        fds_descaled.DeScale();
        for( size_t model_n( 0), n_models( m_Models.GetSize()); model_n < n_models; ++model_n)
        {
          linal::Matrix< float> output( ( *m_Models( model_n))( fds_descaled).GetMatrix());
          ave_results += output;
          fds_descaled.DeScale();
          for( size_t result_n( 0), n_results( GetNumberOutputs()); result_n < n_results; ++result_n)
          {
            for
            (
              size_t feature_id( 0), number_features( FEATURES.GetNumberFeatures());
              feature_id < number_features;
              ++feature_id
            )
            {

              // get the corresponding spline
              const math::CubicSplineDamped &spline
              (
                m_Splines( m_NodesSharedDistanceMetric ? 0 : winning_nodes( feature_id))( result_n)
              );

              // distance between kohonen map node and feature
              const float distance( output( feature_id, result_n));

              // get spline value. Handle fencepost cases too
              const float spline_fract
              (
                distance > spline.GetXValues().Last()
                ? spline.GetYValues().Last()
                : distance < spline.GetXValues().First()
                  ? spline.GetYValues().First()
                  : spline( distance)
              );
              output( feature_id, result_n) = spline_fract;
            }
          }
          all_results( model_n) = output;
          ave_local_ppv += output;
        }
        const linal::Matrix< float> &ave_result( ave_results.GetAverage());
        for( size_t result_n( 0), n_results( GetNumberOutputs()); result_n < n_results; ++result_n)
        {
          for
          (
            size_t feature_id( 0), number_features( FEATURES.GetNumberFeatures());
            feature_id < number_features;
            ++feature_id
          )
          {

            // get the corresponding spline
            const math::CubicSplineDamped &spline
            (
              m_Splines( m_NodesSharedDistanceMetric ? 0 : winning_nodes( feature_id))( result_n)
            );

            // distance between kohonen map node and feature
            const float distance( ave_result( feature_id, result_n));

            // get spline value. Handle fencepost cases too
            const float spline_fract
            (
              distance > spline.GetXValues().Last()
              ? spline.GetYValues().Last()
              : distance < spline.GetXValues().First()
                ? spline.GetYValues().First()
                : spline( distance)
            );

            results( feature_id, result_n) = spline_fract;
          }
        }
      }

      return FeatureDataSet< float>( results);
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> KohonenNetworkApplicabilityDomain::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_RescaleInput);
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KohonenNetworkApplicabilityDomain::Read( std::istream &ISTREAM)
    {
      KohonenNetworkAverage::Read( ISTREAM);
      io::Serialize::Read( m_Splines, ISTREAM);
      io::Serialize::Read( m_NodesSharedDistanceMetric, ISTREAM);
      io::Serialize::Read( m_ModelRetriever, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &KohonenNetworkApplicabilityDomain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      KohonenNetworkAverage::Write( OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Splines, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NodesSharedDistanceMetric, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ModelRetriever, OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief Trains the cubic spline based on distances to NODE
    //! @param MAPPED_POINTS all points mapped to NODE
    math::CubicSplineDamped KohonenNetworkApplicabilityDomain::MapData
    (
      const storage::Vector< linal::VectorConstReference< float> > &MAPPED_POINTS,
      const linal::VectorConstInterface< float> &NODE
    )
    {
      const size_t number_points( MAPPED_POINTS.GetSize());
      // compute all distances
      linal::Vector< double> distances( number_points);
      for( size_t point_id( 0); point_id < number_points; ++point_id)
      {
        distances( point_id) = linal::Distance( MAPPED_POINTS( point_id), NODE);
      }

      // sort the vector
      std::sort( distances.Begin(), distances.End());

      // extract some representative points from the series 0, 1/2, 3/4, 7/8, 15/16, 31/32, 63/64, 127/128, 1
      // interpolate between closest points. These points are chosen with a heavy tail because distances from kohonen
      // nodes have a high slope towards the end of this distribution
      int n_interpolation_points( 0);
      for( size_t x( number_points); x; x >>= 1)
      {
        ++n_interpolation_points;
      }
      n_interpolation_points = std::max( n_interpolation_points, 2) + 5;

      linal::Vector< double> cumul_dist_function( n_interpolation_points, double( 0.0));
      linal::Vector< double> rep_distances( n_interpolation_points, double( 0.0));

      cumul_dist_function( 0) = 0.0;
      cumul_dist_function( 1) = 0.01;
      cumul_dist_function( 2) = 0.125;
      cumul_dist_function( 3) = 0.25;
      cumul_dist_function( 4) = 0.375;
      cumul_dist_function( 5) = 0.5;
      cumul_dist_function( 6) = 0.625;
      for( int i( 7), denominator( 4); i < n_interpolation_points; ++i, denominator <<= 1)
      {
        cumul_dist_function( i) = double( denominator - 1) / double( denominator);
      }
      cumul_dist_function( n_interpolation_points - 1) = 1.0;
      if( !distances.IsEmpty())
      {
        rep_distances( n_interpolation_points - 1) = distances.Last();
      }
      else
      {
        rep_distances( n_interpolation_points - 1) = 0.01;
      }
      for( int i( 1); i < n_interpolation_points - 1; ++i)
      {
        const size_t lower_bound( cumul_dist_function( i) * number_points);
        const double divided_lo( double( lower_bound) / double( number_points));
        rep_distances( i) = distances( lower_bound);
        if( divided_lo == lower_bound || lower_bound + 1 == number_points)
        {
          continue;
        }
        const size_t upper_bound( lower_bound + 1);
        const double dy( 1.0 / double( number_points));
        const double dx( distances( upper_bound) - distances( lower_bound));
        // interpolation based on nearest points
        rep_distances( i) += ( cumul_dist_function( i) - divided_lo) * dx / dy;
      }

      // train the spline
      math::CubicSplineDamped spline;
      spline.Train
      (
        rep_distances,
        cumul_dist_function,
        0,
        0
      );
      return spline;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool KohonenNetworkApplicabilityDomain::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( m_ModelRetriever.IsDefined())
      {
        m_Models = m_ModelRetriever->RetrieveEnsemble();
      }
      return true;
    }

  } // namespace model
} // namespace bcl
