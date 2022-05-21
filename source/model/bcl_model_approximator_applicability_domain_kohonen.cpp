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
#include "model/bcl_model_approximator_applicability_domain_kohonen.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_min_max.h"
#include "model/bcl_model_kohonen_network_applicability_domain.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "model/bcl_model_retrieve_interface.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ApproximatorApplicabilityDomainKohonen::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorApplicabilityDomainKohonen( true))
    );
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ApproximatorApplicabilityDomainKohonen::s_ContingencyInstance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorApplicabilityDomainKohonen( false))
    );

    //! @brief default constructor
    ApproximatorApplicabilityDomainKohonen::ApproximatorApplicabilityDomainKohonen( const bool &USE_DISTANCE_METRIC) :
      ApproximatorKohonenNetwork(),
      m_ShareDistanceMetric( USE_DISTANCE_METRIC),
      m_UseDistanceMetric( USE_DISTANCE_METRIC),
      m_Monotonize( true)
    {
      ApproximatorKohonenNetwork::m_NoTrack = true;
    }

    //! @brief constructor from all necessary parameters
    //! @param MAP_DIMENSIONS dimensions of the map
    //! @param INITIAL_LENGTH how many iterations to train for (if applicable).
    //! @param INITAL_RADIUS the initial neighborhood radius
    //! @param OBJECTIVE_FUNCTION the objective function from the approximator framework
    //! @param UPDATE_EVERY_NTH_FEATURE update the nodes after seeing this many features
    //! @param NEIGHBOR_KERNEL the neighbor kernel type
    ApproximatorApplicabilityDomainKohonen::ApproximatorApplicabilityDomainKohonen
    (
      const linal::Vector< double> &MAP_DIMENSIONS,
      const size_t &INITIAL_LENGTH,
      const float &INITIAL_RADIUS,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
      const size_t &UPDATE_EVERY_NTH_FEATURE,
      const NeighborKernel &NEIGHBOR_KERNEL,
      const util::Implementation< RetrieveInterface> &RETREIVER
    ) :
      ApproximatorKohonenNetwork
      (
        MAP_DIMENSIONS,
        INITIAL_LENGTH,
        INITIAL_RADIUS,
        OBJECTIVE_FUNCTION,
        UPDATE_EVERY_NTH_FEATURE,
        NEIGHBOR_KERNEL
      ),
      m_Retriever( RETREIVER),
      m_ShareDistanceMetric( true),
      m_UseDistanceMetric( true),
      m_Monotonize( true)
    {
      ApproximatorKohonenNetwork::m_NoTrack = true;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorApplicabilityDomainKohonen::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorApplicabilityDomainKohonen::GetAlias() const
    {
      static const std::string s_alias( "ApplicabilityDomainKohonen"), s_con_alias( "ApplicabilityMeasureKohonen");
      return m_UseDistanceMetric ? s_alias : s_con_alias;
    }

    //! @brief set training data set for a specific iterate in approximator framework
    //! @param DATA training data set
    void ApproximatorApplicabilityDomainKohonen::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      m_FinalModel.Reset();
      this->GetTracker() = opti::Tracker< util::ShPtr< Interface>, float>( opti::e_SmallerEqualIsBetter);
      if( m_Retriever.IsDefined())
      {
        BCL_Assert( !m_Retriever->GetAllKeys().IsEmpty(), "No models in: " + m_Retriever->GetString());
        BCL_Assert
        (
          DATA->GetFeatures().GetFeatureLabelSet().IsDefined(),
          "Features are required for model applicability, otherwise predictions could not be made"
        );
        BCL_Assert
        (
          m_Retriever->RetrieveUniqueDescriptor().GetArguments()
          == DATA->GetFeatures().GetFeatureLabelSet()->GetMemberLabels().InternalData(),
          "The applicability domain kohonen map is presently only trainable if it has the  be trained on the same set "
          "of features as the underlying model"
        );
        m_ResultDataset = m_Retriever->ReadMergedIndependentPredictions();
        BCL_Assert
        (
          m_ResultDataset->GetFeatureSize() == DATA->GetResultSize(),
          "Different # of results for model than in given dataset!"
        );
        BCL_Assert
        (
          math::EqualWithinTolerance
          (
            DATA->GetResults().GetMatrix(),
            m_ResultDataset->GetResults().GetMatrix()
          ),
          "Independent set results did not match dataset results!"
        );
      }
      ApproximatorKohonenNetwork::SetTrainingData( DATA);

      // compute average distance between gridpoints to give user idea of how stable the map is
      math::RunningAverage< double> ave_distance_between_gridpoints;
      for
      (
        size_t i( 0), n_gridpoints( ApproximatorKohonenNetwork::m_Network.GetCodeBook().GetSize());
        i < n_gridpoints;
        ++i
      )
      {
        const linal::VectorConstInterface< float> &position_i
        (
          ApproximatorKohonenNetwork::m_Network.GetCodeBook()( i).GetPosition()
        );
        for( size_t j( 0); j < i; ++j)
        {
          ave_distance_between_gridpoints +=
            linal::Distance( position_i, ApproximatorKohonenNetwork::m_Network.GetCodeBook()( j).GetPosition());
        }
        // on-diagonal
        ave_distance_between_gridpoints.AddWeightedObservation( 0.0, 0.5);
      }
      m_AverageDistanceBetweenGridpoints = ave_distance_between_gridpoints.GetAverage();

      if( !this->ShouldContinue())
      {
        MakeFinalModel();
        this->GetTracker().Track( GetCurrentApproximation(), opti::e_Improved);
      }
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorApplicabilityDomainKohonen::GetCurrentModel() const
    {
      return GetCurrentApproximation()->First();
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > ApproximatorApplicabilityDomainKohonen::GetCurrentApproximation() const
    {
      if( m_FinalModel.IsDefined())
      {
        return
          util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
          (
            new storage::Pair< util::ShPtr< Interface>, float>
            (
              m_FinalModel,
              m_AverageDisplacement / m_AverageDistanceBetweenGridpoints
            )
          );
      }
      util::ShPtr< KohonenNetworkApplicabilityDomain> app_domain
      (
        new KohonenNetworkApplicabilityDomain( m_Network, m_ShareDistanceMetric, m_Retriever)
      );
      util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > approximation
      (
        new storage::Pair< util::ShPtr< Interface>, float>
        (
          app_domain,
          m_AverageDisplacement / m_AverageDistanceBetweenGridpoints
        )
      );

      // more iterations ahead; just setup the splines without regards to any metric
      app_domain->SetupSplines
      (
        *m_TrainingData->GetFeaturesPtr(),
        ApproximatorKohonenNetwork::m_LastClosestNodes,
        m_ShareDistanceMetric
      );

      return approximation;
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorApplicabilityDomainKohonen::Next()
    {
      m_AverageDisplacement = 0;
      storage::Vector< size_t> last_round_assigned( ApproximatorKohonenNetwork::m_LastClosestNodes);
      ApproximatorKohonenNetwork::Next();
      math::RunningAverage< double> ave_distance_moved, ave_dist_between_points;
      for( size_t i( 0), n_data( last_round_assigned.GetSize()); i < n_data; ++i)
      {
        const size_t prev_node( last_round_assigned( i));
        const size_t curr_node( ApproximatorKohonenNetwork::m_LastClosestNodes( i));

        if( prev_node != curr_node)
        {
          ave_distance_moved +=
            linal::Distance
            (
              m_Network.GetCodeBook()( prev_node).GetPosition(),
              m_Network.GetCodeBook()( curr_node).GetPosition()
            );
        }
        else
        {
          ave_distance_moved += 0.0;
        }
      }
      m_AverageDisplacement = ave_distance_moved.GetAverage();
      BCL_MessageStd
      (
        util::Format()( m_AverageDisplacement) + " average distance each feature moved on map, "
        + util::Format().FFP( 2)
          (
            100.0 - std::min( 100.0, m_AverageDisplacement / m_AverageDistanceBetweenGridpoints * 100.0)
          )
        + " % solidified"
      );
      this->GetTracker().Track( GetCurrentApproximation());
      if( !this->ShouldContinue())
      {
        MakeFinalModel();
        this->GetTracker().Track( GetCurrentApproximation(), opti::e_Improved);
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorApplicabilityDomainKohonen::GetSerializer() const
    {
      io::Serializer parameters( ApproximatorKohonenNetwork::GetSerializer());
      parameters.SetClassDescription
      (
        "A kohonen-map based implementation to detect whether a point is within the applicability domain of a model. "
        "Inspired by Assessment of applicability domain for multivariate counter-propagation artificial neural network "
        "predictive models by minimum Euclidean distance space analysis: A case study"
      );

      // fix the objective function to Constant
      parameters.SetDefault( "objective function", "Constant(direction=SmallerIsBetter)");

      if( m_UseDistanceMetric)
      {
        parameters.AddInitializer
        (
          "share distance metric",
          "If true, all nodes will use the same spline for computing applicability. This implies an assumption that the "
          "model has the most difficulty predicting things far from any node center, regardless of which node center it "
          "is. If False, all nodes will have their own distance metric, which is valid if the model is capable of "
          "distinguishing between classes of features (e.g. a kohonen map itself).",
          io::Serialization::GetAgent( &m_ShareDistanceMetric),
          "True"
        );
      }
      else
      {
        parameters.AddInitializer
        (
          "model",
          "The model to use to make predictions.",
          io::Serialization::GetAgent( &m_Retriever)
        );
        parameters.AddInitializer
        (
          "monotonic",
          "True if monotonicity should be enforced",
          io::Serialization::GetAgent( &m_Monotonize),
          "True"
        );
        parameters.AddInitializer
        (
          "measure",
          "Contingency matrix measure to compute spline over.",
          io::Serialization::GetAgent( &m_Measure)
        );
      }
      return parameters;
    }

    //! @brief produce the final model
    //! Updates m_FinalModel
    void ApproximatorApplicabilityDomainKohonen::MakeFinalModel()
    {
      m_FinalModel =
        util::ShPtr< KohonenNetworkApplicabilityDomain>
        (
          new KohonenNetworkApplicabilityDomain( m_Network, m_ShareDistanceMetric, m_Retriever)
        );

      if( m_UseDistanceMetric)
      {
        // fix the current approximation
        storage::Vector< sched::Mutex> empty_mutex_vector;
        ApproximatorKohonenNetwork::m_Network.GetWinningNodeIndices
        (
          *ApproximatorKohonenNetwork::m_TrainingData->GetFeaturesPtr(),
          math::Range< size_t>( 0, ApproximatorKohonenNetwork::m_Schedule.GetOrder().GetSize()),
          storage::Vector< size_t>(),
          ApproximatorKohonenNetwork::m_LastClosestNodes,
          empty_mutex_vector
        );
        m_FinalModel->SetupSplines
        (
          *m_TrainingData->GetFeaturesPtr(),
          ApproximatorKohonenNetwork::m_LastClosestNodes,
          m_ShareDistanceMetric
        );
        return;
      }

      BCL_Assert( m_Retriever.IsDefined(), "Must have a defined set of models to compute the applicability metrics for!");

      linal::MatrixConstReference< float> independent_results( m_ResultDataset->GetFeaturesReference());
      linal::MatrixConstReference< float> true_results( m_ResultDataset->GetResultsReference());
      const size_t n_results( independent_results.GetNumberCols());

      // winning nodes selection
      storage::Vector< size_t> last_closest_nodes;
      if( !m_ShareDistanceMetric)
      {
        storage::Vector< sched::Mutex> empty_mutex_vector;
        ApproximatorKohonenNetwork::m_Network.GetWinningNodeIndices
        (
          *m_TrainingData->GetFeaturesPtr(),
          math::Range< size_t>( 0, m_TrainingData->GetSize()),
          storage::Vector< size_t>(),
          last_closest_nodes,
          empty_mutex_vector
        );
      }

      const math::ContingencyMatrixMeasures measure
      (
        math::ContingencyMatrixMeasures::MeasureEnum( m_Measure->GetAlias())
      );

      const size_t n_nodes( this->m_Network.GetCodeBook().GetSize());
      // create ROC curves for each result, optionally per-node
      storage::Vector< storage::Vector< math::CubicSplineDamped> > splines
      (
        m_ShareDistanceMetric ? size_t( 1) : n_nodes,
        storage::Vector< math::CubicSplineDamped>( n_results)
      );
      m_FinalModel->SetupModels( m_Retriever);
      util::Implementation< ObjectiveFunctionInterface> objective
      (
        m_Retriever->RetrieveCommonCVInfo().GetObjective()
      );
      for( size_t result_n( 0); result_n < n_results; ++result_n)
      {
        storage::List< storage::Pair< double, double> > test_results;
        for( size_t i( 0), dataset_n( m_TrainingData->GetSize()); i < dataset_n; ++i)
        {
          test_results.PushBack
          (
            storage::Pair< double, double>( independent_results( i, result_n), true_results( i, result_n))
          );
        }
        storage::Map< double, math::ContingencyMatrix> roc_map
        (
          math::ROCCurve
          (
            test_results,
            objective->GetThreshold(),
            objective->GetRankingParity()
          ).ToMap()
        );

        const bool computing_local_ppv
        (
          m_Measure->GetAlias()
          == math::ContingencyMatrixMeasures::GetMeasureName( math::ContingencyMatrixMeasures::e_LocalPPV)
        );

        math::PiecewiseFunction local_ppv;
        if( computing_local_ppv)
        {
          local_ppv = math::ROCCurve
                      (
                        test_results,
                        objective->GetThreshold(),
                        objective->GetRankingParity()
                      ).GetLocalPPVCurve();
        }

        if( m_ShareDistanceMetric)
        {
          storage::Vector< storage::Pair< double, double> > dataset_roc_points;
          dataset_roc_points.AllocateMemory( m_TrainingData->GetSize());
          for( size_t i( 0), dataset_n( m_TrainingData->GetSize()); i < dataset_n; ++i)
          {
            dataset_roc_points.PushBack
            (
              storage::Pair< double, double>
              (
                independent_results( i, result_n),
                computing_local_ppv
                ? local_ppv( independent_results( i, result_n))
                : measure( roc_map[ independent_results( i, result_n)])
              )
            );
          }
          // sort ROC curve points
          dataset_roc_points.Sort( std::less< storage::Pair< double, double> >());

          math::RunningMinMax< double> min_max_y;
          for( size_t i( 0), dataset_n( m_TrainingData->GetSize()); i < dataset_n; ++i)
          {
            min_max_y += dataset_roc_points( i).Second();
          }
          const double gradation( min_max_y.GetRange() / std::min( size_t( 100), m_TrainingData->GetSize()));
          storage::Vector< double> selected_x, selected_y;
          // only select ROC points that are significantly different in PPV or whatever measure
          if( !m_Monotonize)
          {
            double last_y( -1000.0);
            for( size_t i( 0), dataset_n( m_TrainingData->GetSize()); i < dataset_n; ++i)
            {
              if( math::Absolute( dataset_roc_points( i).Second() - last_y) > gradation)
              {
                selected_x.PushBack( dataset_roc_points( i).First());
                selected_y.PushBack( dataset_roc_points( i).Second());
                last_y = dataset_roc_points( i).Second();
              }
            }
          }
          else if( objective->GetRankingParity() == measure.GetOptimizationParity())
          {
            double last_y( -1000.0);
            for( size_t i( 0), dataset_n( m_TrainingData->GetSize()); i < dataset_n; ++i)
            {
              if( ( dataset_roc_points( i).Second() - last_y) > gradation)
              {
                selected_x.PushBack( dataset_roc_points( i).First());
                selected_y.PushBack( dataset_roc_points( i).Second());
                last_y = dataset_roc_points( i).Second();
              }
            }
          }
          else
          {
            double last_y( 1000.0);
            for( size_t i( 0), dataset_n( m_TrainingData->GetSize()); i < dataset_n; ++i)
            {
              if( ( last_y - dataset_roc_points( i).Second()) > gradation)
              {
                selected_x.PushBack( dataset_roc_points( i).First());
                selected_y.PushBack( dataset_roc_points( i).Second());
                last_y = dataset_roc_points( i).Second();
              }
            }
          }
          splines( 0)( result_n).Train
          (
            linal::Vector< double>( selected_x.Begin(), selected_x.End()),
            linal::Vector< double>( selected_y.Begin(), selected_y.End())
          );
        }
        else
        {
          storage::Vector< storage::Vector< storage::Pair< double, double> > > dataset_roc_points
          (
            n_nodes,
            storage::Vector< storage::Pair< double, double> >
            (
              size_t( 1),
              storage::Pair< double, double>
              (
                roc_map.Begin()->first,
                computing_local_ppv
                ? local_ppv( roc_map.Begin()->first)
                : ( *m_Measure)( roc_map[ roc_map.Begin()->first])
              )
            )
          );
          storage::Vector< storage::Vector< storage::Pair< double, double> > > dataset_points
          (
            n_nodes
          );

          for( size_t i( 0), dataset_n( m_TrainingData->GetSize()); i < dataset_n; ++i)
          {
            dataset_points( last_closest_nodes( i)).PushBack
            (
              storage::Pair< double, double>
              (
                independent_results( i, result_n),
                true_results( i, result_n)
              )
            );
          }
          for( size_t node_n( 0); node_n < n_nodes; ++node_n)
          {
            // sort ROC curve points
            storage::Vector< storage::Pair< double, double> > &node_points( dataset_points( node_n));

            storage::List< storage::Pair< double, double> > node_points_list( node_points.Begin(), node_points.End());
            math::ROCCurve local_roc( node_points_list, objective->GetThreshold(), objective->GetRankingParity());
            math::PiecewiseFunction node_local_ppv( local_roc.GetLocalPPVCurve());
            storage::Map< double, math::ContingencyMatrix> local_roc_map( local_roc.ToMap());

            for( size_t i( 0), n_points( node_points.GetSize()); i < n_points; ++i)
            {
              dataset_roc_points( node_n).PushBack
              (
                storage::Pair< double, double>
                (
                  node_points( i).First(),
                  computing_local_ppv
                  ? node_local_ppv( node_points( i).First())
                  : ( *m_Measure)( local_roc_map[ node_points( i).First()])
                )
              );
            }
          }
          for( size_t node_n( 0); node_n < n_nodes; ++node_n)
          {
            // sort ROC curve points
            storage::Vector< storage::Pair< double, double> > &node_roc_points( dataset_roc_points( node_n));

            node_roc_points.Sort( std::less< storage::Pair< double, double> >());

            // only select ROC points that are significantly different in PPV or whatever measure
            storage::Vector< double> selected_x, selected_y;
            math::RunningMinMax< double> min_max_y;
            for( size_t i( 0), dataset_n( node_roc_points.GetSize()); i < dataset_n; ++i)
            {
              min_max_y += node_roc_points( i).Second();
            }
            const double gradation( ( min_max_y.GetRange() + 1.0e-5) / std::min( size_t( 1000), node_roc_points.GetSize()));
            // only select ROC points that are significantly different in PPV or whatever measure
            if( !m_Monotonize)
            {
              double last_y( -1000.0);
              for( size_t i( 0), dataset_n( node_roc_points.GetSize()); i < dataset_n; ++i)
              {
                if( math::Absolute( node_roc_points( i).Second() - last_y) > gradation)
                {
                  selected_x.PushBack( node_roc_points( i).First());
                  selected_y.PushBack( node_roc_points( i).Second());
                  last_y = node_roc_points( i).Second();
                }
              }
            }
            else if( objective->GetRankingParity() == measure.GetOptimizationParity())
            {
              double last_y( -1000.0);
              for( size_t i( 0), dataset_n( node_roc_points.GetSize()); i < dataset_n; ++i)
              {
                if( ( node_roc_points( i).Second() - last_y) > gradation)
                {
                  selected_x.PushBack( node_roc_points( i).First());
                  selected_y.PushBack( node_roc_points( i).Second());
                  last_y = node_roc_points( i).Second();
                }
              }
            }
            else
            {
              double last_y( 1000.0);
              for( size_t i( 0), dataset_n( node_roc_points.GetSize()); i < dataset_n; ++i)
              {
                if( ( last_y - node_roc_points( i).Second()) > gradation)
                {
                  selected_x.PushBack( node_roc_points( i).First());
                  selected_y.PushBack( node_roc_points( i).Second());
                  last_y = node_roc_points( i).Second();
                }
              }
            }
            splines( node_n)( result_n).Train
            (
              linal::Vector< double>( selected_x.Begin(), selected_x.End()),
              linal::Vector< double>( selected_y.Begin(), selected_y.End())
            );
          }
        }
      }
      m_FinalModel->SetupSplines( splines);
    }

  } // namespace model
} // namespace bcl
