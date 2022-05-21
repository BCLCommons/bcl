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
#include "model/bcl_model_neural_network.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_operations_interface.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "model/bcl_model_feature_data_reference.h"
#include "sched/bcl_sched_binary_function_job_with_data.h"
#include "sched/bcl_sched_scheduler_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! the default input range for neural network transfer functions
    const math::Range< float> NeuralNetwork::s_DefaultInputRange( -1, 1);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetwork::s_Instance
    (
      GetObjectInstances().AddInstance( new NeuralNetwork())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    NeuralNetwork::NeuralNetwork()
    {
    }

    //! @brief construct from all necessary parameters
    //! @param RESCALE_INPUT
    //! @param RESCALE_OUTPUT
    //! @param BIAS
    //! @param WEIGHT
    //! @param TRANSFER_FUNCTION
    NeuralNetwork::NeuralNetwork
    (
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT,
      const util::ShPtr< RescaleFeatureDataSet> &RESCALE_OUTPUT,
      const storage::Vector< linal::Vector< float> > &BIAS,
      const storage::Vector< linal::Matrix< float> > &WEIGHT,
      const util::Implementation< TransferFunctionInterface> &TRANSFER_FUNCTION
    ) :
      m_TransferFunction( TRANSFER_FUNCTION),
      m_RescaleInput( RESCALE_INPUT),
      m_RescaleOutput
      (
        RESCALE_OUTPUT.IsDefined()
        ? RESCALE_OUTPUT
        : util::ShPtr< RescaleFeatureDataSet>
          (
            new RescaleFeatureDataSet
            (
              linal::Matrix< float>( size_t( 1), WEIGHT.LastElement().GetNumberRows()),
              math::Range< float>(),
              RescaleFeatureDataSet::e_None
            )
          )
      ),
      m_Bias( BIAS),
      m_Weight( WEIGHT)
    {
      SetImplicitArchitecture();
    }

    //! @brief copy constructor
    //! @param NETWORK the nn that will be copied
    NeuralNetwork::NeuralNetwork( const NeuralNetwork &NETWORK) :
      m_TransferFunction( NETWORK.m_TransferFunction),
      m_RescaleInput( NETWORK.m_RescaleInput),
      m_RescaleOutput( NETWORK.m_RescaleOutput),
      m_Bias( NETWORK.m_Bias),
      m_Weight( NETWORK.m_Weight),
      m_Architecture(),
      m_Mutex(),
      m_Allocated( NETWORK.m_Allocated),
      m_Available( NETWORK.m_Available)
    {
      m_Mutex.Lock();
      m_Available.Reset();
      m_Allocated.Reset();
      m_Mutex.Unlock();

      SetImplicitArchitecture();
    }

    //! copy constructor
    NeuralNetwork *NeuralNetwork::Clone() const
    {
      return new NeuralNetwork( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetwork::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetwork::GetAlias() const
    {
      static const std::string s_Name( "NeuralNetwork");
      return s_Name;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns architecture of this neural network
    //! @return first entry is the number of input neurons
    //!         middle entries are the number of hidden neurons in each layer
    //!         last number in vector is the number of output neurons
    const storage::Vector< size_t> &NeuralNetwork::GetArchitecture() const
    {
      return m_Architecture;
    }

    //! @brief set the architecture architecture of this neural network
    //! @param ARCHITECTURE first entry is the number of input neurons
    //!         middle entries are the number of hidden neurons in each layer
    //!         last number in vector is the number of output neurons
    void NeuralNetwork::SetArchitecture( const storage::Vector< size_t> &ARCHITECTURE)
    {
      m_Architecture = ARCHITECTURE;
      m_Bias.Reset();
      m_Weight.Reset();

      for( size_t i( 1); i < ARCHITECTURE.GetSize(); ++i)
      {
        m_Bias.PushBack( linal::Vector< float>( ARCHITECTURE( i)));
        m_Weight.PushBack( linal::Matrix< float>( ARCHITECTURE( i), ARCHITECTURE( i - 1)));
      }
    }

    //! get number inputs
    size_t NeuralNetwork::GetNumberInputs() const
    {
      return m_Weight.FirstElement().GetNumberCols();
    }

    //! get number output neurons
    size_t NeuralNetwork::GetNumberOutputs() const
    {
      return m_Weight.LastElement().GetNumberRows();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void NeuralNetwork::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      if( !m_RescaleInput.IsDefined())
      {
        FEATURE.DeScale();
      }
      else if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_RescaleInput);
      }
    }

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> NeuralNetwork::PredictWithoutRescaling
    (
      const FeatureDataSetInterface< float> &FEATURE
    ) const
    {
      const size_t num_outs( m_Bias.LastElement().GetSize());
      const size_t num_ins( FEATURE.GetFeatureSize());
      const size_t num_points( FEATURE.GetNumberFeatures());
      BCL_Assert
      (
        num_ins == m_Weight( 0).GetNumberCols(),
        "Given " + util::Format()( num_ins) + " inputs; model expected " + util::Format()( m_Weight( 0).GetNumberCols())
      );
      linal::Matrix< float> result_matrix( num_points, num_outs);

      BCL_Assert
      (
        num_ins == m_Weight( 0).GetNumberCols(),
        "Given " + util::Format()( num_ins) + " inputs; model expected " + util::Format()( m_Weight( 0).GetNumberCols())
      );

      const size_t number_of_available_cpus( std::min( sched::GetNumberCPUs(), num_points));

      if( num_points < size_t( 12) || number_of_available_cpus == size_t( 1))
      {
        // no point in using threads if there are so few features
        PredictWithThread( FEATURE, result_matrix);

        return FeatureDataSet< float>( result_matrix, *m_RescaleOutput);
      }

      // data interval for jobs
      const size_t interval( num_points / number_of_available_cpus);

      // calculate how many threads will operate on interval + 1 features
      const size_t number_threads_with_interval_plus_one( num_points % number_of_available_cpus);

      // scheduler for managing concurrent jobs
      util::ShPtrVector< sched::JobInterface> schedule;

      // number of jobs
      size_t number_jobs( number_of_available_cpus);

      schedule.AllocateMemory( number_jobs);

      // create a vector to hold the results
      storage::List< linal::MatrixReference< float> > results;
      storage::List< FeatureDataReference< float> > inputs;

      // create ranges
      for
      (
        size_t job_number( 0), end_position( 0);
        job_number < number_jobs;
        ++job_number
      )
      {
        // save the previous end position as the current start position
        const size_t start_position( end_position);

        // add the interval
        end_position += interval;

        // add one if necessary
        if( job_number < number_threads_with_interval_plus_one)
        {
          ++end_position;
        }

        const size_t number_rows( end_position - start_position);
        results.PushBack( linal::MatrixReference< float>( number_rows, num_outs, result_matrix[ start_position]));
        inputs.PushBack(
          FeatureDataReference< float>
          (
            linal::MatrixConstReference< float>
            (
              number_rows,
              num_ins,
              FEATURE.GetMatrix()[ start_position]
            )
          )
        );

        // create and pushback jobs into scheduler
        schedule.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::BinaryFunctionJobWithData
            <
              const FeatureDataSetInterface< float>,
              linal::MatrixInterface< float>,
              void,
              NeuralNetwork
            >
            (
              1,
              *this,
              &NeuralNetwork::PredictWithThread,
              inputs.LastElement(),
              results.LastElement(),
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );

        // submit this particular job to start a thread
        sched::GetScheduler().RunJob( schedule.LastElement());
      }

      // join all jobs ( need to be finished)
      for( size_t job_id( 0); job_id < number_jobs; ++job_id)
      {
        sched::GetScheduler().Join( schedule( job_id));
      }

      return FeatureDataSet< float>( result_matrix, *m_RescaleOutput);
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> NeuralNetwork::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !m_RescaleInput.IsDefined() && !FEATURE.GetScaling().IsDefined())
      {
        // no rescaling to perform
        // handle neural networks that do not require rescaling
        // data is already rescaled
        return PredictWithoutRescaling( FEATURE).DeScale();
      }
      else if( !m_RescaleInput.IsDefined())
      {
        FeatureDataSet< float> feature( FEATURE);
        feature.DeScale();
        return PredictWithoutRescaling( feature).DeScale();
      }
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_RescaleInput);
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

    //! @brief predict result with model and compute input sensitivity using a rescaled feature vector
    //! @param FEATURE feature of interest, MUST be rescaled
    //! @return predicted result and input sensitivity (N_Descriptors X N_Results)
    storage::Pair< linal::Vector< float>, linal::Matrix< float> >
    NeuralNetwork::ComputeResultInputSensitivity( const linal::VectorConstInterface< float> &FEATURE) const
    {
      const size_t num_outs( m_Bias.LastElement().GetSize());
      const size_t num_ins( FEATURE.GetSize());
      linal::Vector< float> result( num_outs);
      linal::Vector< float> tmp_input_sensitivity( num_ins);
      float *itr_result( result.Begin());

      storage::List< storage::Vector< linal::Vector< float> > >::iterator
        itr_hidden( AcquireHiddenVectors()), itr_hidden_input( AcquireHiddenVectors()),
        itr_err( AcquireHiddenVectors());
      storage::Vector< linal::Vector< float> > &hidden( *itr_hidden), &hidden_input( *itr_hidden_input);
      storage::Vector< linal::Vector< float> > &errors( *itr_err);

      const size_t number_hidden_layers( hidden.GetSize());

      // input layer
      // perform hidden_input( 0) = m_Weight( 0) * FEATURE + m_Bias( 0); without creating new vectors
      hidden_input( 0) = m_Bias( 0);
      linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( 0), m_Weight( 0), FEATURE);

      // run hidden_input through the transfer functions to get the 1st hidden layer output
      m_TransferFunction->F( hidden( 0), hidden_input( 0));

      // remaining layers
      for( size_t k( 1); k < number_hidden_layers; ++k)
      {
        // perform hidden_input( k) = m_Weight( k) * hidden(k-1) + m_Bias( k); without creating new vectors
        hidden_input( k) = m_Bias( k);
        linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( k), m_Weight( k), hidden( k - 1));

        // run hidden_input through the transfer functions to get the k+1th hidden layer output
        m_TransferFunction->F( hidden( k), hidden_input( k));
      }

      // copy the results into a results matrix
      std::copy( hidden.LastElement().Begin(), hidden.LastElement().End(), itr_result);

      // rescale the result
      linal::MatrixReference< float> result_ref( size_t( 1), num_outs, result.Begin());
      m_RescaleOutput->DeScaleMatrix( result_ref);

      // now compute input sensitivity
      linal::Matrix< float> input_sensitivity( num_ins, num_outs);

      for( size_t output( 0); output < num_outs; ++output)
      {
        // set all values in the error vector to 0, except for the output currently under investigation
        linal::Vector< float> &last_err( errors.LastElement());

        last_err = float( 0.0);
        last_err( output) = 1.0;

        // the following code performs this operation
        // errors.LastElement() = ( RESULT - hidden.LastElement())
        //                          * m_TransferFunction->dF( hidden_input.LastElement(), hidden.LastElement())
        // without creating any (expensive) temporary vectors
        m_TransferFunction->MultiplyBydF
        (
          last_err,
          hidden_input.LastElement(),
          hidden.LastElement()
        );

        // all other layers
        for( size_t i( hidden.GetSize() - 1); i > 0; --i)
        {
          // errors( i - 1) = errors( i) * m_Weight( i)
          linal::GetDefaultOperations< float>().VectorEqualsVectorTimesMatrix( errors( i - 1), errors( i), m_Weight( i));
          m_TransferFunction->MultiplyBydF
          (
            errors( i - 1),
            hidden_input( i - 1),
            hidden( i - 1)
          );
        }

        linal::GetDefaultOperations< float>().VectorEqualsVectorTimesMatrix( tmp_input_sensitivity, errors( 0), m_Weight( 0));
        for( size_t input( 0); input < num_ins; ++input)
        {
          input_sensitivity( input, output) = tmp_input_sensitivity( input);
        }
      }

//    The following code can be used to verify that the above algorithm works
//
//      BCL_MessageStd( "input sensitivity computed via analytic method: " + util::Format()( input_sensitivity));
//
//      linal::Vector< float> feature_copy( FEATURE);
//      linal::Vector< float> original_result( hidden.LastElement());
//      linal::Matrix< float> numeric_input_sensitivity( num_ins, num_outs);
//
//      for( size_t feature( 0); feature < num_ins; ++feature)
//      {
//        feature_copy( feature) += 0.01;
//        // input layer
//
//        // perform hidden_input( 0) = m_Weight( 0) * FEATURE + m_Bias( 0); without creating new vectors
//        hidden_input( 0) = m_Bias( 0);
//        math::VectorPlusEqualsMatrixTimesVector( hidden_input( 0), m_Weight( 0), feature_copy);
//
//        // run hidden_input through the transfer functions to get the 1st hidden layer output
//        m_TransferFunction->F( hidden( 0), hidden_input( 0));
//
//        // remaining layers
//        for( size_t k( 1); k < number_hidden_layers; ++k)
//        {
//          // perform hidden_input( k) = m_Weight( k) * hidden(k-1) + m_Bias( k); without creating new vectors
//          hidden_input( k) = m_Bias( k);
//          math::VectorPlusEqualsMatrixTimesVector( hidden_input( k), m_Weight( k), hidden( k - 1));
//
//          // run hidden_input through the transfer functions to get the k+1th hidden layer output
//          m_TransferFunction->F( hidden( k), hidden_input( k));
//        }
//
//        feature_copy( feature) -= 0.01;
//        for( size_t out( 0); out < num_outs; ++out)
//        {
//          numeric_input_sensitivity( feature, out) = ( hidden.LastElement()( out) - original_result( out)) / 0.01;
//        }
//      }
//
//      BCL_MessageStd( "input sensitivity computed numerically: " + util::Format()( numeric_input_sensitivity));

      ReleaseHiddenVectors( itr_hidden);
      ReleaseHiddenVectors( itr_hidden_input);
      ReleaseHiddenVectors( itr_err);
      return storage::Pair< linal::Vector< float>, linal::Matrix< float> >( result, input_sensitivity);
    }

    //! @brief compute the result and deviation of that result using test-time dropout
    //! @param DATASET feature set of interest
    //! @param DROPOUT_RATIOS ratios of neurons to dropout from each layer (except output; excess layer ratios ignored)
    //! @param NREPEATS number of times to repeat dropout mask to get an estimate of the actual result
    storage::VectorND< 3, FeatureDataSet< float> > NeuralNetwork::TestWithDropout
    (
      const FeatureDataSetInterface< float> &FEATURE,
      const storage::Vector< double> &DROPOUT_RATIOS,
      const size_t &NREPEATS
    ) const
    {
      FeatureDataSet< float> rescaled( FEATURE);
      Rescale( rescaled);

      const size_t num_outs( m_Bias.LastElement().GetSize());
      const size_t num_ins( FEATURE.GetFeatureSize());
      const size_t num_points( FEATURE.GetNumberFeatures());
      BCL_Assert
      (
         num_ins == m_Weight( 0).GetNumberCols(),
        "Given " + util::Format()( num_ins) + " inputs; model expected " + util::Format()( m_Weight( 0).GetNumberCols())
      );
      linal::Matrix< float> result_matrix( num_points, num_outs);

      linal::Vector< size_t> n_to_drop( m_Architecture.GetSize() - size_t( 1), size_t( 0));
      storage::Vector< storage::Vector< size_t> > dropout_candidates( m_Architecture.GetSize() - size_t( 1));
      storage::Vector< storage::Vector< size_t> > chosen_to_drop( m_Architecture.GetSize() - size_t( 1));
      auto weightadj( m_Weight);
      for( size_t i( 0), sz( n_to_drop.GetSize()); i < sz; ++i)
      {
        n_to_drop( i) = i < DROPOUT_RATIOS.GetSize() ? size_t( DROPOUT_RATIOS( i) * m_Architecture( i)) : size_t( 0);
        dropout_candidates( i) = storage::CreateIndexVector( m_Architecture( i));
        chosen_to_drop( i).Resize( n_to_drop( i));
        weightadj( i) /= float( 1.0 - double( n_to_drop( i)) / double( weightadj( i).GetNumberCols()));
      }
      storage::List< storage::Vector< linal::Vector< float> > >::iterator
        itr_hidden( AcquireHiddenVectors()), itr_hidden_input( AcquireHiddenVectors()),
        itr_err( AcquireHiddenVectors());
      storage::Vector< linal::Vector< float> > &hidden( *itr_hidden), &hidden_input( *itr_hidden_input);

      const size_t number_hidden_layers( hidden.GetSize());

      linal::Vector< float> feature_with_dropout( num_ins, float( 0.0));
      math::RunningAverageSD< linal::Matrix< float> > ave_sd_result;
      for( size_t rep_nr( 0); rep_nr < NREPEATS; ++rep_nr)
      {
        util::GetLogger().LogStatus( "Sampling #" + util::Format()( rep_nr) + " / " + util::Format()( NREPEATS));
        for( size_t i( 0), sz( n_to_drop.GetSize()); i < sz; ++i)
        {
          storage::Vector< size_t> &chosen_array( chosen_to_drop( i));
          // shuffle the first n-indices in the m_NeuronIndices array
          storage::Vector< size_t> &indices_array( dropout_candidates( i));
          const size_t layer_size( indices_array.GetSize());
          for( size_t neuron( 0), d( n_to_drop( i)); neuron < d; ++neuron)
          {
            // TODO check whether this should be -1 or not
            std::swap( indices_array( neuron), indices_array( random::GetGlobalRandom().Random( neuron, layer_size - 1)));
            chosen_array( neuron) = indices_array( neuron);
          }
          chosen_array.Sort( std::less< size_t>());
        }
        for( size_t ftr_nr( 0); ftr_nr < num_points; ++ftr_nr)
        {
          // input layer
          // perform hidden_input( 0) = m_Weight( 0) * FEATURE + m_Bias( 0); without creating new vectors
          hidden_input( 0) = m_Bias( 0);

          if( n_to_drop( 0))
          {
            auto tvec( rescaled( ftr_nr));
            std::copy( tvec.Begin(), tvec.End(), feature_with_dropout.Begin());
            const storage::Vector< size_t> &drop_ids( chosen_to_drop( 0));
            for
            (
              storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
              itr_drop != itr_drop_end;
              ++itr_drop
            )
            {
              feature_with_dropout( *itr_drop) = 0.0;
            }
            //feature_with_dropout /= float( 1.0 - double( drop_ids.GetSize()) / double( num_ins));
            linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( 0), weightadj( 0), feature_with_dropout);
          }
          else
          {
            linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( 0), weightadj( 0), rescaled( ftr_nr));
          }

          // run hidden_input through the transfer functions to get the 1st hidden layer output
          m_TransferFunction->F( hidden( 0), hidden_input( 0));

          // remaining layers
          for( size_t k( 1); k < number_hidden_layers; ++k)
          {
            // perform hidden_input( k) = m_Weight( k) * hidden(k-1) + m_Bias( k); without creating new vectors
            hidden_input( k) = m_Bias( k);
            if( n_to_drop( k))
            {
              linal::Vector< float> &hidden_dropped( hidden( k - 1));
              const storage::Vector< size_t> &drop_ids( chosen_to_drop( k));
              for
              (
                storage::Vector< size_t>::const_iterator itr_dropped( drop_ids.Begin()), itr_dropped_end( drop_ids.End());
                itr_dropped != itr_dropped_end;
                ++itr_dropped
              )
              {
                hidden_dropped( *itr_dropped) = 0.0;
              }
              //hidden_dropped /= float( 1.0 - double( drop_ids.GetSize()) / double( dropout_candidates( k).GetSize()));
            }
            linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( k), weightadj( k), hidden( k - 1));

            // run hidden_input through the transfer functions to get the k+1th hidden layer output
            m_TransferFunction->F( hidden( k), hidden_input( k));
          }
          result_matrix.GetRow( ftr_nr).CopyValues( hidden.LastElement());
        }
        m_RescaleOutput->DeScaleMatrix( result_matrix);
        ave_sd_result += result_matrix;
      }
      FeatureDataSet< float> av( ave_sd_result.GetAverage());
      FeatureDataSet< float> sd( ave_sd_result.GetStandardDeviation());
      FeatureDataSet< float> res( PredictWithoutRescaling( rescaled));
      res.DeScale();
      storage::VectorND< 3, FeatureDataSet< float> > ave_sd_results( av, sd, res);
      return ave_sd_results;
    }

    //! @brief remove output layer and set last hidden layer to output layer
    void NeuralNetwork::RemoveOutputLayer()
    {
      m_Bias.RemoveElements( m_Bias.GetSize() - 1, size_t( 1));
      m_Weight.RemoveElements( m_Weight.GetSize() - 1, size_t( 1));
      m_Architecture.RemoveElements( m_Architecture.GetSize() - 1, size_t( 1));
      m_RescaleOutput = util::ShPtr< RescaleFeatureDataSet>
                        (
                          new RescaleFeatureDataSet
                          (
                            linal::Matrix< float>( size_t( 1), m_Weight.LastElement().GetNumberRows()),
                            math::Range< float>(),
                            RescaleFeatureDataSet::e_None
                          )
                        );
      m_Mutex.Lock();
      m_Available.Reset();
      m_Allocated.Reset();
      m_Mutex.Unlock();

      // determine architecture from weights and bias
      SetImplicitArchitecture();
    }

    //! @brief remove input layer and set last hidden layer to input layer
    void NeuralNetwork::RemoveInputLayer()
    {
      m_Bias.RemoveElements( 0, size_t( 1));
      m_Weight.RemoveElements( 0, size_t( 1));
      m_Architecture.RemoveElements( 0, size_t( 1));
//      m_RescaleOutput = util::ShPtr< RescaleFeatureDataSet>
//                        (
//                          new RescaleFeatureDataSet
//                          (
//                            linal::Matrix< float>( size_t( 1), m_Weight.LastElement().GetNumberRows()),
//                            math::Range<float>(),
//                            RescaleFeatureDataSet::e_None
//                          )
//                        );
      m_Mutex.Lock();
      m_Available.Reset();
      m_Allocated.Reset();
      m_Mutex.Unlock();

      // determine architecture from weights and bias
      SetImplicitArchitecture();
    }

    //! @brief append two networks if layers are compatible
    void NeuralNetwork::Append( const util::ShPtr< NeuralNetwork> &NETWORK)
    {
      // stop if the two networks do not fit together
      BCL_Assert
      (
        m_Architecture.LastElement() == NETWORK->GetArchitecture().FirstElement(),
        "Append: networks are not compatible! " + util::Format()( m_Architecture.LastElement())
        + " to " + util::Format()( NETWORK->GetArchitecture().FirstElement())
      );

      // append all necessary components
      m_Bias.Append( NETWORK->GetBias());
      m_Weight.Append( NETWORK->GetWeight());
      m_RescaleOutput = NETWORK->m_RescaleOutput;

      m_Mutex.Lock();
      m_Available.Reset();
      m_Allocated.Reset();
      m_Mutex.Unlock();

      SetImplicitArchitecture();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read NeuralNetwork from std::istream
    std::istream &NeuralNetwork::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_TransferFunction, ISTREAM);
      io::Serialize::Read( m_RescaleInput, ISTREAM);
      io::Serialize::Read( m_RescaleOutput, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Weight, ISTREAM);

      SetImplicitArchitecture();
      // end
      return ISTREAM;
    }

    //! write NeuralNetwork into std::ostream
    std::ostream &NeuralNetwork::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_TransferFunction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleInput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleOutput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Weight, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetwork::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "see http://en.wikipedia.org/wiki/Neural_network"
      );

      parameters.AddInitializer
      (
        "transfer function",
        "function that translates input from neurons in the prior layer into the output of each hidden layer",
        io::Serialization::GetAgent( &m_TransferFunction),
        "Sigmoid"
      );

      parameters.AddInitializer
      (
        "architecture",
        "# of neurons in each layer (input, hidden layer(s), output), e.g. (100, 4, 4, 1)",
        io::Serialization::GetAgentContainerWithCheck( &m_Architecture, io::Serialization::GetAgentWithMin( size_t( 1)))
      );

      // TODO: Add data members

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool NeuralNetwork::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // now set up the architecture
      SetArchitecture( m_Architecture);
      return true;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determines the architecture already present in the weight matrix/vectors of the neural network
    void NeuralNetwork::SetImplicitArchitecture()
    {
      // architecture
      m_Architecture.Reset();

      if( m_Weight.GetSize())
      {
        m_Architecture.AllocateMemory( GetNumberLayers() + 1);

        // number of input neurons
        m_Architecture.PushBack( GetNumberInputs());

        // iterate over all layers to get their size
        storage::Vector< linal::Vector< float> >::const_iterator itr_bias( m_Bias.Begin());
        for
        (
          storage::Vector< linal::Matrix< float> >::const_iterator
            layer_itr( m_Weight.Begin()), layer_itr_end( m_Weight.End());
          layer_itr != layer_itr_end;
          ++layer_itr, ++itr_bias
        )
        {
          m_Architecture.PushBack( layer_itr->GetNumberRows());
          BCL_Assert( itr_bias->GetSize() == layer_itr->GetNumberRows(), "Non-matching bias to weight dimensions!");
        }
      }
    }

    //! @brief acquire a hidden input or hidden vector set (same size as m_Bias)
    //! @return iterator to the hidden input or hidden vector set
    storage::List< storage::Vector< linal::Vector< float> > >::iterator NeuralNetwork::AcquireHiddenVectors() const
    {
      m_Mutex.Lock();
      // all hidden vectors have been allocated
      if( m_Available.IsEmpty())
      {
        m_Available.PushBack( m_Bias);
      }
      // splice the last hidden vector set from the available list onto the allocated list
      storage::List< storage::Vector< linal::Vector< float> > >::iterator avail_itr( m_Available.Begin());
      m_Allocated.InternalData().splice( m_Allocated.Begin(), m_Available.InternalData(), avail_itr);

      // save the iterator to the now-allocated array set
      avail_itr = m_Allocated.Begin();
      m_Mutex.Unlock();
      return avail_itr;
    }

    //! @brief release a given hidden vector set
    //! @param ITR iterator to the hidden input or hidden vector set
    void NeuralNetwork::ReleaseHiddenVectors
    (
      const storage::List< storage::Vector< linal::Vector< float> > >::iterator &ITR
    ) const
    {
      m_Mutex.Lock();
      // splice the iterator from the allocated back onto available
      m_Available.InternalData().splice( m_Available.Begin(), m_Allocated.InternalData(), ITR);
      m_Mutex.Unlock();
    }

    //! @brief function called by each thread to accomplish threaded forward propagation
    //! @param FEATURES features to train on
    //! @param STORAGE storage for the result
    void NeuralNetwork::PredictWithThread
    (
      const FeatureDataSetInterface< float> &FEATURES,
      linal::MatrixInterface< float> &STORAGE
    ) const
    {
      storage::List< storage::Vector< linal::Vector< float> > >::iterator
        itr_hidden( AcquireHiddenVectors()), itr_hidden_input( AcquireHiddenVectors());
      storage::Vector< linal::Vector< float> > &hidden( *itr_hidden), &hidden_input( *itr_hidden_input);

      const size_t num_points( FEATURES.GetNumberFeatures());
      const size_t num_outs( m_Bias.LastElement().GetSize());
      const size_t number_hidden_layers( hidden.GetSize());

      float *itr_result( STORAGE.Begin());
      const float *itr_hidden_result( hidden.LastElement().Begin());
      const float *itr_hidden_result_end( hidden.LastElement().End());

      for( size_t feature_number( 0); feature_number < num_points; ++feature_number, itr_result += num_outs)
      {
        // input layer
        // perform hidden_input( 0) = m_Weight( 0) * FEATURE + m_Bias( 0); without creating new vectors
        hidden_input( 0) = m_Bias( 0);

        linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( 0), m_Weight( 0), FEATURES( feature_number));

        // run hidden_input through the transfer functions to get the 1st hidden layer output
        m_TransferFunction->F( hidden( 0), hidden_input( 0));

        // remaining layers
        for( size_t k( 1); k < number_hidden_layers; ++k)
        {
          // perform hidden_input( k) = m_Weight( k) * hidden(k-1) + m_Bias( k); without creating new vectors
          hidden_input( k) = m_Bias( k);
          linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( k), m_Weight( k), hidden( k - 1));

          // run hidden_input through the transfer functions to get the k+1th hidden layer output
          m_TransferFunction->F( hidden( k), hidden_input( k));
        }

        // copy the results into a results matrix
        std::copy( itr_hidden_result, itr_hidden_result_end, itr_result);
      }

      ReleaseHiddenVectors( itr_hidden);
      ReleaseHiddenVectors( itr_hidden_input);
    }

  } // namespace model
} // namespace bcl
