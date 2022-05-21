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
#include "model/bcl_model_objective_function_categorical_max.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionCategoricalMax::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionCategoricalMax()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new Evaluator
    ObjectiveFunctionCategoricalMax *ObjectiveFunctionCategoricalMax::Clone() const
    {
      return new ObjectiveFunctionCategoricalMax( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ObjectiveFunctionCategoricalMax::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionCategoricalMax::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      size_t results_size( DATA.GetFeatureSize());
      size_t dataset_size( DATA.GetNumberFeatures());
      linal::MatrixConstReference< float> results_ref( DATA.GetMatrix());
      linal::Vector< size_t> classification_store( results_size, size_t( 0));
      size_t number_classifications( 0);

      util::SiPtr< const RescaleFeatureDataSet> rescale( DATA.GetScaling());
      {
        size_t classification_id( 0);
        const float *itr_result( results_ref[ 0]);
        if( rescale.IsDefined())
        {
          for( size_t result_id( 0); result_id < results_size; ++result_id, ++itr_result)
          {
            if( rescale->DescaleValue( result_id, *itr_result) > 0.5)
            {
              classification_store( classification_id) = result_id;
              ++classification_id;
            }
          }
        }
        else
        {
          for( size_t result_id( 0); result_id < results_size; ++result_id, ++itr_result)
          {
            if( *itr_result > 0.5)
            {
              classification_store( classification_id) = result_id;
              ++classification_id;
            }
          }
        }
        m_ResultsBestIndex = linal::Matrix< size_t>( dataset_size, classification_id);
        classification_store = classification_store.CreateSubVector( classification_id);
        number_classifications = classification_id;
        m_ResultsBestIndex.ReplaceRow( 0, classification_store);
        m_ClassBoundaries = classification_store;
      }
      for( size_t feature( 0); feature < dataset_size; ++feature)
      {
        size_t classification_id( 0);
        const float *itr_result( results_ref[ feature]);
        for( size_t result_id( 0); result_id < results_size; ++result_id, ++itr_result)
        {
          if
          (
            ( rescale.IsDefined() && rescale->DescaleValue( result_id, *itr_result) > 0.5)
            || ( !rescale.IsDefined() && *itr_result > 0.5)
          )
          {
            BCL_Assert
            (
              classification_id < number_classifications,
              "MultiClassification requires a constant # of classifications per feature row, 1st row had "
              + util::Format()( number_classifications) + " but row " + util::Format()( feature) + " had "
              + util::Format()( classification_id + 1) + " or more "
            );
            classification_store( classification_id) = result_id;
            m_ClassBoundaries( classification_id)
              = std::min( m_ClassBoundaries( classification_id), classification_store( classification_id));
            ++classification_id;
          }
        }
        BCL_Assert
        (
          classification_id == number_classifications,
          "MultiClassification requires a constant # of classifications per feature row, 1st row had "
          + util::Format()( number_classifications) + " but row " + util::Format()( feature) + " had "
          + util::Format()( classification_id)
        );
        m_ResultsBestIndex.ReplaceRow( feature, classification_store);
      }
      {
        storage::Vector< size_t> class_bounds( m_ClassBoundaries.Begin() + 1, m_ClassBoundaries.End());
        class_bounds.PushBack( results_size);
        m_ClassBoundaries = linal::Vector< size_t>( class_bounds.Begin(), class_bounds.End());
      }
      m_PredictedClasses = m_ResultsBestIndex;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief classify. Obtain a matrix of '0'/'1' values (as char, for memory) as to whether the objective function is
    //!        satisfied with each result.  Denote results that the objective function ignores by \0
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with '0'/'1'/\0 values indicating Bad/Good/Not considered prediction
    linal::Matrix< char> ObjectiveFunctionCategoricalMax::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // number of data points in dataset
      const size_t data_set_size( EXPERIMENTAL.GetNumberFeatures());
      const size_t result_size( EXPERIMENTAL.GetFeatureSize());

      BCL_Assert
      (
        EXPERIMENTAL.GetNumberFeatures() == PREDICTED.GetNumberFeatures(),
        "Experimental and predicted values do not have the same number of elements!"
      );
      BCL_Assert
      (
        result_size == PREDICTED.GetFeatureSize(),
        "Experimental and predicted values have different result sizes"
      );
      BCL_Assert
      (
        m_ClassBoundaries.GetSize() && result_size == m_ClassBoundaries.Last(),
        "SetData was not called (based on class boundaries: " + util::Format()( m_ClassBoundaries)
        + "; max should be: " + util::Format()( result_size)
      );

      linal::Matrix< char> prediction_classification( data_set_size, result_size, '\0');

      if( &EXPERIMENTAL != &PREDICTED)
      {
        BCL_Assert
        (
          data_set_size == m_ResultsBestIndex.GetNumberRows(),
          "Set data was not called (using different experimental and predicted"
        );
        // iterate over all experimental and predicted values and check for agreement on which index is largest
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          // get the vector containing the original maxes
          linal::VectorConstReference< size_t> results_best_index( m_ResultsBestIndex.GetRow( counter));

          size_t index_start( 0);
          const size_t *itr_best( results_best_index.Begin()), *itr_best_end( results_best_index.End());
          const size_t *itr_category( m_ClassBoundaries.Begin());
          const float *prediction( PREDICTED[ counter]);
          char *pred_classification_row( prediction_classification[ counter]);

          // track overall position on this row
          size_t col( 0);

          for( ; itr_best != itr_best_end; ++itr_best, ++itr_category)
          {
            // get the value from what should be the max column
            const size_t desired_max( *itr_best);

            // get the value out of the prediction
            const float desired_prediction_value( prediction[ desired_max]);

            // initially, set the classification of the max to P (TP); it will be set to n (FN) if it is found not
            // to be the max
            pred_classification_row[ desired_max] = 'P';

            // get the actual maximum
            for
            (
              const float *itr_prediction( prediction + index_start), *itr_prediction_end( prediction + *itr_category);
              itr_prediction != itr_prediction_end;
              ++itr_prediction, ++col
            )
            {
              if( *itr_prediction < desired_prediction_value)
              {
                pred_classification_row[ col] = 'N';
              }
              else if( *itr_prediction > desired_prediction_value || col != desired_max)
              {
                pred_classification_row[ col] = 'p';
                pred_classification_row[ desired_max] = 'n';
              }
            }

            index_start = *itr_category;
          }
        }
      }
      else
      {
        // just go by the real maximum
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          // get the vector containing the original maxes
          linal::VectorConstReference< size_t> results_best_index( m_ResultsBestIndex.GetRow( counter));

          size_t index_start( 0);
          const size_t *itr_category( m_ClassBoundaries.Begin()), *itr_category_end( m_ClassBoundaries.End());
          const float *prediction( PREDICTED[ counter]);
          char *pred_classification_row( prediction_classification[ counter]);

          // track overall position on this row
          size_t col( 0);

          for( ; itr_category != itr_category_end; ++itr_category)
          {
            // get the value from what should be the max column
            const size_t desired_max
            (
              std::max_element( prediction + index_start, prediction + *itr_category) - prediction
            );

            // get the value out of the prediction
            const float desired_prediction_value( prediction[ desired_max]);

            // initially, set the classification of the max to P (TP); it will be set to n (FN) if it is found not
            // to be the max
            pred_classification_row[ desired_max] = 'P';

            // get the actual maximum
            for
            (
              const float *itr_prediction( prediction + index_start), *itr_prediction_end( prediction + *itr_category);
              itr_prediction != itr_prediction_end;
              ++itr_prediction, ++col
            )
            {
              if( *itr_prediction < desired_prediction_value)
              {
                pred_classification_row[ col] = 'N';
              }
              else if( *itr_prediction > desired_prediction_value || col != desired_max)
              {
                pred_classification_row[ col] = 'p';
                pred_classification_row[ desired_max] = 'n';
              }
            }

            index_start = *itr_category;
          }
        }
      }
      return prediction_classification;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    float ObjectiveFunctionCategoricalMax::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // number of data points in dataset
      const size_t data_set_size( EXPERIMENTAL.GetNumberFeatures());

      BCL_Assert
      (
        EXPERIMENTAL.GetNumberFeatures() == PREDICTED.GetNumberFeatures(),
        "Experimental and predicted values do not have the same number of elements!"
      );
      BCL_Assert
      (
        EXPERIMENTAL.GetFeatureSize() == PREDICTED.GetFeatureSize(),
        "Experimental and predicted values have different result sizes"
      );
      BCL_Assert
      (
        data_set_size == m_ResultsBestIndex.GetNumberRows(),
        "Set data was not called"
      );

      // count the number of predictions that were on the same side of the cutoff as the experimental predictions
      linal::Vector< size_t> accurate_prediction_count( m_ClassBoundaries.GetSize(), size_t( 0));
      storage::Vector< linal::Matrix< size_t> > contingency_matrices( m_ClassBoundaries.GetSize());
      size_t last_class_start( 0), class_set_number( 0);
      for
      (
        auto itr_pre( m_ClassBoundaries.Begin()), itr_pre_end( m_ClassBoundaries.End());
        itr_pre != itr_pre_end;
        ++itr_pre, ++class_set_number
      )
      {
        const size_t class_size( *itr_pre - last_class_start);
        contingency_matrices( class_set_number) = linal::Matrix< size_t>( class_size, class_size, size_t( 0));
        last_class_start = *itr_pre;
      }

      // iterate over all experimental and predicted values and check for agreement on which index is largest
      for( size_t counter( 0); counter < data_set_size; ++counter)
      {
        // get the vector containing the original maxes
        linal::VectorConstReference< size_t> results_best_index( m_ResultsBestIndex.GetRow( counter));
        linal::VectorReference< size_t> predicted_best_index( m_PredictedClasses.GetRow( counter));

        size_t index_start( 0);
        const size_t *itr_best( results_best_index.Begin()), *itr_best_end( results_best_index.End());
        size_t *itr_predicted( predicted_best_index.Begin());
        const size_t *itr_category( m_ClassBoundaries.Begin());
        storage::Vector< linal::Matrix< size_t> >::iterator itr_contingency_matrices( contingency_matrices.Begin());
        size_t *itr_accurate( accurate_prediction_count.Begin());
        const float *prediction( PREDICTED[ counter]);

        for( ; itr_best != itr_best_end; ++itr_best, ++itr_category, ++itr_accurate, ++itr_contingency_matrices, ++itr_predicted)
        {
          // get the value from what should be the max column
          const size_t desired_max( *itr_best);

          // get the actual maximum
          const size_t actual_max
          (
            std::max_element( prediction + index_start, prediction + *itr_category) - prediction
          );
          *itr_predicted = actual_max;

          // compare the predicted max with the desired max for this class to determine if it was accurate
          if( actual_max == desired_max)
          {
            ++*itr_accurate;
          }

          // update the contingency matrix
          ++( *itr_contingency_matrices)( desired_max - index_start, actual_max - index_start);

          // update start index for the next category
          index_start = *itr_category;
        }
      }

      if( util::GetMessenger().GetCurrentMessageLevel() >= util::Message::e_Verbose)
      {
        std::ostringstream oss;
        oss << "Accuracies By Class:\t";
        for
        (
          size_t class_boundary( 0), n_class_boundaries( m_ClassBoundaries.GetSize());
          class_boundary < n_class_boundaries;
          ++class_boundary
        )
        {
          oss << class_boundary << ":\t"
              << float( accurate_prediction_count( class_boundary)) / float( data_set_size) << '\t';
        }
        util::GetLogger() << oss.str();
        oss.str( std::string());
        oss << "Contingency matrices by Class:\n";
        for
        (
          size_t class_boundary( 0), n_class_boundaries( m_ClassBoundaries.GetSize());
          class_boundary < n_class_boundaries;
          ++class_boundary
        )
        {
          oss << "\nClass " << class_boundary << " Counts\n";

          // get the current matrix
          const linal::Matrix< size_t> &contingency_matrix( contingency_matrices( class_boundary));
          const size_t subclasses( contingency_matrix.GetNumberRows());

          // write the header
          oss << "vPredicted >Experimental\t";
          for( size_t sub_class_a( 0); sub_class_a < subclasses; ++sub_class_a)
          {
            oss << sub_class_a << '\t';
          }
          oss << '\n';
          for( size_t sub_class_a( 0); sub_class_a < subclasses; ++sub_class_a)
          {
            oss << sub_class_a << '\t';
            for( size_t sub_class_b( 0); sub_class_b < subclasses; ++sub_class_b)
            {
              oss << contingency_matrix( sub_class_a, sub_class_b) << '\t';
            }
            oss << '\n';
          }
        }
        util::GetLogger() << oss.str();
      }

      // return accuracy
      return double( accurate_prediction_count.Sum()) / double( data_set_size * m_ClassBoundaries.GetSize());
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionCategoricalMax::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the fraction of predictions that were correct for results which are series of one or more blocks of "
        "classifications, with the threshold for a true classification set at 0.5.  For example, "
        "(0 0.1 0.4 1) represents an output with a single classification (as only one value is above 0.5, while "
        "(1 0 0 0 0 1 0 0) represents an output with multiple classifications.  The number of classification categories "
        "must remain constant in a given dataset; e.g. (0 1 0 1) cannot belong to the same dataset as (0 0 0 1) for this "
        "objective function.  Detection of classification boundaries is performed automatically, and each classification "
        "block is given the same weight, so for X classification blocks, the output will range from 0 to X.  Message level "
        "can be set to verbose to enable output of contingency matrices for all classes."
      );
      parameters.AddDataMember
      (
        "borders",
        io::Serialization::GetAgent( &m_ClassBoundaries)
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
