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
#include "model/bcl_model_objective_function_segment_overlap.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_data_set_select_columns.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionSegmentOverlap::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionSegmentOverlap()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionSegmentOverlap::ObjectiveFunctionSegmentOverlap() :
      ObjectiveFunctionCategoricalMax(),
      m_PerIdOutput( false),
      m_OutputSubclassOverlaps( false),
      m_WeightPerElement( false),
      m_SequenceIdLabel( "ProteinId"),
      m_ElementIdLabel( "AASeqID"),
      m_MissingIds(),
      m_BoundaryIds(),
      m_ClassCounts(),
      m_SequenceNames(),
      m_SequenceSizes(),
      m_ActualSegments()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Evaluator
    ObjectiveFunctionSegmentOverlap *ObjectiveFunctionSegmentOverlap::Clone() const
    {
      return new ObjectiveFunctionSegmentOverlap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ObjectiveFunctionSegmentOverlap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionSegmentOverlap::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      // reset all data members
      m_MissingIds.Reset();
      m_BoundaryIds.Reset();
      m_SequenceNames.Reset();
      m_SequenceSizes.Reset();
      m_ActualSegments.Reset();

      // set data on the base class
      ObjectiveFunctionCategoricalMax::SetData( DATA, IDS);

      // get the feature labels so that sequences can be identified
      const FeatureLabelSet model_feature_labels( *IDS.GetFeatureLabelSet());
      DataSetSelectColumns sequence_id_selector
      (
        model_feature_labels.GetSize(),
        IDS.GetFeatureLabelSet()->GetPropertyIndices( m_SequenceIdLabel)
      );
      DataSetSelectColumns element_id_selector
      (
        model_feature_labels.GetSize(),
        IDS.GetFeatureLabelSet()->GetPropertyIndices( m_ElementIdLabel)
      );
      BCL_Assert( sequence_id_selector.GetOutputFeatureSize(), "Could not find " + m_SequenceIdLabel.ToString() + " in IDs!");
      BCL_Assert( element_id_selector.GetOutputFeatureSize(), "Could not find " + m_ElementIdLabel.ToString() + " in IDs!");
      const size_t dataset_size( DATA.GetNumberFeatures());
      BCL_Assert( dataset_size == IDS.GetNumberFeatures(), "IDs and Data have different size!");

      // add boundary for the first sequence, unless the dataset is empty
      if( dataset_size)
      {
        m_BoundaryIds.PushBack( 0);
      }
      {
        // locate all the sequence boundaries and ids
        std::string last_sequence_id( sequence_id_selector.GetOutputFeatureSize(), ' ');
        std::string this_sequence_id( last_sequence_id);
        for( size_t i( 0); i < dataset_size; ++i)
        {
          sequence_id_selector( IDS( i), &this_sequence_id[ 0]);
          if( this_sequence_id != last_sequence_id)
          {
            if( i)
            {
              m_SequenceSizes.PushBack( i - m_BoundaryIds.LastElement());
              m_BoundaryIds.PushBack( i);
            }
            last_sequence_id = this_sequence_id;
            m_SequenceNames.PushBack( this_sequence_id);
          }
        }
      }

      // add the final sequence's size and boundary id information
      if( dataset_size)
      {
        m_SequenceSizes.PushBack( dataset_size - m_BoundaryIds.LastElement());
        m_BoundaryIds.PushBack( dataset_size);
      }

      // get the number of sequences
      const size_t number_sequences( m_SequenceSizes.GetSize());
      const size_t number_class_groups( this->GetClassBoundaries().GetSize());
      const size_t number_results( DATA.GetFeatureSize());

      // create the other data members

      // missing ids for each sequence
      m_MissingIds.Reset();
      m_MissingIds.Resize( number_sequences);

      // class counts
      m_ClassCounts = linal::Matrix< size_t>( number_sequences, number_results);

      // base-classes matrix containing the class for each row
      const linal::Matrix< size_t> &actual_classes( this->GetActualClasses());

      // segments
      m_ActualSegments.Reset();
      m_ActualSegments.Resize( number_sequences);
      m_ActualSegments.SetAllElements
      (
        storage::Vector< storage::Vector< storage::Triplet< size_t, size_t, size_t> > >( number_class_groups)
      );

      for( auto itr_perm( m_AllowedPermutations.Begin()), itr_perm_end( m_AllowedPermutations.End()); itr_perm != itr_perm_end; ++itr_perm)
      {
        BCL_Assert( itr_perm->GetSize() == number_results, "All permutations should have the same size as the result output!");
      }
      m_AllowedPermutations.InsertElements( 0, storage::CreateIndexVector( number_results));

      size_t sequence_number( 0);
      // members for missing ids in the sequence
      std::string this_element_id_str( element_id_selector.GetOutputFeatureSize(), ' ');
      for
      (
        storage::Vector< size_t>::const_iterator
          itr_bounds( m_BoundaryIds.Begin()),
          itr_bounds_next( m_BoundaryIds.Begin() + 1),
          itr_bounds_end( m_BoundaryIds.End());
        itr_bounds_next != itr_bounds_end;
        ++itr_bounds, ++itr_bounds_next, ++sequence_number
      )
      {
        const size_t protein_start_row( *itr_bounds);
        const size_t protein_end_row( *itr_bounds_next);
        size_t last_element_id( util::GetUndefinedSize_t());

        // get a reference to the relevant row in the class-bounds matrix
        linal::VectorReference< size_t> class_counts( m_ClassCounts.GetRow( sequence_number));

        // get a reference to the vector where missing ids will be placed
        storage::Vector< size_t> &missing_id_vec( m_MissingIds( sequence_number));

        for( size_t row( protein_start_row); row < protein_end_row; ++row)
        {
          // accumulate class counts
          for( size_t subclass( 0); subclass < number_class_groups; ++subclass)
          {
            ++class_counts( actual_classes( row, subclass));
          }

          // get this row's ids
          element_id_selector( IDS( row), &this_element_id_str[ 0]);

          // get the id as an integer
          const size_t this_id( util::ConvertStringToNumericalValue< size_t>( this_element_id_str));

          // compare the id to the previous id
          if( row == protein_start_row)
          {
            // first row, just accept this id as the previous id
            last_element_id = this_id;
          }
          else
          {
            // after the first row
            if( this_id != last_element_id + 1 && actual_classes.GetRow( row) != actual_classes.GetRow( row - 1))
            {
              missing_id_vec.PushBack( row);
            }
            last_element_id = this_id;
          }
        }
        missing_id_vec.PushBack( protein_end_row);

        // get the segments for each subclass of this one
        for( size_t subclass( 0); subclass < number_class_groups; ++subclass)
        {
          m_ActualSegments( sequence_number)( subclass) = GetSegments( subclass, sequence_number, actual_classes);
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    float ObjectiveFunctionSegmentOverlap::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // compute QX
      const float basic_state_accuracy( ObjectiveFunctionCategoricalMax::operator()( EXPERIMENTAL, PREDICTED));

      BCL_Assert
      (
        EXPERIMENTAL.GetNumberFeatures() == size_t( 0) || !m_SequenceSizes.IsEmpty(),
        "SetData was not called before operator()!"
      );

      // compute average SOVs, one for each subclass
      const size_t number_subclasses( EXPERIMENTAL.GetFeatureSize());
      storage::Vector< math::RunningAverage< float> > segment_overlap_aves( number_subclasses);

      // compute average SOVs, one for each class
      const size_t number_classifications( this->GetClassBoundaries().GetSize());
      storage::Vector< math::RunningAverage< float> > segment_overlap_overall_aves( number_classifications);
      math::RunningAverage< float> combined_ave;

      // create a vector with the class boundaries and then the number of results
      storage::Vector< size_t> all_class_boundaries( size_t( 1), size_t( 0));
      all_class_boundaries.Append
      (
        storage::Vector< size_t>( this->GetClassBoundaries().Begin(), this->GetClassBoundaries().End())
      );
      const bool consider_perms( m_PerIdOutput && m_AllowedPermutations.GetSize() > size_t( 1));

      std::ostringstream output;
      output << "\nQ-overall: " << basic_state_accuracy << '\n';
      if( m_PerIdOutput)
      {
        output << "Sequence\t";
      }
      for( size_t classification_number( 0); classification_number < number_classifications; ++classification_number)
      {
        output << "SOV Overall Class " << classification_number << '\t';
        if( m_OutputSubclassOverlaps)
        {
          const size_t subclass_start( all_class_boundaries( classification_number));
          const size_t subclass_end( all_class_boundaries( classification_number + 1));

          for( size_t subclass_id( subclass_start); subclass_id < subclass_end; ++subclass_id)
          {
            output << "SOV Subclass " << subclass_id << '\t';
          }
        }
        if( consider_perms)
        {
          output << "Permutation #";
        }
      }
      output << '\n';
      // compute SOV for each sequence and subclass type, if desired
      for
      (
        size_t sequence_number( 0), number_sequences( m_SequenceSizes.GetSize());
        sequence_number < number_sequences;
        ++sequence_number
      )
      {
        if( m_PerIdOutput)
        {
          output << m_SequenceNames( sequence_number) << '\t';
        }
        const float seq_size( m_SequenceSizes( sequence_number));

        // for each classification
        for( size_t classification_number( 0); classification_number < number_classifications; ++classification_number)
        {
          // compute overall segment overlaps
          const storage::Pair< float, size_t> sov_this_sequence_classification
          (
            GetSegmentOverlap( classification_number, sequence_number)
          );

          if( !m_WeightPerElement)
          {
            segment_overlap_overall_aves( classification_number) += sov_this_sequence_classification.First();
            combined_ave += sov_this_sequence_classification.First();
          }
          else
          {
            segment_overlap_overall_aves( classification_number).AddWeightedObservation
            (
              sov_this_sequence_classification.First(),
              double( seq_size)
            );
            combined_ave.AddWeightedObservation( sov_this_sequence_classification.First(), double( seq_size));
          }
          if( m_PerIdOutput)
          {
            output << sov_this_sequence_classification.First() << '\t';
          }

          // compute segment overlaps for every subclass of this sequence, if desired
          if( m_OutputSubclassOverlaps)
          {
            const size_t subclass_start( all_class_boundaries( classification_number));
            const size_t subclass_end( all_class_boundaries( classification_number + 1));
            for( size_t subclass_id( subclass_start); subclass_id < subclass_end; ++subclass_id)
            {
              const float sov_this_sequence_subclassification
              (
                GetSegmentOverlap( classification_number, sequence_number, subclass_id, sov_this_sequence_classification.Second()).First()
              );

              if( !m_WeightPerElement)
              {
                segment_overlap_aves( subclass_id) += sov_this_sequence_subclassification;
              }
              else
              {
                segment_overlap_aves( subclass_id).AddWeightedObservation
                (
                  sov_this_sequence_subclassification,
                  seq_size
                );
              }

              if( m_PerIdOutput)
              {
                output << sov_this_sequence_subclassification << '\t';
              }
            }
          }
          if( consider_perms)
          {
            output << sov_this_sequence_classification.Second() << '\t';
          }
        }
        if( m_PerIdOutput)
        {
          output << '\n';
        }
      }
      if( m_PerIdOutput)
      {
        output << "Overall\t";
      }
      for( size_t classification_number( 0); classification_number < number_classifications; ++classification_number)
      {
        output << segment_overlap_overall_aves( classification_number).GetAverage() << '\t';
        if( m_OutputSubclassOverlaps)
        {
          const size_t subclass_start( all_class_boundaries( classification_number));
          const size_t subclass_end( all_class_boundaries( classification_number + 1));

          for( size_t subclass_id( subclass_start); subclass_id < subclass_end; ++subclass_id)
          {
            output << segment_overlap_aves( subclass_id).GetAverage() << '\t';
          }
        }
      }
      BCL_MessageStd( output.str());
      // return accuracy
      return combined_ave.GetAverage();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionSegmentOverlap::GetSerializer() const
    {
      io::Serializer parameters( ObjectiveFunctionCategoricalMax::GetSerializer());
      parameters.SetClassDescription
      (
        "Segment overlap is a common method for benchmarking secondary structure prediction methods. "
        "The version implemented here is the updated definition from Zemla et al. - PROTEINS: Structure, Function, and "
        "Genetics, 34, 1999, pp. 220-223"
      );
      parameters.AddInitializer
      (
        "output sequence info",
        "whether to output information for every sequence member",
        io::Serialization::GetAgent( &m_PerIdOutput),
        "False"
      );
      parameters.AddInitializer
      (
        "output subclass overlaps",
        "Whether to output SOV for every subclass; for secondary structure predictions, these are usually termed "
        "SOV-Helix, SOV-Strand, SOV-Coil",
        io::Serialization::GetAgent( &m_OutputSubclassOverlaps),
        "False"
      );
      parameters.AddInitializer
      (
        "element weight",
        "Whether to weight each sequence's contribution to the SOV by its number of elements (e.g. residues)",
        io::Serialization::GetAgent( &m_WeightPerElement),
        "False"
      );
      parameters.AddInitializer
      (
        "sequence id",
        "label from the ids that indicates the sequence",
        io::Serialization::GetAgent( &m_SequenceIdLabel),
        "ProteinId"
      );
      parameters.AddInitializer
      (
        "element id",
        "label from the ids that indicates the element's position within the sequence",
        io::Serialization::GetAgent( &m_ElementIdLabel),
        "AASeqID"
      );
      parameters.AddOptionalInitializer
      (
        "equivalences",
        "If, for a given sequence, a permutation of the labels is equally valid, it can be indicated here. For example, "
        "for membrane proteins if we have the labels Inside/Outside/Membrane, for most purposes it would be fine to flip "
        "the labels for inside and outside, so long as it is done consistently across the protein. Which sequences used "
        "a different permutation will be indicated in the output. For multiple class groups, the indexing will not start over in "
        "each, so a permutation label like (0,1,2,3,4,5) would be the identity permutation regardless of how many internal "
        "class groups there happen to be",
        io::Serialization::GetAgent( &m_AllowedPermutations)
      );
      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool ObjectiveFunctionSegmentOverlap::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      m_ElementIdLabel.SetName( "", true);
      m_SequenceIdLabel.SetName( "", true);
      return true;
    }

    //! @brief get the segments for a given subclass and sequence
    //! @param SUBCLASS subclass index; if just computing secondary structure using on method, this is always 0; but
    //!        when using multiple methods, it may range from 0-number classifications
    //! @param SEQUENCE sequence index
    //! @param CHOSEN_SUBCLASS_MATRIX matrix of chosen sub-classes
    //! @return vector of segments for the given sequence for the given subclass
    ObjectiveFunctionSegmentOverlap::SegmentContainer
      ObjectiveFunctionSegmentOverlap::GetSegments
      (
        const size_t &SUBCLASS,
        const size_t &SEQUENCE,
        const linal::MatrixConstInterface< size_t> &CHOSEN_CLASSES,
        const size_t &PERMUTATION
      ) const
    {
      SegmentContainer segments;
      // get the bounds for the sequence
      const size_t start_row( m_BoundaryIds( SEQUENCE)), end_row( m_BoundaryIds( SEQUENCE + 1));

      // detect empty sequences
      if( start_row == end_row)
      {
        return segments;
      }

      // get the chosen class at the first position
      size_t prev_class( CHOSEN_CLASSES( start_row, SUBCLASS));

      // insert it into the segments
      segments.PushBack( storage::Triplet< size_t, size_t, size_t>( prev_class, start_row, start_row));

      // get an iterator to the missing ids for this sequence
      storage::Vector< size_t>::const_iterator itr_missing( m_MissingIds( SEQUENCE).Begin());

      // get the permutation
      const storage::Vector< size_t> &perm( m_AllowedPermutations( PERMUTATION));

      for( size_t pos( start_row + 1); pos < end_row; ++pos)
      {
        // get the chosen class at this position
        const size_t pos_class( perm( CHOSEN_CLASSES( pos, SUBCLASS)));

        // test whether there was a gap between the previous row and this one
        if( *itr_missing == pos)
        {
          // end the previous segment, start a new one
          segments.PushBack( storage::Triplet< size_t, size_t, size_t>( pos_class, pos, pos));
          ++itr_missing;
          prev_class = pos_class;
        }
        else if( pos_class != prev_class)
        {
          // start a new segment
          segments.PushBack( storage::Triplet< size_t, size_t, size_t>( pos_class, pos, pos));
          prev_class = pos_class;
        }
        else
        {
          // extension to the previous segment
          ++segments.LastElement().Third();
        }
      }
      return segments;
    }

    //! @brief compute segment overlap from predicted to experimental
    //! @param SUBCLASS subclass index; if just computing secondary structure using on method, this is always 0; but
    //!        when using multiple methods, it may range from 0-number classifications
    //! @param SEQUENCE_ID sequence index
    //! @param RESTRICT_CLASS optional value; if set, only consider sub
    //! @param RESTRICT_PERMUTATION optional value; if set, only consider the indicated permutation
    //! @return the segment overlap and the permutation index used
    storage::Pair< float, size_t> ObjectiveFunctionSegmentOverlap::GetSegmentOverlap
    (
      const size_t &SUBCLASS,
      const size_t &SEQUENCE,
      const size_t &RESTRICT_CLASS,
      const size_t &RESTRICT_PERMUTATION
    ) const
    {
      // handle differential permutations
      if( RESTRICT_PERMUTATION == util::GetUndefined< size_t>())
      {
        storage::Pair< float, size_t> res_perm
        (
          GetSegmentOverlap( SUBCLASS, SEQUENCE, util::GetUndefined< size_t>(), size_t( 0))
        );
        for( size_t perm_id( 1), n_perms( m_AllowedPermutations.GetSize()); perm_id < n_perms; ++perm_id)
        {
          storage::Pair< float, size_t> res_perm_tmp
          (
            GetSegmentOverlap( SUBCLASS, SEQUENCE, util::GetUndefined< size_t>(), perm_id)
          );
          if( res_perm_tmp.First() > res_perm.First())
          {
            res_perm = res_perm_tmp;
          }
        }
        if( !util::IsDefined( RESTRICT_CLASS))
        {
          return res_perm;
        }
        return GetSegmentOverlap( SUBCLASS, SEQUENCE, RESTRICT_CLASS, res_perm.Second());
      }

      // get the predicted segments for this sequence
      SegmentContainer predicted_segments( GetSegments( SUBCLASS, SEQUENCE, this->GetPredictedClasses(), RESTRICT_PERMUTATION));

      // get the actual segments for the sequence
      const SegmentContainer &actual_segments( m_ActualSegments( SEQUENCE)( SUBCLASS));

      // test whether a single-class segment overlap was requested
      const bool single_class( util::IsDefined( RESTRICT_CLASS));

      // normalization factor; == sum of actual counts in all classes considered
      double normalization( 0.0);

      double segment_overlap( 0.0);

      // iterate over the actual segments
      for
      (
        SegmentContainer::const_iterator
          itr_actual_seg( actual_segments.Begin()), itr_actual_seg_end( actual_segments.End());
        itr_actual_seg != itr_actual_seg_end;
        ++itr_actual_seg
      )
      {
        bool found_one( false);
        const size_t actual_state( itr_actual_seg->First());

        // get the start and end of this segment
        const size_t &actual_start( itr_actual_seg->Second()), &actual_end( itr_actual_seg->Third());

        // compute the real length of this segment
        const size_t actual_length( actual_end - actual_start + 1);

        // iterate over the predicted segments
        for
        (
          SegmentContainer::const_iterator
            itr_pred_seg( predicted_segments.Begin()), itr_pred_seg_end( predicted_segments.End());
          itr_pred_seg != itr_pred_seg_end;
          ++itr_pred_seg
        )
        {
          // get the predicted segments state
          const size_t pred_state( itr_pred_seg->First());

          // continue if only considering one state and this is not the correct state
          if( single_class && pred_state != RESTRICT_CLASS && actual_state != RESTRICT_CLASS)
          {
            continue;
          }

          // get the start and end of this predicted segment
          const size_t &pred_start( itr_pred_seg->Second()), &pred_end( itr_pred_seg->Third());

          // skip non-overlapping ranges
          if( pred_end < actual_start || pred_start > actual_end)
          {
            continue;
          }

          // compute the predicted length of this segment
          const size_t pred_length( pred_end - pred_start + 1);

          // only consider matching states
          if( actual_state != pred_state)
          {
            continue;
          }

          found_one = true;

          // update normalization factor
          normalization += actual_length;

          // compute the intersection and union of the predicted and actual ranges
          size_t intersection_start( 0), intersection_end( 0), union_start( 0), union_end( 0);
          if( actual_end > pred_end)
          {
            union_end = actual_end;
            intersection_end = pred_end;
          }
          else
          {
            union_end = pred_end;
            intersection_end = actual_end;
          }
          if( actual_start < pred_start)
          {
            union_start = actual_start;
            intersection_start = pred_start;
          }
          else
          {
            union_start = pred_start;
            intersection_start = actual_start;
          }

          // determine the actual length of the intersection and the union
          const size_t intersection_length( intersection_end - intersection_start + 1);
          const size_t union_length( union_end - union_start + 1);

          // compute the delta, which is the number of residues that could be wrong due to trivial factors like an extra
          // residue or two at the ends of an SSE.  The delta length is restricted such that no segment overlap will be
          // > 1
          const size_t delta_length
          (
            std::min
            (
              size_t( union_length - intersection_length), // prevents segment overlap from going larger than 1
              std::min
              (
                intersection_length,
                size_t( std::min( pred_length, actual_length) / 2) // half of the predicted/actual segment length
              )
            )
          );

          segment_overlap += float( intersection_length + delta_length) * float( actual_length) / float( union_length);
//          if( !SEQUENCE)
//          {
//            BCL_Debug( pred_start);
//            BCL_Debug( pred_end);
//            BCL_Debug( actual_start);
//            BCL_Debug( actual_end);
//
//            BCL_Debug( intersection_length);
//            BCL_Debug( delta_length);
//            BCL_Debug( actual_length);
//            BCL_Debug( union_length);
//            BCL_Debug( float( intersection_length + delta_length) * float( actual_length) / float( union_length));
//            BCL_Debug( segment_overlap);
//          }
        }
        if( !found_one && ( actual_state == RESTRICT_CLASS || RESTRICT_CLASS == util::GetUndefined< size_t>()))
        {
          normalization += actual_length;
        }
      }

      if( normalization)
      {
        segment_overlap /= normalization;
      }
      else if( !m_ClassCounts( SEQUENCE, SUBCLASS) && predicted_segments.IsEmpty())
      {
        segment_overlap = 1.0;
      }
      else
      {
        segment_overlap = 0.0;
      }
      //if( !SEQUENCE)
      //{
      //  BCL_Debug( m_SequenceNames( SEQUENCE));
      //  BCL_Debug( SUBCLASS);
      //  BCL_Debug( this->GetSegmentString( SUBCLASS, SEQUENCE));
      //  BCL_Debug( predicted_segments);
      //  BCL_Debug( actual_segments);
      //  BCL_Debug( RESTRICT_CLASS);
      //  BCL_Debug( segment_overlap);
      //  BCL_Debug( normalization);
      //}
      return storage::Pair< float, size_t>( segment_overlap, RESTRICT_PERMUTATION);
    }

    //! @brief create the segment string, useful for checking results
    //! @param SUBCLASS subclass index; if just computing secondary structure using on method, this is always 0; but
    //!        when using multiple methods, it may range from 0-number classifications
    //! @param SEQUENCE_ID sequence index
    //! @return string with format: XXXXX ID,  PSEC,  OSEC
    std::string ObjectiveFunctionSegmentOverlap::GetSegmentString( const size_t &SUBCLASS, const size_t &SEQUENCE) const
    {
      std::ostringstream output;

      // get the actual segments for the sequence
      const SegmentContainer &actual_segments( m_ActualSegments( SEQUENCE)( SUBCLASS));
      SegmentContainer::const_iterator itr_actual_segments( actual_segments.Begin());

      // get the bounds for the sequence
      const size_t start_row( m_BoundaryIds( SEQUENCE)), end_row( m_BoundaryIds( SEQUENCE + 1));

      // get an iterator to the missing ids for this sequence
      storage::Vector< size_t>::const_iterator itr_missing( m_MissingIds( SEQUENCE).Begin());

      size_t missing_so_far( 0);
      for( size_t pos( start_row); pos < end_row; ++pos)
      {
        // test whether there was a gap between the previous row and this one
        if( *itr_missing == pos)
        {
          output << "XXXXX " << pos + missing_so_far << ",-1,-1\n";
          ++missing_so_far;
          ++itr_missing;
          ++itr_actual_segments;
        }
        if( itr_actual_segments->Third() < pos)
        {
          ++itr_actual_segments;
        }
        output << "XXXXX " << pos + missing_so_far << ','
               << this->GetPredictedClasses()( pos, SUBCLASS) << ','
               << itr_actual_segments->First() << '\n';
      }
      return output.str();
    }

  } // namespace model
} // namespace bcl
