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

#ifndef BCL_DESCRIPTOR_DATASET_BUILDER_HPP_
#define BCL_DESCRIPTOR_DATASET_BUILDER_HPP_

// include the header of this class
#include "bcl_descriptor_dataset_builder.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_dataset.h"
#include "bcl_descriptor_iterator.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param FEATURE descriptor for the features
    //! @param RESULT descriptor for the results
    //! @param ID descriptor for the ID
    template< typename t_DataType>
    DatasetBuilder< t_DataType>::DatasetBuilder()
    {
    }

    //! @brief constructor from info
    //! @param FEATURE descriptor for the features
    //! @param RESULT descriptor for the results
    //! @param ID descriptor for the ID
    template< typename t_DataType>
    DatasetBuilder< t_DataType>::DatasetBuilder
    (
      const util::ObjectDataLabel &FEATURE,
      const util::ObjectDataLabel &RESULT,
      const util::ObjectDataLabel &ID
    ) :
      m_Feature(),
      m_Result(),
      m_Identification()
    {
      SetFeatureCode( FEATURE);
      SetResultsCode( RESULT);
      SetIdCode( ID);
    }

    //! @brief Clone function
    //! @return pointer to new DatasetBuilder
    template< typename t_DataType>
    DatasetBuilder< t_DataType> *DatasetBuilder< t_DataType>::Clone() const
    {
      return new DatasetBuilder< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &DatasetBuilder< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the size of each feature in the dataset
    //! @return return the size of each feature in the dataset
    template< typename t_DataType>
    size_t DatasetBuilder< t_DataType>::GetFeatureSize() const
    {
      return m_Feature.GetSizeOfFeatures();
    }

    //! @brief return the size of each id (# characters) in the dataset
    //! @return return the size of each id (# characters) in the dataset
    template< typename t_DataType>
    size_t DatasetBuilder< t_DataType>::GetIdSize() const
    {
      return m_Identification.GetSizeOfFeatures();
    }

    //! @brief return the size of each result in the dataset
    //! @return return the size of each result in the dataset
    template< typename t_DataType>
    size_t DatasetBuilder< t_DataType>::GetResultSize() const
    {
      return m_Result.GetSizeOfFeatures();
    }

    //! @brief set the feature code
    //! @param FEATURES the desired feature code
    template< typename t_DataType>
    void DatasetBuilder< t_DataType>::SetFeatureCode( const util::ObjectDataLabel &FEATURES)
    {
      if( FEATURES.IsEmpty())
      {
        return;
      }
      util::OwnPtr< const util::ObjectDataLabel> features_ptr( &FEATURES, false);
      if( FEATURES.GetValue() != m_Feature.GetAlias())
      {
        features_ptr =
          util::OwnPtr< const util::ObjectDataLabel>
          (
            new util::ObjectDataLabel( "", m_Feature.GetAlias(), storage::Vector< util::ObjectDataLabel>( 1, FEATURES)),
            true
          );
      }
      io::ValidationResult val( m_Feature.ValidateRead( *features_ptr, util::GetLogger()));
      if( !val.IsAllowed())
      {
        if( val.IsHelp())
        {
          BCL_ExitWithoutCallstack( "", 0);
        }
        BCL_Exit( "Feature label could not be read", -1);
      }
      Type type( m_Type);
      type.GeneralizeToHandle( m_Feature.GetType());
      if( type.GetDimension() != m_Type.GetDimension())
      {
        m_Type = type;
        m_Result.SetDimension( m_Type.GetDimension());
        m_Identification.SetDimension( m_Type.GetDimension());
      }
      else
      {
        m_Feature.SetDimension( type.GetDimension());
      }
      m_FeatureLabels = util::ShPtr< model::FeatureLabelSet>( m_Feature.GetLabelsWithSizes().Clone());
    }

    //! @brief set the results code
    //! @param RESULTS the desired results code
    template< typename t_DataType>
    void DatasetBuilder< t_DataType>::SetResultsCode( const util::ObjectDataLabel &RESULTS)
    {
      if( RESULTS.IsEmpty())
      {
        return;
      }
      util::OwnPtr< const util::ObjectDataLabel> results_ptr( &RESULTS, false);
      if( RESULTS.GetValue() != m_Result.GetAlias())
      {
        results_ptr =
          util::OwnPtr< const util::ObjectDataLabel>
          (
            new util::ObjectDataLabel( "", m_Result.GetAlias(), storage::Vector< util::ObjectDataLabel>( 1, RESULTS)),
            true
          );
      }
      io::ValidationResult val( m_Result.ValidateRead( *results_ptr, util::GetLogger()));
      if( !val.IsAllowed())
      {
        if( val.IsHelp())
        {
          BCL_ExitWithoutCallstack( "", 0);
        }
        BCL_Exit( "Result label could not be read", -1);
      }
      Type type( m_Type);
      type.GeneralizeToHandle( m_Result.GetType());
      if( type.GetDimension() != m_Type.GetDimension())
      {
        m_Type = type;
        m_Feature.SetDimension( m_Type.GetDimension());
        m_Identification.SetDimension( m_Type.GetDimension());
      }
      else
      {
        m_Result.SetDimension( type.GetDimension());
      }
      m_ResultLabels = util::ShPtr< model::FeatureLabelSet>( m_Result.GetLabelsWithSizes().Clone());
    }

    //! @brief set the id code
    //! @param IDS the desired id code
    template< typename t_DataType>
    void DatasetBuilder< t_DataType>::SetIdCode( const util::ObjectDataLabel &IDS)
    {
      if( IDS.IsEmpty())
      {
        return;
      }
      util::OwnPtr< const util::ObjectDataLabel> ids_ptr( &IDS, false);
      if( IDS.GetValue() != m_Result.GetAlias())
      {
        ids_ptr =
          util::OwnPtr< const util::ObjectDataLabel>
          (
            new util::ObjectDataLabel( "", m_Result.GetAlias(), storage::Vector< util::ObjectDataLabel>( 1, IDS)),
            true
          );
      }
      io::ValidationResult val( m_Identification.ValidateRead( *ids_ptr, util::GetLogger()));
      if( !val.IsAllowed())
      {
        if( val.IsHelp())
        {
          BCL_ExitWithoutCallstack( "", 0);
        }
        BCL_Exit( "Result label could not be read", -1);
      }
      Type type( m_Type);
      type.GeneralizeToHandle( m_Identification.GetType());
      if( type.GetDimension() != m_Type.GetDimension())
      {
        m_Type = type;
        m_Feature.SetDimension( m_Type.GetDimension());
        m_Result.SetDimension( m_Type.GetDimension());
      }
      else
      {
        m_Identification.SetDimension( type.GetDimension());
      }
      m_IdLabels = util::ShPtr< model::FeatureLabelSet>( m_Identification.GetLabelsWithSizes().Clone());
    }

    //! @brief return the feature code
    //! @return return the feature code
    template< typename t_DataType>
    const Combine< t_DataType, float> &DatasetBuilder< t_DataType>::GetFeatureCode() const
    {
      return m_Feature;
    }

    //! @brief return the result code
    //! @return return the result code
    template< typename t_DataType>
    const Combine< t_DataType, float> &DatasetBuilder< t_DataType>::GetResultCode() const
    {
      return m_Result;
    }

    //! @brief return the id code
    //! @return return the id code
    template< typename t_DataType>
    const Combine< t_DataType, char> &DatasetBuilder< t_DataType>::GetIdCode() const
    {
      return m_Identification;
    }

    //! @brief Type of dataset to build (e.g. atom-based pairwise)
    //! @return the type of dataset that will be built
    template< typename t_DataType>
    Type DatasetBuilder< t_DataType>::GetType() const
    {
      return m_Type;
    }

    //! @brief Override the computed type of descriptor
    //! @param TYPE the actual type of descriptors desired
    template< typename t_DataType>
    void DatasetBuilder< t_DataType>::SetType( const Type &TYPE)
    {
      if( TYPE != m_Type)
      {
        m_Type = TYPE;
        m_Feature.SetDimension( m_Type.GetDimension());
        m_Result.SetDimension( m_Type.GetDimension());
        m_Identification.SetDimension( m_Type.GetDimension());
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief create a dataset from the given sequence interface
    //! @param SEQUENCE sequence of interest
    template< typename t_DataType>
    Dataset DatasetBuilder< t_DataType>::operator()( const SequenceInterface< t_DataType> &SEQUENCE)
    {
      Iterator< t_DataType> itr( m_Type, SEQUENCE);
      Dataset dataset( itr.GetSize(), GetFeatureSize(), GetResultSize(), GetIdSize());
      dataset.GetFeatures().SetFeatureLabelSet( m_FeatureLabels);
      dataset.GetResults().SetFeatureLabelSet( m_ResultLabels);
      dataset.GetIds().SetFeatureLabelSet( m_IdLabels);
      m_Feature.SetObject( SEQUENCE);
      m_Result.SetObject( SEQUENCE);
      m_Identification.SetObject( SEQUENCE);
      for( ; itr.NotAtEnd(); ++itr)
      {
        dataset.AddData( itr.GetPosition(), m_Feature( itr), m_Result( itr), m_Identification( itr));
      }
      dataset.RemoveUndefinedExamples();
      return dataset;
    }

    //! @brief create a dataset from the given sequence interfaces
    //! @param ITERATOR iterator to the sequence interfaces
    //! @param START start feature inside the given iterator
    //! @param MAX_FEATURES maximum number of features to generate
    template< typename t_DataType>
    Dataset DatasetBuilder< t_DataType>::operator()
    (
      const iterate::Generic< const SequenceInterface< t_DataType> > &ITERATOR,
      const size_t &START,
      const size_t &MAX_FEATURES
    )
    {
      // begin by determining the total size of the dataset
      size_t n_features( 0), features_before_start( 0);

      // create a dataset to hold all the features
      for
      (
        iterate::Generic< const SequenceInterface< t_DataType> > itr_seq( ITERATOR);
        itr_seq.NotAtEnd() && n_features < MAX_FEATURES;
        ++itr_seq
      )
      {
        size_t sequence_features( m_Type.GetNumberFeatures( itr_seq->GetSize()));
        if( features_before_start < START)
        {
          if( features_before_start + sequence_features <= START)
          {
            features_before_start += sequence_features;
            sequence_features = 0;
          }
          else
          {
            sequence_features = features_before_start + sequence_features - START;
            features_before_start = START;
          }
        }
        n_features += sequence_features;
      }
      n_features = std::min( n_features, MAX_FEATURES);

      // create a dataset to hold all the features
      Dataset dataset( n_features, GetFeatureSize(), GetResultSize(), GetIdSize());
      dataset.GetFeatures().SetFeatureLabelSet( m_FeatureLabels);
      dataset.GetResults().SetFeatureLabelSet( m_ResultLabels);
      dataset.GetIds().SetFeatureLabelSet( m_IdLabels);

      // create a reusable iterator
      Iterator< t_DataType> itr( m_Type);

      // keep track of the number of features seen so far
      size_t overall_pos( 0), main_itr_pos( 0);
      util::SiPtr< const SequenceInterface< t_DataType> > previous_seq_ptr;
      for( iterate::Generic< const SequenceInterface< t_DataType> > itr_seq( ITERATOR); itr_seq.NotAtEnd(); ++itr_seq)
      {
        // call set object on all relevant objects
        itr.SetObject( *itr_seq);
        if( itr.GetSize() + main_itr_pos < START)
        {
          main_itr_pos += itr.GetSize();
          continue;
        }

        // move forward to the first descriptor to generate
        while( main_itr_pos < START)
        {
          ++main_itr_pos;
          ++itr;
        }

        m_Feature.SetObject( *itr_seq);
        m_Result.SetObject( *itr_seq);
        m_Identification.SetObject( *itr_seq);
        // generate the descriptors for this object
        for( ; itr.NotAtEnd() && overall_pos < MAX_FEATURES; ++itr, ++overall_pos)
        {
          dataset.AddData( overall_pos, m_Feature( itr), m_Result( itr), m_Identification( itr));
        }
        if( previous_seq_ptr.IsDefined())
        {
          previous_seq_ptr->ResetCache();
        }
        previous_seq_ptr = util::ToSiPtr( *itr_seq);
      }

      const size_t num_initial_features( dataset.GetSize());
      const size_t num_ommitted( dataset.RemoveUndefinedExamples());
      // write a message if any features were omitted
      BCL_Message
      (
        num_ommitted ? util::Message::e_Standard : util::Message::e_Verbose,
        "Generated " + util::Format()( num_initial_features) + " features, removed " + util::Format()( num_ommitted)
        + " due to presence of undefined values, " + util::Format()( dataset.GetSize()) + " features remain"
      );
      return dataset;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &DatasetBuilder< t_DataType>::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Feature, ISTREAM);
      io::Serialize::Read( m_Identification, ISTREAM);
      io::Serialize::Read( m_Result, ISTREAM);
      io::Serialize::Read( m_Type, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_DataType>
    std::ostream &DatasetBuilder< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Feature, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Identification, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Result, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Type, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_DATASET_BUILDER_H_
