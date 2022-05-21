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
#include "descriptor/bcl_descriptor_dataset.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Dataset::s_Instance
    (
      GetObjectInstances().AddInstance( new Dataset( 0, 0, 0, 0))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Dataset::Dataset() :
      m_Features( new model::FeatureDataSet< float>()),
      m_Results( new model::FeatureDataSet< float>()),
      m_Identification( new model::FeatureDataSet< char>())
    {
    }

    //! @brief constructor from sizes
    //! @param NR_EXAMPLES number of training examples the dataset will contain
    //! @param FEATURE_SIZE the size of each feature in the dataset
    //! @param RESULT_SIZE the size of each result in the dataset
    //! @param ID_SIZE the size of each id in the dataset
    Dataset::Dataset
    (
      const size_t &NR_EXAMPLES,
      const size_t &FEATURE_SIZE,
      const size_t &RESULT_SIZE,
      const size_t &ID_SIZE
    ) :
      m_Features( new model::FeatureDataSet< float>( NR_EXAMPLES, FEATURE_SIZE, util::GetUndefined< float>())),
      m_Results( new model::FeatureDataSet< float>( NR_EXAMPLES, RESULT_SIZE, util::GetUndefined< float>())),
      m_Identification( new model::FeatureDataSet< char>( NR_EXAMPLES, ID_SIZE, ' '))
    {
    }

    //! @brief constructor from labels
    //! @param NR_EXAMPLES number of training examples the dataset will contain
    //! @param FEATURE_LABELS labels for all features in the dataset
    //! @param RESULT_LABELS labels for each result in the dataset
    //! @param ID_LABELS labels for each id in the dataset
    Dataset::Dataset
    (
      const size_t &NR_EXAMPLES,
      const model::FeatureLabelSet &FEATURE_LABELS,
      const model::FeatureLabelSet &RESULT_LABELS,
      const model::FeatureLabelSet &ID_LABELS
    ) :
      m_Features( new model::FeatureDataSet< float>( NR_EXAMPLES, FEATURE_LABELS.GetSize(), util::GetUndefined< float>())),
      m_Results( new model::FeatureDataSet< float>( NR_EXAMPLES, RESULT_LABELS.GetSize(), util::GetUndefined< float>())),
      m_Identification( new model::FeatureDataSet< char>( NR_EXAMPLES, ID_LABELS.GetSize(), ' '))
    {
      m_Features->SetFeatureLabelSet( FEATURE_LABELS);
      m_Results->SetFeatureLabelSet( RESULT_LABELS);
      m_Identification->SetFeatureLabelSet( ID_LABELS);
    }

    //! @brief constructor
    //! @param FEATURES the features
    //! @param RESULTS the results
    Dataset::Dataset
    (
      const util::ShPtr< model::FeatureDataSet< float> > &FEATURES,
      const util::ShPtr< model::FeatureDataSet< float> > &RESULTS
    ) :
      m_Features( FEATURES),
      m_Results( RESULTS),
      m_Identification( new model::FeatureDataSet< char>())
    {
      // make sure that features and results have the same size and are either both defined or both undefined
      BCL_Assert
      (
        FEATURES->GetNumberFeatures() == RESULTS->GetNumberFeatures(),
        GetClassIdentifier() + " requires features and results to have the same size, but received "
        + util::Format()( FEATURES->GetNumberFeatures()) + " features and "
        + util::Format()( RESULTS->GetNumberFeatures()) + " results"
      );
    }

    //! @brief constructor
    //! @param FEATURES the features
    //! @param RESULTS the results
    Dataset::Dataset
    (
      const linal::MatrixConstInterface< float> &FEATURES,
      const linal::MatrixConstInterface< float> &RESULTS
    ) :
      m_Features( new model::FeatureDataSet< float>( FEATURES)),
      m_Results( new model::FeatureDataSet< float>( RESULTS)),
      m_Identification( new model::FeatureDataSet< char>())
    {
      // make sure that features and results have the same size and are either both defined or both undefined
      BCL_Assert
      (
        FEATURES.GetNumberRows() == RESULTS.GetNumberRows(),
        GetClassIdentifier() + " requires features and results to have the same size, but received "
        + util::Format()( FEATURES.GetNumberRows()) + " features and "
        + util::Format()( RESULTS.GetNumberRows()) + " results"
      );
    }

    //! @brief constructor
    //! @param FEATURES the features
    //! @param RESULTS the results
    //! @param IDS the ids
    Dataset::Dataset
    (
      const linal::MatrixConstInterface< float> &FEATURES,
      const linal::MatrixConstInterface< float> &RESULTS,
      const linal::MatrixConstInterface< char> &IDS
    ) :
      m_Features( new model::FeatureDataSet< float>( FEATURES)),
      m_Results( new model::FeatureDataSet< float>( RESULTS)),
      m_Identification( new model::FeatureDataSet< char>( IDS))
    {
      // make sure that features and results have the same size and are either both defined or both undefined
      BCL_Assert
      (
        FEATURES.GetNumberRows() == RESULTS.GetNumberRows(),
        GetClassIdentifier() + " requires features and results to have the same size, but received "
        + util::Format()( FEATURES.GetNumberRows()) + " features and "
        + util::Format()( RESULTS.GetNumberRows()) + " results"
      );
      BCL_Assert
      (
        FEATURES.GetNumberRows() == IDS.GetNumberRows(),
        GetClassIdentifier() + " requires features and ids to have the same size, but received "
        + util::Format()( FEATURES.GetNumberRows()) + " features and "
        + util::Format()( IDS.GetNumberRows()) + " ids"
      );
    }

    //! @brief constructor from vector of features and results
    //! @param DATA_SET the vector of features and results
    //! @note this constructor is DEPRECATED : use matrices in newly written code
    Dataset::Dataset
    (
      const storage::Vector< storage::VectorND< 2, linal::Vector< float> > > &DATA_SET
    ) :
      m_Features( new model::FeatureDataSet< float>()),
      m_Results( new model::FeatureDataSet< float>()),
      m_Identification( new model::FeatureDataSet< char>())
    {
      if( DATA_SET.IsEmpty())
      {
        return;
      }
      const size_t n_features( DATA_SET.GetSize());
      const size_t feature_size( DATA_SET.FirstElement().First().GetSize());
      const size_t result_size( DATA_SET.FirstElement().Second().GetSize());

      linal::Matrix< float> features( n_features, feature_size);
      linal::Matrix< float> results( n_features, result_size);

      for( size_t i( 0); i < n_features; ++i)
      {
        const linal::Vector< float> &feature( DATA_SET( i).First());
        const linal::Vector< float> &result( DATA_SET( i).Second());

        BCL_Assert
        (
          feature.GetSize() == feature_size,
          "Expected all features to be the same size as the first feature ( size=" + util::Format()( feature_size)
          + ") but feature #" + util::Format()( i) + " had size " + util::Format()( feature.GetSize())
        );
        BCL_Assert
        (
          result.GetSize() == result_size,
          "Expected all results to be the same size as the first result ( size=" + util::Format()( result_size)
          + ") but result #" + util::Format()( i) + " had size " + util::Format()( result.GetSize())
        );
        features.ReplaceRow( i, linal::VectorConstReference< float>( feature_size, feature.Begin()));
        results.ReplaceRow( i, linal::VectorConstReference< float>( result_size, result.Begin()));
      }
      m_Features = util::ShPtr< model::FeatureDataSetInterface< float> >( new model::FeatureDataSet< float>( features));
      m_Results = util::ShPtr< model::FeatureDataSetInterface< float> >( new model::FeatureDataSet< float>( results));
    }

    //! @brief Clone function
    //! @return pointer to new Dataset
    Dataset *Dataset::Clone() const
    {
      return new Dataset( *this);
    }

    //! @brief Hard copy -> copies internal data as well
    //! @return pointer to new Dataset
    Dataset *Dataset::HardCopy() const
    {
      Dataset *n_dataset( new Dataset( m_Features.HardCopy(), m_Results.HardCopy()));
      n_dataset->SetIds( m_Identification.HardCopy());
      return n_dataset;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Dataset::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return number of examples in the data set
    //! @return returns number of examples in the Dataset
    size_t Dataset::GetSize() const
    {
      return m_Features->GetNumberFeatures();
    }

    //! @brief return the size of each feature in the dataset
    //! @return return the size of each feature in the dataset
    size_t Dataset::GetFeatureSize() const
    {
      return m_Features->GetFeatureSize();
    }

    //! @brief return the size of each id (# characters) in the dataset
    //! @return return the size of each id (# characters) in the dataset
    size_t Dataset::GetIdSize() const
    {
      return m_Identification->GetFeatureSize();
    }

    //! @brief return the size of each result in the dataset
    //! @return return the size of each result in the dataset
    size_t Dataset::GetResultSize() const
    {
      return m_Results->GetFeatureSize();
    }

    //! @brief check whether features or results contain a number of zero rows
    //! @return true if features or results have zero rows in matrix
    bool Dataset::IsEmpty() const
    {
      return GetSize() == 0;
    }

    //! @brief returns the features pointer
    //! @return features pointer
    const util::ShPtr< model::FeatureDataSet< float> > &Dataset::GetFeaturesPtr() const
    {
      return m_Features;
    }

    //! @brief returns the results pointer
    //! @return the results pointer
    const util::ShPtr< model::FeatureDataSet< float> > &Dataset::GetResultsPtr() const
    {
      return m_Results;
    }

    //! @brief returns the ids pointer
    //! @return the ids pointer
    const util::ShPtr< model::FeatureDataSet< char> > &Dataset::GetIdsPtr() const
    {
      return m_Identification;
    }

    //! @brief returns the features pointer
    //! @return features pointer
    model::FeatureDataSet< float> &Dataset::GetFeatures()
    {
      return *m_Features;
    }

    //! @brief returns the results pointer
    //! @return the results pointer
    model::FeatureDataSet< float> &Dataset::GetResults()
    {
      return *m_Results;
    }

    //! @brief returns the ids dataset
    //! @return the ids dataset
    model::FeatureDataSet< char> &Dataset::GetIds()
    {
      return *m_Identification;
    }

    //! @brief access to the features matrix
    //! @return a features matrix reference
    linal::MatrixReference< float> Dataset::GetFeaturesReference()
    {
      return m_Features->GetMatrix();
    }

    //! @brief access to the features matrix
    //! @return a features matrix reference
    linal::MatrixReference< float> Dataset::GetResultsReference()
    {
      return m_Results->GetMatrix();
    }

    //! @brief access to the features matrix
    //! @return a features matrix reference
    linal::MatrixReference< char> Dataset::GetIdsReference()
    {
      return m_Identification->GetMatrix();
    }

    //! @brief const access to the features matrix
    //! @return a features matrix reference
    linal::MatrixConstReference< float> Dataset::GetFeaturesReference() const
    {
      return m_Features->GetMatrix();
    }

    //! @brief const access to the features matrix
    //! @return a features matrix reference
    linal::MatrixConstReference< float> Dataset::GetResultsReference() const
    {
      return m_Results->GetMatrix();
    }

    //! @brief const access to the features matrix
    //! @return an id matrix reference
    linal::MatrixConstReference< char> Dataset::GetIdsReference() const
    {
      return m_Identification->GetMatrix();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief remove all undefined examples
    //! @return the number of examples actually removed
    size_t Dataset::RemoveUndefinedExamples()
    {
      storage::Vector< size_t> undefined_row_ids( linal::FindUndefinedRows( GetFeaturesReference()));
      undefined_row_ids.Append( linal::FindCompletelyUndefinedRows( GetResultsReference()));

      const size_t old_size( GetSize());
      // remove the undefined rows from all three matrices
      m_Features->GetRawMatrix().RemoveRows( undefined_row_ids);
      m_Results->GetRawMatrix().RemoveRows( undefined_row_ids);
      if( m_Identification->GetFeatureSize())
      {
        const linal::Matrix< char> &ids( m_Identification->GetRawMatrix());
        for
        (
          storage::Vector< size_t>::const_iterator itr( undefined_row_ids.Begin()), itr_end( undefined_row_ids.End());
          itr != itr_end;
          ++itr
        )
        {
          BCL_MessageVrb
          (
            "Removing row with undefined values: "
            + std::string( ids.GetRow( *itr).Begin(), m_Identification->GetFeatureSize())
          );
        }
      }
      m_Identification->GetRawMatrix().RemoveRows( undefined_row_ids);
      return old_size - GetSize();
    }

    //! @brief append another data set onto this one
    //! @param DATASET the dataset to append
    void Dataset::Append( const Dataset &DATASET)
    {
      m_Features->GetRawMatrix().Append( DATASET.GetFeaturesReference());
      m_Results->GetRawMatrix().Append( DATASET.GetResultsReference());
      m_Identification->GetRawMatrix().Append( DATASET.GetIdsReference());
    }

    //! @brief add a specified number of extra rows to this dataset
    //! @param N_ROWS the number of rows to add
    void Dataset::AddRows( const size_t &N_ROWS)
    {
      if( N_ROWS == size_t( 0))
      {
        return;
      }

      m_Features->GetRawMatrix().AddNRows( N_ROWS, util::GetUndefined< float>());
      m_Results->GetRawMatrix().AddNRows( N_ROWS, util::GetUndefined< float>());
      m_Identification->GetRawMatrix().AddNRows( N_ROWS, ' ');
    }

    //! @brief shrink rows (remove unused rows at the end of the dataset)
    //! @param NEW_SIZE the new size of the matrix
    void Dataset::ShrinkRows( const size_t &NEW_SIZE)
    {
      if( NEW_SIZE < GetSize())
      {
        m_Features->GetRawMatrix().ShrinkRows( NEW_SIZE);
        m_Results->GetRawMatrix().ShrinkRows( NEW_SIZE);
        m_Identification->GetRawMatrix().ShrinkRows( NEW_SIZE);
      }
    }

    //! @brief keep only the specified rows
    //! @param KEEPERS indices of rows to keep
    void Dataset::KeepRows( const storage::Vector< size_t> &KEEPERS)
    {
      m_Features->GetRawMatrix().KeepRows( KEEPERS);
      m_Results->GetRawMatrix().KeepRows( KEEPERS);
      m_Identification->GetRawMatrix().KeepRows( KEEPERS);
    }

    //! @brief get a new dataset consisting of only the specified rows
    //! @param ROWS the rows desired
    util::ShPtr< Dataset> Dataset::GetRows( const storage::Vector< size_t> &KEEPERS) const
    {
      util::ShPtr< Dataset> new_d
      (
        new Dataset
        (
          KEEPERS.GetSize(),
          *m_Features->GetFeatureLabelSet(),
          *m_Results->GetFeatureLabelSet(),
          m_Identification->GetFeatureLabelSet().IsDefined() ? *m_Identification->GetFeatureLabelSet() : model::FeatureLabelSet()
        )
      );
      for( size_t i( 0), sz( KEEPERS.GetSize()); i < sz; ++i)
      {
        new_d->GetFeaturesReference().GetRow( i).CopyValues( m_Features->GetMatrix().GetRow( KEEPERS( i)));
        new_d->GetResultsReference().GetRow( i).CopyValues( m_Results->GetMatrix().GetRow( KEEPERS( i)));
        if( m_Identification.IsDefined() && m_Identification->GetFeatureSize())
        {
          new_d->GetIdsReference().GetRow( i).CopyValues( m_Identification->GetMatrix().GetRow( KEEPERS( i)));
        }
      }
      return new_d;
    }

    //! @brief Shuffle the rows randomly
    void Dataset::Shuffle()
    {
      storage::Vector< size_t> ordering( GetSize());
      for( size_t i( 0), sz( GetSize()); i < sz; ++i)
      {
        ordering( i) = i;
      }
      ordering.Shuffle();
      m_Features->GetRawMatrix().ReorderRows( ordering);
      m_Results->GetRawMatrix().ReorderRows( ordering);
      m_Identification->GetRawMatrix().ReorderRows( ordering);
    }

    //! @brief Shuffle only the results (y-scrambling)
    void Dataset::YScramble()
    {
      storage::Vector< size_t> ordering( GetSize());
      for( size_t i( 0), sz( GetSize()); i < sz; ++i)
      {
        ordering( i) = i;
      }
      ordering.Shuffle();
      m_Results->GetRawMatrix().ReorderRows( ordering);
    }

    //! @brief add data for a specific row
    //! @param ROW the row to add data to
    //! @param FEATURE the feature to add
    //! @param RESULT the result to add
    //! @param ID the id information to add
    void Dataset::AddData
    (
      const size_t &ROW,
      const linal::VectorConstInterface< float> &FEATURE,
      const linal::VectorConstInterface< float> &RESULT,
      const linal::VectorConstInterface< char> &ID
    )
    {
      m_Features->GetRawMatrix().ReplaceRow( ROW, FEATURE);
      m_Results->GetRawMatrix().ReplaceRow( ROW, RESULT);
      m_Identification->GetRawMatrix().ReplaceRow( ROW, ID);
    }

    //! @brief set features
    //! @param FEATURES new model::FeatureDataSet
    void Dataset::SetIds( const util::ShPtr< model::FeatureDataSet< char> > &IDS)
    {
      BCL_Assert
      (
        IDS->GetNumberFeatures() == m_Results->GetNumberFeatures(),
        "Number of new ids has to be the same as existing ids!"
      );

      m_Identification = IDS;
    }

    //! @brief set features
    //! @param FEATURES new model::FeatureDataSet
    void Dataset::SetFeatures( const util::ShPtr< model::FeatureDataSet< float> > &FEATURES)
    {
      BCL_Assert
      (
        FEATURES->GetNumberFeatures() == m_Results->GetNumberFeatures(),
        "Number of new features has to be the same as existing features!"
      );

      m_Features = FEATURES;
    }

    //! @brief set results
    //! @param RESULTS new model::FeatureDataSet
    void Dataset::SetResults( const util::ShPtr< model::FeatureDataSet< float> > &RESULTS)
    {
      BCL_Assert
      (
        RESULTS->GetNumberFeatures() == m_Features->GetNumberFeatures(),
        "Number of new results has to be the same as existing results!"
      );

      m_Results = RESULTS;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Dataset::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Features, ISTREAM);
      io::Serialize::Read( m_Results, ISTREAM);
      io::Serialize::Read( m_Identification, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Dataset::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Features, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Results, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Identification, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace descriptor
} // namespace bcl
