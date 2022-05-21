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
#include "model/bcl_model_cross_validation_info.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_retrieve_data_set_from_delimited_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CrossValidationInfo::s_Instance
    (
      GetObjectInstances().AddInstance( new CrossValidationInfo())
    );

    //! @brief Name for each prediction file format
    //! @param FORMAT the format to convert to a string
    //! @return string representing FORMAT
    const std::string &CrossValidationInfo::GetPredictionOutputFormatString( const PredictionOutputFormat &FORMAT)
    {
      static const std::string s_names[] =
      {
        "",
        "Matrix",
        "CSV",
        GetStaticClassName< PredictionOutputFormat>()
      };
      return s_names[ size_t( FORMAT)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CrossValidationInfo::CrossValidationInfo() :
      m_Result( util::GetUndefined< float>()),
      m_ImprovementType(),
      m_Independent(),
      m_Monitoring(),
      m_Training(),
      m_NumberOutputs( 0),
      m_PredictionOutputFormat( e_ToBeDetermined),
      m_Filename()
    {
    }

    //! @brief default constructor
    CrossValidationInfo::CrossValidationInfo( const std::string &FILENAME) :
      m_Result( util::GetUndefined< float>()),
      m_ImprovementType(),
      m_Independent(),
      m_Monitoring(),
      m_Training(),
      m_NumberOutputs( 0),
      m_PredictionOutputFormat( e_ToBeDetermined),
      m_Filename( FILENAME)
    {
      io::IFStream input;
      io::File::MustOpenIFStream( input, FILENAME);
      // read this object from an object data label
      io::Serialize::Read( *this, input);
      io::File::CloseClearFStream( input);

      // commonly, if the models are used later on and are moved somewhere, the predictions files are moved
      // into the same directory as the models
      if( !m_IndependentPredictions.empty() && !io::DirectoryEntry( m_IndependentPredictions).DoesExist())
      {
        const std::string basename( io::File::SplitToPathAndFileName( m_IndependentPredictions).Second());
        const std::string initial_path( io::File::SplitToPathAndFileName( m_Filename).First());
        if( io::DirectoryEntry( initial_path + "/" + basename).DoesExist())
        {
          BCL_MessageVrb
          (
            "Could not locate independent predictions at: " + m_IndependentPredictions
            + "; using file from " + initial_path + " with same name"
          );
          m_IndependentPredictions = initial_path + "/" + basename;
        }
      }
    }

    //! @brief constructor from members
    CrossValidationInfo::CrossValidationInfo
    (
      const float &RESULT,
      const opti::ImprovementTypeEnum &IMPROVEMENT_TYPE,
      const util::ObjectDataLabel &INDEPENDENT,
      const util::ObjectDataLabel &MONITORING,
      const util::ObjectDataLabel &TRAINING,
      const util::ObjectDataLabel &OBJECTIVE,
      const util::ObjectDataLabel &ITERATE,
      const std::string &INDEPENDENT_PREDICTIONS,
      const FeatureLabelSet &ID_LABELS,
      const size_t &NUMBER_OUTPUTS
    ) :
      m_Result( RESULT),
      m_ImprovementType( IMPROVEMENT_TYPE),
      m_Independent( INDEPENDENT),
      m_Monitoring( MONITORING),
      m_Training( TRAINING),
      m_Objective( OBJECTIVE),
      m_Iterate( ITERATE),
      m_IndependentPredictions( INDEPENDENT_PREDICTIONS),
      m_NumberOutputs( NUMBER_OUTPUTS),
      m_IdLabels( ID_LABELS.GetLabel()),
      m_CharsPerId( ID_LABELS.GetPropertySizes()),
      m_PredictionOutputFormat( e_Delimited) // all new files will be written out with delimited format
    {
      if( !INDEPENDENT_PREDICTIONS.empty())
      {
        m_IndependentPredictions = io::File::MakeAbsolutePath( m_IndependentPredictions);
      }
    }

    //! @brief Clone function
    //! @return pointer to new CrossValidationInfo
    CrossValidationInfo *CrossValidationInfo::Clone() const
    {
      return new CrossValidationInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CrossValidationInfo::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &CrossValidationInfo::GetAlias() const
    {
      static const std::string s_name( "CVMetadata");
      return s_name;
    }

    //! @brief get the id feature label set
    //! @return the id feature label set
    FeatureLabelSet CrossValidationInfo::GetIDsFeatureLabelSet() const
    {
      return FeatureLabelSet( "Combine", m_IdLabels.GetArguments(), m_CharsPerId);
    }

    //! @brief keep only the information in this cv-info that is held in common with OTHER
    //! @param OTHER the other cross-validation info object
    void CrossValidationInfo::KeepSharedInfo( const CrossValidationInfo &OTHER)
    {
      m_Result = m_Result == OTHER.m_Result ? m_Result : util::GetUndefined< float>();
      BCL_Assert( m_ImprovementType == OTHER.m_ImprovementType, "Cannot combine CV-ids with different improvement types");
      m_Independent = m_Independent == OTHER.m_Independent ? m_Independent : util::ObjectDataLabel();
      m_Monitoring = m_Monitoring == OTHER.m_Monitoring ? m_Monitoring : util::ObjectDataLabel();
      m_Training = m_Training == OTHER.m_Training ? m_Training : util::ObjectDataLabel();
      m_Objective = m_Objective == OTHER.m_Objective ? m_Objective : util::ObjectDataLabel();
      m_Iterate = m_Independent == OTHER.m_Iterate ? m_Iterate : util::ObjectDataLabel();
      m_IndependentPredictions = m_IndependentPredictions == OTHER.m_IndependentPredictions ? m_IndependentPredictions : std::string();
      m_NumberOutputs = m_NumberOutputs == OTHER.m_NumberOutputs ? m_NumberOutputs : util::GetUndefined< size_t>();
      m_IdLabels = m_IdLabels == OTHER.m_IdLabels ? m_IdLabels : util::ObjectDataLabel();
      m_CharsPerId = m_CharsPerId == OTHER.m_CharsPerId ? m_CharsPerId : storage::Vector< size_t>();
      m_PredictionOutputFormat = m_PredictionOutputFormat == OTHER.m_PredictionOutputFormat ? m_PredictionOutputFormat.GetEnum() : e_Delimited;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CrossValidationInfo::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Cross-validation related data for a model");
      parameters.AddInitializer
      (
        "independent",
        "independent dataset retriever",
        io::Serialization::GetAgent( &m_Independent)
      );
      parameters.AddInitializer
      (
        "monitoring",
        "monitoring dataset retriever",
        io::Serialization::GetAgent( &m_Monitoring)
      );
      parameters.AddInitializer
      (
        "training",
        "training dataset retriever",
        io::Serialization::GetAgent( &m_Training),
        ""
      );
      parameters.AddInitializer
      (
        "result",
        "Objective function performance on the independent dataset",
        io::Serialization::GetAgent( &m_Result)
      );
      parameters.AddInitializer
      (
        "improvement type",
        "Directionality of results e.g. LargerIsBetter indicates that larger result values are better",
        io::Serialization::GetAgent( &m_ImprovementType)
      );
      parameters.AddInitializer
      (
        "objective",
        "Final objective function that was trained on",
        io::Serialization::GetAgent( &m_Objective),
        ""
      );
      parameters.AddInitializer
      (
        "iterate",
        "Actual optimization method",
        io::Serialization::GetAgent( &m_Iterate),
        ""
      );
      // Formerly, GetAgentInputFilename was used for validation, but when models are moved around the independent predictions will
      // necessarily do the same so we perform a limited search for them in ReadInitializerSuccessHook if the given path doesn't exist
      parameters.AddOptionalInitializer
      (
        "independent predictions",
        "File containing the independent predictions for this training",
        io::Serialization::GetAgent( &m_IndependentPredictions)
      );
      parameters.AddInitializer
      (
        "number outputs",
        "# of outputs (0 if unknown)",
        io::Serialization::GetAgent( &m_NumberOutputs),
        "0"
      );
      parameters.AddInitializer
      (
        "ids",
        "ID labels specified for the training",
        io::Serialization::GetAgent( &m_IdLabels),
        ""
      );
      parameters.AddInitializer
      (
        "id sizes",
        "Sizes for each id label",
        io::Serialization::GetAgentWithSizeLimits( &m_CharsPerId, 0, 10000),
        ""
      );
      parameters.AddInitializer
      (
        "prediction output format",
        "Format for the prediction output file",
        io::Serialization::GetAgent( &m_PredictionOutputFormat),
        GetPredictionOutputFormatString( e_ToBeDetermined)
      );
      return parameters;
    }

    //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool CrossValidationInfo::ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      // commonly, if the models are used later on and are moved somewhere, the predictions files are moved
      // into the same directory as the models
      if( !m_IndependentPredictions.empty() && !io::DirectoryEntry( m_IndependentPredictions).DoesExist())
      {
        const std::string basename( io::File::SplitToPathAndFileName( m_IndependentPredictions).Second());
        const std::string initial_path( io::File::SplitToPathAndFileName( m_Filename).First());
        if( io::DirectoryEntry( initial_path + "/" + basename).DoesExist())
        {
          BCL_MessageVrb
          (
            "Could not locate independent predictions at: " + m_IndependentPredictions
            + "; using file from " + initial_path + " with same name"
          );
          m_IndependentPredictions = initial_path + "/" + basename;
        }
        else
        {
          ERR_STREAM << "Could not location predictions file at: "
                     << m_IndependentPredictions << " or " << initial_path + "/" + basename << '\n';
          return false;
        }
      }
      return true;
    }

    //! @brief read predictions from any file
    //! @param FILENAME the filename
    util::ShPtr< descriptor::Dataset> CrossValidationInfo::ReadPredictions( const std::string &FILENAME)
    {
      const std::string info_filename( FILENAME + ".info");
      CrossValidationInfo info;
      if( io::DirectoryEntry( info_filename).DoesExist())
      {
        // read the info file
        io::IFStream input;
        io::File::MustOpenIFStream( input, info_filename);
        input >> info;
        io::File::CloseClearFStream( input);
      }
      info.SetIndependentPredictionsFilename( FILENAME);
      return info.ReadIndependentPredictions();
    }

    //! @brief read the independent dataset predictions into a dataset object, whose features are the predicted values
    //!        results are the same result values, and ids are the given ids
    util::ShPtr< descriptor::Dataset> CrossValidationInfo::ReadIndependentPredictions()
    {
      util::ShPtr< descriptor::Dataset> sp_dataset;
      io::DirectoryEntry directory_entry( m_IndependentPredictions);
      if
      (
        m_IndependentPredictions.empty()
        || !directory_entry.DoesExist()
        || directory_entry.IsType( io::Directory::e_Dir)
      )
      {
        // return if no output location was given
        return sp_dataset;
      }
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_IndependentPredictions);
      if( m_PredictionOutputFormat == e_ToBeDetermined)
      {
        // determine the input format
        std::string first_line;
        std::getline( input, first_line);
        if( util::TrimString( first_line) == linal::Matrix< float>().GetClassIdentifier())
        {
          m_PredictionOutputFormat = e_Matrix;
        }
        else
        {
          m_PredictionOutputFormat = e_Delimited;
          const size_t number_commas( std::count( first_line.begin(), first_line.end(), ','));
          BCL_Assert( number_commas, "Bad prediction file!");
          if( number_commas % 2)
          {
            m_IdLabels = "Combine";
            m_CharsPerId.Reset();
            m_NumberOutputs = ( number_commas + 1) / 2;
          }
          else
          {
            m_NumberOutputs = number_commas / 2;
            m_IdLabels = "Combine(UnknownIdChar)";
            m_CharsPerId.Reset();
            m_CharsPerId.PushBack( first_line.find( ','));
          }
        }

        // close and reopen file
        io::File::CloseClearFStream( input);
        io::File::MustOpenIFStream( input, m_IndependentPredictions);
      }

      if( m_PredictionOutputFormat == e_Matrix)
      {
        // matrix format
        linal::Matrix< float> experimental_predictions;
        input >> experimental_predictions;
        const size_t feature_size( experimental_predictions.GetNumberCols() / 2);
        const size_t n_training_examples( experimental_predictions.GetNumberRows());
        sp_dataset =
          util::ShPtr< descriptor::Dataset>
          (
            new descriptor::Dataset( n_training_examples, feature_size, feature_size, size_t( 0))
          );
        linal::MatrixReference< float> predictions_ref( sp_dataset->GetFeaturesReference());
        linal::MatrixReference< float> experimental_ref( sp_dataset->GetResultsReference());
        for( size_t row_id( 0); row_id < n_training_examples; ++row_id)
        {
          linal::VectorConstReference< float> ep( experimental_predictions.GetRow( row_id));
          experimental_ref.GetRow( row_id).CopyValues( ep.Reference( 0, feature_size));
          predictions_ref.GetRow( row_id).CopyValues( ep.Reference( feature_size, feature_size));
        }
      }
      else
      {
        // retrieve as a delimited file
        RetrieveDataSetFromDelimitedFile retriever
        (
          m_IndependentPredictions,
          m_NumberOutputs,
          GetIDsFeatureLabelSet().GetSize()
        );
        sp_dataset = retriever.GenerateDataSet();
        sp_dataset->GetIds().SetFeatureLabelSet( GetIDsFeatureLabelSet());
      }
      io::File::CloseClearFStream( input);
      return sp_dataset;
    }

    //! @brief write independent predictions using the current format, to the independent prediction file
    //! @param EXPERIMENTAL matrix with experimental data; non-const to allow descaling
    //! @param PREDICTIONS predicted results; non-const to allow descaling
    //! @param IDS ids matrix
    void CrossValidationInfo::WritePredictions
    (
      FeatureDataSet< float> &EXPERIMENTAL,
      FeatureDataSet< float> &PREDICTIONS,
      const FeatureDataSetInterface< char> &IDS,
      const std::string &FILENAME,
      const PredictionOutputFormat &FORMAT
    )
    {
      // descale predictions and results, just in case they are still scaled
      EXPERIMENTAL.DeScale();
      PREDICTIONS.DeScale();
      CrossValidationInfo::WritePredictions
      (
        EXPERIMENTAL.GetMatrix(),
        PREDICTIONS.GetMatrix(),
        IDS.GetMatrix(),
        FILENAME,
        FORMAT
      );
    }

    //! @brief write independent predictions using the current format, to the independent prediction file
    //! @param EXPERIMENTAL matrix with experimental data; non-const to allow descaling
    //! @param PREDICTIONS predicted results; non-const to allow descaling
    //! @param IDS ids matrix
    //! @param FILENAME file to write the predictions out to
    //! @param FORMAT the desired format to write the predictions in
    void CrossValidationInfo::WritePredictions
    (
      const linal::MatrixConstInterface< float> &EXPERIMENTAL,
      const linal::MatrixConstInterface< float> &PREDICTIONS,
      const linal::MatrixConstInterface< char> &IDS,
      const std::string &FILENAME,
      const PredictionOutputFormat &FORMAT
    )
    {
      // write experimental and predicted data
      io::OFStream output;
      io::File::MustOpenOFStream( output, FILENAME);

      if( FORMAT == e_Matrix)
      {
        // create matrix with twice the number of cols of the result vector
        linal::Matrix< float> matrix( EXPERIMENTAL.GetNumberRows(), size_t( 2) * PREDICTIONS.GetNumberCols(), float( 0.0));

        // off set to access indices of matrix for result values
        const size_t offset( PREDICTIONS.GetNumberCols());

        for( size_t data_id( 0), data_set_size( matrix.GetNumberRows()); data_id < data_set_size; ++data_id)
        {
          for( size_t row_element_id( 0), row_size( offset); row_element_id < row_size; ++row_element_id)
          {
            // get the actual result
            matrix[ data_id][ row_element_id] = EXPERIMENTAL[ data_id][ row_element_id];
            matrix[ data_id][ row_element_id + offset] = PREDICTIONS[ data_id][ row_element_id];
          }
        }
        output << matrix;
      }
      else // if( FORMAT == e_Delimited || FORMAT == e_ToBeDetermined)
      {
        const size_t id_size( IDS.GetNumberCols()), feature_size( PREDICTIONS.GetNumberCols());
        std::string rowid( id_size, ' ');
        const char delimiter( ',');
        const size_t feature_size_m1( feature_size - 1);
        for( size_t data_id( 0), data_set_size( PREDICTIONS.GetNumberRows()); data_id < data_set_size; ++data_id)
        {
          if( id_size)
          {
            // write id for this row
            linal::VectorConstReference< char> ids_row( IDS.GetRow( data_id));
            std::copy( ids_row.Begin(), ids_row.End(), rowid.begin());
            std::replace( rowid.begin(), rowid.end(), '\0', ' ');
            output << rowid << delimiter;
          }
          linal::VectorConstReference< float> prediction( PREDICTIONS.GetRow( data_id));
          linal::VectorConstReference< float> results( EXPERIMENTAL.GetRow( data_id));
          for( size_t row_element_id( 0); row_element_id < feature_size; ++row_element_id)
          {
            output << prediction( row_element_id) << delimiter;
          }
          for( size_t row_element_id( 0); row_element_id < feature_size_m1; ++row_element_id)
          {
            output << results( row_element_id) << delimiter;
          }
          output << results( feature_size_m1);
          output << '\n';
        }
      }

      io::File::CloseClearFStream( output);
    }

  } // namespace model
} // namespace bcl
