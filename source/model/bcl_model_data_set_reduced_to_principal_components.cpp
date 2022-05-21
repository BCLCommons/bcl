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

// include forward header of this class
#include "model/bcl_model_data_set_reduced_to_principal_components.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_principal_component_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DataSetReducedToPrincipalComponents::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new DataSetReducedToPrincipalComponents())
    );

    //! @brief constructor from parameters
    //! @param FRACTION fraction of total sum to keep
    //! @param MAX_COMPONENTS_TO_KEEP maximum components to keep
    DataSetReducedToPrincipalComponents::DataSetReducedToPrincipalComponents
    (
      const double &FRACTION,
      const size_t &MAX_COMPONENTS_TO_KEEP
    ) :
      m_Fraction( FRACTION),
      m_MaxToKeep( MAX_COMPONENTS_TO_KEEP),
      m_NumberKept( 0)
    {
      BCL_Assert( m_Fraction >= 0.0 && m_Fraction <= 1.0, "Fraction must be between 0.0 and 1.0");
    }

    //! @brief Clone function
    //! @return pointer to new DataSetReducedToPrincipalComponents
    DataSetReducedToPrincipalComponents *DataSetReducedToPrincipalComponents::Clone() const
    {
      return new DataSetReducedToPrincipalComponents( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DataSetReducedToPrincipalComponents::GetClassIdentifier() const
    {
      return GetStaticClassName< DataSetReducedToPrincipalComponents>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &DataSetReducedToPrincipalComponents::GetAlias() const
    {
      static const std::string s_Name( "PCA");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    void DataSetReducedToPrincipalComponents::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Retriever->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @param CODE the new code
    void DataSetReducedToPrincipalComponents::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Retriever->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void DataSetReducedToPrincipalComponents::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Retriever->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToPrincipalComponents::GetFeatureLabelsWithSizes() const
    {
      FeatureLabelSet labels;
      labels.PushBack( util::ObjectDataLabel( "PCA"), m_NumberKept);
      return labels;
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToPrincipalComponents::GetResultCodeWithSizes() const
    {
      return m_Retriever->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToPrincipalComponents::GetIdCodeWithSizes() const
    {
      return m_Retriever->GetIdCodeWithSizes();
    }

    //! @brief set retriever
    //! @param RETRIEVER the retriever used to get the data set
    void DataSetReducedToPrincipalComponents::SetRetriever( const util::Implementation< RetrieveDataSetBase> &RETRIEVER)
    {
      m_Retriever = RETRIEVER;
    }

    //! @brief sets the filename of the rescaler, eigenvectors, and eigenvalues
    //! @param FILENAME the filename
    void DataSetReducedToPrincipalComponents::SetFilename( const std::string &FILENAME)
    {
      m_Filename = FILENAME;
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > DataSetReducedToPrincipalComponents::GetNumberPartitionsAndIds() const
    {
      return m_Retriever->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief generate dataset, reduced to the desired # of principal components
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    DataSetReducedToPrincipalComponents::GenerateDataSet()
    {
      // load in the info from the eigenvector matrix file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);
      RescaleFeatureDataSet rescaler;
      io::Serialize::Read( rescaler, input);
      storage::Pair< linal::Matrix< float>, linal::Vector< float> > eigenvectors_values;
      io::Serialize::Read( eigenvectors_values, input);
      io::File::CloseClearFStream( input);

      // get the base data set
      util::ShPtr< descriptor::Dataset> feature_result( m_Retriever->GenerateDataSet());

      // get a reference to the features matrix
      linal::MatrixReference< float> features( feature_result->GetFeaturesReference());
      rescaler.RescaleMatrix( features);

      // reduce input matrix
      linal::Matrix< float> new_matrix;
      linal::PrincipalComponentAnalysis< float>::ReduceInputMatrix
      (
        features,
        new_matrix,
        eigenvectors_values.First(),
        eigenvectors_values.Second(),
        m_Fraction,
        m_MaxToKeep
      );

      // copy the feature matrix, now that it is reduced
      util::ShPtr< FeatureDataSet< float> > features_sp( feature_result->GetFeaturesPtr());
      features_sp->GetRawMatrix() = new_matrix;
      features_sp->SetFeatureLabelSet( GetFeatureLabelsWithSizes());

      // return a new data set initialized with the feature results from the pca
      return feature_result;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DataSetReducedToPrincipalComponents::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Reduces dataset to principal component");

      member_data.AddInitializer
      (
        "dataset",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Retriever)
      );

      member_data.AddInitializer
      (
        "fraction",
        "fraction of total eigenvalue sum to keep",
        io::Serialization::GetAgentWithRange( &m_Fraction, 0.0, 1.0)
      );

      member_data.AddInitializer
      (
        "components",
        "max number of components to allow in the output (independent of the fraction of total component)",
        io::Serialization::GetAgent( &m_MaxToKeep),
        "1000000"
      );

      member_data.AddInitializer
      (
        "filename",
        "file where rescaling, eigenvectors and values are stored",
        io::Serialization::GetAgentInputFilename( &m_Filename)
      );

      return member_data;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t DataSetReducedToPrincipalComponents::GetNominalSize() const
    {
      return m_Retriever->GetNominalSize();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool DataSetReducedToPrincipalComponents::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // load in the info from the eigenvector matrix file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);
      RescaleFeatureDataSet rescaler;
      io::Serialize::Read( rescaler, input);
      storage::Pair< linal::Matrix< float>, linal::Vector< float> > eigenvectors_values;
      io::Serialize::Read( eigenvectors_values, input);
      io::File::CloseClearFStream( input);
      linal::Matrix< float> features( 1, m_Retriever->GetFeatureLabelsWithSizes().GetSize(), 0.0);
      linal::Matrix< float> new_matrix;
      linal::PrincipalComponentAnalysis< float>::ReduceInputMatrix
      (
        features,
        new_matrix,
        eigenvectors_values.First(),
        eigenvectors_values.Second(),
        m_Fraction,
        m_MaxToKeep
      );
      m_NumberKept = new_matrix.GetNumberCols();
      return true;
    }

  } // namespace model
} // namespace bcl
