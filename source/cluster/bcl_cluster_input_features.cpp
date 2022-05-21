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
#include "cluster/bcl_cluster_input_features.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_distances_euclidean.h"
#include "model/bcl_model_retrieve_data_set_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> InputFeatures::s_Instance
    (
      GetObjectInstances().AddInstance( new InputFeatures())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    InputFeatures::InputFeatures() :
      m_DatasetRetriever(),
      m_RetrievedDataset(),
      m_Scale( false)
    {
    }

    //! @brief constructor taking members
    //! @param RETRIEVER method for retrieving feature vectors
    //! @param SCALE bool indicating if feature vectors should be scaled or not
    InputFeatures::InputFeatures( const util::Implementation< model::RetrieveDataSetBase> &RETRIEVER, bool SCALE) :
      m_DatasetRetriever( RETRIEVER),
      m_RetrievedDataset(),
      m_Scale( SCALE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new InputFeatures
    InputFeatures *InputFeatures::Clone() const
    {
      return new InputFeatures( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &InputFeatures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief HandleInput gets data from input; puts it into the data construct with pointers to the objects
    //         in m_Objects
    //! @param IFSTREAM is the stream from which the input will be read
    //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
    util::ShPtr
    <
      math::FunctionInterface
      <
        storage::VectorND< 2, util::SiPtr< const linal::Vector< float> > >, float
      >
    >
    InputFeatures::HandleInput( io::IFStream &IFSTREAM)
    {
      // retrieve the data set of features and results and set the member variable
      m_RetrievedDataset = m_DatasetRetriever->GenerateDataSet();

      // true if the data set needs to be scaled
      if( m_Scale)
      {
        // scale
        m_RetrievedDataset->GetFeaturesReference().CopyValues( Rescale( m_RetrievedDataset->GetFeaturesReference()));
      }

      linal::MatrixConstReference< float> features( m_RetrievedDataset->GetFeaturesReference());

      // make list that will hold the objects to be clustered
      util::ShPtr< storage::List< linal::Vector< float> > > cluster_objects( new storage::List< linal::Vector< float> >());

      // get the size of each feature
      const size_t feature_size( features.GetNumberCols());

      for
      (
        size_t feature_number( 0), number_features( features.GetNumberRows());
        feature_number < number_features;
        ++feature_number
      )
      {
        cluster_objects->PushBack( linal::Vector< float>( feature_size, features[ feature_number]));
      }

      SetInputObjects( cluster_objects);

      return util::ShPtr
        <
          math::FunctionInterface
          <
            storage::VectorND< 2, util::SiPtr< const linal::Vector< float> > >, float
          >
        >( new DistancesEuclidean< float>());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &InputFeatures::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DatasetRetriever, ISTREAM);
      io::Serialize::Read( m_RetrievedDataset, ISTREAM);
      io::Serialize::Read( m_Scale, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &InputFeatures::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DatasetRetriever, OSTREAM, INDENT);
      io::Serialize::Write( m_RetrievedDataset, OSTREAM, INDENT);
      io::Serialize::Write( m_Scale, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief rescales a matrix
    //! @param DATA the matrix to rescale
    //! @return matrix that has been rescaled
    linal::Matrix< float> InputFeatures::Rescale( const linal::MatrixConstInterface< float> &DATA)
    {
      // store rescale information here
      storage::Vector< math::Range< float> > rescale_ranges;
      const size_t feature_size( DATA.GetNumberCols());
      rescale_ranges.AllocateMemory( feature_size);

      // scale all features and results between 0 and 1, if possible
      math::Range< float> autoscale_range( 0.0, 1.0);

      // creating vector of ranges
      float tmp( 0);
      for( size_t col( 0); col < feature_size; ++col)
      {
        // set initial min and max
        float min( std::numeric_limits< float>::infinity());
        float max( -std::numeric_limits< float>::infinity());

        // iterate through all values in that feature col
        for( size_t row( 0), num_rows( DATA.GetNumberRows()); row < num_rows; ++row)
        {
          // set as min or max if qualifies
          tmp = DATA( row, col);
          if( tmp < min)
          {
            min = tmp;
          }
          else if( tmp > max)
          {
            max = tmp;
          }
        }
        // handle the case where min == max
        if( min == max)
        {
          min -= 0.5;
          max += 0.5;
        }
        // add to range vector
        rescale_ranges.PushBack( math::Range< float>( min, max));
      }

      // copy data so it can be rescaled
      linal::Matrix< float> rescaled_matrix( DATA);

      for( size_t row( 0), num_rows( DATA.GetNumberRows()); row < num_rows; ++row)
      {
        float *feature_row( rescaled_matrix[ row]);
        // rescale each feature
        for( size_t col( 0); col < feature_size; ++col)
        {
          feature_row[ col] = autoscale_range( DATA( row, col), rescale_ranges( col));
        }
      }

      return rescaled_matrix;
    }

  } // namespace cluster

} // namespace bcl
