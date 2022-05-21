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

#ifndef BCL_MODEL_KOHONEN_NETWORK_APPLICABILITY_DOMAIN_H_
#define BCL_MODEL_KOHONEN_NETWORK_APPLICABILITY_DOMAIN_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_kohonen_network_average.h"
#include "bcl_model_retrieve_interface.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KohonenNetworkApplicabilityDomain
    //! @brief takes the average output of the training data won by each node.
    //!
    //! @see @link example_model_kohonen_network_applicability_domain.cpp @endlink
    //! @author mendenjl
    //! @date Sep 19, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KohonenNetworkApplicabilityDomain :
      public KohonenNetworkAverage
    {

    //////////
    // data //
    //////////

      //! Splines associated with each kohonen node and each result (inner)
      storage::Vector< storage::Vector< math::CubicSplineDamped> > m_Splines;

      //! Bool: whether to only use a single spline for all network distances
      bool m_NodesSharedDistanceMetric;

      //! Optional implementation of a model retrieval used to compute the input values.
      //! If given, when this model is called to make a prediction, it will use the underlying model to compute the
      //! prediction, which will then be transformed in some manner (specified by the approximator that created this
      //! model). A common application is transforming predicted values into PPV values
      util::Implementation< RetrieveInterface> m_ModelRetriever;

      //! Models for all predictions
      mutable storage::Vector< util::OwnPtr< Interface> > m_Models;

    public:

      //! add self framework enums
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      KohonenNetworkApplicabilityDomain() :
        m_NodesSharedDistanceMetric( true)
      {
      }

      //! @brief constructor from parameters
      //! @param NETWORK the network to copy
      //! @param SHARE_DISTANCE_METRIC Whether to share the distance metric across all nodes
      //! @param MODEL_RETRIEVER retriever for the models
      KohonenNetworkApplicabilityDomain
      (
        const KohonenNetworkAverage &NETWORK,
        const bool &SHARE_DISTANCE_METRIC = true,
        const util::Implementation< RetrieveInterface> &MODEL_RETRIEVER = util::Implementation< RetrieveInterface>()
      );

      //! @brief Clone function
      //! @return pointer to new KohonenNetworkApplicabilityDomain
      KohonenNetworkApplicabilityDomain *Clone() const
      {
        return new KohonenNetworkApplicabilityDomain( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief gets the splines
      //! @return return constant reference to map distance function
      const storage::Vector< storage::Vector< math::CubicSplineDamped> > &GetSplines() const
      {
        return m_Splines;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief setup all internal splines
      //! @param NORMALIZED_DATA all the training data, already normalized
      //! @param PREVIOUS_WINNERS previous rounds node winners, will be updated
      //! @param SHARE_DISTANCE_METRIC Whether to share the distance metric across all nodes
      void SetupSplines
      (
        const FeatureDataSetInterface< float> &NORMALIZED_DATA,
        const storage::Vector< size_t> &PREVIOUS_WINNERS,
        const bool &SHARE_DISTANCE_METRIC = true
      );

      //! @brief set the model that this object should call to obtain the un-transformed prediction
      //! @param MODELS the models to call to compute the underlying prediction
      void SetupModels( const util::Implementation< RetrieveInterface> &MODEL_RETRIEVER);

      //! @brief setup all internal splines
      //! @param SPLINES splines, one for each node, or just one (for whole network)
      void SetupSplines( const storage::Vector< storage::Vector< math::CubicSplineDamped> > &SPLINES);

    ///////////////
    // operators //
    ///////////////

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling
      (
        const FeatureDataSetInterface< float> &FEATURE
      ) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief Trains the cubic spline based on distances to NODE
      //! @param MAPPED_POINTS all points mapped to NODE
      static math::CubicSplineDamped MapData
      (
        const storage::Vector< linal::VectorConstReference< float> > &MAPPED_POINTS,
        const linal::VectorConstInterface< float> &NODE
      );

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class KohonenNetworkApplicabilityDomain

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_KOHONEN_NETWORK_APPLICABILITY_DOMAIN_H_
