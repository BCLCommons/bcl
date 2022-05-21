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

#ifndef BCL_MODEL_APPROXIMATOR_APPLICABILITY_DOMAIN_KOHONEN_H_
#define BCL_MODEL_APPROXIMATOR_APPLICABILITY_DOMAIN_KOHONEN_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_approximator_kohonen_network.h"
#include "bcl_model_kohonen_network_applicability_domain.h"
#include "bcl_model_kohonen_network_average.h"
#include "bcl_model_objective_function_wrapper.h"
#include "bcl_model_retrieve_data_set_base.h"
#include "bcl_model_retrieve_interface.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorApplicabilityDomainKohonen
    //! @brief A kohonen-map based implementation to detect whether a point is within the applicability domain of a model
    //! Inspired by
    //! @see @link http://www.sciencedirect.com/science/article/pii/S000326701201642X @endlink
    //!
    //! @see @link example_model_approximator_applicability_domain_kohonen.cpp @endlink
    //! @author mendenjl
    //! @date Sep 17, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorApplicabilityDomainKohonen :
      public ApproximatorKohonenNetwork
    {

    //////////
    // data //
    //////////

    public:

      //! static instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_ContingencyInstance;

      //! object data label for the model storage from which to load the models
      //! May be left undefined, in which case this can be trained as a traditional model
      util::Implementation< RetrieveInterface> m_Retriever;

      //! Average movement in grid space for molecule assignments last iteration
      double m_AverageDisplacement;

      //! Average distance between gridpoints
      double m_AverageDistanceBetweenGridpoints;

      //! Whether to share the distance metric between all nodes
      bool m_ShareDistanceMetric;

      //! Bool. True if using a distance metric, false if a measure should be specified
      bool m_UseDistanceMetric;

      //! Monotonize. True if using a contingency matrix metric that is normally expected to be monotonic
      bool m_Monotonize;

      //! Measure to use for spline computation. If specified, this replaces the distance metric
      util::Implementation< util::FunctionInterfaceSerializable< math::ContingencyMatrix, double> > m_Measure;

      //! Result dataset, used for contingency matrix-based metric
      util::ShPtr< descriptor::Dataset> m_ResultDataset;

      //! Final model, produced only when the tracker indicates training is complete
      util::ShPtr< KohonenNetworkApplicabilityDomain> m_FinalModel;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorApplicabilityDomainKohonen( const bool &USE_DISTANCE_METRIC = true);

      //! @brief constructor from all necessary parameters
      //! @param MAP_DIMENSIONS dimensions of the map
      //! @param INITIAL_LENGTH how many iterations to train for (if applicable).
      //! @param INITAL_RADIUS the initial neighborhood radius
      //! @param OBJECTIVE_FUNCTION the objective function from the approximator framework
      //! @param UPDATE_EVERY_NTH_FEATURE update the nodes after seeing this many features
      //! @param NEIGHBOR_KERNEL the neighbor kernel type
      ApproximatorApplicabilityDomainKohonen
      (
        const linal::Vector< double> &MAP_DIMENSIONS,
        const size_t &INITIAL_LENGTH,
        const float &INITAL_RADIUS,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
        const size_t &UPDATE_EVERY_NTH_FEATURE,
        const NeighborKernel &NEIGHBOR_KERNEL,
        const util::Implementation< RetrieveInterface> &RETREIVER =
            util::Implementation< RetrieveInterface>()
      );

      //! @brief clone function
      //! @return pointer to new ApproximatorApplicabilityDomainKohonen
      ApproximatorApplicabilityDomainKohonen *Clone() const
      {
        return new ApproximatorApplicabilityDomainKohonen( *this);
      }

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief set training data set for a specific iterate in approximator framework
      //! @param DATA training data set
      void SetTrainingData
      (
        util::ShPtr< descriptor::Dataset> &DATA
      );

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief produce the final model
      //! Updates m_FinalModel
      void MakeFinalModel();

    }; // class ApproximatorApplicabilityDomainKohonen

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_APPLICABILITY_DOMAIN_KOHONEN_H_

