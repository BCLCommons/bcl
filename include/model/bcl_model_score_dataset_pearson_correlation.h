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

#ifndef BCL_MODEL_SCORE_DATASET_PEARSON_CORRELATION_H
#define BCL_MODEL_SCORE_DATASET_PEARSON_CORRELATION_H

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_score_dataset_interface.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreDatasetPearsonCorrelation
    //! @brief calculate pearson correlation values for every column of a dataset
    //!
    //! f(i)= Ave((X_avg - X_avg) * (Y_avg - Y_avg))
    //!
    //!
    //! @author mendenjl
    //! @see @link example_model_score_dataset_pearson_correlation.cpp @endlink
    //! @date Apr 02, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDatasetPearsonCorrelation :
      public ScoreDatasetInterface
    {
    private:

      bool m_AutoCorrelation; //!< True to compute max correlation of each descriptor against any other descriptor

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_AutoInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from m_AutoCorrelation parameter
      ScoreDatasetPearsonCorrelation( const bool &AUTO_CORRELATION = false) :
        m_AutoCorrelation( AUTO_CORRELATION)
      {
      }

      //! @brief Clone function
      //! @return pointer to new ScoreDatasetPearsonCorrelation
      ScoreDatasetPearsonCorrelation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief score a given dataset
      //! @param DATASET dataset of interest
      //! @return scores of the dataset
      linal::Vector< float> Score( const descriptor::Dataset &DATASET) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ScoreDatasetPearsonCorrelation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DATASET_PEARSON_CORRELATION_H

