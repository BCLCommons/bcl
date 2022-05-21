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

#ifndef BCL_MODEL_SCORE_DATASET_F_SCORE_H
#define BCL_MODEL_SCORE_DATASET_F_SCORE_H

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
    //! @class ScoreDatasetFScore
    //! @brief calculate f-score values for every column in combined active and inactive dataset
    //!
    //!                                     (X_avg_pos - X_avg)^2 + (X_avg_neg - X_avg)^2
    //! f(i)= -------------------------------------------------------------------------------------------------------
    //!       1/(n_pos - 1) SUM k=1 to n_pos (X_k_pos - X_avg)^2 + 1/(n_neg - 1) SUM k=1 to n_neg (X_k_neg - X_avg)^2
    //!
    //!
    //! @author mendenjl
    //! @see @link example_model_score_dataset_f_score.cpp @endlink
    //! @date Apr 06, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDatasetFScore :
      public ScoreDatasetInterface
    {
    private:

    //////////
    // data //
    //////////

      float m_Cutoff; //!< Cutoff between actives/inactives

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ScoreDatasetFScore
      ScoreDatasetFScore *Clone() const;

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

    }; // class ScoreDatasetFScore

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DATASET_F_SCORE_H

