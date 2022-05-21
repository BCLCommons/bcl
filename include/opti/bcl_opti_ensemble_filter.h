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

#ifndef BCL_OPTI_ENSEMBLE_FILTER_H_
#define BCL_OPTI_ENSEMBLE_FILTER_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_optimization_interface.h"
#include "assemble/bcl_assemble_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "score/bcl_score_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnsembleFilter
    //! @brief Filters an ensemble of protein models using a provided scoring function.
    //!
    //! @see @link example_opti_ensemble_filter.cpp @endlink
    //! @author fischea
    //! @date Apr 19, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API EnsembleFilter :
      public OptimizationInterface< assemble::Ensemble< assemble::ProteinModel> >
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! scoring function used to evaluate the sampled structures
      util::Implementation< score::ProteinModel> m_ScoreFunction;

      //! percentage of models to keep
      double m_KeepPercentage;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      EnsembleFilter();

      //! @brief construct from members
      //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
      //! @param KEEP percentage of model to keep
      EnsembleFilter( const score::ProteinModel &SCORE_FUNCTION, double KEEP);

      //! @brief clone function
      //! @return pointer to a new EnsembleFilter
      EnsembleFilter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief filters the provided ensemble
      //! @param ENSEMBLE ensemble to be filtered
      void Optimize( assemble::Ensemble< assemble::ProteinModel> &ENSEMBLE) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class OptimizationMCM

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_ENSEMBLE_FILTER_H_
