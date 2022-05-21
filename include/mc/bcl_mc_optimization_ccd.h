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

#ifndef BCL_MC_OPTIMIZATION_CCD_H_
#define BCL_MC_OPTIMIZATION_CCD_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_optimization_mcm.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OptimizationCCD
    //! @brief Constructs loop regions using cyclic coordinate descent for a given protein model
    //!
    //! @see @link example_mc_optimization_ccd.cpp @endlink
    //! @author fischea
    //! @date Oct 27, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API OptimizationCCD :
      public OptimizationMCM
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

      //! @brief default constructor
      OptimizationCCD();

      //! @brief construct from members
      //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
      //! @param MUTATES mutates to sample models
      //! @param CRITERION termination criterion
      //! @param METROPOLIS metropolis criterion to decide whether mutates are accepted
      OptimizationCCD
      (
        const score::ProteinModelScoreSum &SCORE_FUNCTION,
        const math::MutateInterface< assemble::ProteinModel> &MUTATES,
        const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
        const Metropolis< double> &METROPOLIS
      );

      //! @brief clone function
      //! @return pointer to a new OptimizationCCD
      OptimizationCCD *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief performs pre-processing of the argument before the optimization is performed
      //! @detail adds initial coordinates for missing loop regions and adds a loop domain locator to the protein model data
      //! @param ARGUMENT data to be pre-processed
      //! @return pre-processed data
      void PreProcess( assemble::ProteinModel &ARGUMENT) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class OptimizationCCD

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_OPTIMIZATION_CCD_H_
