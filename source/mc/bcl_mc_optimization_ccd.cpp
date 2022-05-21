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
#include "mc/bcl_mc_optimization_ccd.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_handler_locator_loop_domain.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_protocol_loop_close.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> OptimizationCCD::s_Instance
    (
      util::Enumerated< opti::OptimizationInterface< assemble::ProteinModel> >::AddInstance( new OptimizationCCD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    OptimizationCCD::OptimizationCCD()
    {
    }

    //! @brief construct from members
    //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
    //! @param MUTATES mutates to sample models
    //! @param CRITERION termination criterion
    //! @param METROPOLIS metropolis criterion to decide whether mutates are accepted
    OptimizationCCD::OptimizationCCD
    (
      const score::ProteinModelScoreSum &SCORE_FUNCTION,
      const math::MutateInterface< assemble::ProteinModel> &MUTATES,
      const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
      const Metropolis< double> &METROPOLIS
    ) :
      OptimizationMCM( SCORE_FUNCTION, MUTATES, CRITERION, METROPOLIS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new OptimizationCCD
    OptimizationCCD *OptimizationCCD::Clone() const
    {
      return new OptimizationCCD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &OptimizationCCD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &OptimizationCCD::GetAlias() const
    {
      static const std::string s_alias( "CCDOptimizer");
      return s_alias;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief performs pre-processing of the argument before the optimization is performed
    //! @detail adds initial coordinates for missing loop regions and adds a loop domain locator to the protein model data
    //! @param ARGUMENT data to be pre-processed
    //! @return pre-processed data
    void OptimizationCCD::PreProcess( assemble::ProteinModel &ARGUMENT) const
    {
      const fold::ProtocolLoopClose &protocol( fold::ProtocolLoopClose::GetInstance());

      // complete the model with coordinates in missing places
      protocol.AddLoopCoordinates( ARGUMENT);

      // split the coils, where they are not peptide bonded
      protocol.SplitCoilsAtNonPetideBond( ARGUMENT);

      // add backbone hydrogens
      ARGUMENT = *protocol.AddNitrogenHydrogens( ARGUMENT);

      // add a loop locator to the protein model data
      const util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > sp_locator
      (
        fold::HandlerLocatorLoopDomain().CreateBidirectionalLocatorsForInteriorCoil( ARGUMENT)
      );
      util::ShPtr< assemble::ProteinModelData> sp_model_data( ARGUMENT.GetProteinModelData());
      if( !sp_model_data->GetData( assemble::ProteinModelData::e_LoopDomainLocators).IsDefined())
      {
        sp_model_data->Insert( assemble::ProteinModelData::e_LoopDomainLocators, sp_locator);
      }
      else
      {
        sp_model_data->Replace( assemble::ProteinModelData::e_LoopDomainLocators, sp_locator);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace mc
} // namespace bcl
