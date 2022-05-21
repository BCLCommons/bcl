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
#include "fold/bcl_fold_mutate_protein_model_loop_domain_ccd.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_collector_loop_domain.h"
#include "fold/bcl_fold_collector_loop_domain_random.h"
#include "fold/bcl_fold_loop_domain.h"
#include "fold/bcl_fold_mutate_loop_domain_dihedral.h"
#include "fold/bcl_fold_phi_psi_generator_ccd.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelLoopDomainCCD::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelLoopDomainCCD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelLoopDomainCCD::MutateProteinModelLoopDomainCCD() :
      m_DomainCollector(),
      m_DomainLocator(),
      m_RandomNumberGenerator( random::GetGlobalRandom()),
      m_RandomFraction( 1.0, 1.0)
    {
    }

    //! @brief constructor taking member parameters
    //! @param DOMAIN_COLLECTOR the collector to collect possible loop domains for mutating
    //! @param RANDOM_NUMBER_GENERATOR the random number generator that should be used
    //! @param RANDOM_FRACTION_RANGE range a random fraction is drawn from and multiplied with suggested rotation
    MutateProteinModelLoopDomainCCD::MutateProteinModelLoopDomainCCD
    (
      const util::ShPtr< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> > &DOMAIN_COLLECTOR,
      const random::DistributionInterface &RANDOM_NUMBER_GENERATOR,
      const math::Range< double> &RANDOM_FRACTION_RANGE
    ) :
      m_DomainCollector( DOMAIN_COLLECTOR),
      m_DomainLocator(),
      m_RandomNumberGenerator( RANDOM_NUMBER_GENERATOR),
      m_RandomFraction( RANDOM_FRACTION_RANGE)
    {
    }

    //! @brief constructor taking member parameters
    //! @param DOMAIN_LOCATOR the locator to locate loop domain for mutating
    //! @param RANDOM_NUMBER_GENERATOR the random number generator that should be used
    //! @param RANDOM_FRACTION bool indicating if random fraction of the optimal rotation should be taken
    MutateProteinModelLoopDomainCCD::MutateProteinModelLoopDomainCCD
    (
      const util::ShPtr< find::LocatorInterface< util::ShPtr< LoopDomain>, assemble::DomainInterface> > &DOMAIN_LOCATOR,
      const random::DistributionInterface &RANDOM_NUMBER_GENERATOR,
      const math::Range< double> &RANDOM_FRACTION_RANGE
    ) :
      m_DomainCollector(),
      m_DomainLocator( DOMAIN_LOCATOR),
      m_RandomNumberGenerator( RANDOM_NUMBER_GENERATOR),
      m_RandomFraction( RANDOM_FRACTION_RANGE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelLoopDomainCCD
    MutateProteinModelLoopDomainCCD *MutateProteinModelLoopDomainCCD::Clone() const
    {
      return new MutateProteinModelLoopDomainCCD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelLoopDomainCCD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModelLoopDomainCCD::GetAlias() const
    {
      static const std::string s_name( "MutateProteinModelLoopDomainCCD");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateProteinModelLoopDomainCCD::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Changes dihedral angles of residues in loop regions.");
      serializer.AddInitializer
      (
        "random fraction range",
        "range of the random multiplier to apply to the suggested rotation angle",
        io::Serialization::GetAgent( &m_RandomFraction),
        "(1.0, 1.0)"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool MutateProteinModelLoopDomainCCD::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM
    )
    {
      // create the domain locator
      util::ShPtr< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> > sp_unclosed_loop_collector
      (
        new CollectorLoopDomain( true, 0.08)
      );
      m_DomainCollector = sp_unclosed_loop_collector;

      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param PROTEIN_MODEL Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< assemble::ProteinModel> MutateProteinModelLoopDomainCCD::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static variable to hold undefined model pointer
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // initialize empty loop domain
      util::ShPtr< LoopDomain> loop_domain;

      // if locator is defined
      if( m_DomainLocator.IsDefined())
      {
        loop_domain = m_DomainLocator->Locate( PROTEIN_MODEL);
      }
      // otherwise use collector and pick one randomly
      else
      {
        // collect loop domains
        util::ShPtrList< LoopDomain> loop_domains( m_DomainCollector->Collect( PROTEIN_MODEL));

        // if no loop domains
        if( loop_domains.IsEmpty())
        {
          BCL_MessageCrt( "No loop domains were found in the model");
          return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
        }

        // pick one randomly
        loop_domain = *random::GetGlobalRandom().Iterator( loop_domains.Begin(), loop_domains.End(), loop_domains.GetSize());
      }

      // create a method for collecting residues to mutate
      const util::ShPtr< CollectorLoopDomainRandom> mutation_residue_collector
      (
        // initialize with CollectorLoopDomainRandom which will randomly select one of the residues in "loop_domain"
        new CollectorLoopDomainRandom( 1, m_RandomNumberGenerator)
      );

      // create list of TargetAndMovingPointPair objects
      const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> target_and_moving_points
      (
        loop_domain->TargetAndMovingPointsForCCD( PROTEIN_MODEL)
      );

      // need points
      if( target_and_moving_points.IsEmpty())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // create a ShPtr to a PhiPsiGeneratorCCD which will be used to generate phi and psi angles
      const util::ShPtr< PhiPsiGeneratorCCD> phi_psi_generator
      (
        new PhiPsiGeneratorCCD
        (
          target_and_moving_points,
          m_RandomNumberGenerator,
          loop_domain->GetSequenceDirection(),
          m_RandomFraction
        )
      );

      // create a MutateLoopDomainDihedral which will be used to mutate "loop_domain"
      // and initialize it with "mutation_residue_collector" and "phi_psi_generator"
      const MutateLoopDomainDihedral domain_mutator( mutation_residue_collector, phi_psi_generator);

      // create a mutated loop domain "mutated_loop" as provided by "domain_mutator"
      util::ShPtr< LoopDomain> mutated_loop( domain_mutator( *loop_domain).GetArgument());

      // make sure that "mutated_loop" is defined
      if( !mutated_loop.IsDefined())
      {
        BCL_MessageCrt( "mutated_loop_domain is not defined");
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

//      // make sure the atoms of the anchor atom is set correctly
//      mutated_loop->SetAnchorResidueOxygen( PROTEIN_MODEL);

      // create a ShPtr to the protein model where its sses have been replaced with those of "mutated_loop"
      util::ShPtr< assemble::ProteinModel> replaced_sse_model
      (
        // updated the protein mdoel according to the mutated loop
        mutated_loop->UpdateProteinModel( PROTEIN_MODEL)
      );

      // create mutate result "mutate_result" out of "replaced_sse_model"
      return math::MutateResult< assemble::ProteinModel>( replaced_sse_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
