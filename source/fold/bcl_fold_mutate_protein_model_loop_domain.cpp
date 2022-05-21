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
#include "fold/bcl_fold_mutate_protein_model_loop_domain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_loop_domain.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelLoopDomain::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelLoopDomain())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelLoopDomain::MutateProteinModelLoopDomain() :
      m_DomainCollector(),
      m_DomainMutate()
    {
    }

    //! @brief constructor taking member variable types
    //! @param LOOP_DOMAIN_COLLECTOR the method for collecting loop domains in a protein model
    //! @param DOMAIN_MUTATE the method for mutating a loop domain
    MutateProteinModelLoopDomain::MutateProteinModelLoopDomain
    (
      const util::ShPtr< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> > &LOOP_DOMAIN_COLLECTOR,
      const util::ShPtr< math::MutateInterface< LoopDomain> > &DOMAIN_MUTATE
    ) :
      m_DomainCollector( *LOOP_DOMAIN_COLLECTOR),
      m_DomainMutate( *DOMAIN_MUTATE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelLoopDomainDihedral
    MutateProteinModelLoopDomain *MutateProteinModelLoopDomain::Clone() const
    {
      return new MutateProteinModelLoopDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelLoopDomain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModelLoopDomain::GetAlias() const
    {
      static const std::string s_name( "MutateProteinModelLoopDomain");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateProteinModelLoopDomain::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Changes the loop domain of a protein model.");
      serializer.AddInitializer
      (
        "collector",
        "method of collecting a loop domain in a protein model",
        io::Serialization::GetAgent( &m_DomainCollector)
      );
      serializer.AddInitializer
      (
        "mutate",
        "mutate for changing the loop domain",
        io::Serialization::GetAgent( &m_DomainMutate)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param PROTEIN_MODEL Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< assemble::ProteinModel> MutateProteinModelLoopDomain::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static variable to hold undefined model pointer
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // collect loop domains
      util::ShPtrList< LoopDomain> loop_domains( m_DomainCollector->Collect( PROTEIN_MODEL));

      // if no loop domains
      if( loop_domains.IsEmpty())
      {
        BCL_MessageCrt( "No loop domains were found in the model");
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // pick one randomly
      util::ShPtr< LoopDomain> selected_domain
      (
        *random::GetGlobalRandom().Iterator( loop_domains.Begin(), loop_domains.End(), loop_domains.GetSize())
      );

      // mutate selected domain
      util::ShPtr< LoopDomain> mutated_loop( m_DomainMutate->operator()( *selected_domain).GetArgument());

      // make sure "mutated_loop" is defined
      BCL_Assert( mutated_loop.IsDefined(), "mutated_loop is not defined");

      // get ShPtr to protein model where its sses have been replaced with sses of mutated_loop
      util::ShPtr< assemble::ProteinModel> replaced_sse_model
      (
        // iterate through the loop domain and replace the sses in the protein model
        mutated_loop->UpdateProteinModel( PROTEIN_MODEL)
      );

      // create a mutate result from "replaced_sse_model"
      math::MutateResult< assemble::ProteinModel> mutate_result( replaced_sse_model, *this);

      // return "mutate_result"
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
