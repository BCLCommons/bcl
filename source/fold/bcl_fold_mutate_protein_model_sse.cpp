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
#include "fold/bcl_fold_mutate_protein_model_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "find/bcl_find_locator_interface.h"
#include "io/bcl_io_serialization.h"
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

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSE::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelSSE())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSE::MutateProteinModelSSE() :
      m_Locator(),
      m_Mutate(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSE>())
    {
    }

    //! @brief constructor from a SSE locator, a SSE Mutate and a scheme
    //! @param LOCATOR function that chooses the sse
    //! @param MUTATE function that performs the mutate on the sse
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSE::MutateProteinModelSSE
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &LOCATOR,
      const util::ShPtr< math::MutateInterface< assemble::SSE> > &MUTATE,
      const std::string &SCHEME
    ) :
      m_Locator( *LOCATOR),
      m_Mutate( *MUTATE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateProteinModelSSE *MutateProteinModelSSE::Clone() const
    {
      return new MutateProteinModelSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateProteinModelSSE::GetAlias() const
    {
      static const std::string s_alias( "MutateProteinModelSSE");
      return s_alias;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateProteinModelSSE::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Changes a protein model using the provided mutate.");
      serializer.AddInitializer
      (
        "locator",
        "locates SSEs to mutate",
        io::Serialization::GetAgent( &m_Locator)
      );
      serializer.AddInitializer
      (
        "mutate",
        "mutates the protein model",
        io::Serialization::GetAgent( &m_Mutate)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSE::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // choose a random sse and make copy
      const util::SiPtr< const assemble::SSE> located_sse( m_Locator->Locate( PROTEIN_MODEL));

      // if there was no sse found, return an empty result
      if( !located_sse.IsDefined())
      {
        BCL_MessageVrb( "could not find sse for mutating: " + util::Format()( *m_Locator));

        // return result with no model
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // create a copy of the model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // report selected sse from protein model to be moved
      BCL_MessageVrb( "selected sse to be mutated: " + located_sse->GetIdentification());

      // call the mutate and get the result
      math::MutateResult< assemble::SSE> result_sse( m_Mutate->operator ()( *located_sse));

      // if not successful
      if( !result_sse.GetArgument().IsDefined())
      {
        BCL_MessageVrb( "could not mutate sse : " + located_sse->GetIdentification());

        // return result with no model
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // replace the sse with the mutated copy of the same sse
      new_model->Replace( result_sse.GetArgument());

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
