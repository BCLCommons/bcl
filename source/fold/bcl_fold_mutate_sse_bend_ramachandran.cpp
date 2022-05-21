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
#include "fold/bcl_fold_mutate_sse_bend_ramachandran.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "biol/bcl_biol_ramachandran.h"
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

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateSSEBendRamachandran::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::SSE> >::AddInstance( new MutateSSEBendRamachandran())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateSSEBendRamachandran::MutateSSEBendRamachandran() :
      m_NrChanges( math::Range< size_t>( 1, 1)),
      m_BendingDirection( biol::AASequenceFlexibility::e_Bidirectional),
      m_Scheme( GetStaticClassName< MutateSSEBendRamachandran>())
    {
    }

    //! @brief constructor from a range of number of residues to change and a scheme
    //! @param NR_RESIDUES_TO_CHANGE_RANGE ange of number of residues to change
    //! @param BENDING_DIRECTION direction the phi/psi changes should be propagated towards
    //! @param SCHEME Scheme to be used
    MutateSSEBendRamachandran::MutateSSEBendRamachandran
    (
      const math::Range< size_t> &NR_RESIDUES_TO_CHANGE_RANGE,
      const biol::AASequenceFlexibility::SequenceDirection BENDING_DIRECTION,
      const std::string &SCHEME
    ) :
      m_NrChanges( NR_RESIDUES_TO_CHANGE_RANGE),
      m_BendingDirection( BENDING_DIRECTION),
      m_Scheme( SCHEME)
    {
      // make sure the min range is one
      BCL_Assert
      (
        m_NrChanges.GetMin() >= 1,
        "The minimum number of changes has to be at least one" + util::Format()( m_NrChanges)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateSSEBendRamachandran
    MutateSSEBendRamachandran *MutateSSEBendRamachandran::Clone() const
    {
      return new MutateSSEBendRamachandran( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateSSEBendRamachandran::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateSSEBendRamachandran::GetAlias() const
    {
      static const std::string s_alias( "MutateSSEBendRamachandran");
      return s_alias;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateSSEBendRamachandran::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Applies bendings to an SSE based on the Ramachandran distribution.");
      serializer.AddInitializer
      (
        "num changes",
        "number of residues to change simultaneously",
        io::Serialization::GetAgent( &m_NrChanges)
      );
      serializer.AddInitializer
      (
        "direction",
        "sequence direction in which to propagate the change",
        io::Serialization::GetAgent( &m_BendingDirection)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes an SSE and bends amino acids and returns the new SSE
    //! @param THIS_SSE SSE of interest to bend
    //! @return math::MutateResult that has a new SSE with one or more amino acids with new phi/psi values
    math::MutateResult< assemble::SSE> MutateSSEBendRamachandran::operator()( const assemble::SSE &THIS_SSE) const
    {
      return this->operator()( THIS_SSE, util::SiPtr< const biol::Membrane>());
    }

    //! @brief operator that takes an SSE and bends amino acids and returns the new SSE
    //! @param THIS_SSE SSE of interest to bend
    //! @return math::MutateResult that has a new SSE with one or more amino acids with new phi/psi values
    math::MutateResult< assemble::SSE> MutateSSEBendRamachandran::operator()
    (
      const assemble::SSE &THIS_SSE,
      const util::SiPtr< const biol::Membrane> &MEMBRANE
    ) const
    {
      // make a copy of the selected SSE
      util::ShPtr< assemble::SSE> new_sse( THIS_SSE.Clone());

      // determine the number of changes
      const size_t total_nr_changes
      (
        random::GetGlobalRandom().SizeT
        (
          math::Range< size_t>
          (
            std::min( m_NrChanges.GetMin(), THIS_SSE.GetSize()), std::min( m_NrChanges.GetMax(), THIS_SSE.GetSize())
          )
        )
      );

      BCL_MessageVrb
      (
        "The number of amino acids to set phi/psi " + util::Format()( total_nr_changes)
      );

      // get the amino acids
      util::SiPtrVector< const biol::AABase> amino_acids( THIS_SSE.GetMembers());

      // iterate until requested nr_changes
      for( size_t nr_change( 0); nr_change != total_nr_changes; ++nr_change)
      {
        // get a random amino acid
        const util::SiPtr< const biol::AABase> this_aa( amino_acids.RemoveRandomElement());

        // skip unknown amino acids since we don't have a probability distribution for them
        if( this_aa->GetType() == biol::GetAATypes().UNK)
        {
          BCL_MessageVrb
          (
            "amino acid " + this_aa->GetIdentification() + " skipped"
          );
          continue;
        }

        // get environment type in which the current amino acid can be found
        const biol::EnvironmentType current_environment_type
        (
          MEMBRANE.IsDefined() && MEMBRANE->IsDefined() ?
          MEMBRANE->DetermineEnvironmentType( this_aa->GetFirstSidechainAtom().GetCoordinates()) :
          biol::GetEnvironmentTypes().e_Solution
        );

        // skip undefined environment
        if( !current_environment_type.IsDefined())
        {
          continue;
        }

        // get a random phi psi using ramachandran to set for this amino acid
        storage::VectorND< 2, double> new_phi_psi
        (
          biol::Ramachandran::GetDefaultInstance().GetRandomPhiPsi( this_aa->GetType(), THIS_SSE.GetType(), current_environment_type)
        );

        // message
        BCL_MessageVrb
        (
          "amino acid " + this_aa->GetIdentification() +
          " new phi: " + util::Format()( new_phi_psi.First()) + " psi: " + util::Format()( new_phi_psi.Second())
        );

        // use AASequence flexibility to set phi/psi
        biol::AASequenceFlexibility::SetPhiPsi( *new_sse, this_aa->GetSeqID(), new_phi_psi, m_BendingDirection);
      }

      // set the geometries of the sse
      new_sse->SetGeometry();

      // return the result
      return math::MutateResult< assemble::SSE>( new_sse, *this);
    }

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateSSEBendRamachandran::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      auto sses( PROTEIN_MODEL.GetSSEs());

      // get random iterator on the list
      util::SiPtrVector< const assemble::SSE>::const_iterator random_iterator
      (
        random::GetGlobalRandom().Iterator( sses.Begin(), sses.End(), sses.GetSize())
      );

      // if there was no sse found, return an empty result
      if( random_iterator == sses.End())
      {
        // return result with no model
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // randomly pick an sse
      util::SiPtr< const assemble::SSE> random_sse( *random_iterator);

      // get membrane for current protein model
      const util::ShPtr< biol::Membrane> sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // create a copy of the model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // report selected sse from protein model to be moved
      BCL_MessageVrb( "selected sse to be mutated: " + random_sse->GetIdentification());

      // call the mutate and get the result
      math::MutateResult< assemble::SSE> result_sse( this->operator()( *random_sse, sp_membrane));

      // if not successful
      if( !result_sse.GetArgument().IsDefined())
      {
        BCL_MessageVrb( "could not mutate sse : " + random_sse->GetIdentification());

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
