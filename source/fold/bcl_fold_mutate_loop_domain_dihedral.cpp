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
#include "fold/bcl_fold_mutate_loop_domain_dihedral.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_loop_domain.h"
#include "fold/bcl_fold_mutation_residue.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateLoopDomainDihedral::s_Instance
    (
      util::Enumerated< math::MutateInterface< LoopDomain> >::AddInstance( new MutateLoopDomainDihedral())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateLoopDomainDihedral::MutateLoopDomainDihedral() :
      m_ResiduesToChangeCollector(),
      m_DihedralGenerator()
    {
    }

    //! @brief constructor taking parameters
    //! @param RESIDUES_TO_CHANGE_COLLECTOR the method for collecting residues to change their dihedral angles
    //! @param DIHEDRAL_GENERATOR the method that determines what the dihedral angles should be set to
    MutateLoopDomainDihedral::MutateLoopDomainDihedral
    (
      const util::ShPtr
      <
        find::CollectorInterface< storage::List< MutationResidue>, LoopDomain>
      > &RESIDUES_TO_CHANGE_COLLECTOR,
      const util::ShPtr
      <
        math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
      > &DIHEDRAL_GENERATOR
    ) :
      m_ResiduesToChangeCollector( *RESIDUES_TO_CHANGE_COLLECTOR),
      m_DihedralGenerator( *DIHEDRAL_GENERATOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateLoopDomainDihedral
    MutateLoopDomainDihedral *MutateLoopDomainDihedral::Clone() const
    {
      return new MutateLoopDomainDihedral( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateLoopDomainDihedral::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MutateLoopDomainDihedral::GetAlias() const
    {
      static const std::string s_name( "MutateLoopDomainDihedral");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateLoopDomainDihedral::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Changes the dihedral angles of a loop domain in a protein model.");
      serializer.AddInitializer
      (
        "collector",
        "method of collecting a loop domain in a protein model",
        io::Serialization::GetAgent( &m_ResiduesToChangeCollector)
      );
      serializer.AddInitializer
      (
        "generator",
        "method to determine the new dihedral angles",
        io::Serialization::GetAgent( &m_DihedralGenerator)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
    //! @param LOOP_DOMAIN Argument of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< LoopDomain> MutateLoopDomainDihedral::operator()( const LoopDomain &LOOP_DOMAIN) const
    {
      // make a clone copy into a ShPtr of "LOOP_DOMAIN" which will be mutated
      util::ShPtr< LoopDomain> new_loop_domain( LOOP_DOMAIN.Clone());

      // collet all of the residues that will be mutated in "new_loop_domain"
      storage::List< MutationResidue> aas_to_change( m_ResiduesToChangeCollector->Collect( *new_loop_domain));

      // iterate through the list of residues to mutate and mutate them
      for
      (
        storage::List< MutationResidue>::iterator
          aa_to_change_itr( aas_to_change.Begin()), aa_to_change_itr_end( aas_to_change.End());
        aa_to_change_itr != aa_to_change_itr_end;
        ++aa_to_change_itr
      )
      {
        // get a set of phi and psi angles which will be given to the residue currently denoted by "aa_to_change_itr"
        const storage::VectorND< 2, double> phi_psi( m_DihedralGenerator->operator()( *aa_to_change_itr));

        // make sure the mutation residue is defined
        BCL_Assert( aa_to_change_itr->GetMutationResidue().IsDefined(), "mutation residue is not defined");

        // true if phi is defined
        if( util::IsDefined( phi_psi.First()))
        {
          BCL_MessageDbg
          (
            "setting phi of resi " + aa_to_change_itr->GetMutationResidue()->GetIdentification()
            + " to " + util::Format()( phi_psi.First())
          );
          new_loop_domain->SetPhi( *aa_to_change_itr->GetMutationResidue(), phi_psi.First());
        }
        // true if psi is defined
        if( util::IsDefined( phi_psi.Second()))
        {
          BCL_MessageDbg
          (
            "setting psi of resi " + aa_to_change_itr->GetMutationResidue()->GetIdentification()
            + " to " + util::Format()( phi_psi.Second())
          );
          new_loop_domain->SetPsi( *aa_to_change_itr->GetMutationResidue(), phi_psi.Second());
        }
      } //< loop through mutation residues

      // iterate through the domain and set the geometry of the non rigid sses
      for
      (
        storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator
          segment_itr( new_loop_domain->GetSegments().Begin()), segment_itr_end( new_loop_domain->GetSegments().End());
        segment_itr != segment_itr_end;
        ++segment_itr
      )
      {
        // if not rigid need to set geometry
        if( !segment_itr->IsRigid())
        {
          // set geometry
          util::ShPtr< assemble::SSE> new_sse( segment_itr->GetSSE()->Clone());
          new_sse->SetGeometry();
          std::pair< storage::Set< LoopSegment, LoopSegmentSequenceOrder>::const_iterator, bool> replace_status
          (
            new_loop_domain->ReplaceSegment( LoopSegment( new_sse, segment_itr->IsRigid()))
          );
          segment_itr = replace_status.first;
        }
      }

      // create mutate result "mutate_result" out of "new_loop_domain"
      math::MutateResult< LoopDomain> mutate_result( new_loop_domain, *this);

      // return "mutate_result"
      return mutate_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl
