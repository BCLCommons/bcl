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
#include "fold/bcl_fold_collector_loop_domain_all_non_rigid.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_loop_domain.h"
#include "fold/bcl_fold_mutation_residue.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorLoopDomainAllNonRigid::s_Instance
    (
      GetObjectInstances().AddInstance( new CollectorLoopDomainAllNonRigid())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorLoopDomainAllNonRigid::CollectorLoopDomainAllNonRigid()
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateLoopDomainAllNonRigid
    CollectorLoopDomainAllNonRigid *CollectorLoopDomainAllNonRigid::Clone() const
    {
      return new CollectorLoopDomainAllNonRigid( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CollectorLoopDomainAllNonRigid::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual Collect taking an ARGUMENT and returning a list of MutationResidues
    //! @param LOOP_DOMAIN LoopDomain of interest
    //! @return storage::List< MutationResidue> that results from collecting from "LOOP_DOMAIN"
    storage::List< MutationResidue> CollectorLoopDomainAllNonRigid::Collect( const LoopDomain &LOOP_DOMAIN) const
    {
      // create list which will be used to hold mutation residues for all of the non rigid residues
      storage::List< MutationResidue> mutation_residues;

      storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > residues( LOOP_DOMAIN.GetResidues());

      // iterate through the residues
      for
      (
        storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> >::const_iterator
          resi_itr( residues.Begin()), resi_itr_end( residues.End());
        resi_itr != resi_itr_end;
        ++resi_itr
      )
      {
        // true if the segment is rigid, meaning it cannot be changed
        if( resi_itr->Second())
        {
          // go to next residue
          continue;
        }

        BCL_MessageDbg( "currently making mutation residue for " + resi_itr->First()->GetIdentification());
        // create mutation residue
        const MutationResidue current_mutation_residue( resi_itr, residues);

        mutation_residues.PushBack( current_mutation_residue);
      }

      // return
      return mutation_residues;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
