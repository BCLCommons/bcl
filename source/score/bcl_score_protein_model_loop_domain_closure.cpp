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
#include "score/bcl_score_protein_model_loop_domain_closure.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_loop_domain.h"
#include "fold/bcl_fold_mutate_protein_model_loop_domain_ccd.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelLoopDomainClosure::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelLoopDomainClosure())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelLoopDomainClosure::ProteinModelLoopDomainClosure()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelLoopDomainClosure
    ProteinModelLoopDomainClosure *ProteinModelLoopDomainClosure::Clone() const
    {
      return new ProteinModelLoopDomainClosure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelLoopDomainClosure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &ProteinModelLoopDomainClosure::GetAlias() const
    {
      static const std::string s_name( "ProteinModelLoopDomainClosure");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelLoopDomainClosure::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Scores how close the a loop is to being closed according to the psuedo residue, target residue, and atoms"
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score after replicating the protein model using symmetry
    //! @param PROTEIN_MODEL protein model to be scored
    //! @return the score after replicating the protein model using symmetry
    double ProteinModelLoopDomainClosure::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // get the loop domain locator from model data
      util::ShPtr< util::ShPtrList< fold::LocatorLoopDomain> > sp_loop_domains(
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_LoopDomainLocators));

      // initialize score
      double score( 0.0);

      // make sure it's defined
      if( !sp_loop_domains.IsDefined())
      {
        BCL_MessageDbg( "Loop domain locators are not stored with the given model");
        return score;
      }

      // iterate through the list of loop domains
      for
      (
        auto loop_itr( sp_loop_domains->Begin()), loop_itr_end( sp_loop_domains->End());
        loop_itr != loop_itr_end;
        ++loop_itr
      )
      {
        util::ShPtr< fold::LoopDomain> sp_loop_domain( ( *loop_itr)->Locate( PROTEIN_MODEL));
        if( !sp_loop_domain.IsDefined())
        {
          continue;
        }

        // create list of TargetAndMovingPointPair objects using the current model
        const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> target_and_moving_points
        (
          sp_loop_domain->TargetAndMovingPointsForCCD( PROTEIN_MODEL)
        );

        if( target_and_moving_points.IsEmpty())
        {
          // score does not change if no target and moving points. this unfortunately means emptiness is favored
          // but protein model should be complete at anyways. score entropy scores inverse (i.e. empty model)
          BCL_MessageStd(
            "There are no target and moving points for loop domain " + util::Format()( ( *loop_itr)->GetIdentification()));
          continue;
        }

        double rmsd( fold::LocatorLoopDomain::CalculateSquareDistanceSum( target_and_moving_points));

        BCL_MessageDbg(
          " loop domain score sum_square_distance for: " + ( *loop_itr)->GetLoopSegments().Begin()->GetIdentification() + " is: " + util::Format()( rmsd));

        rmsd /= target_and_moving_points.GetSize();
        rmsd = math::Sqrt( rmsd);

        static const double s_peptide_bond_length( biol::GetAtomTypes().C->GetBondLength( biol::GetAtomTypes().N));

        // sum up the score
        score += std::log( rmsd / s_peptide_bond_length);
      }

      // return the score
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace score
} // namespace bcl
