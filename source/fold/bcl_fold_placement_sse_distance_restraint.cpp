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
#include "fold/bcl_fold_placement_sse_distance_restraint.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_placement_sse_next_to_sse.h"
#include "restraint/bcl_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PlacementSSEDistanceRestraint::s_Instance
    (
      GetObjectInstances().AddInstance( new PlacementSSEDistanceRestraint())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PlacementSSEDistanceRestraint::PlacementSSEDistanceRestraint()
    {
    }

    //! @brief construct from restraints
    //! @param RESTRAINTS restraints to be used
    PlacementSSEDistanceRestraint::PlacementSSEDistanceRestraint
    (
      const util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > &RESTRAINTS
    ) :
      m_Restraints( RESTRAINTS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PlacementSSENOE
    PlacementSSEDistanceRestraint *PlacementSSEDistanceRestraint::Clone() const
    {
      return new PlacementSSEDistanceRestraint( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PlacementSSEDistanceRestraint::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate placement for the given SSE at a random orientation wrt to a located SSE from the given model
    //! @param SELECTED_SSE SiPtr to SSE to be placed
    //! @param PROTEIN_MODEL to which the SSE is going to be added
    storage::Pair< math::TransformationMatrix3D, bool> PlacementSSEDistanceRestraint::Place
    (
      const assemble::SSE &SELECTED_SSE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      static const math::TransformationMatrix3D s_def_transformation;

      // if no restraints were found
      if( !m_Restraints.IsDefined() || m_Restraints->IsEmpty())
      {
        return storage::Pair< math::TransformationMatrix3D, bool>( s_def_transformation, false);
      }

      // get the protein model sses
      const util::SiPtrVector< const assemble::SSE> model_sses( PROTEIN_MODEL.GetSSEs());

      // if the model is empty
      if( model_sses.IsEmpty())
      {
        // place the sse at the origin
        return storage::Pair< math::TransformationMatrix3D, bool>( s_def_transformation, true);
      }

      // create map to hold # restraints between SSEs
      storage::Map< util::SiPtr< const assemble::SSE>, size_t> restraint_map;

      // iterate over the SSEs
      size_t number_restraints( 0);
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( model_sses.Begin()),
          sse_itr_end( model_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // iterate through the restraints and track the total number
        for
        (
          util::ShPtrVector< restraint::AtomDistance>::const_iterator restraint_itr( m_Restraints->Begin()),
            restraint_itr_end( m_Restraints->End());
          restraint_itr != restraint_itr_end; ++restraint_itr
        )
        {
          // if the restraint is between the SSEs
          if( ContainsRestraint( SELECTED_SSE, **sse_itr, **restraint_itr))
          {
            ++number_restraints;
            // if the SSE is in the map
            if( restraint_map.Has( *sse_itr))
            {
              // increment the count
              ++restraint_map[ *sse_itr];
            }
            // the SSE is not in the map
            else
            {
              // set the count to 1
              restraint_map[ *sse_itr] = 1;
            }
          }
        }
      }

      // randomly select an SSE weighted by its number of restraints
      int weight_index( random::GetGlobalRandom().Random< int>( 0, number_restraints));

      // iterate over the map
      util::SiPtr< const assemble::SSE> sp_random_sse;
      for
      (
        storage::Map< util::SiPtr< const assemble::SSE>, size_t>::const_iterator map_itr( restraint_map.Begin()),
          map_itr_end( restraint_map.End());
         map_itr != map_itr_end; ++map_itr
      )
      {
        weight_index -= map_itr->second;
        if( weight_index <= 0)
        {
          sp_random_sse = map_itr->first;
        }
      }

      // return a placement using PlacementSSENextToSSE
      return PlacementSSENextToSSE::Place( SELECTED_SSE, *sp_random_sse);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PlacementSSEDistanceRestraint::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PlacementSSEDistanceRestraint::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief checks if the restraint is between 2 SSEs
    //! @param SSE_A first SSE
    //! @param SSE_B second SSE
    //! @param RESTRAINT distance restraint
    //! @return whether the restraint is between 2 SSEs
    bool PlacementSSEDistanceRestraint::ContainsRestraint
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B,
      const restraint::AtomDistance &RESTRAINT
    )
    {
      return
        (
          HasResidue( SSE_A, RESTRAINT.GetData().First()) &&
          HasResidue( SSE_B, RESTRAINT.GetData().Second())
        ) ||
        (
          HasResidue( SSE_B, RESTRAINT.GetData().First()) &&
          HasResidue( SSE_A, RESTRAINT.GetData().Second())
        );
    }

    //! @brief checks if the SSE has the residue in the locator
    //! @param SELECTED_SSE SSE to be checked
    //! @param LOCATOR locator from the restraint
    //! @return whether the SSE has the residue in the locator
    bool PlacementSSEDistanceRestraint::HasResidue
    (
      const assemble::SSE &SELECTED_SSE,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR
    )
    {
      // if the chain id doesn't match
      if( !LOCATOR.IsDefined() || SELECTED_SSE.GetChainID() != LOCATOR->GetChainID() || SELECTED_SSE.GetSize() == 0)
      {
        return false;
      }

      // construct a range
      const math::Range< int> seq_id_range
      (
        SELECTED_SSE.GetFirstAA()->GetSeqID(), SELECTED_SSE.GetLastAA()->GetSeqID()
      );

      // return whether the seq id is in the range
      return seq_id_range.IsWithin( LOCATOR->GetSeqID());
    }

  } // namespace fold
} // namespace bcl
