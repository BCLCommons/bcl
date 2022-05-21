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

#ifndef BCL_FOLD_PLACEMENT_SSE_DISTANCE_RESTRAINT_H_
#define BCL_FOLD_PLACEMENT_SSE_DISTANCE_RESTRAINT_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "assemble/bcl_assemble_protein_model_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementSSEDistanceRestraint
    //! @brief Places an SSE using NOE restraints
    //! @details Randomly choses an SSE in the protein model, weighted by the number of restraints between it and the
    //!          selected SSE.  The SSE is then placed using PlacementSSENextToSSE
    //!
    //! @see @link example_fold_placement_sse_distance_restraint.cpp @endlink
    //! @author weinerbe
    //! @date Jan 11, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementSSEDistanceRestraint :
      public PlacementInterface< assemble::SSE, assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! the atom distance restraint data
      util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > m_Restraints;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PlacementSSEDistanceRestraint();

      //! @brief construct from restraints
      //! @param RESTRAINTS restraints to be used
      PlacementSSEDistanceRestraint( const util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > &RESTRAINTS);

      //! @brief Clone function
      //! @return pointer to new PlacementSSEDistanceRestraint
      PlacementSSEDistanceRestraint *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief generate placement for the given SSE at a random orientation wrt to a located SSE from the given model
      //! @param SELECTED_SSE SiPtr to SSE to be placed
      //! @param PROTEIN_MODEL to which the SSE is going to be added
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SELECTED_SSE,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief checks if the restraint is between 2 SSEs
      //! @param SSE_A first SSE
      //! @param SSE_B second SSE
      //! @param RESTRAINT distance restraint
      //! @return whether the restraint is between 2 SSEs
      static bool ContainsRestraint
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B,
        const restraint::AtomDistance &RESTRAINT
      );

      //! @brief checks if the SSE has the residue in the locator
      //! @param SELECTED_SSE SSE to be checked
      //! @param LOCATOR locator from the restraint
      //! @return whether the SSE has the residue in the locator
      static bool HasResidue
      (
        const assemble::SSE &SELECTED_SSE,
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR
      );

    }; // class PlacementSSEDistanceRestraint

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PLACEMENT_SSE_DISTANCE_RESTRAINT_H_
