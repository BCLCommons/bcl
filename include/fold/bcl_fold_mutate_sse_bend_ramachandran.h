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

#ifndef BCL_FOLD_MUTATE_SSE_BEND_RAMACHANDRAN_H_
#define BCL_FOLD_MUTATE_SSE_BEND_RAMACHANDRAN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSSEBendRamachandran
    //! @brief This class sets the phi/psi of a residue in the given SSE to a random phi/psi from Ramachandran plot
    //! @details This class uses biol::Ramachandran class random number generator to get a random phi/psi value
    //! for a randomly selected residue in the given SSE and sets the phi/psi values of the residue accordingly
    //!
    //! @see @link example_fold_mutate_sse_bend_ramachandran.cpp @endlink
    //! @author karakam
    //! @date Jan 25, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MutateSSEBendRamachandran :
      public math::MutateInterface< assemble::SSE>,
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! range for number of residues to set phi/psi for
      math::Range< size_t> m_NrChanges;

      //! propagation direction of the phi/psi changes
      biol::AASequenceFlexibility::DirectionEnum m_BendingDirection;

      //! scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateSSEBendRamachandran();

      //! @brief constructor from a range of number of residues to change and a scheme
      //! @param NR_RESIDUES_TO_CHANGE_RANGE ange of number of residues to change
      //! @param BENDING_DIRECTION direction the phi/psi changes should be propagated towards
      //! @param SCHEME Scheme to be used
      MutateSSEBendRamachandran
      (
        const math::Range< size_t> &NR_RESIDUES_TO_CHANGE_RANGE,
        const biol::AASequenceFlexibility::SequenceDirection BENDING_DIRECTION = biol::AASequenceFlexibility::e_Bidirectional,
        const std::string &SCHEME = GetStaticClassName< MutateSSEBendRamachandran>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateSSEBendRamachandran
      MutateSSEBendRamachandran *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return range of number of residues to change
      //! @return range of number of residues to change
      const math::Range< size_t> &GetNrResiduesChangeRange() const
      {
        return m_NrChanges;
      }

      //! @brief return bending direction
      //! @return bending direction
      biol::AASequenceFlexibility::SequenceDirection GetBendingDirection() const
      {
        return m_BendingDirection;
      }

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes an SSE and bends amino acids and returns the new SSE
      //! @param THIS_SSE SSE of interest to bend
      //! @return math::MutateResult that has a new SSE with one or more amino acids with new phi/psi values
      math::MutateResult< assemble::SSE> operator()( const assemble::SSE &THIS_SSE) const;

      //! @brief operator that takes an SSE and bends amino acids and returns the new SSE
      //! @param THIS_SSE SSE of interest to bend
      //! @return math::MutateResult that has a new SSE with one or more amino acids with new phi/psi values
      math::MutateResult< assemble::SSE> operator()
      (
        const assemble::SSE &THIS_SSE,
        const util::SiPtr< const biol::Membrane> &MEMBRANE
      ) const;

      //! @brief operator that takes an ProteinModel and bends amino acids in one random SSE and returns the new ProteinModel
      //! @param PROTEIN_MODEL the protein model of in
      //! @return math::MutateResult that has a new SSE with one or more amino acids with new phi/psi values
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateSSEBendRamachandran

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SSE_BEND_RAMACHANDRAN_H_
