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

#ifndef BCL_FOLD_MUTATE_SSE_BEND_RANDOM_H_
#define BCL_FOLD_MUTATE_SSE_BEND_RANDOM_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSSEBendRandom
    //! @brief bends a given SSE randomly within given phi and psi ranges
    //! @details This class selects one random amino acid from given SSE and changes its phi and psi angles
    //! randomly within the ranges defined by the user.
    //!
    //! @see @link example_fold_mutate_sse_bend_random.cpp @endlink
    //! @author karakam
    //! @date 06/08/09
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSSEBendRandom :
      public math::MutateInterface< assemble::SSE>
    {
    private:

    //////////
    // data //
    //////////

      //! range for phi change
      math::Range< double> m_PhiChangeRange;

      //! range for psi change
      math::Range< double> m_PsiChangeRange;

      //! propagation direction of the phi/psi changes
      biol::AASequenceFlexibility::DirectionEnum m_BendingDirection;

      //! scheme
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateSSEBendRandom();

      //! @brief constructor phi/psi change ranges and a scheme
      //! @param PHI_CHANGE_RANGE range for phi changes
      //! @param PSI_CHANGE_RANGE range for psi changes
      //! @param BENDING_DIRECTION direction the phi/psi changes should be propagated towards
      //! @param SCHEME Scheme to be used
      MutateSSEBendRandom
      (
        const math::Range< double> &PHI_CHANGE_RANGE,
        const math::Range< double> &PSI_CHANGE_RANGE,
        const biol::AASequenceFlexibility::SequenceDirection BENDING_DIRECTION = biol::AASequenceFlexibility::e_Bidirectional,
        const std::string &SCHEME = GetStaticClassName< MutateSSEBendRandom>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateSSEBendRandom
      MutateSSEBendRandom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get phi change range
      //! @return phi change range
      const math::Range< double> &GetPhiChangeRange() const;

      //! @brief get phi change range
      //! @return phi change range
      const math::Range< double> &GetPsiChangeRange() const;

      //! @brief return bending direction
      //! @return bending direction
      biol::AASequenceFlexibility::SequenceDirection GetBendingDirection() const;

      //! @brief get bending direction
      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that changes phi/psi angles of a random amino acid in the given SSE
      //! @param THIS_SSE SSE interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::SSE> operator()( const assemble::SSE &THIS_SSE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateSSEBendRandom

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_MUTATE_SSE_BEND_RANDOM_H_
