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

#ifndef BCL_ASSEMBLE_SSE_POOL_MOVE_AA_H_
#define BCL_ASSEMBLE_SSE_POOL_MOVE_AA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPoolMoveAA
    //! @brief moves residues between two adjacent secondary structure elements
    //!
    //! @see @link example_assemble_sse_pool_move_aa.cpp @endlink
    //! @author woetzen
    //! @date Jun 18, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPoolMoveAA :
      public math::MutateInterface< SSEPool>
    {

    private:

    //////////
    // data //
    //////////

      //! max number of residues to move
      math::Range< size_t> m_ResdiuesToMoveRange;

      //! scheme
      std::string m_Scheme;

    public:

      //! @brief the default scheme for this class
      static const std::string &GetDefaultScheme();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from max number of residues to move
      //! @param RESIDUES_TO_MOVE_RANGE range of number of residues to move
      //! @param SCHEME the scheme
      SSEPoolMoveAA
      (
        const math::Range< size_t> &RESIDUES_TO_MOVE_RANGE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new SSEPoolMoveAA
      SSEPoolMoveAA *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that mutates a pool by shifting AAs between SSEs
      //! @param SSE_POOL
      //! @return MutateResult that results from mutating to the SSE_POOL
      math::MutateResult< SSEPool> operator()( const SSEPool &SSE_POOL) const;

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

    }; // class SSEPoolMoveAA

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_POOL_MOVE_AA_H_ 
