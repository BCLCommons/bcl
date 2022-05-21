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

#ifndef BCL_ASSEMBLE_SSE_POOL_INSERT_COIL_INTO_SSE_H_
#define BCL_ASSEMBLE_SSE_POOL_INSERT_COIL_INTO_SSE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "sspred/bcl_sspred_methods.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPoolInsertCoilIntoSSE
    //! @brief mutates a pool by picking an sse and splitting it into three SSEs, inserting a coil of lowest average probability
    //!
    //! @see @link example_assemble_sse_pool_insert_coil_into_sse.cpp @endlink
    //! @author woetzen
    //! @date Jun 18, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPoolInsertCoilIntoSSE :
      public math::MutateInterface< SSEPool>
    {

    private:

    //////////
    // data //
    //////////

      //! range for how many residues from the middle should be split
      sspred::Method m_Method;

      //! locator of an sse in a domain
      util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > m_SSELocator;

      //! length of coil to insert
      size_t m_CoilLength;

      //! scheme for this mutate
      std::string m_Scheme;

    public:

      //! @brief the default scheme for this class
      static const std::string &GetDefaultScheme();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from locator and mutate
      //! @param SS_METHOD the method to use to locate the largest drop in the ss prediction for the located sse
      //! @param SP_LOCATE_SSE picker of sse from domain
      //! @param COIL_LENGTH length of the coil to insert
      //! @param SCHEME the scheme
      SSEPoolInsertCoilIntoSSE
      (
        const sspred::Method &SS_METHOD,
        const util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > &SP_LOCATE_SSE,
        const size_t COIL_LENGTH,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new SSEPoolInsertCoilIntoSSE
      SSEPoolInsertCoilIntoSSE *Clone() const;

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

      //! @brief operator that mutates a pool by splitting a single sse
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

    }; // class SSEPoolInsertCoilIntoSSE

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_POOL_INSERT_COIL_INTO_SSE_H_ 
