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

#ifndef BCL_FOLD_MUTATE_SSE_TYPE_H_
#define BCL_FOLD_MUTATE_SSE_TYPE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSSEType
    //! @brief changes an SSEs type and conformation accordingly
    //! @details from a defined set of ss types, an sstype is chose that the sse is not, this will be set and the
    //!          conformation will be changed accordingly, preserving the current position
    //!
    //! @see @link example_fold_mutate_sse_type.cpp @endlink
    //! @author woetzen
    //! @date Jun 16, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSSEType :
      public math::MutateInterface< assemble::SSE>
    {

    private:

    //////////
    // data //
    //////////

      //! the type that the sse can be set to
      storage::Set< biol::SSType> m_PossibleTypes;

      //! scheme of mutate
      std::string m_Scheme;

    public:

      //! @brief the default scheme
      //! @return reference to default scheme as string
      static const std::string &GetDefaultScheme();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param POSSIBLE_TYPES set of types that sses can be set to
      //! @param SCHEME the scheme of the mutate
      MutateSSEType
      (
        const storage::Set< biol::SSType> &POSSIBLE_TYPES,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new MutateSSEType
      MutateSSEType *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief change the sse type and the conformation accordingly
      //! @param ELEMENT the sse to change
      //! @return MutateResult that results from mutating to the SSE
      math::MutateResult< assemble::SSE> operator()( const assemble::SSE &ELEMENT) const;

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

    }; // class MutateSSEType

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SSE_TYPE_H_ 
