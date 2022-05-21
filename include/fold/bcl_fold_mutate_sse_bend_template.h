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

#ifndef BCL_FOLD_MUTATE_SSE_BEND_TEMPLATE_H_
#define BCL_FOLD_MUTATE_SSE_BEND_TEMPLATE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSSEBendTemplate
    //! @brief bends an SSE according to a randomly selected geometry
    //! @details selects a random geometry from a fold template that satisfies the comparison member function and then
    //!          will bend the given SSE to match it
    //!
    //! @see @link example_fold_mutate_sse_bend_template.cpp @endlink
    //! @author weinerbe
    //! @date Jul 1, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSSEBendTemplate :
      public math::MutateInterface< assemble::SSE>
    {

    private:

    //////////
    // data //
    //////////

      //! comparison function for deciding if template geometry can be used with passed SSE
      util::ShPtr
      <
        math::BinaryFunctionInterface< assemble::SSE, assemble::SSEGeometryPhiPsi, bool>
      > m_SSEGeometryCompare;

      //! scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateSSEBendTemplate();

      //! @brief construct from SSE geometry comparison function
      //! @param SSE_GEOMETRY_COMPARE SSE geometry comparison function
      //! @param SCHEME Scheme to be used
      MutateSSEBendTemplate
      (
        const math::BinaryFunctionInterface< assemble::SSE, assemble::SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE,
        const std::string &SCHEME = GetStaticClassName< MutateSSEBendTemplate>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateSSEBendTemplate
      MutateSSEBendTemplate *Clone() const;

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

      //! @brief operator that takes an SSE and bends amino acids and returns the new SSE
      //! @param THIS_SSE SSE of interest to bend
      //! @return math::MutateResult that has a new SSE bent according to the random geometry from a fold template
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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateSSEBendTemplate

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SSE_BEND_TEMPLATE_H_ 
