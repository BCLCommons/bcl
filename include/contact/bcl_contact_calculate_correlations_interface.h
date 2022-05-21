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

#ifndef BCL_CONTACT_CALCULATE_CORRELATIONS_INTERFACE_H_
#define BCL_CONTACT_CALCULATE_CORRELATIONS_INTERFACE_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "align/bcl_align.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_correlation_matrix.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CalculateCorrelationsInterface
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @tparam t_Member TODO: add a description for the template parameter
    //!
    //! @see @link example_contact_CalculateCorrelationsInterface.cpp @endlink
    //! @author teixeipl
    //! @date Jul 31, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! TODO: remove this and derive straight from function interface

    template< typename t_Member>
    class CalculateCorrelationsInterface :
      public util::FunctionInterface< align::AlignmentInterface< t_Member>, CorrelationMatrix>
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new CalculateCorrelationsInterface< t_Member>
      virtual CalculateCorrelationsInterface< t_Member> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief all derived classes should calculate a correlation matrix if given an AlignmentInterface of type t_Member
      //! @param ALIGNMENT_INTERFACE is the MSA representation from which the CM is calculated
      //! @return Returns a correlation matrix of dimensions N by N where N is the length of the MSA, containing doubles
      virtual CorrelationMatrix operator()( const align::AlignmentInterface< t_Member> &ALIGNMENT_INTERFACE) const = 0;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // template class CalculateCorrelationsInterface

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_CALCULATE_CORRELATIONS_INTERFACE_H_
