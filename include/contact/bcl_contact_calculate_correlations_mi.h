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

#ifndef BCL_CONTACT_CALCULATE_CORRELATIONS_MI_H_
#define BCL_CONTACT_CALCULATE_CORRELATIONS_MI_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "align/bcl_align.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_calculate_correlations_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CalculateCorrelationsMI
    //! @brief Calculates a correlation matrix based on entropy/mutual information algorithm
    //! @details Creates a correlation matrix based on the mutual information available for all pairs of positions
    //!            in columns i and j where i != j of the sum of P(xi,ji) log_b( P(xi,yj)/((P(xi)P(xj)))) with b set to 20
    //!
    //! @see @link example_contact_calculate_correlations_mi.cpp @endlink
    //! @author teixeipl
    //! @date Aug 4, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CalculateCorrelationsMI :
      public CalculateCorrelationsInterface< biol::AABase>
    {

    private:
      // Base for log
      // TODO: Check reference and switch to 2 for mutual information output
      // TODO: Add entropy to statistics
      static const size_t s_LogBase;

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CalculateCorrelationsMI();

      //! @brief Clone function
      //! @return pointer to new CalculateCorrelationsMI
      CalculateCorrelationsMI *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // TODO: ADD COMMENT
      CorrelationMatrix operator ()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT_INTERFACE) const;

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

    private:

    }; // class CalculateCorrelationsMI

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_CALCULATE_CORRELATIONS_MI_H_
