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

#ifndef BCL_SCORE_AA_ASSIGNMENT_IDENTITY_H_
#define BCL_SCORE_AA_ASSIGNMENT_IDENTITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "function/bcl_function_binary_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAAssignmentIdentity
    //! @brief This is a Function derived class for comparing pairs of AA
    //!
    //! @see @link example_score_aa_assignment_identity.cpp @endlink
    //! @author meilerj
    //! @date 21.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAAssignmentIdentity :
      public function::BinaryInterface< const biol::AABase, const biol::AABase, double>
    {
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
      AAAssignmentIdentity();

      //! @brief virtual copy constructor
      //! @return pointer to a new AAAssignmentIdentity copied from this one
      AAAssignmentIdentity *Clone() const;

      //! @brief destructor
      ~AAAssignmentIdentity();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the score between two assigned members
      //! @param MEMBER_A amino acid A that is compared
      //! @param MEMBER_B amino acid A that is compared
      //! @return 1.0 if they are identical - 0 if they are not
      double operator()( const biol::AABase &MEMBER_A, const biol::AABase &MEMBER_B) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

    }; // class AAAssignmentIdentity

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_ASSIGNMENT_IDENTITY_H_
