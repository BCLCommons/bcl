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

#ifndef BCL_RESTRAINT_EPR_DISTANCE_ASSIGNMENT_H_
#define BCL_RESTRAINT_EPR_DISTANCE_ASSIGNMENT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EPRDistanceAssignment
    //! @brief This class stores the exposure, Ca and Cb atoms and assigned experimental distance for an amino acid pair
    //!
    //! @see @link example_restraint_epr_distance_assignment.cpp @endlink
    //! @author fischea
    //! @date May 20, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EPRDistanceAssignment :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! Ca atom of the first residue
      biol::Atom m_FirstAtomCa;

      //! Cb atom of the first residue
      biol::Atom m_FirstAtomCb;

      //! Ca atom of the second residue
      biol::Atom m_SecondAtomCa;

      //! Cb atom of the second residue
      biol::Atom m_SecondAtomCb;

      //! experimentally determined distance
      double m_ExpDistance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from atoms and experimental distance
      //! @param FIRST_ATOM_A Ca atom of the first residue
      //! @param FIRST_ATOM_B Cb atom of the first residue
      //! @param SECOND_ATOM_A Ca atom of the second residue
      //! @param SECOND_ATOM_B Cb atom of the second residue
      EPRDistanceAssignment
      (
        const biol::Atom &FIRST_ATOM_A,
        const biol::Atom &FIRST_ATOM_B,
        const biol::Atom &SECOND_ATOM_A,
        const biol::Atom &SECOND_ATOM_B,
        const double DISTANCE
      );

      //! @brief returns a pointer to a new EPRDistanceAssignment
      //! @return pointer to a new EPRDistanceAssignment
      EPRDistanceAssignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads in members from stream
      //! @param ISTREAM stream to read members from
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write members into a stream
      //! @param OSTREAM stream to write members into
      //! @INDENT number of indentations to use
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class EPRDistanceAssignment

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_EPR_DISTANCE_ASSIGNMENT_H_
