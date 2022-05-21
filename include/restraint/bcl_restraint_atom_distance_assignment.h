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

#ifndef BCL_RESTRAINT_ATOM_DISTANCE_ASSIGNMENT_H_
#define BCL_RESTRAINT_ATOM_DISTANCE_ASSIGNMENT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_distance.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomDistanceAssignment
    //! @brief Stores two located atoms and the experimental distance between them.
    //! @details Contains two atoms, typically generated using the GenerateAssignment function in AtomDistance.  The
    //!          calculated distance can be simply calculated from the two atom coordinates.  Collections of these
    //!          assignments can then be scored by scoring classes.
    //!
    //! @see @link example_restraint_atom_distance_assignment.cpp @endlink
    //! @author weinerbe
    //! @date Mar 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomDistanceAssignment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! first atom
      biol::Atom m_AtomA;

      //! second atom
      biol::Atom m_AtomB;

      //! experimental distance
      util::ShPtr< Distance> m_Distance;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AtomDistanceAssignment();

      //! @brief construct from atoms and distance
      //! @param ATOM_A first atom
      //! @param ATOM_B second atom
      //! @param DISTANCE experimental distance
      AtomDistanceAssignment
      (
        const biol::Atom &ATOM_A,
        const biol::Atom &ATOM_B,
        const util::ShPtr< Distance> &DISTANCE
      );

      //! @brief Clone function
      //! @return pointer to new AtomDistanceAssignment
      AtomDistanceAssignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets atom A
      //! @return atom A
      const biol::Atom &GetAtomA() const
      {
        return m_AtomA;
      }

      //! @brief gets atom B
      //! @return atom B
      const biol::Atom &GetAtomB() const
      {
        return m_AtomB;
      }

      //! @brief gets the distance object
      //! @return the distance object
      const Distance &GetDistanceObject() const
      {
        return *m_Distance;
      }

      //! @brief gets the experimental distance
      //! @return the experimental distance
      double GetDistance() const
      {
        return m_Distance->GetDistance();
      }

      //! @brief gets the experimental upper bound
      //! @return the experimental upper bound
      double GetUpperBound() const
      {
        return m_Distance->UpperBound();
      }

      //! @brief gets the experimental lower bound
      //! @return the experimental lower bound
      double GetLowerBound() const
      {
        return m_Distance->LowerBound();
      }

      //! @brief calculates the distance between the member atoms
      //! @return the distance between the member atoms
      double CalculateAtomDistance() const;

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

    }; // class AtomDistanceAssignment

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ATOM_DISTANCE_ASSIGNMENT_H_ 
