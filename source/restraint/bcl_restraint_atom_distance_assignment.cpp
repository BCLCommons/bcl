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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "restraint/bcl_restraint_atom_distance_assignment.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AtomDistanceAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomDistanceAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomDistanceAssignment::AtomDistanceAssignment() :
      m_AtomA(),
      m_AtomB(),
      m_Distance()
    {
    }

    //! @brief construct from atoms and distance
    //! @param ATOM_A first atom
    //! @param ATOM_B second atom
    //! @param DISTANCE experimental distance
    AtomDistanceAssignment::AtomDistanceAssignment
    (
      const biol::Atom &ATOM_A,
      const biol::Atom &ATOM_B,
      const util::ShPtr< Distance> &DISTANCE
    ) :
      m_AtomA( ATOM_A),
      m_AtomB( ATOM_B),
      m_Distance( DISTANCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AtomDistanceAssignment
    AtomDistanceAssignment *AtomDistanceAssignment::Clone() const
    {
      return new AtomDistanceAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomDistanceAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief calculates the distance between the member atoms
    //! @return the distance between the member atoms
    double AtomDistanceAssignment::CalculateAtomDistance() const
    {
      return linal::Distance( m_AtomA.GetCoordinates(), m_AtomB.GetCoordinates());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomDistanceAssignment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AtomA, ISTREAM);
      io::Serialize::Read( m_AtomB, ISTREAM);
      io::Serialize::Read( m_Distance, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AtomDistanceAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AtomA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomB, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Distance, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
  
} // namespace bcl
