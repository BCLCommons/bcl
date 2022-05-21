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

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_pick_atom_random.h"

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PickAtomRandom::s_Instance
    (
      GetObjectInstances().AddInstance( new PickAtomRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    PickAtomRandom *PickAtomRandom::Clone() const
    {
      return new PickAtomRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickAtomRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Picks a random atom
    //! @param ATOMS is the SiPtrList which provides the pool of Atom to pick from
    //! @return returns SiPtr to a random Atom in the list
    util::SiPtr< const AtomConformationalInterface>
      PickAtomRandom::Pick( const util::SiPtrList< const AtomConformationalInterface> &ATOMS) const
    {
      // if not found, return empty ShPtr< Atom>
      BCL_Assert( !ATOMS.IsEmpty(), "atom list is empty!!!")

      // return
      return *random::GetGlobalRandom().Iterator( ATOMS.Begin(), ATOMS.End(), ATOMS.GetSize());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PickAtomRandom::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @return ostream which was read from
    std::ostream &PickAtomRandom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
