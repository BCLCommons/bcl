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
#include "chemistry/bcl_chemistry_collector_valence.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! clone constructor
    CollectorValence *CollectorValence::Clone() const
    {
      return new CollectorValence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorValence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &CollectorValence::GetAlias() const
    {
      static const std::string s_name( "CollectorValence");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorValence::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects atoms containing a free valence from a small molecule.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate translates the CollectorValence denoting an atom with its valence electrons into the actual atom
    //! @param FRAGMENT fragment of interest
    //! @return pointer to atom that was located with valences
    util::SiPtrList< const AtomConformationalInterface>
      CollectorValence::Collect( const FragmentComplete &FRAGMENT) const
    {
      // ShPtrList of atoms with open valences
      util::SiPtrList< const AtomConformationalInterface> valences;

      // iterate over all atoms in MOLECULE
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( FRAGMENT.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        // determine the # of implict hydrogens
        const size_t number_valences( itr_atoms->GetNumberValenceBonds());
        if( number_valences && util::IsDefined( number_valences))
        {
          // if free valences exist, return atom containing those valences
          valences.Append( util::SiPtrList< const AtomConformationalInterface>( number_valences, *itr_atoms));
        }
      } // end itr_atoms

      // return
      return valences;
    }

    //! @brief read from std::istream
    std::istream &CollectorValence::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &CollectorValence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
