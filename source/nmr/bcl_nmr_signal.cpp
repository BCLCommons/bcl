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
#include "nmr/bcl_nmr_signal.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Signal::s_Instance
    (
      GetObjectInstances().AddInstance( new Signal())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief check if given atom is involved in that signal
    //! @param ATOM the atom of interest
    //! @return true, if any Signal1D has the given atom
    bool Signal::ContainsAtom( const chemistry::AtomConformationalInterface &ATOM) const
    {
      const util::SiPtr< const chemistry::AtomConformationalInterface> si_ptr_atom( &ATOM);
      for
      (
        util::ShPtrVector< Signal1D>::const_iterator itr( m_Signals1D.Begin()), itr_end( m_Signals1D.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetAtomInvolvedInSignal() == si_ptr_atom)
        {
          return true;
        }
      }

      // no involved atoms agree
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Signal from std::istream
    //! @param ISTREAM input stream that contains Signal object
    std::istream &Signal::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_Signals1D, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write Signal into std::ostream
    //! @param OSTREAM output stream into that the Signal object is written
    //! @param INDENT indentation
    std::ostream &Signal::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_Signals1D, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace nmr
} // namespace bcl

