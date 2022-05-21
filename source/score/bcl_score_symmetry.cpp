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
#include "score/bcl_score_symmetry.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    //! @brief construct from file that contains movable locators
    //! this function has to be specialized for each t_ArgumentType within the cpp
    //! @param FILENAME filename that contains some kind of movable locators
    //! @param SCHEME
    template<>
    Symmetry< assemble::ProteinModel>::Symmetry
    (
      const std::string &FILENAME,
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME),
      m_MovableLocators()
    {
      // create stream to file
      io::IFStream read;
      io::File::MustOpenIFStream( read, FILENAME);
      BCL_MessageStd( "read symmetry file: " + FILENAME);

      // create "symmetry_definitions" and initialize with the information from "read"
      storage::Vector< storage::Vector< std::string> > symmetry_definitions
      (
        util::SplittedStringLineListFromIStream( read)
      );

      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          line_iter( symmetry_definitions.Begin()), line_iter_end( symmetry_definitions.End());
        line_iter != line_iter_end;
        ++line_iter
      )
      {
        util::ShPtrVector< find::LocatorInterface< linal::Vector3D, assemble::ProteinModel> > symmetry_unit;
        for
        (
          storage::Vector< std::string>::const_iterator
            point_iter( line_iter->Begin()), point_iter_end( line_iter->End());
          point_iter != point_iter_end;
          ++point_iter
        )
        {
          const char chain_char( point_iter->operator[]( 0));
          BCL_MessageDbg( "chain is " + util::Format()( chain_char));

          // go to next (it will be residue number)
          ++point_iter;
          const int residue( util::ConvertStringToNumericalValue< int>( *point_iter));
          BCL_MessageDbg( "residue is " + util::Format()( residue));

          // go to next (it will be atom name)
          ++point_iter;
          const biol::AtomType atom_type( biol::GetAtomTypes().TypeFromPDBAtomName( *point_iter));
          BCL_MessageDbg( "atom name is " + util::Format()( atom_type->GetName()));

          const util::ShPtr< find::LocatorInterface< linal::Vector3D, assemble::ProteinModel> > coordinate_locator
          (
            new assemble::LocatorAtom( chain_char, residue, atom_type)
          );

          symmetry_unit.PushBack( coordinate_locator);
        }
        BCL_MessageDbg( "Number of points in symmetry unit is " + util::Format()( symmetry_unit.GetSize()));
        m_MovableLocators.PushBack( symmetry_unit);
      }

      BCL_MessageDbg( "Number of symmetry units is " + util::Format()( m_MovableLocators.GetSize()));
    }

  } // namespace score
} // namespace bcl
