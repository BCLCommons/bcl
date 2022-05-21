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
#include "chemistry/bcl_chemistry_molecule_ensemble.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct from an input stream
    //! @param ISTREAM input stream
    //! @param RANGE the range of small molecules to load from the stream
    MoleculeEnsemble::MoleculeEnsemble
    (
      std::istream &ISTREAM,
      const math::Range< size_t> &RANGE
    )
    {
      ReadMoreFromMdl( ISTREAM, RANGE);
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeEnsemble
    MoleculeEnsemble *MoleculeEnsemble::Clone() const
    {
      return new MoleculeEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief read additional molecules into the ensemble from an mdl stream
    //! @param ISTREAM input stream, reads in SDF format
    //! @param RANGE a range of small molecules to load, by default loads all small molecules in the stream
    void MoleculeEnsemble::ReadMoreFromMdl
    (
      std::istream &ISTREAM,
      const math::Range< size_t> &RANGE
    )
    {
      // get the original size
      const size_t original_size( GetSize());

      // close the borders on the range for ease of use
      math::Range< size_t> closed_range( RANGE.CloseBorders());

      // intended # to read; max is necessary in case the range is the range of size_t
      const size_t n_to_read( std::max( closed_range.GetWidth(), closed_range.GetWidth() + size_t( 1)));

      // read in from the stream until we reach the end of the file or the last index
      for( FragmentFeed feed( ISTREAM, sdf::e_Saturate, n_to_read, closed_range.GetMin()); feed.NotAtEnd(); ++feed)
      {
        m_MoleculeEnsemble.PushBack( *feed);
      }

      BCL_MessageStd( "finished reading ensemble with " + util::Format()( GetSize() - original_size) + " molecules.");
    }

  //////////////////////
  // input and output //
  //////////////////////

    std::istream &MoleculeEnsemble::Read( std::istream &ISTREAM)
    {
      //read
      ReadMoreFromMdl( ISTREAM);
      return ISTREAM;
    }

    std::ostream &MoleculeEnsemble::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      WriteMDL( OSTREAM);
      return OSTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &MoleculeEnsemble::WriteMDL( std::ostream &OSTREAM) const
    {
      // iterate through all molecules in the given ensemble
      for
      (
        storage::List< MoleculeComplete>::const_iterator
          itr_mols( m_MoleculeEnsemble.Begin()),
          itr_mols_end( m_MoleculeEnsemble.End());
        itr_mols != itr_mols_end;
        ++itr_mols
      )
      {
        itr_mols->WriteMDL( OSTREAM);
      }
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
