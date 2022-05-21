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
#include "chemistry/bcl_chemistry_reaction_ensemble.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_rxn_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct from an input stream
    //! @param ISTREAM input stream
    //! @param RANGE the range of reactions to load from the stream
    ReactionEnsemble::ReactionEnsemble
    (
      std::istream &ISTREAM,
      const math::Range< size_t> &RANGE
    )
    {
      ReadMoreFromRXN( ISTREAM, RANGE);
    }

    //! @brief Clone function
    //! @return pointer to new ReactionEnsemble
    ReactionEnsemble *ReactionEnsemble::Clone() const
    {
      return new ReactionEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ReactionEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief read additional molecules into the ensemble
    //! @param ISTREAM input stream, reads in RXN format
    //! @param RANGE a range of reactions to load, by default loads all reactions in the stream
    void ReactionEnsemble::ReadMoreFromRXN
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
      sdf::RXNHandler handler;
      for( size_t rxn_no( 0); !ISTREAM.eof() && rxn_no < n_to_read; ++rxn_no)
      {
        handler.ReadFromRXN( ISTREAM);
        
        if( rxn_no >= RANGE.GetMin())
        {
          m_ReactionEnsemble.PushBack
          ( 
            sdf::RXNFactory::MakeReactionComplete( handler)
          );
        }
      }

      BCL_MessageVrb( "successfully read ensemble of " + util::Format()( GetSize() - original_size) + " reactions.");
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ReactionEnsemble::WriteRXN( std::ostream &OSTREAM) const
    {
      // iterate through all rxnecules in the given ensemble
      for
      (
        storage::List< ReactionComplete>::const_iterator
          itr_rxns( m_ReactionEnsemble.Begin()),
          itr_rxns_end( m_ReactionEnsemble.End());
        itr_rxns != itr_rxns_end;
        ++itr_rxns
      )
      {
        itr_rxns->WriteRXN( OSTREAM);
      }
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ReactionEnsemble::Read( std::istream &ISTREAM)
    {
      ReadMoreFromRXN( ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ReactionEnsemble::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      WriteRXN( OSTREAM);
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
