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
#include "fold/bcl_fold_locator_missing_coordinates.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> LocatorMissingCoordinates::s_Instance
    (
      util::Enumerated< LocatorMissingCoordinates>::AddInstance( new LocatorMissingCoordinates())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from members
    //! @param IGNORE_TERMS ignore terminal loops
    LocatorMissingCoordinates::LocatorMissingCoordinates( bool IGNORE_TERMS) :
      m_IgnoreTerms( IGNORE_TERMS)
    {
    }

    //! @brief clone function
    //! @return pointer to a new LocatorMissingCoordinates
    LocatorMissingCoordinates *LocatorMissingCoordinates::Clone() const
    {
      return new LocatorMissingCoordinates( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &LocatorMissingCoordinates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &LocatorMissingCoordinates::GetAlias() const
    {
      static const std::string s_alias( "LocatorMissingCoordinates");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorMissingCoordinates::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Locates sequence spans with missing backbone coordinates in a protein model.");
      serializer.AddInitializer
      (
        "ignore termini",
        "ignore terminal loop regions",
        io::Serialization::GetAgent( &m_IgnoreTerms)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return regions with undefined backbone coordinates in the given protein model
    //! @param MODEL protein model for which to return regions with undefined backbone coordinates
    //! @return regions with undefined backbone coordinates
    storage::Vector< LocatorMissingCoordinates::Span> LocatorMissingCoordinates::Locate
    (
      const assemble::ProteinModel &MODEL
    ) const
    {
      storage::Vector< Span> spans;
      const util::SiPtrVector< const biol::AASequence> chains( MODEL.GetSequences());
      for( auto chain_it( chains.Begin()), chain_it_end( chains.End()); chain_it != chain_it_end; ++chain_it)
      {
        spans.Append( Locate( **chain_it));
      }

      return spans;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return regions with undefined backbone coordinates in the given sequence
    //! @param CHAIN sequence for which to return regions with undefined backbone coordinates
    //! @return regions with undefined backbone coordinates
    storage::Vector< LocatorMissingCoordinates::Span> LocatorMissingCoordinates::Locate( const biol::AASequence &CHAIN) const
    {
      storage::Vector< Span> spans;
      const util::ShPtrVector< biol::AABase> aas( CHAIN.GetData());
      util::ShPtr< biol::AABase> def_n;
      util::ShPtr< biol::AABase> def_c;
      for( auto res_it( aas.Begin()), res_it_end( aas.End()); res_it != res_it_end; ++res_it)
      {
        const biol::AABase &current_aa( **res_it);
        const bool current_defined( current_aa.HasDefinedCoordinates());
        if( !current_defined)
        {
          // find the next residue with defined backbone coordinates
          for( auto res_find_it( res_it + 1); res_find_it != res_it_end; ++res_find_it)
          {
            if( ( **res_find_it).HasDefinedCoordinates())
            {
              def_c = *res_find_it;
              res_it = res_find_it + 1;
              break;
            }
          }

          // determine the parameters of the span
          int n_id( CHAIN.GetFirstAA()->GetSeqID());
          int c_id( CHAIN.GetLastAA()->GetSeqID());
          const char chain_id( CHAIN.GetChainID());
          if( def_n.IsDefined())
          {
            n_id = def_n->GetSeqID() + 1;
          }
          if( def_c.IsDefined())
          {
            c_id = def_c->GetSeqID() - 1;
          }

          if( !m_IgnoreTerms || ( n_id != CHAIN.GetFirstAA()->GetSeqID() && c_id != CHAIN.GetLastAA()->GetSeqID()))
          {
            const Span span( n_id, c_id, chain_id);
            spans.PushBack( span);
          }

          def_n = def_c;
          def_c.Reset();
        }
        else
        {
          def_n = *res_it;
        }
      }

      return spans;
    }

  } // namespace fold
} // namespace bcl
