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
#include "restraint/bcl_restraint_handler_epr_decay.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_hash_map.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> HandlerEPRDecay::s_Instance
    (
      util::Enumerated< HandlerBase< storage::Vector< EPRDecay> > >::AddInstance( new HandlerEPRDecay())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from defaults
    HandlerEPRDecay::HandlerEPRDecay()
    {
    }

    //! @brief Clone function
    //! @return pointer to new HandlerEPRDecay
    HandlerEPRDecay *HandlerEPRDecay::Clone() const
    {
      return new HandlerEPRDecay( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &HandlerEPRDecay::GetAlias() const
    {
      static const std::string s_name( "EPRDecay");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer HandlerEPRDecay::GetSerializer() const
    {
      io::Serializer serializer( HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads files listing the results of EPR decay measurements. Lines that begin with # or ! are automatically ignored."
      );
      return serializer;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerEPRDecay::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads atom distance restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    storage::Vector< EPRDecay> HandlerEPRDecay::ReadRestraints( std::istream &ISTREAM) const
    {
      // create all lines
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // create a mapping from the spin-labeling site to the EPR decay data
      storage::HashMap< std::string, EPRDecay> epr_decay_map;

      // iterate over all lines and add each valid restraint to the corresponding EPRDecay instance
      for( auto line_itr( lines.Begin()), line_itr_end( lines.End()); line_itr != line_itr_end; ++line_itr)
      {
        // skip empty lines or comment lines indicated by a leading # or !
        if( line_itr->empty() || ( *line_itr)[ 0] == '!' || ( *line_itr)[ 0] == '#')
        {
          continue;
        }

        // split the line and convert the values
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, " \t"));
        char chain_id_first, chain_id_second;
        int seq_id_first, seq_id_second;
        double time, decay;
        util::TryConvertFromString( chain_id_first, *split_line[ 0], util::GetLogger());
        util::TryConvertFromString( seq_id_first, *split_line[ 1], util::GetLogger());
        util::TryConvertFromString( chain_id_second, *split_line[ 2], util::GetLogger());
        util::TryConvertFromString( seq_id_second, *split_line[ 3], util::GetLogger());
        util::TryConvertFromString( time, *split_line[ 4], util::GetLogger());
        util::TryConvertFromString( decay, *split_line[ 5], util::GetLogger());

        // compute the unique key for this spin-labeling pair
        const std::string key( *split_line[ 0] + *split_line[ 1] + *split_line[ 2] + *split_line[ 3]);

        // if EPR decay data has already been added for this spin-labeling pair, add it to the same set
        if( epr_decay_map.Find( key) == epr_decay_map.End())
        {
          EPRDecay epr_decay_set( chain_id_first, seq_id_first, chain_id_second, seq_id_second);
          epr_decay_set.AddMeasurement( time, decay);
          epr_decay_map.Insert( storage::Pair< std::string, EPRDecay>( key, epr_decay_set));
        }
        else
        {
          epr_decay_map[ key].AddMeasurement( time, decay);
        }
      }

      // create the final list containing all EPR decay data sets
      storage::Vector< EPRDecay> decay_data;
      for( auto map_it( epr_decay_map.Begin()), map_it_end( epr_decay_map.End()); map_it != map_it_end; ++map_it)
      {
        decay_data.PushBack( map_it->second);
      }

      return decay_data;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace restraint
} // namespace bcl
