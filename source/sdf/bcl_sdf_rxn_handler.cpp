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
#include "sdf/bcl_sdf_rxn_handler.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "sdf/bcl_sdf_mdl_handler.h"
#include "sdf/bcl_sdf_mdl_header.h"
#include "sdf/bcl_sdf_mdl_line_types.h"
#include "sdf/bcl_sdf_mdl_property.h"
#include "storage/bcl_storage_vector.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    const size_t RXNHandler::s_NumberDescriptionLines = 3;

    //! @brief standard constructor
    RXNHandler::RXNHandler() :
      m_Description(),
      m_NumberReactants( 0),
      m_NumberProducts( 0),
      m_ReactantMolfiles(),
      m_ProductMolfiles(),
      m_IsValid( false)
    {
    }

    //! @brief constructor with initial input
    RXNHandler::RXNHandler( std::istream &ISTREAM) :
      m_Description(),
      m_NumberReactants( 0),
      m_NumberProducts( 0),
      m_ReactantMolfiles(),
      m_ProductMolfiles(),
      m_IsValid( false)
    {
      ReadFromRXN( ISTREAM);
    }

    //! @brief virtual copy constructor
    RXNHandler *RXNHandler::Clone() const
    {
      return new RXNHandler( *this);
    }

    //! @brief read from std::istream
    //! @note THIS READS TO THE END OF THE STREAM.  RXN files are designed to contain only one
    //!       reaction, which is an assumption made here.  All of the lines from the stream may not
    //!       be processed but the whole stream WILL be processed
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RXNHandler::ReadFromRXN( std::istream &ISTREAM)
    {
      // RXN files are only supposed to contain one reaction
      // assume ISTREAM contains only this and buffer everything appropriately

      // check that istream is good
      if( !ISTREAM.good() || ISTREAM.eof())
      {
        if( !ISTREAM.good())
        {
          BCL_MessageCrt( "Error reading RXN: passed bad stream");
        }
        return ISTREAM;
      }

      storage::List< std::string> rxn_buff;
      std::string buf;
      while( std::getline( ISTREAM, buf))
      {
        rxn_buff.PushBack( buf);
      }

      // read information from the buffer
      ReadFromRXN( rxn_buff.Begin(), rxn_buff.End());
      
      return ISTREAM;
    }

    std::ostream &RXNHandler::WriteToRXN
    (
      std::ostream &OSTREAM,
      const std::string &DESCRIPTION,
      const storage::Vector< MolfileHandler> &REACTANT_MOLFILES,
      const storage::Vector< MolfileHandler> &PRODUCT_MOLFILES
    )
    {
      size_t n_reactants( REACTANT_MOLFILES.GetSize());
      size_t n_products( REACTANT_MOLFILES.GetSize());
      if( n_reactants < 1000 && n_products < 1000)
      {
        std::stringstream oss;
        
        // Add the $RXN label and description
        oss << GetDefaultLine( e_RXNStartLine) << std::endl;
        oss << MolfileHandler::StandardizeDescription( DESCRIPTION) << std::endl;

        // Add the reaction header, i.e. (n reactants)(n products)
        std::string rxn_header( GetDefaultLine( e_RXNHeaderLine));
        GetMdlEntryTypes().RXNHeader_NumberReactantLines->Set( rxn_header, n_reactants);
        GetMdlEntryTypes().RXNHeader_NumberProductLines->Set( rxn_header, n_products);
        oss << rxn_header << std::endl;

        size_t mol_no( 0);
        for
        (
          storage::Vector< MolfileHandler>::const_iterator itr_mf( REACTANT_MOLFILES.Begin()),
            itr_mf_end( REACTANT_MOLFILES.End());
          itr_mf != itr_mf_end;
          ++itr_mf, ++mol_no
        )
        {
          if( !itr_mf->IsValid())
          {
            BCL_MessageCrt
            ( 
              "Cannot write RXN file: invalid molfile was provided for reactant " + util::Format()( mol_no)
              + " with description \"" + itr_mf->GetDescription() + "\""
            );
            return OSTREAM;
          }

          // write $MOL
          oss << GetDefaultLine( e_RXNMolStartLine) << std::endl;

          // write out the molfile
          itr_mf->WriteMolfile( oss);
        }

        mol_no = 0;
        for
        (
          storage::Vector< MolfileHandler>::const_iterator itr_mf( PRODUCT_MOLFILES.Begin()),
            itr_mf_end( PRODUCT_MOLFILES.End());
          itr_mf != itr_mf_end;
          ++itr_mf, ++mol_no
        )
        {
          if( !itr_mf->IsValid())
          {
            BCL_MessageCrt
            ( 
              "Cannot write RXN file: invalid molfile was provided for product " + util::Format()( mol_no)
              + " with description \"" + itr_mf->GetDescription() + "\""
            );
            return OSTREAM;
          }

          // write $MOL
          oss << GetDefaultLine( e_RXNMolStartLine) << std::endl;

          // write the molfile out
          itr_mf->WriteMolfile( oss);
        }

        OSTREAM << oss.str() << "\n";
      }
      else
      {
        BCL_MessageCrt
        ( 
          "Cannot write RXN file with >999 reactants or products (" + util::Format()( n_reactants)
          + " and " + util::Format()( n_products) + " were given, respectively)"
        );
        BCL_MessageCrt( "Unable to write RXN with description \"" + DESCRIPTION + "\"");
      }
      return OSTREAM;
    }

    std::ostream &RXNHandler::WriteToRXN
    (
      std::ostream &OSTREAM, 
      const std::string &DESCRIPTION,
      const storage::Vector< chemistry::FragmentComplete> &REACTANTS,
      const storage::Vector< chemistry::FragmentComplete> &PRODUCTS,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTANT_ATOM_MAP,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &PRODUCT_ATOM_MAP
    )
    {
      storage::Vector< std::string> reactant_descriptions;
      reactant_descriptions.AllocateMemory( REACTANTS.GetSize());

      storage::Vector< std::string> product_descriptions;
      product_descriptions.AllocateMemory( PRODUCTS.GetSize());

      storage::Vector< storage::Vector< AtomInfo> > reactant_atom_info;
      reactant_atom_info.AllocateMemory( REACTANTS.GetSize());
      storage::Vector< storage::Vector< BondInfo> > reactant_bond_info;
      reactant_bond_info.AllocateMemory( REACTANTS.GetSize());
      storage::Vector< storage::Vector< AtomInfo> > product_atom_info;
      product_atom_info.AllocateMemory( REACTANTS.GetSize());
      storage::Vector< storage::Vector< BondInfo> > product_bond_info;
      product_bond_info.AllocateMemory( REACTANTS.GetSize());

      for( size_t rno( 0), end_rno( REACTANTS.GetSize()); rno < end_rno; ++rno)
      {
        reactant_descriptions.PushBack( REACTANTS( rno).GetName());
        reactant_atom_info.PushBack( REACTANTS( rno).GetAtomInfo());
        reactant_bond_info.PushBack( REACTANTS( rno).GetBondInfo());
      }

      for( size_t pno( 0), end_pno( PRODUCTS.GetSize()); pno < end_pno; ++pno)
      {
        product_descriptions.PushBack( PRODUCTS( pno).GetName());
        product_atom_info.PushBack( PRODUCTS( pno).GetAtomInfo());
        product_bond_info.PushBack( PRODUCTS( pno).GetBondInfo());
      }

      return WriteToRXN
             ( 
               OSTREAM, DESCRIPTION, reactant_atom_info, 
               reactant_bond_info, product_atom_info, product_bond_info, 
               REACTANT_ATOM_MAP, PRODUCT_ATOM_MAP,
               reactant_descriptions, product_descriptions
             );
    } 

    //! @brief write a RXN-formatted reaction to an output stream 
    //! @param OSTREAM the output stream
    //! @param DESCRIPTION the description of the reaction
    //! @param REACTANTS the reactants to write
    //! @param PRODUCTS the products to write
    //! @param REACTANT_ATOM_MAP reactive atom map for reactants
    //! @param PRODUCT_ATOM_MAP reactive atom map for products
    //! @param TERMINATION_LINE the line to terminate the RXN with
    //! @return the output stream that was written
    std::ostream &RXNHandler::WriteToRXN
    (
      std::ostream &OSTREAM, 
      const std::string &DESCRIPTION,
      const storage::Vector< storage::Vector< AtomInfo> > &REACTANT_ATOM_INFO,
      const storage::Vector< storage::Vector< BondInfo> > &REACTANT_BOND_INFO,
      const storage::Vector< storage::Vector< AtomInfo> > &PRODUCT_ATOM_INFO,
      const storage::Vector< storage::Vector< BondInfo> > &PRODUCT_BOND_INFO,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &REACTANT_ATOM_MAP,
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &PRODUCT_ATOM_MAP,
      const storage::Vector< std::string> &REACTANT_DESCRIPTIONS,
      const storage::Vector< std::string> &PRODUCT_DESCRIPTIONS
    ) 
    {
      size_t n_reactants( REACTANT_ATOM_INFO.GetSize());
      size_t n_products( PRODUCT_ATOM_INFO.GetSize());
      
      if( REACTANT_BOND_INFO.GetSize() != n_reactants)
      {
        BCL_MessageStd
        ( 
          "Cannot write reaction, atom infos were specified for " + util::Format()( n_reactants) + " reactant(s), "
          "but bond infos were specified for " + util::Format()( REACTANT_BOND_INFO.GetSize()) + " reactant(s).  "
          "These should be equivalent"
        );
        return OSTREAM;
      }

      if( PRODUCT_BOND_INFO.GetSize() != n_reactants)
      {
        BCL_MessageStd
        ( 
          "Cannot write reaction, atom infos were specified for " + util::Format()( n_reactants) + " product(s), "
          "but bond infos were specified for " + util::Format()( PRODUCT_BOND_INFO.GetSize()) + " product(s).  "
          "These should be equivalent"
        );
        return OSTREAM;
      }
      
      if( n_reactants > 999 || n_products > 999)
      {
        BCL_MessageCrt
        ( 
          "Ignoring request to write a reaction with " + util::Format()( n_reactants) 
          + " reactants and " + util::Format()( n_products) + " products"
        );
        BCL_MessageCrt( "The MDL specification does not allow for reactions with more than 999 molecules");
        return OSTREAM;
      }

      // translate the reactive atom maps to something that maps atom index (key) to mapped value (value)
      // for each reactant/product (vector index)
      storage::Vector< storage::Map< size_t, size_t> > reactive_atoms_reactants( n_reactants);
      for
      (
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator itr_map( REACTANT_ATOM_MAP.Begin()),
          itr_map_end( REACTANT_ATOM_MAP.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->second.First() >= n_reactants)
        {
          BCL_MessageStd
          ( 
            "Cannot output reaction: malformed atom mapping.  Tried to add a map value of " + 
            util::Format()( itr_map->first) + " to reactant " + util::Format()( itr_map->second.First()) +
            " but maximum reactant index is " + util::Format()( n_reactants - 1)
          ); 
          return OSTREAM;
        }
        reactive_atoms_reactants( itr_map->second.First())[ itr_map->second.Second()] = itr_map->first;
      }

      storage::Vector< storage::Map< size_t, size_t> > reactive_atoms_products( n_products);
      for
      (
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator itr_map( PRODUCT_ATOM_MAP.Begin()),
          itr_map_end( PRODUCT_ATOM_MAP.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->second.First() >= n_products)
        {
          BCL_MessageStd
          ( 
            "Cannot output reaction: malformed atom mapping.  Tried to add a map value of " + 
            util::Format()( itr_map->first) + " to product " + util::Format()( itr_map->second.First()) +
            " but maximum product index is " + util::Format()( n_products - 1)
          ); 
          return OSTREAM;
        }
        reactive_atoms_products( itr_map->second.First())[ itr_map->second.Second()] = itr_map->first;
      }

      storage::Vector< MolfileHandler> reactant_molfiles;
      reactant_molfiles.AllocateMemory( n_reactants);
      storage::Vector< MolfileHandler> product_molfiles;
      product_molfiles.AllocateMemory( n_products);

      bool use_reactant_descriptions( REACTANT_DESCRIPTIONS.GetSize() == n_reactants);
      bool use_product_descriptions( PRODUCT_DESCRIPTIONS.GetSize() == n_products);

      // generate reactant molfile handlers
      for( size_t rno( 0); rno < n_reactants; ++rno)
      {
        std::string desc;
        if( use_reactant_descriptions)
        {
          desc = REACTANT_DESCRIPTIONS( rno);
        }
        else
        {
          desc = std::string( "Reactant " + util::Format()( rno) + "\n\n\n");
        }
        desc = MolfileHandler::StandardizeDescription( desc);
        reactant_molfiles.PushBack( MolfileHandler( desc, REACTANT_ATOM_INFO( rno), REACTANT_BOND_INFO( rno)));

        MolfileHandler &r_molfile( reactant_molfiles.LastElement());

        // set the atom mapping
        for
        (
          storage::Map< size_t, size_t>::const_iterator itr_map( reactive_atoms_reactants( rno).Begin()), 
            itr_map_end( reactive_atoms_reactants( rno).End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          r_molfile.SetAtomMapping( itr_map->first, itr_map->second);
        }

        if( !r_molfile.IsValid())
        {
          BCL_MessageStd( "Cannot write reaction, reactant " + util::Format()( rno) + " could not be formatted properly");
          return OSTREAM;
        }
      }

      // generate product molfile handlers
      for( size_t pno( 0); pno < n_products; ++pno)
      {
        std::string desc;
        if( use_product_descriptions)
        {
          desc = PRODUCT_DESCRIPTIONS( pno);
        }
        else
        {
          desc = std::string( "Product " + util::Format()( pno) + "\n\n\n");
        }
        desc = StandardizeDescription( desc);
        product_molfiles.PushBack( MolfileHandler( desc, PRODUCT_ATOM_INFO( pno), PRODUCT_BOND_INFO( pno)));

        MolfileHandler &p_molfile( product_molfiles.LastElement());

        // set the atom mapping
        for
        (
          storage::Map< size_t, size_t>::const_iterator itr_map( reactive_atoms_products( pno).Begin()), 
            itr_map_end( reactive_atoms_products( pno).End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          p_molfile.SetAtomMapping( itr_map->first, itr_map->second);
        }

        if( !p_molfile.IsValid())
        {
          BCL_MessageStd( "Cannot write reaction, product " + util::Format()( pno) + " could not be formatted properly");
          return OSTREAM;
        }
      }

      return WriteToRXN( OSTREAM, DESCRIPTION, reactant_molfiles, product_molfiles);
    }

    //! @brief gets the class name
    //! @return the class name
    const std::string &RXNHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get all description lines stored in handler
    //! @return list of mdl description lines stored in handler
    const std::string &RXNHandler::GetDescription() const
    {
      return m_Description;
    }

    //! @brief check if the handler is valid
    //! @return true if the mdl handler is in a valid state
    bool RXNHandler::IsValid() const
    {
      return m_IsValid;
    }

    //! @brief Get the number of Reactants
    const size_t &RXNHandler::GetNumberReactants() const
    {
      return m_NumberReactants;
    }

    //! @brief Get the number of Reactants
    const size_t &RXNHandler::GetNumberProducts() const
    {
      return m_NumberProducts;
    }

    //! @brief Get reacting atom mapping on the reactants
    const storage::Map< size_t, storage::Pair< size_t, size_t> > &RXNHandler::GetReactiveAtomsInReactants() const
    {
      return m_ReactiveAtomsReactants;
    }

    //! @brief Get reacting atom mapping on the products
    const storage::Map< size_t, storage::Pair< size_t, size_t> > &RXNHandler::GetReactiveAtomsInProducts() const
    {
      return m_ReactiveAtomsProducts;
    }

    //! @brief read RXN from a file buffer as iterators to strings
    //! @param LINE_BEGIN iterator to first line to begin reading (should be "$RXN" but may be blank/whitespace only)
    //! @param LINE_END iterator to one past the last line that should be read
    //! @return iterator to the first line that was not read.  If all were read, should be LINE_END
    storage::List< std::string>::const_iterator RXNHandler::ReadFromRXN
    (
      const storage::List< std::string>::const_iterator &LINE_BEGIN,
      const storage::List< std::string>::const_iterator &LINE_END
    )
    {
      // reset all members
      m_Description.erase();
      m_NumberReactants = m_NumberProducts = 0;
      m_ReactiveAtomsReactants.Reset();
      m_ReactiveAtomsProducts.Reset();
      m_IsValid = false;

      if( LINE_BEGIN == LINE_END)
      {
        return LINE_BEGIN;
      }

      storage::List< std::string>::const_iterator itr_line( LINE_BEGIN);

      // Skip lines until the header line is reached, then skip the header line
      size_t n_skipped_lines( 0);
      for( ; itr_line != LINE_END && !IsRXNDelimiter( *itr_line); ++n_skipped_lines, ++itr_line);
      ++itr_line; // move past $RXN line

      if( n_skipped_lines)
      {
        BCL_MessageCrt
        ( 
          "Warning: there were " + util::Format()( n_skipped_lines) + 
          " blank lines before $RXN.  This is non-standard (should be 0); continuing anyway"
        );
      }

      // Read 3 lines in as the description
      for( size_t d( 0); d < s_NumberDescriptionLines && itr_line != LINE_END; ++d, ++itr_line)
      {
        m_Description += *itr_line + "\n";
      }
      m_Description = MolfileHandler::StandardizeDescription( m_Description);

      // Next two lines should be the header line and the first $MOL line
      if
      ( 
        !GetMdlEntryTypes().RXNHeader_NumberReactantLines->IsUnsignedInt( *itr_line) ||
        !GetMdlEntryTypes().RXNHeader_NumberProductLines->IsUnsignedInt( *itr_line)
      )
      {
        BCL_MessageCrt
        ( 
          "Line \"" + *itr_line + "\" should have been an RXN header line, but it was not parseable as one.  "
          "Not reading remaining information for RXN with description \"" + m_Description + "\""
        );
        return ++itr_line;
      }

      // read in reactant and product information
      if( itr_line != LINE_END)
      {
        m_NumberReactants = GetMdlEntryTypes().RXNHeader_NumberReactantLines->GetUnsignedInt( *itr_line);
        m_NumberProducts = GetMdlEntryTypes().RXNHeader_NumberProductLines->GetUnsignedInt( *itr_line);
      }
      else
      {
        BCL_MessageCrt( "Unexpected end of RXN file: no header was found for RXN with description \"" + m_Description + "\"");
        return itr_line;
      }
      ++itr_line;

      // if there are no reactants or products then we have technically read in a legitimate RXN file
      if( !m_NumberReactants && !m_NumberProducts)
      {
        BCL_MessageStd( "Note: RXN with description " + m_Description + " contains no molecules");
        m_IsValid = true;
        return itr_line;
      }

      while( itr_line != LINE_END && m_ReactantMolfiles.GetSize() < m_NumberReactants)
      {
        // find the $MOL
        // TODO more error checking here
        n_skipped_lines = 0;
        for( ; itr_line != LINE_END && !IsRXNMolDelimiter( *itr_line); ++n_skipped_lines, ++itr_line);
        ++itr_line; // skip past the $MOL
        
        // if we're at the end then bail out
        if( itr_line == LINE_END)
        {
          BCL_MessageCrt( "Unexpected end of RXN: reactants were not all read");
          return itr_line;
        }
        
        // post a warning about skipped lines, files really shouldn't have any blanks but we'll tolerate it
        if( n_skipped_lines)
        {
          BCL_MessageStd( "Warning: reactant block of RXN: there were " + util::Format()( n_skipped_lines) + " blank lines before $MOL, this is non-standard");
        }

        // read the reactant molfile
        m_ReactantMolfiles.PushBack( MolfileHandler());
        itr_line = m_ReactantMolfiles.LastElement().ReadMolfile( itr_line, LINE_END);

        // if there was a problem, stop reading
        if( !m_ReactantMolfiles.LastElement().IsValid())
        {
          BCL_MessageCrt( "Error reading RXN file: could not read reactant #" + util::Format()( m_ReactantMolfiles.GetSize()) + "; exitting");
          return itr_line;
        }
      }

      // make sure we read all reactants
      if( m_ReactantMolfiles.GetSize() != m_NumberReactants)
      {
        BCL_MessageCrt
        ( 
          "Unexpected end of RXN file: could not read the required number of reactants.  "
          "Read " + util::Format()( m_ReactantMolfiles.GetSize()) + " but file specified " + util::Format()( m_NumberReactants)
        );
        return itr_line;
      }

      // now do the same thing for products
      while( itr_line != LINE_END && m_ProductMolfiles.GetSize() < m_NumberProducts)
      {
        // find the $MOL
        // TODO more error checking here
        n_skipped_lines = 0;
        for( ; itr_line != LINE_END && !IsRXNMolDelimiter( *itr_line); ++n_skipped_lines, ++itr_line);
        ++itr_line; // skip past the $MOL
        
        // if we're at the end then bail out
        if( itr_line == LINE_END)
        {
          BCL_MessageCrt( "Unexpected end of RXN: products were not all read");
          return itr_line;
        }
        
        // post a warning about skipped lines, files really shouldn't have any blanks but we'll tolerate it
        if( n_skipped_lines)
        {
          BCL_MessageStd( "Warning: product block in RXN: there were " + util::Format()( n_skipped_lines) + " blank lines before $MOL, this is non-standard");
        }

        // read the reactant molfile
        m_ProductMolfiles.PushBack( MolfileHandler());
        itr_line = m_ProductMolfiles.LastElement().ReadMolfile( itr_line, LINE_END);

        // if there was a problem, stop reading
        if( !m_ProductMolfiles.LastElement().IsValid())
        {
          BCL_MessageCrt( "Error reading RXN file: could not read product #" + util::Format()( m_ProductMolfiles.GetSize()) + "; exitting");
          return itr_line;
        }
      }

      // make sure we read all reactants
      if( m_ProductMolfiles.GetSize() != m_NumberProducts)
      {
        BCL_MessageCrt
        ( 
          "Unexpected end of RXN file: could not read the required number of products.  "
          "Read " + util::Format()( m_ReactantMolfiles.GetSize()) + " but file specified " + util::Format()( m_NumberProducts)
        );
        return itr_line;
      }

      // everything is good, itr_line should now point to either LINE_END or whatever one line after the last 'M  END' line is
      // finish off the parsing and set up necessary internal structures
      m_IsValid = ValidateInput(); 

      return itr_line;
    }

    //! @brief finalize parsing; translates m_Stream into internal members, sets m_Parsed to true, cleans m_Stream
    bool RXNHandler::ValidateInput() const
    {
      bool did_succeed( false);

      // set up mappings for reactants
      for( size_t r_no( 0); r_no < m_NumberReactants; ++r_no)
      {
        storage::Vector< size_t> r_mapping( m_ReactantMolfiles( r_no).GetAtomMapping());
        if( r_mapping.GetSize() != m_ReactantMolfiles( r_no).GetAtomInfo().GetSize())
        {
          BCL_MessageCrt( "Error parsing RXN file: reactant " + util::Format()( r_no) + " atoms did not have a proper mapping");
          return did_succeed;
        }
        
        // if atom mappings are non-zero then add them to the reactant atom map
        for( size_t a_no( 0), end_a( r_mapping.GetSize()); a_no < end_a; ++a_no)
        {
          const size_t &map_val( r_mapping( a_no));

          // if the mapping value is 0 then skip it
          if( !map_val)
          {
            continue;
          }

          if( m_ReactiveAtomsReactants.Has( map_val))
          {
            const storage::Pair< size_t, size_t> &old_val( m_ReactiveAtomsReactants.GetValue( map_val));
            BCL_MessageCrt( "Error parsing RXN file: multiple atoms were mapped to value " + util::Format()( map_val));
            BCL_MessageCrt
            ( 
              "  Reactant " + util::Format()( r_no) + " atom index " + util::Format()( a_no) + " and reactant "
              + util::Format()( old_val.First()) + " atom index " + util::Format()( old_val.Second()) + " both "
              "contained this value"
            );
            return did_succeed;
          }

          m_ReactiveAtomsReactants[ map_val].First() = r_no;
          m_ReactiveAtomsReactants[ map_val].Second() = a_no;
        }
      }

      // set up mappings for products
      for( size_t p_no( 0); p_no < m_NumberProducts; ++p_no)
      {
        storage::Vector< size_t> p_mapping( m_ProductMolfiles( p_no).GetAtomMapping());
        if( p_mapping.GetSize() != m_ProductMolfiles( p_no).GetAtomInfo().GetSize())
        {
          BCL_MessageCrt( "Error parsing RXN file: product " + util::Format()( p_no) + " atoms did not have a proper mapping");
          return did_succeed;
        }
        
        // if atom mappings are non-zero then add them to the product atom map
        for( size_t a_no( 0), end_a( p_mapping.GetSize()); a_no < end_a; ++a_no)
        {
          const size_t &map_val( p_mapping( a_no));

          // if the mapping value is 0 then skip it
          if( !map_val)
          {
            continue;
          }

          if( m_ReactiveAtomsProducts.Has( map_val))
          {
            const storage::Pair< size_t, size_t> &old_val( m_ReactiveAtomsProducts.GetValue( map_val));
            BCL_MessageCrt( "Error parsing RXN file: multiple atoms were mapped to value " + util::Format()( map_val));
            BCL_MessageCrt
            ( 
              "  Product " + util::Format()( p_no) + " atom index " + util::Format()( a_no) + " and product "
              + util::Format()( old_val.First()) + " atom index " + util::Format()( old_val.Second()) + " both "
              "contained this value"
            );
            return did_succeed;
          }

          m_ReactiveAtomsProducts[ map_val].First() = p_no;
          m_ReactiveAtomsProducts[ map_val].Second() = a_no;
        }
      }

      did_succeed = true;
      return did_succeed;
    }

    //! @brief read RXNHandler object from std::istream
    //! @param ISTREAM istream that contains RXNHandler object
    //! @return istream after RXNHandler object was extracted
    std::istream &RXNHandler::Read( std::istream &ISTREAM)
    {
      return ReadFromRXN( ISTREAM);
    }

    //! @brief write RXNHandler into std::ostream
    //! @param OSTREAM ostream that gets RXNHandler object
    //! @param INDENT indentation
    //! @return ostream after RXNHandler object was inserted
    std::ostream &RXNHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the description; checks # of newlines, ensures that it is exactly s_NumberDescriptionLines
    //! If the # of newlines is < s_NumberDescriptionLines - 1, adds new lines as necessary
    //! Any newlines >= s_NumberDescriptionLines are replaced with spaces
    std::string RXNHandler::StandardizeDescription( const std::string &DESCRIPTION)
    {
      const size_t last_non_space( DESCRIPTION.find_last_not_of( " \n\t\r"));
      const size_t description_size( last_non_space + 1);
      std::string description;
      description.reserve( description_size);

      // keep track of the number of new lines seen so far
      size_t number_new_lines( 0);

      for( size_t i( 0); i < description_size; ++i)
      {
        if( DESCRIPTION[ i] == '\n' && ++number_new_lines >= s_NumberDescriptionLines)
        {
          description += ' ';
          continue;
        }
        // skip carriage returns (occur when reading sdfs from windows on non-windows machine)
        else if( DESCRIPTION[ i] == '\r')
        {
          continue;
        }
        description += DESCRIPTION[ i];
      }
      // add new lines until there are s_NumberDescriptionLines - 1 of them
      while( ++number_new_lines < s_NumberDescriptionLines)
      {
        description += '\n';
      }
      return description;
    }

    //! @brief test whether a line contains only spaces or is otherwise empty
    //! @param STRING the string to test
    bool RXNHandler::ContainsNonspaceCharacters( const std::string &STRING)
    {
      for( std::string::const_iterator itr( STRING.begin()), itr_end( STRING.end()); itr != itr_end; ++itr)
      {
        if( !isspace( *itr))
        {
          return true;
        }
      }
      return false;
    }

    //! @brief write mdl lines into std::ostream
    //! @param OSTREAM ostream that gets RXNHandler object
    void RXNHandler::WriteMiscProperties
    (
      std::ostream &OSTREAM,
      const storage::Map< std::string, std::string> &MISC_PROPERTIES
    )
    {
      static const std::string pre_property_name_str( "> <");  // goes before each property's name
      static const std::string post_property_name_str( ">\n"); // goes after each property's name

      // iterate over all MDL misc property lines
      for
      (
        storage::Map< std::string, std::string>::const_iterator
          itr_map( MISC_PROPERTIES.Begin()),
          itr_map_end( MISC_PROPERTIES.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        if( itr_map->first.size() != 0 && itr_map->second.size() != 0) // the property has a name and value
        {
          OSTREAM << pre_property_name_str        // property name start delimiter
                  << itr_map->first               // property name
                  << post_property_name_str;      // property name deliminater

          const std::string &value( itr_map->second);

          // write the value of the given property line
          OSTREAM << value;

          // if the last character of the string was not a new-line, add one here
          if( value[ value.size() - 1] != '\n')
          {
            OSTREAM << '\n'; // followed by a newline
          }
          OSTREAM << '\n';     // put a blank line follows the last line of the property value
        }
      }
    }

    //! @brief check if a line is the terminating line
    //! @return true if the line just contains $RXN
    bool RXNHandler::IsRXNDelimiter( const std::string &LINE)
    {
      static const std::string s_default_terminal_line( GetDefaultLine( e_RXNStartLine));
      return util::StartsWith( LINE, s_default_terminal_line);
    }

    //! @brief check if a line is the terminating line
    //! @return true if the line just contains $MOL
    bool RXNHandler::IsRXNMolDelimiter( const std::string &LINE)
    {
      static const std::string s_default_terminal_line( GetDefaultLine( e_RXNMolStartLine));
      return util::StartsWith( LINE, s_default_terminal_line);
    }

    //! @brief return the datalable, if line is a datalabel line
    //! @param LINE line from mdl section
    //! @return string that constains datalable, string will be empty for non-data lable lines
    std::string RXNHandler::GetMDLDataLabel( const std::string &LINE)
    {
      // data label delimiter left
      static const char s_misc_property_delimiter_left( '<');

      // data label delimiter right
      static const char s_misc_property_delimiter_right( '>');

      // handle empty lines and lines that do not start with <
      if( LINE.empty() || LINE[ 0] != s_misc_property_delimiter_right)
      {
        return std::string();
      }

      // find label start and end
      const std::string::size_type label_start( LINE.find( s_misc_property_delimiter_left, 1));
      const std::string::size_type label_end( LINE.rfind( s_misc_property_delimiter_right));

      // return an empty string for non data label line
      if( label_start == std::string::npos || label_end == std::string::npos)
      {
        return std::string();
      }

      // misc property name
      return LINE.substr( label_start + 1, label_end - label_start - 1);
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RXNHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new RXNHandler())
    );

  } // namespace sdf
} // namespace bcl

