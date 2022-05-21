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
#include "sdf/bcl_sdf_ctab_handler.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "sdf/bcl_sdf_mdl_header.h"
#include "sdf/bcl_sdf_mdl_line_types.h"
#include "sdf/bcl_sdf_mdl_property.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief standard constructor
    CTabHandler::CTabHandler() :
      m_IsValid( false),
      m_Header(),
      m_AtomInfos(),
      m_BondInfos(),
      m_AtomMapping(),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
    }

    //! @brief constructor from an input stream
    CTabHandler::CTabHandler( std::istream &ISTREAM) :
      m_IsValid( false),
      m_Header(),
      m_AtomInfos(),
      m_BondInfos(),
      m_AtomMapping(),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
      ReadCTab( ISTREAM);
    }

    //! @brief constructor from a pre-read set of lines
    CTabHandler::CTabHandler( const storage::List< std::string> &LINES) :
      m_IsValid( false),
      m_Header(),
      m_AtomInfos(),
      m_BondInfos(),
      m_AtomMapping(),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
      ReadCTab( LINES.Begin(), LINES.End());
    }

    //! @brief constructor from a list of AtomInfo and BondInfos
    CTabHandler::CTabHandler
    (
      const storage::Vector< AtomInfo> &ATOM_INFOS,
      const storage::Vector< BondInfo> &BOND_INFOS
    ) :
      m_IsValid( true),
      m_Header( ATOM_INFOS.GetSize(), BOND_INFOS.GetSize()),
      m_AtomInfos( ATOM_INFOS),
      m_BondInfos( BOND_INFOS),
      m_AtomMapping( ATOM_INFOS.GetSize(), size_t( 0)),
      m_AtomTypesRead( false),
      m_BondTypesRead( false),
      m_DoubleBondIsometryWasRead( false),
      m_ChiralityWasRead( false)
    {
      size_t n_atoms( ATOM_INFOS.GetSize());

      // check sanity of the bond data
      for
      (
        storage::Vector< BondInfo>::const_iterator itr_bond( m_BondInfos.Begin()),
          itr_bond_end( m_BondInfos.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        if( itr_bond->GetAtomIndexLow() >= n_atoms || itr_bond->GetAtomIndexHigh() >= n_atoms)
        {
          BCL_MessageCrt
          (
            "Could not make a connection table, invalid data given.  Bond tried to join atoms "
            + util::Format()( itr_bond->GetAtomIndexLow()) + " and " + util::Format()( itr_bond->GetAtomIndexHigh()) +
            " but maximum index is " + util::Format()( n_atoms - 1)
          );
          m_AtomInfos.Reset();
          m_BondInfos.Reset();
          m_AtomMapping.Reset(),
          m_Header = MdlHeader();
          m_IsValid = false;
          break;
        }
      }
    }

    //! @brief virtual copy constructor
    CTabHandler *CTabHandler::Clone() const
    {
      return new CTabHandler( *this);
    }

    //! @brief Read the CTab using iterators to strings
    //! @param LINE_BEGIN a line that represents a header/counts line
    //! @param LINE_END one-past-end of possible lines
    //! @details if the CTab ends before LINE_END, not all lines will be read
    //! @return an iterator to the first line that was not read
    storage::List< std::string>::const_iterator CTabHandler::ReadCTab
    (
      const storage::List< std::string>::const_iterator &LINE_BEGIN,
      const storage::List< std::string>::const_iterator &LINE_END
    )
    {

      Reset();

      storage::List< std::string>::const_iterator itr( LINE_BEGIN);
      if( LINE_BEGIN == LINE_END)
      {
        BCL_MessageStd( "CTabHandler: nothing to read");
        return LINE_BEGIN;
      }

      // Read the header, should be the first line
      m_Header.SetFromMdlLine( *itr, 0);

      if( !m_Header.IsValid())
      {
        BCL_MessageStd( "CTab header \"" + *itr + "\" is invalid, not reading any further");
        return ++itr;
      }
      ++itr;

      const size_t n_atoms( m_Header.GetNumberAtoms());
      const size_t n_bonds( m_Header.GetNumberBonds());

      m_AtomInfos.AllocateMemory( n_atoms);
      m_BondInfos.AllocateMemory( n_bonds);
      m_AtomMapping.AllocateMemory( n_atoms);

      // read atom information
      size_t found_atoms = 0;
      for( ; itr != LINE_END && found_atoms < n_atoms; ++itr)
      {
        if( !ContainsNonspaceCharacters( *itr))
        {
          BCL_MessageCrt( "Should not find an empty line in the atom line section! - skipping line");
          continue;
        }

        if( !AtomInfo::FormattedAsMdlAtomLine( *itr))
        {
          BCL_MessageCrt( "Line \"" + *itr + "\" was not formatted as an MDL atom line; not reading any further in this molecule");
          break;
        }

        // get the basic atom information and mapping
        m_AtomInfos.PushBack( AtomInfo().ExtractMdlAtomLineInfo( *itr));
        if( GetMdlEntryTypes().Atom_AtomMappingNumber->IsUnsignedInt( *itr))
        {
          m_AtomMapping.PushBack( GetMdlEntryTypes().Atom_AtomMappingNumber->GetUnsignedInt( *itr));
        }
        else
        {
          m_AtomMapping.PushBack( 0);
        }
        ++found_atoms;
      }

      // check we read enough atoms
      if( found_atoms != n_atoms)
      {
        BCL_MessageCrt
        (
          "Did not find the appropriate number of atoms!  The MDL header line specified " +
          util::Format()( n_atoms) + " atom lines but only found " + util::Format()( found_atoms)
        );
        return itr; // this will be one after the last atom line that was read
      }

      // read bonds
      size_t found_bonds = 0;
      for( ; itr != LINE_END && found_bonds < n_bonds; ++itr)
      {
        if( !ContainsNonspaceCharacters( *itr))
        {
          BCL_MessageCrt( "Should not find an empty line in the bond line section! - skipping line");
          continue;
        }

        // check that formatting is correct
        if( !BondInfo::FormattedAsMdlBondLine( *itr))
        {
          BCL_MessageCrt( "Line \"" + *itr + "\" was not formatted as an MDL bond line; not reading any further in this molecule");
          break;
        }

        // create the bond line
        BondInfo bond;
        bond.ExtractMdlBondLineInfo( *itr);

        // check the high atom.  If < n_atoms we're fine, if > n_atoms then at least this one is incorrect
        if( bond.GetAtomIndexHigh() >= n_atoms)
        {
          BCL_MessageCrt
          (
            "Warning: line \"" + *itr + "\" references an invalid atom (" +
            util::Format()( bond.GetAtomIndexHigh() + 1) + ")."
            "  Maximum atom number for this molecule is " + util::Format()( n_atoms) + ". Skipping this bond"
          );
          continue;
        }

        // check that the bond type is recognized, otherwise skip it
        if( !bond.GetConfigurationalBondType().IsDefined())
        {
          BCL_MessageCrt
          (
            "Warning: bond line \"" + *itr + "\" contains an unrecognized bond order; skipping this bond"
          );

          // increment this because we still read everything properly, we just don't know what to do with it
          ++found_bonds;
          continue;
        }
        m_BondInfos.PushBack( bond);
        ++found_bonds;
      }

      // check that enough bonds were read
      if( found_bonds != n_bonds)
      {
        BCL_MessageCrt
        (
          "Did not find the appropriate number of bonds!  The MDL header specified " +
          util::Format()( n_bonds) + " bond lines but only found " + util::Format()( found_bonds)
        );
        return itr;
      }

      // read in MDL properties
      for( ; itr != LINE_END; ++itr)
      {
        MdlProperty prop( *itr);

        //if( !MdlProperty::IsBCLPropertyLine( *itr))
        if( !prop.IsValid())
        {
          BCL_MessageStd( "Unrecognized MDL property line \"" + *itr + "\"; skipping line.");
          continue;
        }

        // found the end of the block, so everything was read correctly
        if( prop.GetProperty() == MdlProperty::e_BlockTerminator)
        {
          ++itr; // go to the next line, i.e. last unread line
          m_IsValid = true; // we are only valid if we actually read a terminal line
          break;
        }

        // some, but not all, properties need to be stored with the associated ctab.  check for this
        if( prop.ShouldCache())
        {
          m_MDLProperties[ prop.GetLabel()] = prop.GetPropertyStrings();
        }

        ApplyMDLProperty( prop);
      }
      return itr;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CTabHandler::ReadCTab( std::istream &ISTREAM)
    {
      // buffer for the stream
      storage::List< std::string> lines;

      // check that istream is good
      if( !ISTREAM.good() || ISTREAM.eof())
      {
        BCL_MessageCrt( "Error reading CTab: passed bad istream");
        return ISTREAM;
      }

      while( !ISTREAM.eof())
      {
        lines.PushBack( std::string());
        std::string &last_line( lines.LastElement());
        std::getline( ISTREAM, lines.LastElement());
        if( IsTerminalLine( last_line))
        {
          ReadCTab( lines.Begin(), lines.End());
          break;
        }
      }
      return ISTREAM;
    }

    //! @brief write to std::ostream in mdl format
    //! @param OSTREAM the stream to write to
    //! @param FORCE_WRITE_ATOM_TYPES whether to force writing BCL atom types (and other info)
    std::ostream &CTabHandler::WriteCTab
    (
      std::ostream &OSTREAM,
      const bool &FORCE_WRITE_ATOM_TYPES
    ) const
    {

      // don't write anything if our CTab isn't valid
      if( !IsValid())
      {
        return OSTREAM;
      }

      // write to a stringstream first and output at the end of writing
      // if anything bad happens (can't convert int etc) then we don't want to write a fragmented file
      std::stringstream tmp_ostream;

      // v2000 format allots 3 characters for # of atoms / bonds, so it cannot handle more than 999 atoms or bonds
      if( m_AtomInfos.GetSize() > size_t( 999) || m_BondInfos.GetSize() > size_t( 999))
      {
        BCL_MessageCrt
        (
          "Ignoring request to write a connection table (CTab);  "
          "CTab contains " + util::Format()( m_AtomInfos.GetSize()) + " atoms and "
          + util::Format()( m_BondInfos.GetSize()) + " bonds, but the V2000 format can only"
          "support up to 999 of each"
        );
        return OSTREAM;
      }
      //std::ostringstream tmp_ostream;

      // write out the header
      tmp_ostream << m_Header.ToMdlLine() << '\n';

      // write out all atom info
      size_t atom_no( 0);
      for
      (
        storage::Vector< AtomInfo>::const_iterator
          itr_atom( m_AtomInfos.Begin()), itr_atom_end( m_AtomInfos.End());
        itr_atom != itr_atom_end;
        ++itr_atom, ++atom_no
      )
      {
        // basic atom information
        std::string line( itr_atom->ToMdlAtomLine());

        // add mapping
        GetMdlEntryTypes().Atom_AtomMappingNumber->Set( line, m_AtomMapping( atom_no));

        tmp_ostream << line << '\n';
      }

      // write out all bond lines
      for
      (
        storage::Vector< BondInfo>::const_iterator
          itr_bond( m_BondInfos.Begin()), itr_bond_end( m_BondInfos.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        tmp_ostream << itr_bond->ToMdlBondLine() << '\n';
      }

      // determine which MDL properties need to be added
      storage::Vector< MdlProperty::PropertyEnum> properties_to_add( GetNecessaryMDLProperties( FORCE_WRITE_ATOM_TYPES));

      for
      (
        storage::Vector< MdlProperty::PropertyEnum>::const_iterator
          itr_prop( properties_to_add.Begin()), itr_prop_end( properties_to_add.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        // create the string from the property
        const std::string property_str( MdlProperty( *itr_prop, m_AtomInfos, m_BondInfos).GetString());

        if( !property_str.empty())
        {
          // write the given property
          tmp_ostream << property_str << '\n';
        }
        else
        {
          BCL_MessageVrb( "Warning: Property \"" + MdlProperty::GetPropertyName( *itr_prop) + "\" was needed but no info was written");
        }
      }

      // write properties block terminator
      tmp_ostream << MdlProperty( MdlProperty::e_BlockTerminator).GetString() << "\n";

      OSTREAM << tmp_ostream.str();
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if a line is the terminating line
    //! @return true if the line just contains 'M  END'
    bool CTabHandler::IsTerminalLine( const std::string &LINE)
    {
      return LINE == MdlProperty( MdlProperty::e_BlockTerminator).GetString();
    }

    //! @brief test whether a line contains only spaces or is otherwise empty
    //! @param STRING the string to test
    bool CTabHandler::ContainsNonspaceCharacters( const std::string &STRING)
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

    //! @brief apply instructions from an MDL property to this class
    //! @param PROP the property to apply
    void CTabHandler::ApplyMDLProperty( const MdlProperty &PROP)
    {
      if( !PROP.IsValid())
      {
        BCL_MessageStd( "Warning: Property \"" + PROP.GetString() + "\" is not a valid property, not reading it");
      }
      if( PROP.GetProperty() < MdlProperty::s_NumberProperties && PROP.IsValid())
      {
        PROP.ApplyProperty( m_AtomInfos, m_BondInfos);
        if( PROP.GetProperty() == MdlProperty::e_BclAtomType)
        {
          m_AtomTypesRead = true;
        }
        else if( PROP.GetProperty() == MdlProperty::e_BclBondType)
        {
          m_BondTypesRead = true;
        }
        else if( PROP.GetProperty() == MdlProperty::e_BclChirality)
        {
          m_ChiralityWasRead = true;
        }
        else if( PROP.GetProperty() == MdlProperty::e_BclDoubleBondIsometry)
        {
          m_DoubleBondIsometryWasRead = true;
        }
      }
    }

    //! @brief determine the MDL properties that must be used when writing CTab information
    //! @return a vector of MDL PropertyEnums that denote which properties must be computed
    storage::Vector< MdlProperty::PropertyEnum> CTabHandler::GetNecessaryMDLProperties( const bool &FORCE_ATOM_TYPES) const
    {
      storage::Vector< MdlProperty::PropertyEnum> properties_to_add;

      // determine whether more than 1 coordinate is undefined / 0, in which case all
      // isometry information must be written as well
      // If the positions are all defined, and no more than 1 atom is at the origin, then
      // the isometry can be determined later on
      bool positions_invalid( false);
      {
        size_t zero_position_count( 0);
        for
        (
          storage::Vector< AtomInfo>::const_iterator
            itr_atom( m_AtomInfos.Begin()), itr_atom_end( m_AtomInfos.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          // keep track of whether the coordinates are all defined and that no more than one of them is at the origin
          if( !itr_atom->GetCoordinates().IsDefined())
          {
            positions_invalid = true;
            break;
          }
          else if( math::EqualWithinAbsoluteTolerance( itr_atom->GetCoordinates().Norm(), double( 0.0), 1.0e-4))
          {
            if( ++zero_position_count == 2)
            {
              positions_invalid = true;
              break;
            }
          }
        }
      }

      // check whether charges must be written out as a separate property
      bool charges_outside_of_mdl_range( false);
      {
        for
        (
          storage::Vector< AtomInfo>::const_iterator
            itr_atom( m_AtomInfos.Begin()), itr_atom_end( m_AtomInfos.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          // if the absolute value of a charge is greater than 3 it must be put in the "M  CHG' field
          if( std::abs( itr_atom->GetAtomType()->GetFormalCharge()) > short( 3))
          {
            charges_outside_of_mdl_range = true;
            break;
          }
        }
      }

      // if valences are present, the bcl atom type must be written out, otherwise
      // the atom type may be determined incorrectly when the file is reloaded
      if( FORCE_ATOM_TYPES || TestMustWriteAtomTypes( m_AtomInfos, m_BondInfos))
      {
        properties_to_add.PushBack( MdlProperty::e_BclAtomType);
        properties_to_add.PushBack( MdlProperty::e_BclBondType);
      }

      // if positions were invalid, or atom types are being written, isometry info must be written as well
      if( positions_invalid || properties_to_add.GetSize())
      {
        properties_to_add.PushBack( MdlProperty::e_BclChirality);
        properties_to_add.PushBack( MdlProperty::e_BclDoubleBondIsometry);
      }

      // if any charges were outside the normal mdl range, add them as a property
      if( charges_outside_of_mdl_range)
      {
        properties_to_add.PushBack( MdlProperty::e_Charge);
      }

      return properties_to_add;
    }

    //! @brief from a set of atoms and bonds, determine whether atom types must be written
    //! @param ATOMS vector of atom info
    //! @param BONDS vector of bond info
    //! @return true if atom types must be written to ensure that they will be the same upon rereading the molecule
    bool CTabHandler::TestMustWriteAtomTypes
    (
      const storage::Vector< AtomInfo> &ATOMS,
      const storage::Vector< BondInfo> &BONDS
    )
    {
      storage::Vector< size_t> bond_counts( ATOMS.GetSize(), size_t( 0));
      // add up the bond counts for each atom
      for( storage::Vector< BondInfo>::const_iterator itr( BONDS.Begin()), itr_end( BONDS.End()); itr != itr_end; ++itr)
      {
        ++bond_counts( itr->GetAtomIndexLow());
        ++bond_counts( itr->GetAtomIndexHigh());
      }

      // test whether any atom has a valence
      bool had_valence( false);
      size_t atom_index( 0);
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End());
        itr != itr_end;
        ++itr, ++atom_index
      )
      {
        // only consider well-defined atom types
        if( itr->GetAtomType().IsDefined() && itr->GetAtomType()->IsGasteigerAtomType())
        {
          if( bond_counts( atom_index) < itr->GetAtomType()->GetNumberBonds())
          {
            had_valence = true;
            break;
          }
        }
      }

      // return true if there were any valences
      return had_valence;
    }

    //! @brief read CTabHandler object from std::istream
    //! @param ISTREAM istream that contains CTabHandler object
    //! @return istream after CTabHandler object was extracted
    std::istream &CTabHandler::Read( std::istream &ISTREAM)
    {
      return ReadCTab( ISTREAM);
    }

    //! @brief write CTabHandler into std::ostream
    //! @param OSTREAM ostream that gets CTabHandler object
    //! @param INDENT indentation
    //! @return ostream after CTabHandler object was inserted
    std::ostream &CTabHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CTabHandler::s_Instance
    (
      GetObjectInstances().AddInstance( new CTabHandler())
    );

  } // namespace sdf
} // namespace bcl

