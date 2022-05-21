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
#include "sdf/bcl_sdf_mdl_property.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    //! @brief property as string
    //! @param PROPERTY the property desired
    //! @return the property as string
    const std::string &MdlProperty::GetPropertyName( const MdlProperty::Property &PROPERTY)
    {
      static const std::string s_names[] =
      {
        "CHG",
        "RBC",
        "ALS",
        "BCL ATM",
        "BCL BND",
        "BCL CHI",
        "BCL DBI",
        "BCL ARO", // atom aromaticity
        "BCL RXN",
        "END",
        GetStaticClassName< Property>()
      };
      return s_names[ PROPERTY];
    }

    //! @brief get the possible property prefixes
    const storage::Vector< std::string> &MdlProperty::GetAllPropertyPrefixes()
    {
      static const storage::Vector< std::string> s_prefixes = storage::Vector< std::string>::Create( "M  ", "G  ", "A  ", "V  ");
      return s_prefixes;
    }

    //! @brief property as string, with M  prefix
    //! @param PROPERTY the property desired
    //! @return the property as string with the M  prefix
    const std::string &MdlProperty::GetPropertyNameWithPrefix( const MdlProperty::Property &PROPERTY)
    {
      static const std::string s_names[] =
      {
        "M  " + GetPropertyName( e_Charge),
        "M  " + GetPropertyName( e_RingBondCount),
        "M  " + GetPropertyName( e_AtomList),
        "M  " + GetPropertyName( e_BclAtomType),
        "M  " + GetPropertyName( e_BclBondType),
        "M  " + GetPropertyName( e_BclChirality),
        "M  " + GetPropertyName( e_BclDoubleBondIsometry),
        "M  " + GetPropertyName( e_BclAtomAromaticity),
        "M  " + GetPropertyName( e_BclRXNProperty),
        "M  " + GetPropertyName( e_BlockTerminator),
        GetStaticClassName< Property>()
      };
      return s_names[ PROPERTY];
    }

    //! @brief fixed width of property
    //! @param PROPERTY the property
    //! @return the properties fixed width (undefined if property is not fixed width)
    const size_t &MdlProperty::GetFixedWidthSize( const MdlProperty::PropertyEnum &PROPERTY)
    {
      static const size_t s_fixed_widths[] =
      {
        3, // CHG
        3, // RBC
        util::GetUndefined< size_t>(), // ALS
        util::GetUndefined< size_t>(), // BCL ATM
        util::GetUndefined< size_t>(), // BCL BND
        3, // BCL CHI
        1, // BCL DBI
        util::GetUndefined< size_t>(), // atom aromaticity
        util::GetUndefined< size_t>(), // BCL RXN
        util::GetUndefined< size_t>() // END
      };
      return s_fixed_widths[ PROPERTY];
    }

    //! @brief get the multiplicity of the property (e.g. the number of values per property)
    //! @param PROPERTY the property desired
    //! @return the multiplicity
    const size_t &MdlProperty::GetMultiplicity( const MdlProperty::PropertyEnum &PROPERTY)
    {
      static const size_t s_multiplities[] =
      {
        2, // atom index (1 offset) and charge
        2, // atom index (1 offset) and number of ring bonds
        1, // atom index (1 offset) followed by list of element types
        1, // atom type
        1, // bond type
        2, // atom index (1 offset) and chirality
        1, // just isomorphism string
        2, // list of aromatic atoms
        1, // rxn properties
        0, // no optional strings
        0
      };
      return s_multiplities[ PROPERTY];
    }

    //! @brief determine whether the # of values should be the first number on the line
    //! @param PROPERTY the property desired
    //! @return whether to print the size in the first field
    bool MdlProperty::FirstFieldIsSize( const PropertyEnum &PROPERTY)
    {
      static const bool s_first_field_is_size[] =
      {
        true,  // e_Charge
        true,  // e_RingBondCount
        false, // e_AtomList
        false, // e_BclAtomType
        false, // e_BclBondType
        false, // e_BclChirality
        false, // e_BclDoubleBondIsometry
        false, // e_BclAtomAromaticity
        false, // e_BclRXNProperty
        false, // e_BlockTerminator
        false  // s_NumberProperties
      };
      return s_first_field_is_size[ PROPERTY];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MdlProperty::MdlProperty( const PropertyEnum &PROPERTY) :
      m_Property( PROPERTY)
    {
      SetLabel( PROPERTY);
    }

    //! @brief constructor given atom lines, bond lines, and the property desired
    MdlProperty::MdlProperty
    (
      const PropertyEnum &PROPERTY,
      const storage::Vector< AtomInfo> &ATOM_LINES,
      const storage::Vector< BondInfo> &BOND_LINES
    ) :
      m_Property( PROPERTY)
    {
      SetLabel( PROPERTY);

      // choose the function based on the property
      switch( m_Property)
      {
        case e_Charge:
          RetrieveCharges( ATOM_LINES);
          break;
        case e_RingBondCount:
          break;
        case e_AtomList:
          break;
        case e_BclAtomType:
          RetrieveAtomTypes( ATOM_LINES);
          break;
        case e_BclBondType:
          RetrieveBondTypes( BOND_LINES);
          break;
        case e_BclChirality:
          RetrieveChiralities( ATOM_LINES);
          break;
        case e_BclDoubleBondIsometry:
          RetrieveDoubleBondIsometries( BOND_LINES);
          break;
        case e_BclAtomAromaticity:
          break; // read-only from a Molfile now, can't do anything with atom/bond lines
        case e_BclRXNProperty:
          //RetrieveRXNProperties( ATOM_LINES, BOND_LINES);
          //break;
        case e_BlockTerminator:
        case s_NumberProperties:
        default:
          break;
      }
    }

    //! @brief construct from string containing property
    MdlProperty::MdlProperty( const std::string &STRING) :
      m_Property( s_NumberProperties)
    {
      SetFromString( STRING);
    }

    //! @brief Clone function
    //! @return pointer to new MdlProperty
    MdlProperty *MdlProperty::Clone() const
    {
      return new MdlProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MdlProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the string of the property
    //! @return the string of the property
    std::string MdlProperty::GetString() const
    {

      // unknown/generic property
      if( m_Property == s_NumberProperties)
      {
        //return std::string();

        // write out the properties separated by a space
        std::stringstream strm;
        strm << m_PropertyLabel;
        if( !m_PropertyLabel.empty() && !m_PropertyStrings.IsEmpty())
        {
          strm << " ";
          for( size_t p( 0), end_p( m_PropertyStrings.GetSize() - 1); p < end_p; ++p)
          {
            strm << m_PropertyStrings( p) << " ";
          }
          strm << m_PropertyStrings( m_PropertyStrings.GetSize() - 1);
        }
        return strm.str();
      }

      // if the property string requires values (e.g. multiplicity != 0), but none are available, then
      // just return
      const size_t multiplicity( GetMultiplicity( m_Property));
      if( multiplicity != size_t( 0) && m_PropertyStrings.IsEmpty())
      {
        return std::string();
      }

      std::string property_string( GetPropertyNameWithPrefix( m_Property));

      const size_t fixed_width( GetFixedWidthSize( m_Property));

      util::Format string_formatter;
      if( util::IsDefined( fixed_width))
      {
        string_formatter.W( fixed_width);
      }

      if( FirstFieldIsSize( m_Property))
      {
        property_string += string_formatter( m_PropertyStrings.GetSize() / multiplicity);
      }
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
        itr != itr_end;
        ++itr
      )
      {
        property_string += ' ';
        property_string += string_formatter( *itr);
      }
      return property_string;
    }

    //! @brief report if the property is one recognized by the BCL
    //! @return true if it is a property explicitly handled by the BCL
    bool MdlProperty::IsBCLProperty() const
    {
      return m_Property < s_NumberProperties;
    }

    //! @brief determine if an MDL property is one recognized by the BCL
    //! @param LINE the line to examine
    //! @return true if the property line contains information the BCL can explicitly recognize
    bool MdlProperty::IsBCLPropertyLine( const std::string &LINE)
    {
      bool is_bcl_prop( false);
      if( util::StartsWith( LINE, "M  ")) // all BCL properties start with 'M  ' currently
      {

        // determine whether the property is one of the known properties
        Property test_property;
        for( size_t property_id( 0); property_id < s_NumberProperties; ++property_id)
        {
          test_property = Property( property_id);
          if( util::StartsWith( LINE, GetPropertyNameWithPrefix( test_property)))
          {
            // set variable if the property matches one known to the BCL
            is_bcl_prop = true;
            break;
          }
        }
      }

      return is_bcl_prop;
    }

    //! @brief helper function to determine whether a line is an MDL property line without copying or fully parsing it
    //! @param LINE the line of interest
    //! @return true if the line appears to be an MDL property line
    bool MdlProperty::IsMDLPropertyLine( const std::string &LINE)
    {
      // ensure that the property line looks like a normal starting property line
      const storage::Vector< std::string> &prefixes( GetAllPropertyPrefixes());
      bool has_good_prefix( false);
      for( size_t p( 0), end_p( prefixes.GetSize()); p < end_p; ++p)
      {
        //if( !util::StartsWith( LINE, "M  "))
        if( util::StartsWith( LINE, prefixes( p).c_str()))
        {
          has_good_prefix = true;
          break;
        }
      }
      return has_good_prefix;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief apply the property to a vector of atom and bond lines
    //! @param ATOM_LINES Mdl Atoms on which to apply the property
    //! @param BOND_LINES Mdl Bonds on which to apply the property
    void MdlProperty::ApplyProperty( storage::Vector< AtomInfo> &ATOM_LINES, storage::Vector< BondInfo> &BOND_LINES) const
    {
      // choose the function based on the property
      switch( m_Property)
      {
        case e_Charge:
          ApplyCharges( ATOM_LINES);
          break;
        case e_BclAtomType:
          ApplyAtomTypes( ATOM_LINES);
          break;
        case e_BclBondType:
          ApplyBondTypes( BOND_LINES);
          break;
        case e_BclChirality:
          ApplyChiralities( ATOM_LINES);
          break;
        case e_BclDoubleBondIsometry:
          ApplyDoubleBondIsometries( BOND_LINES);
          break;
        case e_BclRXNProperty:
        case e_BlockTerminator:
        case s_NumberProperties:
        default:
          break;
      }
    }

    //! @brief read MdlProperty object from std::istream
    //! @param ISTREAM istream that contains MdlProperty object t
    //! @return istream after MdlProperty object was extracted
    std::istream &MdlProperty::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Property, ISTREAM);
      io::Serialize::Read( m_PropertyStrings, ISTREAM);
      // return
      return ISTREAM;
    }

    //! @brief write MdlProperty into std::ostream
    //! @param OSTREAM ostream that gets MdlProperty object
    //! @param INDENT indentation
    //! @return ostream after MdlProperty object was inserted
    std::ostream &MdlProperty::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Property, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PropertyStrings, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the label according to the property number
    //! @details looks up the full label if a known property, e.g. "M  END", or leaves it blank
    void MdlProperty::SetLabel( const PropertyEnum &PROPERTY)
    {
      if( PROPERTY < s_NumberProperties)
      {
        m_PropertyLabel = GetPropertyNameWithPrefix( PROPERTY);
      }
    }

    //! @brief return the string of the header
    //! @return the string of the header
    void MdlProperty::SetFromString( const std::string &MDL_STRING)
    {
      m_Property = s_NumberProperties;
      m_PropertyLabel = std::string();
      m_PropertyStrings.Reset();

      // ensure that the property line looks like a normal starting property line
      // test all possible prefixes to see if we match any of them
      // all prefixes are 3 characters long, e.g. 'M  '
      const storage::Vector< std::string> &possible_prefixes( GetAllPropertyPrefixes());
      std::string str_prefix( MDL_STRING.substr( 0, 3));
      bool found_prefix( false);
      for( size_t pre_no( 0), end_no( possible_prefixes.GetSize()); pre_no < end_no; ++pre_no)
      {
        if( str_prefix == possible_prefixes( pre_no))
        {
          found_prefix = true;
          break;
        }
      }
      if( !found_prefix)
      {
        BCL_MessageStd( "String \"" + MDL_STRING + "\" does not contain a valid MDL property prefix; ignoring");
        return;
      }

      const size_t property_name_start( 3);
      size_t property_values_start( 0);

      // determine whether the property is one of the known properties
      bool is_known_property( false);
      Property test_property;
      for( size_t property_id( 0); property_id < s_NumberProperties; ++property_id)
      {
        test_property = Property( property_id);
        const std::string &property_string( GetPropertyName( test_property));
        if( MDL_STRING.substr( property_name_start, property_string.size()) == property_string)
        {
          property_values_start = property_name_start + property_string.size();
          m_Property = test_property;
          is_known_property = true;
          break;
        }
      }

      if( is_known_property)
      {

        // store the label, including the 'M  ' portion
        m_PropertyLabel = MDL_STRING.substr( 0, property_values_start);

        // if the multiplicity is 0, then there should be nothing else to read
        if( GetMultiplicity( m_Property) == size_t( 0))
        {
          return;
        }

        // determine width of this type
        const size_t fixed_width( GetFixedWidthSize( m_Property));

        if( !FirstFieldIsSize( m_Property)) // non-size fields have a space between property name and first value
        {
          ++property_values_start;
        }

        if( util::IsDefined( fixed_width))
        {
          for
          (
            size_t token_start( property_values_start), string_size( MDL_STRING.size());
            token_start < string_size;
            token_start += fixed_width + 1
          )
          {
            m_PropertyStrings.PushBack( util::TrimString( MDL_STRING.substr( token_start, fixed_width)));
          }
        }
        else
        {
          m_PropertyStrings = util::SplitString( MDL_STRING.substr( property_values_start));
        }

        if( FirstFieldIsSize( m_Property))
        {
          m_PropertyStrings.Remove( m_PropertyStrings.Begin());
        }

        // check multiplicity
        BCL_Assert
        (
          !( m_PropertyStrings.GetSize() % GetMultiplicity( m_Property)),
          util::Format()( m_PropertyStrings.GetSize()) + "values given for property: " + GetPropertyName( m_Property)
          + "; should have been a multiple of " + util::Format()( GetMultiplicity( m_Property))
        );
      }
      else
      {
        // unknown property, but we can keep the information anyway.
        // Assume the property label is space-separated from information, i.e. of the format 'X  YYY ..."

        // find the first space character in the string
        size_t label_pos( 3);
        size_t str_len( MDL_STRING.length());

        // skip immediate blank characters
        for( ; label_pos < str_len && isspace( MDL_STRING[ label_pos]); ++label_pos);
        size_t prop_name_start( label_pos);

        // find the first blank space character
        for( ; label_pos < str_len && !isspace( MDL_STRING[ label_pos]); ++label_pos);

        // if label is empty, stop now
        if( label_pos == prop_name_start)
        {
          return;
        }

        // set the property label and then split the remaining part of the line on spaces for the property strings
        m_PropertyLabel = MDL_STRING.substr( 0, label_pos);
        m_PropertyStrings = util::SplitString( MDL_STRING.substr( label_pos));
      }
    }

    //! @brief apply charges to atom lines
    //! @param ATOM_LINES mdl atoms to apply charges to
    void MdlProperty::ApplyCharges( storage::Vector< AtomInfo> &ATOM_LINES) const
    {
      // set the charges on each atom
      for
      (
        size_t charge_number( 0), number_charges( m_PropertyStrings.GetSize());
        charge_number < number_charges;
        charge_number += 2
      )
      {
        const size_t atom_index( util::ConvertStringToNumericalValue< size_t>( m_PropertyStrings( charge_number)) - 1);
        const short charge( util::ConvertStringToNumericalValue< short>( m_PropertyStrings( charge_number + 1)));
        BCL_Assert
        (
          atom_index < ATOM_LINES.GetSize(),
          "Tried to set charge on atom index " + util::Format()( atom_index) + " beyond end"
        );

        ATOM_LINES( atom_index).SetCharge( charge);
      }
    }

    //! @brief apply atom types to atom lines
    //! @param ATOM_LINES mdl atoms to apply atom types to
    void MdlProperty::ApplyAtomTypes( storage::Vector< AtomInfo> &ATOM_LINES) const
    {
      if( ATOM_LINES.GetSize() != m_PropertyStrings.GetSize())
      {
        BCL_MessageCrt
        (
          "Incorrect # of atom types specified!  Needed " + util::Format()( ATOM_LINES.GetSize()) + " but only " +
          util::Format()( m_PropertyStrings.GetSize()) + " were provided. Will leave atom types undefined for this molecule"
        );
        return;
      }
      size_t counter( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        chemistry::AtomType type( chemistry::GetAtomTypes().GetAtomType( *itr));
        BCL_Assert( type.IsDefined(), "Undefined atom type in atom type property: " + *itr);
        BCL_Assert
        (
          type->GetElementType() == ATOM_LINES( counter).GetAtomType()->GetElementType(),
          "Tried to change element type when setting atom type!"
        );
        ATOM_LINES( counter).SetAtomType( type);
      }
    }

    //! @brief apply bond types to bond lines
    //! @param BOND_LINES mdl bond to apply bond types to
    void MdlProperty::ApplyBondTypes( storage::Vector< BondInfo> &BOND_LINES) const
    {
      if( BOND_LINES.GetSize() != m_PropertyStrings.GetSize())
      {
        BCL_MessageCrt
        (
          "Incorrect # of bond types! Needed " + util::Format()( BOND_LINES.GetSize()) + " but " + util::Format()( m_PropertyStrings.GetSize()) +
          " were provided. Will leave bond types on this molecule undefined"
        );
        return;
      }
      size_t counter( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        chemistry::ConstitutionalBondType type( *itr);
        BCL_Assert( type.IsDefined(), "Undefined bond type in bond type property: " + *itr);
        BondInfo old( BOND_LINES( counter));
        BOND_LINES( counter) = BondInfo( old.GetAtomIndexLow(), old.GetAtomIndexHigh(), type);
        BOND_LINES( counter).SetIsometry( old.GetConfigurationalBondType()->GetIsometry());
      }
    }

    //! @brief apply chiralities to mdl atoms
    //! @param ATOM_LINES mdl atoms to apply chirality property to
    void MdlProperty::ApplyChiralities( storage::Vector< AtomInfo> &ATOM_LINES) const
    {
      // reset all chiralities to non-chiral
      for
      (
        storage::Vector< AtomInfo>::iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->SetChirality( chemistry::e_NonChiral);
      }
      // apply chiralities to each atom
      for
      (
        size_t chirality_number( 0), number_chiralities( m_PropertyStrings.GetSize());
        chirality_number < number_chiralities;
        chirality_number += 2
      )
      {
        // written index is 1-offset (for consistency with MDL format), so convert to 0-offset
        const size_t atom_index
        (
          util::ConvertStringToNumericalValue< size_t>( m_PropertyStrings( chirality_number)) - 1
        );
        const chemistry::Chirality chirality
        (
          chemistry::ChiralityEnum( m_PropertyStrings( chirality_number + 1))
        );
        BCL_Assert
        (
          atom_index < ATOM_LINES.GetSize(),
          "Tried to set chirality on atom index " + util::Format()( atom_index) + " beyond end"
        );
        BCL_Assert
        (
          chirality != chemistry::s_NumberChiralities,
          "Bad chirality name " + m_PropertyStrings( chirality_number + 1)
        );
        if( chirality == chemistry::e_UnknownChirality)
        {
          if( ATOM_LINES( atom_index).GetAtomType()->GetNumberBonds() == size_t( 4)
              || !ATOM_LINES( atom_index).GetAtomType()->IsGasteigerAtomType())
          {
            ATOM_LINES( atom_index).SetChirality( chirality);
          }
            continue;
        }
        BCL_Assert
        (
          ATOM_LINES( atom_index).GetAtomType()->GetNumberBonds() == size_t( 4)
          || !ATOM_LINES( atom_index).GetAtomType()->IsGasteigerAtomType(),
          "Tried to set chirality on atom that does not have 4 bonds! "
        );

        ATOM_LINES( atom_index).SetChirality( chirality);
      }
    }

    //! @brief apply double bond isometries to MdlBond with type ConfigurationalBondTypes().e_DoubleBond
    //! @param BOND_LINES mdl bonds to apply isometry property to
    void MdlProperty::ApplyDoubleBondIsometries( storage::Vector< BondInfo> &BOND_LINES) const
    {
      storage::Vector< std::string>::const_iterator itr( m_PropertyStrings.Begin()), itr_end( m_PropertyStrings.End());
      storage::Vector< BondInfo>::iterator itr_bond( BOND_LINES.Begin()), itr_bond_end( BOND_LINES.End());

      for( ; itr != itr_end && itr_bond != itr_bond_end; ++itr, ++itr_bond)
      {
        // get the next bond isometry
        const chemistry::BondIsometry isometry = chemistry::BondIsometryEnum( *itr);

        // ensure it is defined
        BCL_Assert( isometry != chemistry::s_NumberOfIsometries, "Undefined double bond isometry: " + *itr);

        // skip over non-double bonds
        while( itr_bond != itr_bond_end && itr_bond->GetConfigurationalBondType()->GetNumberOfElectrons() != size_t( 4))
        {
          ++itr_bond;
        }

        // stop if there are no further double bonds
        if( itr_bond == itr_bond_end)
        {
          break;
        }

        itr_bond->SetIsometry( isometry);
      }

      if( itr != itr_end)
      {
        BCL_MessageCrt( "More double bond isometries given than double bonds in molecule! File corruption likely.");
        return;
      }

      // skip any remaining non-double bonds
      while( itr_bond != itr_bond_end && itr_bond->GetConfigurationalBondType()->GetNumberOfElectrons() != size_t( 4))
      {
        ++itr_bond;
      }

      if( itr != itr_end)
      {
        BCL_MessageCrt( "Double bond isometries not given for all double bonds! File corruption likely.");
        return;
      }
    }

    //! @brief retrieve charges from atom lines
    //! @param ATOM_LINES mdl atoms to retrieve charges from
    void MdlProperty::RetrieveCharges( const storage::Vector< AtomInfo> &ATOM_LINES)
    {
      m_PropertyStrings.Reset();

      size_t index( 1); // charge indices are 1-offset
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr, ++index
      )
      {
        if( itr->GetAtomType()->GetFormalCharge()) // only write non-zero charges
        {
          m_PropertyStrings.PushBack( util::Format()( index));
          m_PropertyStrings.PushBack( util::Format()( itr->GetAtomType()->GetFormalCharge()));
        }
      }
    }

    //! @brief retrieve atom types to atom lines
    //! @param ATOM_LINES mdl atoms to retrieve atom types from
    void MdlProperty::RetrieveAtomTypes( const storage::Vector< AtomInfo> &ATOM_LINES)
    {
      m_PropertyStrings.Reset();
      m_PropertyStrings.AllocateMemory( ATOM_LINES.GetSize());
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr
      )
      {
        m_PropertyStrings.PushBack( itr->GetAtomType().GetName());
      }
    }

    //! @brief retrieve bond types to bond lines
    //! @param BOND_LINES mdl bonds to retrieve bond types from
    void MdlProperty::RetrieveBondTypes( const storage::Vector< BondInfo> &BOND_LINES)
    {
      m_PropertyStrings.Reset();
      m_PropertyStrings.AllocateMemory( BOND_LINES.GetSize());
      for
      (
        storage::Vector< BondInfo>::const_iterator itr( BOND_LINES.Begin()), itr_end( BOND_LINES.End());
        itr != itr_end;
        ++itr
      )
      {
        m_PropertyStrings.PushBack( itr->GetConstitutionalBondType().GetName());
      }
    }

    //! @brief retrieve chiralities to mdl atoms
    //! @param ATOM_LINES mdl atoms to retrieve chirality property from
    void MdlProperty::RetrieveChiralities( const storage::Vector< AtomInfo> &ATOM_LINES)
    {
      m_PropertyStrings.Reset();

      size_t index( 1); // chirality indices are 1-offset
      for
      (
        storage::Vector< AtomInfo>::const_iterator itr( ATOM_LINES.Begin()), itr_end( ATOM_LINES.End());
        itr != itr_end;
        ++itr, ++index
      )
      {
        if( itr->GetChirality() != chemistry::e_NonChiral && itr->GetChirality() != chemistry::e_UnknownChirality) // only write out chiral centers
        {
          m_PropertyStrings.PushBack( util::Format()( index));
          m_PropertyStrings.PushBack( chemistry::ChiralityEnum( itr->GetChirality()));
        }
      }
    }

    //! @brief retrieve double bond isometries from MdlBonds
    //! @param BOND_LINES mdl bonds to apply isometry property from
    void MdlProperty::RetrieveDoubleBondIsometries( const storage::Vector< BondInfo> &BOND_LINES)
    {
      m_PropertyStrings.Reset();

      for
      (
        storage::Vector< BondInfo>::const_iterator itr_bond( BOND_LINES.Begin()), itr_bond_end( BOND_LINES.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        // skip bonds that could not have isometry
        if( itr_bond->GetConfigurationalBondType()->GetNumberOfElectrons() == size_t( 4))
        {
          // push back the next bond isometry
          m_PropertyStrings.PushBack
          (
            chemistry::BondIsometryEnum( itr_bond->GetConfigurationalBondType()->GetIsometry())
          );
        }
      }
    }
  } // namespace sdf
} // namespace bcl
