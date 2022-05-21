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

#ifndef BCL_SDF_MDL_PROPERTY_H_
#define BCL_SDF_MDL_PROPERTY_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sdf_atom_info.h"
#include "bcl_sdf_bond_info.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MdlProperty
    //! @brief TODO: add an general comment to this class
    //!
    //! @see @link example_sdf_mdl_property.cpp @endlink
    //! @author mendenjl
    //! @date Feb 24, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MdlProperty :
      public util::ObjectInterface
    {
    public:
    //////////
    // enum //
    //////////

      enum Property
      {
        e_Charge,
        e_RingBondCount,
        e_AtomList,
        e_BclAtomType,
        e_BclBondType,
        e_BclChirality,
        e_BclDoubleBondIsometry,
        e_BclAtomAromaticity,
        e_BclRXNProperty,
        e_BlockTerminator,
        s_NumberProperties
      };

      enum PropertyClass
      {
        e_Atom, // applies only to atoms
        e_Bond, // applies only to bonds
        e_Molecule, // applies to whole molecules
        s_NumberPropertyClasses
      };

      //! @brief property as string
      //! @param PROPERTY the property desired
      //! @return the property as string
      static const std::string &GetPropertyName( const Property &PROPERTY);

      //! PropertyEnum simplifies the usage of the Property enum of this class
      typedef util::WrapperEnum< Property, &GetPropertyName, s_NumberProperties> PropertyEnum;

    private:

      //! @brief get the possible property prefixes
      static const storage::Vector< std::string> &GetAllPropertyPrefixes();

      //! @brief property as string, with M  prefix
      //! @param PROPERTY the property desired
      //! @return the property as string with the M  prefix
      static const std::string &GetPropertyNameWithPrefix( const Property &PROPERTY);

      //! @brief fixed width of property
      //! @param PROPERTY the property
      //! @return the properties fixed width (undefined if property is not fixed width)
      static const size_t &GetFixedWidthSize( const PropertyEnum &PROPERTY);

      //! @brief get the multiplicity of the property (e.g. the number of values per property)
      //! @param PROPERTY the property desired
      //! @return the multiplicity
      static const size_t &GetMultiplicity( const PropertyEnum &PROPERTY);

      //! @brief determine whether the # of values should be the first number on the line
      //! @param PROPERTY the property desired
      //! @return whether to print the size in the first field
      static bool FirstFieldIsSize( const PropertyEnum &PROPERTY);

      PropertyEnum                  m_Property;       //! the property value, if a BCL-known property
      std::string                   m_PropertyLabel;  //! the property label, e.g. CHG
      storage::Vector< std::string> m_PropertyStrings; //! the fields associated with this property

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MdlProperty( const PropertyEnum &PROPERTY = PropertyEnum( s_NumberProperties));

      //! @brief constructor given atom lines, bond lines, and the property desired
      MdlProperty
      (
        const PropertyEnum &PROPERTY,
        const storage::Vector< AtomInfo> &ATOM_LINES,
        const storage::Vector< BondInfo> &BOND_LINES
      );

      //! @brief construct from string containing property
      explicit MdlProperty( const std::string &STRING);

      //! @brief Clone function
      //! @return pointer to new MdlProperty
      MdlProperty *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! returns the property type
      const PropertyEnum &GetProperty() const
      {
        return m_Property;
      }

      //! @brief return the string of the property
      //! @return the string of the property
      std::string GetString() const;

      //! @brief check if the handler is valid
      //! @return true if the mdl handler is in a valid state
      bool IsValid() const
      {
        return m_Property != s_NumberProperties;
      }

      //! @brief report if the property is one recognized by the BCL
      //! @return true if it is a property explicitly handled by the BCL
      bool IsBCLProperty() const;

      //! @brief determine if an MDL property is one recognized by the BCL
      //! @param LINE the line to examine
      //! @return true if the property line contains information the BCL can explicitly recognize
      static bool IsBCLPropertyLine( const std::string &LINE);

      //! @brief helper function to determine whether a line is an MDL property line without copying or fully parsing it
      //! @param LINE the line of interest
      //! @return true if the line appears to be an MDL property line
      static bool IsMDLPropertyLine( const std::string &LINE);

      //! @brief get the label for this property, e.g. 'M  CHG'
      //! @return the full label
      const std::string &GetLabel() const
      {
        return m_PropertyLabel;
      }

      const storage::Vector< std::string> &GetPropertyStrings() const
      {
        return m_PropertyStrings;
      }

      //! determine if the property should be cached
      bool ShouldCache() const
      {
        return m_Property == e_BclRXNProperty || m_Property == e_BclAtomAromaticity || m_Property == s_NumberProperties;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief apply the property to a vector of atom and bond lines
      //! @param ATOM_LINES Mdl Atoms on which to apply the property
      //! @param BOND_LINES Mdl Bonds on which to apply the property
      void ApplyProperty( storage::Vector< AtomInfo> &ATOM_LINES, storage::Vector< BondInfo> &BOND_LINES) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief set the label according to the property number
      //! @details looks up the full label if a known property, e.g. "M  END", or leaves it blank
      void SetLabel( const PropertyEnum &PROPERTY);

      //! @brief return the string of the header
      //! @return the string of the header
      void SetFromString( const std::string &MDL_HEADER);

      //! @brief apply charges to atom lines
      //! @param ATOM_LINES mdl atoms to apply charges to
      void ApplyCharges( storage::Vector< AtomInfo> &ATOM_LINES) const;

      //! @brief apply atom types to atom lines
      //! @param ATOM_LINES mdl atoms to apply atom types to
      void ApplyAtomTypes( storage::Vector< AtomInfo> &ATOM_LINES) const;

      //! @brief apply bond types to bond lines
      //! @param BOND_LINES mdl bond to apply bond types to
      void ApplyBondTypes( storage::Vector< BondInfo> &BOND_LINES) const;

      //! @brief apply chiralities to mdl atoms
      //! @param ATOM_LINES mdl atoms to apply chirality property to
      void ApplyChiralities( storage::Vector< AtomInfo> &ATOM_LINES) const;

      //! @brief apply double bond isometries to MdlBonds
      //! @param BOND_LINES mdl bonds to apply isometry property to
      void ApplyDoubleBondIsometries( storage::Vector< BondInfo> &BOND_LINES) const;

      /*
      //! @brief apply an RXN property to atoms or bonds
      //! @param ATOM_LINES the AtomInfos to add any RXN property to
      //! @param BOND_LINES the BondInfos to add any RXN property to
      void ApplyRXNProperty
      ( 
        storage::Vector< AtomInfo> &ATOM_LINES, 
        storage::Vector< BondInfo> &BOND_LINES
      ) const;
      */

      //! @brief retrieve charges from atom lines
      //! @param ATOM_LINES mdl atoms to retrieve charges from
      void RetrieveCharges( const storage::Vector< AtomInfo> &ATOM_LINES);

      //! @brief retrieve atom types to atom lines
      //! @param ATOM_LINES mdl atoms to retrieve atom types from
      void RetrieveAtomTypes( const storage::Vector< AtomInfo> &ATOM_LINES);

      //! @brief retrieve bond types to bond lines
      //! @param BOND_LINES mdl bonds to retrieve bond types from
      void RetrieveBondTypes( const storage::Vector< BondInfo> &BOND_LINES);

      //! @brief retrieve chiralities to mdl atoms
      //! @param ATOM_LINES mdl atoms to retrieve chirality property from
      void RetrieveChiralities( const storage::Vector< AtomInfo> &ATOM_LINES);

      //! @brief retrieve double bond isometries from MdlBonds
      //! @param BOND_LINES mdl bonds to apply isometry property from
      void RetrieveDoubleBondIsometries( const storage::Vector< BondInfo> &BOND_LINES);
      
      /*
      //! @brief get special RXN properties from MdlAtoms and MdlBonds
      //! @note bonds are currently unused
      void RetrieveRXNProperties
      ( 
        const storage::Vector< AtomInfo> &ATOM_LINES, 
        const storage::Vector< BondInfo> &BOND_LINES // ignored for now
      );
      */

    }; // class MdlProperty

  } // namespace sdf
} // namespace bcl

#endif //BCL_SDF_MDL_PROPERTY_H_
