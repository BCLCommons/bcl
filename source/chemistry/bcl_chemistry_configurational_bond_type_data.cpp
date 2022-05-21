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
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitutional_bond_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  ///////////
  // Enums //
  ///////////

    //! @brief Data as string
    //! @param DATA the data whose name is desired
    //! @return the name as string
    const std::string &ConfigurationalBondTypeData::GetDataName( const Data &DATA)
    {
      static const std::string s_Names[ size_t( s_NumberOfData) + 1] =
      {
        "Identity",
        "BondOrder",
        "NumberOfElectrons",
        "Conjugation",
        "IsConjugated",
        "IsAromatic",
        "IsAmide",
        "IsInRing",
        "BondOrderInRingOrAromatic",
        "BondOrderOrAromatic",
        "BondOrderAmideOrAromatic",
        "BondOrderOrAromaticWithRingness",
        "BondOrderAmideOrAromaticWithRingness",
        "FuzzyBondOrderWithIsometryOrAromaticWithRingness",
        "FuzzyBondOrderAmideOrAromaticWithRingness",
        "ConstitutionalBondType",
        "BondOrderWithIsometry",
        "Isometry",
        "IsIsometric",
        "BondOrderWithIsometryOrAromatic",
        "BondOrderAmideWithIsometryOrAromaticWithRingness",
        "ConfigurationalBondType",
        GetStaticClassName< Data>()
      };
      return s_Names[ DATA];
    }

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined bond type
    ConfigurationalBondTypeData::ConfigurationalBondTypeData() :
      m_Isometry( e_NonIsometric)
    {
      // set all the data to undefined
      for( size_t data( 0), number_data( s_NumberOfData); data < number_data; ++data)
      {
        m_Data[ data] = util::GetUndefined< size_t>();
      }
    }

    //! @brief construct bond type from all its data
    //! @param NUMBER_ELECTRONS The total # of electrons involved in the bond
    //! @param SD_FILE_ID the bond id in an sdf file
    //! @param SD_FILE_ALT_ID for conjugated or aromatic types, the other # that can be used for the bond type
    //! @param PUBCHEM_ANNOTATION_ID the pubchem annotation of the bond in an SD-File (describes aromaticity, resonance, and stereochemistry)
    //! @param CONJUGATION the type of conjugation the bond is involved with
    //! @param IS_IN_RING whether the bond is in a ring
    //! @param ISOMETRY isometry of the bond type
    ConfigurationalBondTypeData::ConfigurationalBondTypeData
    (
      const ConstitutionalBondType &CONSTITUTIONAL_BOND_TYPE,
      const BondIsometry &ISOMETRY
    ) :
      m_Isometry( ISOMETRY),
      m_BaseBondType( CONSTITUTIONAL_BOND_TYPE)
    {
      m_Data[ e_Identity                       ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_Identity);
      m_Data[ e_BondOrder                      ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_BondOrder);
      m_Data[ e_NumberOfElectrons              ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_NumberOfElectrons);
      m_Data[ e_Conjugation                    ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_Conjugation);
      m_Data[ e_IsConjugated                   ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_IsConjugated);
      m_Data[ e_IsAromatic                     ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_IsAromatic);
      m_Data[ e_IsAmide                        ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_IsAmide);
      m_Data[ e_IsInRing                       ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_IsInRing);
      m_Data[ e_BondOrderInRingOrAromatic      ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_BondOrderInRingOrAromatic);
      m_Data[ e_BondOrderOrAromatic            ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_BondOrderOrAromatic);
      m_Data[ e_BondOrderAmideOrAromatic       ] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_BondOrderAmideOrAromatic);
      m_Data[ e_BondOrderOrAromaticWithRingness] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_BondOrderOrAromaticWithRingness);
      m_Data[ e_BondOrderAmideOrAromaticWithRingness] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_BondOrderAmideOrAromaticWithRingness);
      m_Data[ e_FuzzyBondOrderOrAromaticWithRingness] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_FuzzyBondOrderOrAromaticWithRingness);
      m_Data[ e_FuzzyBondOrderAmideOrAromaticWithRingness] = m_BaseBondType->GetBondData( ConstitutionalBondTypeData::e_FuzzyBondOrderAmideOrAromaticWithRingness);
      m_Data[ e_ConstitutionalBondType]          = m_BaseBondType.GetIndex();

      m_Data[ e_Isometry] = size_t( m_Isometry);    // Isometry, size_t of the isometry value of the bond (Undef,0,1,2)
      m_Data[ e_IsIsometric] = m_Isometry != e_NonIsometric && m_Isometry != e_UnknownIsometry;

      // like bond order, but distinguishes double bonds with isometry by labeling them 4 (E), or 5 (Z)
      m_Data[ e_BondOrderWithIsometry] = m_Data[ e_BondOrder];
      if( m_Data[ e_IsIsometric])
      {
        m_Data[ e_BondOrderWithIsometry] = size_t( 3) + m_Isometry;
      }

      // e_BondOrderWithIsometry, except aromatic bonds labeled 6
      m_Data[ e_BondOrderWithIsometryOrAromatic] = m_Data[ e_IsAromatic] ? size_t( 6) : m_Data[ e_BondOrderWithIsometry];

      m_Data[ e_BondOrderAmideWithIsometryOrAromaticWithRingness]
              = m_Data[ e_IsAromatic] || m_Data[ e_IsAmide]
                ? ( m_Data[ e_IsAromatic] ? size_t( 7) : size_t( 6))
                : m_Data[ e_IsInRing] * size_t( 7) + m_Data[ e_BondOrderWithIsometry];

      // index of this bond type data in ConfigurationalBondTypes; initially undefined
      m_Data[ e_ConfigurationalBondType] = util::GetUndefined< size_t>();
    }

    //! @brief virtual copy constructor
    ConfigurationalBondTypeData *ConfigurationalBondTypeData::Clone() const
    {
      return new ConfigurationalBondTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConfigurationalBondTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the data member specified by DATA_TO_RETRIEVE
    //! @param DATA_TO_RETRIEVE the data to get from the bond
    //! @return the desired data
    size_t ConfigurationalBondTypeData::GetBondData( const ConfigurationalBondTypeData::Data &DATA_TO_RETRIEVE) const
    {
      return m_Data[ DATA_TO_RETRIEVE];
    }

    //! @brief get the corresponding bond type with same properties except different conjugation
    //! @param CONJUGATION the alternative conjugation
    //! @return the corresponding bond type with same properties except different conjugation
    const ConfigurationalBondType &ConfigurationalBondTypeData::WithConjugation
    (
      const ConstitutionalBondTypeData::Conjugation &CONJUGATION
    ) const
    {
      return *m_ConjugatedTypes[ CONJUGATION];
    }

    //! @brief get the corresponding bond type inside a ring
    //! @return the corresponding bond type inside a ring (this type if it is already in a ring)
    const ConfigurationalBondType &ConfigurationalBondTypeData::WithInRing() const
    {
      return *m_BondTypeInRing;
    }

    //! @brief get the corresponding bond type with a different order
    //! @param ORDER the desired order
    //! @return the corresponding bond type with the desired order
    const ConfigurationalBondType &ConfigurationalBondTypeData::WithOrder( const size_t &ORDER) const
    {
      BCL_Assert
      (
        ORDER <= size_t( 3) && ORDER != size_t( 0),
        "requested bond type must have order between 1-3; requested " + util::Format()( ORDER)
      );
      return *m_AlternateOrderBondType[ ORDER - 1];
    }

    //! @brief get the corresponding bond type with a different isometry
    //! @param ISOMETRY the desired isometry
    //! @return the corresponding bond type with the desired isometry
    const ConfigurationalBondType &ConfigurationalBondTypeData::WithIsometry( const BondIsometry &ISOMETRY) const
    {
      return m_BaseBondType->WithIsometry( ISOMETRY);
    }

    //! @brief get the corresponding bond type with e_NonIsometric; much shorter than calling
    //! Shorthand for WithIsometry( e_NonIsometric)
    //! @return the corresponding bond type with e_NonIsometric
    const ConfigurationalBondType &ConfigurationalBondTypeData::WithoutIsometry() const
    {
      return m_BaseBondType->WithIsometry( e_NonIsometric);
    }

    //! @brief get the bond type, except for aromatic bonds, return the unordered type
    //! @return the bond type, except for aromatic bonds, return the unordered type
    const ConfigurationalBondType &ConfigurationalBondTypeData::WithoutAromaticOrder() const
    {
      return *m_BondTypeSansAromaticOrder;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConfigurationalBondTypeData::Read( std::istream &ISTREAM)
    {
      // read members (use the standard stream operator here because
      io::Serialize::Read( m_BaseBondType, ISTREAM);
      io::Serialize::Read( m_Isometry, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &ConfigurationalBondTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members on one line to make them easier to read
      io::Serialize::Write( m_BaseBondType, OSTREAM, INDENT) << ' ' << m_Isometry;

      // end
      return OSTREAM;
    }

    //! @brief setup connections between this bond type and related types
    //! @param THIS_TYPE the corresponding configurational bond type
    //! This function can only be called by ConfigurationalBondTypes
    void ConfigurationalBondTypeData::AttachToRelatedBondTypes( const ConfigurationalBondType &THIS_TYPE)
    {
      m_Data[ e_ConfigurationalBondType] = THIS_TYPE.GetIndex();

      if( !m_BaseBondType.IsDefined())
      {
        // set up the related types to be all undefined
        for( size_t i( 0); i < ConstitutionalBondTypeData::s_NumberOfConjugations; ++i)
        {
          m_ConjugatedTypes[ i] = &THIS_TYPE;
        }
        for( size_t i( 0), number_orders( 3); i < number_orders; ++i)
        {
          m_AlternateOrderBondType[ i] = &THIS_TYPE;
        }
        m_BondTypeInRing = &THIS_TYPE;
        m_BondTypeSansAromaticOrder = &THIS_TYPE;
        return;
      }

      // compute the related bond types based on the constitution's related bond types
      m_BondTypeInRing = &m_BaseBondType->WithInRing()->WithIsometry( m_Isometry);
      for( size_t i( 0); i < ConstitutionalBondTypeData::s_NumberOfConjugations; ++i)
      {
        m_ConjugatedTypes[ i] =
          &m_BaseBondType->WithConjugation( ConstitutionalBondTypeData::Conjugation( i))->WithIsometry( e_NonIsometric);
      }
      for( size_t bond_order( 1), max_bond_order( 3); bond_order <= max_bond_order; ++bond_order)
      {
        m_AlternateOrderBondType[ bond_order - 1] = &m_BaseBondType->WithOrder( bond_order)->WithIsometry( e_NonIsometric);
      }

      if( m_Isometry != e_NonIsometric)
      {
        m_ConjugatedTypes[ GetConjugation()] = &THIS_TYPE;
        m_AlternateOrderBondType[ m_Data[ e_BondOrder] - 1] = &THIS_TYPE;
        m_BondTypeSansAromaticOrder = &THIS_TYPE;
      }
      else
      {
        m_BondTypeSansAromaticOrder = &m_BaseBondType->WithoutAromaticOrder()->WithIsometry( m_Isometry);
      }

    }

  } // namespace chemistry

} // namespace bcl

