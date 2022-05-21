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
#include "chemistry/bcl_chemistry_constitutional_bond_type_data.h"

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

    //! @brief Conjugation as string
    //! @param CONJUGATION the conjugation whose name is desired
    //! @return the name as string
    const std::string &ConstitutionalBondTypeData::GetConjugationName( const Conjugation &CONJUGATION)
    {
      static const std::string s_Names[ size_t( s_NumberOfConjugations) + 1] =
      {
        "Nonconjugated",
        "Conjugated",
        "Aromatic",
        "Any",
        "Amide",
        GetStaticClassName< Conjugation>()
      };
      return s_Names[ CONJUGATION];
    }

    //! @brief Data as string
    //! @param DATA the data whose name is desired
    //! @return the name as string
    const std::string &ConstitutionalBondTypeData::GetDataName( const Data &DATA)
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
        "FuzzyBondOrderOrAromaticWithRingness",
        "FuzzyBondOrderAmideOrAromaticWithRingness",
        "ConstitutionalBondType",
        GetStaticClassName< Data>()
      };
      return s_Names[ DATA];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined bond type
    ConstitutionalBondTypeData::ConstitutionalBondTypeData() :
      m_SDFileID( 0),
      m_SDFileAltID( util::GetUndefined< size_t>()),
      m_Conjugation( ConstitutionalBondTypeData::e_Nonconjugated),
      m_BondTypeInRing()
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
    //! @param CONJUGATION the type of conjugation the bond is involved with
    //! @param IS_IN_RING whether the bond is in a ring
    ConstitutionalBondTypeData::ConstitutionalBondTypeData
    (
      const size_t &NUMBER_ELECTRONS,
      const size_t &SD_FILE_ID,
      const size_t &SD_FILE_ALT_ID,
      const Conjugation &CONJUGATION,
      const bool &IS_IN_RING
    ) :
      m_SDFileID( SD_FILE_ID),
      m_SDFileAltID( SD_FILE_ALT_ID),
      m_Conjugation( CONJUGATION),
      m_BondTypeInRing()
    {
      m_Data[ e_Identity] = size_t( 1);                // identity = 1 for any bond
      m_Data[ e_NumberOfElectrons] = NUMBER_ELECTRONS; // # e-

      m_Data[ e_Conjugation] = size_t( m_Conjugation); // Conjugation, size_t of the conjugation of the bond (0-3)
      m_Data[ e_IsAromatic] = size_t( m_Conjugation == e_Aromatic); // 1 for any bond that is aromatic, 0 otherwise
      m_Data[ e_IsAmide] = size_t( m_Conjugation == e_Amide); // 1 for any bond that is amide (single bond), 0 otherwise
      m_Data[ e_IsConjugated] = size_t( m_Conjugation == e_Conjugated || m_Conjugation == e_Aromatic || m_Conjugation == e_Amide);

      m_Data[ e_IsInRing] = size_t( IS_IN_RING);

      // the bond order is known if and only if the number of electrons is even
      const bool bond_order_known( !( NUMBER_ELECTRONS & size_t( 1)));

      // bond order is number of electrons / 2, if the bond order is well-defined
      m_Data[ e_BondOrder] = bond_order_known ? NUMBER_ELECTRONS / size_t( 2) : size_t( 0);

      // bond order unless in an aromatic ring, then 4 for aromatic, 0 for all bonds outside rings
      m_Data[ e_BondOrderOrAromatic] = m_Data[ e_IsAromatic] ? size_t( 4) : m_Data[ e_BondOrder];

      m_Data[ e_BondOrderAmideOrAromatic] = m_Data[ e_IsAromatic] || m_Data[ e_IsAmide]
                                            ? ( m_Data[ e_IsAromatic] ? size_t( 5) : size_t( 4))
                                            : m_Data[ e_BondOrder];

      // e_BondOrderOrAromatic, except 0 for all bonds outside rings
      m_Data[ e_BondOrderInRingOrAromatic] = IS_IN_RING ? m_Data[ e_BondOrderOrAromatic] : size_t( 0);

      // e_BondOrderOrAromatic, but distinguishes non-aromatic bonds in rings that have defined bond orders
      m_Data[ e_BondOrderOrAromaticWithRingness]
              = m_Data[ e_IsAromatic]
                ? size_t( 4)
                : bond_order_known
                  ? m_Data[ e_IsInRing] * size_t( 4) + m_Data[ e_BondOrder]
                  : size_t( 0);
      m_Data[ e_BondOrderAmideOrAromaticWithRingness]
              = m_Data[ e_IsAromatic] || m_Data[ e_IsAmide]
                ? ( m_Data[ e_IsAromatic] ? size_t( 5) : size_t( 4))
                : bond_order_known
                  ? m_Data[ e_IsInRing] * size_t( 5) + m_Data[ e_BondOrder]
                  : size_t( 0);

      // e_BondOrderOrAromaticWithRingness is same as e_BondOrderOrAromaticWithRingness, with 2 mapped to 1
      m_Data[ e_FuzzyBondOrderOrAromaticWithRingness] =
        m_Data[ e_BondOrderOrAromaticWithRingness] - !!( m_Data[ e_BondOrderOrAromaticWithRingness] > 1);
      m_Data[ e_FuzzyBondOrderAmideOrAromaticWithRingness] =
        m_Data[ e_BondOrderAmideOrAromaticWithRingness] - !!( m_Data[ e_BondOrderAmideOrAromaticWithRingness] > 1);
      m_Data[ e_ConstitutionalBondType] = util::GetUndefined< size_t>();
    }

    //! @brief virtual copy constructor
    ConstitutionalBondTypeData *ConstitutionalBondTypeData::Clone() const
    {
      return new ConstitutionalBondTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConstitutionalBondTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the data member specified by DATA_TO_RETRIEVE
    //! @param DATA_TO_RETRIEVE the data to get from the bond
    //! @return the desired data
    size_t ConstitutionalBondTypeData::GetBondData( const ConstitutionalBondTypeData::Data &DATA_TO_RETRIEVE) const
    {
      return m_Data[ DATA_TO_RETRIEVE];
    }

    //! @brief get the corresponding bond type with same properties except different conjugation
    //! @param CONJUGATION the alternative conjugation
    //! @return the corresponding bond type with same properties except different conjugation
    const ConstitutionalBondType &ConstitutionalBondTypeData::WithConjugation
    (
      const ConstitutionalBondTypeData::Conjugation &CONJUGATION
    ) const
    {
      return *m_ConjugatedTypes[ CONJUGATION];
    }

    //! @brief get the corresponding bond type inside a ring
    //! @return the corresponding bond type inside a ring (this type if it is already in a ring)
    const ConstitutionalBondType &ConstitutionalBondTypeData::WithInRing() const
    {
      return *m_BondTypeInRing;
    }

    //! @brief get the corresponding bond type with a different order
    //! @param ORDER the desired order
    //! @return the corresponding bond type with the desired order
    const ConstitutionalBondType &ConstitutionalBondTypeData::WithOrder( const size_t &ORDER) const
    {
      BCL_Assert( ORDER <= size_t( 3) && ORDER != size_t( 0), "requested bond type must have order between 1-3");
      return *m_AlternateOrderBondType[ ORDER - 1];
    }

    //! @brief get the corresponding configurational bond type with a certain isometry
    //! @param ISOMETRY the desired isometry
    //! @return the corresponding configurational bond type with a certain isometry
    const ConfigurationalBondType &ConstitutionalBondTypeData::WithIsometry( const BondIsometry &ISOMETRY) const
    {
      // if no one has ever called GetConfigurationalBondTypes before then m_ConfigurationalType is undefined,
      if( !m_IsometryTypes[ ISOMETRY].IsDefined())
      {
        GetConfigurationalBondTypes();
      }

      return *m_IsometryTypes[ ISOMETRY];
    }

    //! @brief get the bond type, except for aromatic bonds, return the unordered type
    //! @return the bond type, except for aromatic bonds, return the unordered type
    const ConstitutionalBondType &ConstitutionalBondTypeData::WithoutAromaticOrder() const
    {
      return *m_BondTypeSansAromaticOrder;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConstitutionalBondTypeData::Read( std::istream &ISTREAM)
    {
      // read members (use the standard stream operator here because
      io::Serialize::Read( m_SDFileID, ISTREAM);
      io::Serialize::Read( m_SDFileAltID, ISTREAM);
      io::Serialize::Read( m_Conjugation, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &ConstitutionalBondTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members on one line to make them easier to read
      io::Serialize::Write( m_SDFileID, OSTREAM, INDENT)
          << ' ' << m_SDFileAltID
          << ' ' << m_Conjugation;

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the index of this bond in the constitutional bond types, and setup related types
    //! @param TYPES the enum; needed to prevent a circular dependency
    //! This function can only be called by ConstitutionalBondTypes
    void ConstitutionalBondTypeData::AttachToRelatedBondTypes( const ConstitutionalBondTypes &TYPES)
    {
      // set all types to undefined
      m_BondTypeInRing = &TYPES.e_Undefined;
      for( size_t i( 0); i < s_NumberOfConjugations; ++i)
      {
        m_ConjugatedTypes[ i] = &TYPES.e_Undefined;
      }
      for( size_t i( 0), number_orders( 3); i < number_orders; ++i)
      {
        m_AlternateOrderBondType[ i] = &TYPES.e_Undefined;
      }
      m_BondTypeSansAromaticOrder = &TYPES.e_Undefined;
      m_Data[ e_ConstitutionalBondType] = TYPES.e_Undefined.GetIndex();

      if( !m_SDFileID) // undefined type, just return
      {
        return;
      }

      for( ConstitutionalBondTypes::const_iterator itr( TYPES.Begin()), itr_end( TYPES.End()); itr != itr_end; ++itr)
      {
        const ConstitutionalBondTypeData &other_bond( **itr);

        // determine whether the bonds are both inside or outside of a ring
        const bool ringness_matches( other_bond.IsBondInRing() == IsBondInRing());

        // determine whether the conjugation of the bonds matches
        const bool conjugation_matches( other_bond.GetConjugation() == GetConjugation());

        // determine how the other bond is related to this one, setup the related bond types accordingly
        if( other_bond.GetNumberOfElectrons() == GetNumberOfElectrons())
        {
          // same bond order
          if( conjugation_matches)
          {
            // same conjugation
            if( other_bond.IsBondInRing() == true)
            {
              // same bond type in a ring
              m_BondTypeInRing = &*itr;
            }
            // same ringness - all these together imply that it's the same bond type
            if( ringness_matches)
            {
              // exact same bond, setup the constitutional bond type index
              m_Data[ e_ConstitutionalBondType] = itr->GetIndex();

              // also setup the link to the same type via conjugation and bond order (if known)
              m_ConjugatedTypes[ other_bond.GetConjugation()] = &*itr;
              if( IsBondOrderKnown())
              {
                m_AlternateOrderBondType[ other_bond.GetBondData( e_BondOrder) - 1] = &*itr;
                if( m_Conjugation != e_Aromatic)
                {
                  // for all non-aromatic types; the bond types sans aromatic order is itself
                  m_BondTypeSansAromaticOrder = &*itr;
                }
              }
              else if( m_Conjugation == e_Aromatic)
              {
                // special case for the aromatic bond type without defined bond orders
                m_BondTypeSansAromaticOrder = &*itr;
              }
            }
          }
          else if( ringness_matches || other_bond.GetConjugation() == e_Aromatic)
          {
            // same # electrons and ring membership, only different conjugation, so set the type related by conjugation
            m_ConjugatedTypes[ other_bond.GetConjugation()] = &*itr;
          }
        }
        else if( ringness_matches && conjugation_matches)
        {
          if( other_bond.IsBondOrderKnown())
          {
            // same conjugation and ring-membership, but differing number of electrons
            // setup as the bond type of alternate order
            m_AlternateOrderBondType[ other_bond.GetBondData( e_BondOrder) - 1] = &*itr;
          }
          else if( m_Conjugation == e_Aromatic)
          {
            // other bond is the aromatic undefined bond type
            m_BondTypeSansAromaticOrder = &*itr;
          }
        }
        else if
        (
          ringness_matches
          && ( m_Conjugation == e_Amide || m_Conjugation == e_Nonconjugated)
          && other_bond.GetConjugation() == e_Conjugated
          && other_bond.IsBondOrderKnown()
        )
        {
          m_AlternateOrderBondType[ other_bond.GetBondData( e_BondOrder) - 1] = &*itr;
        }
      }

      // This code is for testing that the linkages between the conjugations are set up properly
//      std::ostringstream out;
//      out << TYPES.GetEnumFromIndex( m_Data[ e_ConstitutionalBondType]).GetName()
//          << " InRing: " << m_BondTypeInRing->GetName()
//          << " AltOrders: " << m_AlternateOrderBondType[ 0]->GetName() << ' '
//          << m_AlternateOrderBondType[ 1]->GetName() << ' '
//          << m_AlternateOrderBondType[ 2]->GetName() << ' '
//          << " AltConj: "
//          << m_ConjugatedTypes[ e_Nonconjugated]->GetName() << ' '
//          << m_ConjugatedTypes[ e_Conjugated]->GetName() << ' '
//          << m_ConjugatedTypes[ e_Aromatic]->GetName() << ' '
//          << m_ConjugatedTypes[ e_Amide]->GetName() << ' '
//          << " SansAromaticOrder: " <<  m_BondTypeSansAromaticOrder->GetName()
//          << '\n';
//      BCL_MessageStd( "Bond Info Check: " + out.str());
    }

    //! @brief set the related configurational bond type
    //! @param TYPE the corresponding configurational bond type, without isometry
    //! This function can only be called by ConfigurationalBondTypes
    void ConstitutionalBondTypeData::SetConfigurationalBondType( const ConfigurationalBondType &TYPE)
    {
      m_IsometryTypes[ TYPE->GetIsometry()] = &TYPE;
      if( !m_IsometryTypes[ e_UnknownIsometry]->IsDefined())
      {
        m_IsometryTypes[ e_UnknownIsometry] = &TYPE;
      }
    }

    //! @brief initialize all isometries with a particular type
    //! @param TYPE the undefined configurational bond type
    //! This function can only be called by ConfigurationalBondTypes
    void ConstitutionalBondTypeData::InitializeIsometries( const ConfigurationalBondType &TYPE)
    {
      for( size_t i( 0); i < s_NumberOfIsometries; ++i)
      {
        m_IsometryTypes[ i] = &TYPE;
      }
    }

    //! @brief Output operator for chemistry::ConstitutionalBondTypeData::Conjugation
    //! @param OSTREAM stream to write output for
    //! @param CONJUGATION enum to write
    std::ostream &operator <<
    (
      std::ostream &OSTREAM,
      const ConstitutionalBondTypeData::Conjugation &CONJUGATION
    )
    {
      return io::Serialize::Write( CONJUGATION, OSTREAM);
    }

    //! @brief Input operator for chemistry::ConstitutionalBondTypeData::Conjugation
    //! @param ISTREAM stream to read in from the stream
    //! @param CONJUGATION enum to read
    std::istream &operator >>
    (
      std::istream &ISTREAM,
      ConstitutionalBondTypeData::Conjugation &CONJUGATION
    )
    {
      return io::Serialize::Read( CONJUGATION, ISTREAM);
    }

  } // namespace chemistry
} // namespace bcl

