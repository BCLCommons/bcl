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

#ifndef BCL_CHEMISTRY_CONFIGURATIONAL_BOND_TYPE_DATA_H_
#define BCL_CHEMISTRY_CONFIGURATIONAL_BOND_TYPE_DATA_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_bond_isometry.h"
#include "bcl_chemistry_constitutional_bond_type_data.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConfigurationalBondTypeData
    //! @brief stores bond properties
    //! @details This class stores bond properties that are independent of bonded element types
    //!
    //! @see @link example_chemistry_configurational_bond_type_data.cpp @endlink
    //! @author mendenjl
    //! @date Dec 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConfigurationalBondTypeData :
      public util::ObjectInterface
    {

    public:

    /////////////
    // friends //
    /////////////

      friend class ConfigurationalBondTypes;

    //////////
    // Enum //
    //////////

      //! Data that can be retrieved from any bond type
      enum Data
      {
        e_Identity,                                    //!< 1 for any bond
        e_BondOrder,                                   //!< Nominal bond order = 1, 2, or 3
        e_NumberOfElectrons,                           //!< # of electrons in bond (2,4,6)
        e_Conjugation,                                 //!< Conjugation, size_t of the conjugation value of the bond (0-3)
        e_IsConjugated,                                //!< 1 for any bond that is conjugated (or aromatic), 0 otherwise
        e_IsAromatic,                                  //!< 1 for any bond that is aromatic, 0 otherwise
        e_IsAmide,                                     //!< 1 for amide bonds, 0 otherwise
        e_IsInRing,                                    //!< 1 for any bond that is in a ring, 0 otherwise
        e_BondOrderInRingOrAromatic,                   //!< 1-3 for bonds in rings, 4 for aromatic, 0 for all bonds outside rings
        e_BondOrderOrAromatic,                         //!< Nominal bond order (1-3), except aromatic types (4)
        e_BondOrderAmideOrAromatic,                    //!< Nominal bond order (1-3), amide (4), aromatic (5)
        e_BondOrderOrAromaticWithRingness,             //!< same as e_BondOrderOrAromatic, except +3 if the bond is in a ring
        e_BondOrderAmideOrAromaticWithRingness,        //!< same as e_BondOrderAmideOrAromatic, except +3 if the bond is in a ring
        e_FuzzyBondOrderOrAromaticWithRingness,        //!< same as e_BondOrderOrAromaticWithRingness, except double bonds in chains are equiv to single bonds in chains to allow comparison of resonance structures
        e_FuzzyBondOrderAmideOrAromaticWithRingness, //!< same as e_BondOrderOrAromaticWithRingnessAmide, except double bonds in chains are equiv to single bonds in chains to allow comparison of resonance structures
        e_ConstitutionalBondType,                      //!< index of this bond type in the constitutional bond types
        e_BondOrderWithIsometry,                       //!< Nominal bond order = 1, 2, or 3, 4 for E isometry, 5 for Z isometry
        e_Isometry,                                    //!< Isometry, size_t of the isometry value of the bond (Undef,0,1,2)
        e_IsIsometric,                                 //!< 1 for any bond that has E, Z, or unspecified isometry
        e_BondOrderWithIsometryOrAromatic,             //!< Nominal bond order (1-3), aromatic (4), E (5), Z (6)
        e_BondOrderAmideWithIsometryOrAromaticWithRingness, //!< Nominal bond order (1-3), amide (4) aromatic (5), E (6), Z (7)
        e_ConfigurationalBondType,                     //!< index of this bond type in the configurational bond types
        s_NumberOfData
      };

      //! @brief Data as string
      //! @param DATA the data whose name is desired
      //! @return the name as string
      static const std::string &GetDataName( const Data &DATA);

      //! DataEnum simplifies the usage of the Data enum of this class
      typedef util::WrapperEnum< Data, &GetDataName, s_NumberOfData> DataEnum;

    private:

    //////////
    // data //
    //////////

      BondIsometryEnum       m_Isometry;              //!< Isometry of the bond type
      ConstitutionalBondType m_BaseBondType;
      size_t                 m_Data[ s_NumberOfData]; //!< Calculated values for all data

      //! Corresponding bond type with same properties except different conjugation
      util::SiPtr< const ConfigurationalBondType> m_ConjugatedTypes[ ConstitutionalBondTypeData::s_NumberOfConjugations];

      //! Corresponding bond type inside a ring
      util::SiPtr< const ConfigurationalBondType> m_BondTypeInRing;

      //! Corresponding bond type with different bond order
      util::SiPtr< const ConfigurationalBondType> m_AlternateOrderBondType[ 3];

      //! Corresponding bond type without order, if it was an aromatic bond type
      util::SiPtr< const ConfigurationalBondType> m_BondTypeSansAromaticOrder;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined bond type
      ConfigurationalBondTypeData();

      //! @brief construct bond type from all its data
      //! @param ISOMETRY isometry of the bond type
      ConfigurationalBondTypeData
      (
        const ConstitutionalBondType &CONSTITUTIONAL_BOND_TYPE,
        const BondIsometry &ISOMETRY
      );

      //! @brief virtual copy constructor
      ConfigurationalBondTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the number of electrons in the bond
      //! @return the number of electrons in the bond
      size_t GetNumberOfElectrons() const
      {
        return m_BaseBondType->GetNumberOfElectrons();
      }

      //! @return the id of the bond type in an sd-file
      size_t GetSDFileID() const
      {
        return m_BaseBondType->GetSDFileID();
      }

      //! @return the alternative id of the bond type in an sd-file
      size_t GetSDAltFileID() const
      {
        return m_BaseBondType->GetSDAltFileID();
      }

      //! @brief get the conjugation
      //! @return Conjugation
      ConstitutionalBondTypeData::Conjugation GetConjugation() const
      {
        return m_BaseBondType->GetConjugation();
      }

      //! @brief test whether the # of electrons in the bond is known
      //! @return true if # of electrons in the bond is known, false otherwise
      bool IsBondOrderKnown() const
      {
        return m_BaseBondType->IsBondOrderKnown();
      }

      //! @brief test whether this bond is in a ring
      //! @return true if this bond is in a ring
      bool IsBondInRing() const
      {
        return m_BaseBondType->IsBondInRing();
      }

      //! @brief get the isometry
      //! @return Isometry
      BondIsometry GetIsometry() const
      {
        return m_Isometry;
      }

      //! @brief get the ConstitutionalBondType
      //! @return ConstitutionalBondType
      const ConstitutionalBondType &GetConstitutionalBondType() const
      {
        return m_BaseBondType;
      }

      //! @brief get the corresponding bond type with same properties except different conjugation
      //! @param CONJUGATION the alternative conjugation
      //! @return the corresponding bond type with same properties except different conjugation
      const ConfigurationalBondType &WithConjugation( const ConstitutionalBondTypeData::Conjugation &CONJUGATION) const;

      //! @brief get the corresponding bond type inside a ring
      //! @return the corresponding bond type inside a ring (this type if it is already in a ring)
      const ConfigurationalBondType &WithInRing() const;

      //! @brief get the corresponding bond type with a different order
      //! @param ORDER the desired order
      //! @return the corresponding bond type with the desired order
      const ConfigurationalBondType &WithOrder( const size_t &ORDER) const;

      //! @brief get the corresponding bond type with a different isometry
      //! @param ISOMETRY the desired isometry
      //! @return the corresponding bond type with the desired isometry
      const ConfigurationalBondType &WithIsometry( const BondIsometry &ISOMETRY) const;

      //! @brief get the corresponding bond type with e_NonIsometric; much shorter than calling
      //! Shorthand for WithIsometry( ConstitutionalBondTypeData::e_NonIsometric)
      //! @return the corresponding bond type with e_NonIsometric
      const ConfigurationalBondType &WithoutIsometry() const;

      //! @brief get the bond type, except for aromatic bonds, return the unordered type
      //! @return the bond type, except for aromatic bonds, return the unordered type
      const ConfigurationalBondType &WithoutAromaticOrder() const;

      //! @brief get the data member specified by DATA_TO_RETRIEVE
      //! @param DATA_TO_RETRIEVE the data to get from the bond
      //! @return the desired data
      size_t GetBondData( const Data &DATA_TO_RETRIEVE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief setup connections between this bond type and related types
      //! @param THIS_TYPE the corresponding configurational bond type
      //! This function can only be called by ConfigurationalBondTypes
      void AttachToRelatedBondTypes( const ConfigurationalBondType &THIS_TYPE);

    }; // class ConfigurationalBondTypeData

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_CONFIGURATIONAL_BOND_TYPE_DATA_H_

