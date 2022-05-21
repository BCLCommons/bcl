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

#ifndef BCL_CHEMISTRY_CONSTITUTIONAL_BOND_TYPE_DATA_H_
#define BCL_CHEMISTRY_CONSTITUTIONAL_BOND_TYPE_DATA_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_bond_isometry.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstitutionalBondTypeData
    //! @brief stores bond properties
    //! @details This class stores bond properties that are independent of bonded element types
    //!
    //! @see @link example_chemistry_constitutional_bond_type_data.cpp @endlink
    //! @author mendenjl, kothiwsk
    //! @date Dec 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConstitutionalBondTypeData :
      public util::ObjectInterface
    {

    public:

    /////////////
    // friends //
    /////////////

      friend class ConstitutionalBondTypes;  // calls AttachToRelatedBondTypes post-construction
      friend class ConfigurationalBondTypes; // calls SetConfigurationalType when configurational types are constructed

    //////////
    // Enum //
    //////////

      //! Conjugation of the bond type
      enum Conjugation
      {
        e_Nonconjugated,       //!< Not a conjugated bond
        e_Conjugated,          //!< A conjugated bond that is not aromatic
        e_Aromatic,            //!< A bond which is aromatic
        e_Amide,               //!< Partially-resonant conjugation (amides, sulfonamides, etc. where the resonance is shared unequally)
        e_Any,                 //!< A bond that might have any conjugation, used only in queries
        s_NumberOfConjugations
      };

      //! Data that can be retrieved from any bond type
      enum Data
      {
        e_Identity,                        //!< 1 for any bond
        e_BondOrder,                       //!< Nominal bond order = 1, 2, or 3
        e_NumberOfElectrons,               //!< # of electrons in bond (2,4,6)
        e_Conjugation,                     //!< Conjugation, size_t of the conjugation value of the bond (0-4)
        e_IsConjugated,                    //!< 1 for any bond that is conjugated (or aromatic/amide), 0 otherwise
        e_IsAromatic,                      //!< 1 for any bond that is aromatic, 0 otherwise
        e_IsAmide,                         //!< 1 for any bond that is aromatic, 0 otherwise
        e_IsInRing,                        //!< 1 for any bond that is in a ring, 0 otherwise
        e_BondOrderInRingOrAromatic,       //!< 1-3 for bonds in rings, 4 for aromatic, 0 for all bonds outside rings
        e_BondOrderOrAromatic,             //!< Nominal bond order (1-3), except aromatic types (4)
        e_BondOrderAmideOrAromatic,        //!< Nominal bond order (1-3), amide (4), aromatic (5)
        e_BondOrderOrAromaticWithRingness, //!< same as e_BondOrderOrAromatic, except +3 if the bond is in a ring
        e_BondOrderAmideOrAromaticWithRingness, //!< same as e_BondOrderAmideOrAromatic, except +3 if the bond is in a ring
        e_FuzzyBondOrderOrAromaticWithRingness, //!< same as e_BondOrderOrAromaticWithRingness, except double bonds in chains are equiv to single bonds in chains to allow comparison of resonance structures
        e_FuzzyBondOrderAmideOrAromaticWithRingness, //!< same as e_BondOrderOrAromaticWithRingnessAmide, except double bonds in chains are equiv to single bonds in chains to allow comparison of resonance structures
        e_ConstitutionalBondType,          //!< index of this bond type in the constitutional bond types
        s_NumberOfData
      };

      //! @brief Conjugation as string
      //! @param CONJUGATION the conjugation whose name is desired
      //! @return the name as string
      static const std::string &GetConjugationName( const Conjugation &CONJUGATION);

      //! ConjugationEnum simplifies the usage of the Conjugation enum of this class
      typedef util::WrapperEnum< Conjugation, &GetConjugationName, s_NumberOfConjugations> ConjugationEnum;

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

      size_t          m_SDFileID;              //!< The id of the bond type in an SD-File, 0 is the undefined value
      size_t          m_SDFileAltID;           //!< Alt id of the bond type in an SD-File, 0 is the undefined value
      ConjugationEnum m_Conjugation;           //!< Conjugation of the bond
      size_t          m_Data[ s_NumberOfData]; //!< Calculated values for all data

      //! Corresponding bond type with same properties except different conjugation
      util::SiPtr< const ConstitutionalBondType> m_ConjugatedTypes[ s_NumberOfConjugations];

      //! Corresponding bond type inside a ring
      util::SiPtr< const ConstitutionalBondType> m_BondTypeInRing;

      //! Corresponding bond type with different bond order
      util::SiPtr< const ConstitutionalBondType> m_AlternateOrderBondType[ 3];

      //! Corresponding configurational bond type with different isometries
      util::SiPtr< const ConfigurationalBondType> m_IsometryTypes[ s_NumberOfIsometries];

      //! Corresponding bond type without order, if it was an aromatic bond type
      util::SiPtr< const ConstitutionalBondType> m_BondTypeSansAromaticOrder;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined bond type
      ConstitutionalBondTypeData();

      //! @brief construct bond type from all its data
      //! @param NUMBER_ELECTRONS The total # of electrons involved in the bond
      //! @param SD_FILE_ID the bond id in an sdf file
      //! @param SD_FILE_ALT_ID for conjugated or aromatic types, the other # that can be used for the bond type
      //! @param CONJUGATION the type of conjugation the bond is involved with
      //! @param IS_IN_RING whether the bond is in a ring
      ConstitutionalBondTypeData
      (
        const size_t &NUMBER_ELECTRONS,
        const size_t &SD_FILE_ID,
        const size_t &SD_FILE_ALT_ID,
        const Conjugation &CONJUGATION,
        const bool &IS_IN_RING
      );

      //! @brief virtual copy constructor
      ConstitutionalBondTypeData *Clone() const;

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
        return GetBondData( e_NumberOfElectrons);
      }

      //! @return the id of the bond type in an sd-file
      size_t GetSDFileID() const
      {
        return m_SDFileID;
      }

      //! @return the alternative id of the bond type in an sd-file
      size_t GetSDAltFileID() const
      {
        return m_SDFileAltID;
      }

      //! @brief get the conjugation
      //! @return Conjugation
      Conjugation GetConjugation() const
      {
        return m_Conjugation;
      }

      //! @brief test whether the # of electrons in the bond is known
      //! @return true if # of electrons in the bond is known, false otherwise
      bool IsBondOrderKnown() const
      {
        return GetBondData( e_NumberOfElectrons) != 3;
      }

      //! @brief test whether this bond is in a ring
      //! @return true if this bond is in a ring
      bool IsBondInRing() const
      {
        return GetBondData( e_IsInRing);
      }

      //! @brief get the data member specified by DATA_TO_RETRIEVE
      //! @param DATA_TO_RETRIEVE the data to get from the bond
      //! @return the desired data
      size_t GetBondData( const Data &DATA_TO_RETRIEVE) const;

      //! @brief get the corresponding bond type with same properties except different conjugation
      //! @param CONJUGATION the alternative conjugation
      //! @return the corresponding bond type with same properties except different conjugation
      const ConstitutionalBondType &WithConjugation( const Conjugation &CONJUGATION) const;

      //! @brief get the corresponding bond type inside a ring
      //! @return the corresponding bond type inside a ring (this type if it is already in a ring)
      const ConstitutionalBondType &WithInRing() const;

      //! @brief get the corresponding bond type with a different order
      //! @param ORDER the desired order
      //! @return the corresponding bond type with the desired order
      const ConstitutionalBondType &WithOrder( const size_t &ORDER) const;

      //! @brief get the corresponding configurational bond type with a certain isometry
      //! @param ISOMETRY the desired isometry
      //! @return the corresponding configurational bond type with a certain isometry
      const ConfigurationalBondType &WithIsometry( const BondIsometry &ISOMETRY) const;

      //! @brief get the bond type, except for aromatic bonds, return the unordered type
      //! @return the bond type, except for aromatic bonds, return the unordered type
      const ConstitutionalBondType &WithoutAromaticOrder() const;

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

      //! @brief set the index of this bond in the constitutional bond types, and setup related types
      //! @param TYPES the enum; needed to prevent a circular dependency
      //! This function can only be called by ConstitutionalBondTypes
      void AttachToRelatedBondTypes( const ConstitutionalBondTypes &TYPES);

      //! @brief set the related configurational bond type
      //! @param TYPE the corresponding configurational bond type, without isometry
      //! This function can only be called by ConfigurationalBondTypes
      void SetConfigurationalBondType( const ConfigurationalBondType &TYPE);

      //! @brief initialize all isometries with a particular type
      //! @param TYPE the undefined configurational bond type
      //! This function can only be called by ConfigurationalBondTypes
      void InitializeIsometries( const ConfigurationalBondType &TYPE);

    }; // class ConstitutionalBondTypeData

    //! @brief Output operator for chemistry::ConstitutionalBondTypeData::Conjugation
    //! @param OSTREAM stream to write output for
    //! @param CONJUGATION property to write
    //! @return ostream which was written to
    BCL_API std::ostream &operator <<
    (
      std::ostream &OSTREAM,
      const ConstitutionalBondTypeData::Conjugation &CONJUGATION
    );

    //! @brief Input operator for chemistry::ConstitutionalBondTypeData::Conjugation
    //! @param ISTREAM stream to read in from the stream
    //! @param CONJUGATION what property to read
    //! @return istream which was read from
    BCL_API std::istream &operator >>
    (
      std::istream &ISTREAM,
      ConstitutionalBondTypeData::Conjugation &CONJUGATION
    );

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_CONSTITUTIONAL_BOND_TYPE_DATA_H_

