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

#ifndef BCL_CHEMISTRY_ATOM_TYPE_DATA_H_
#define BCL_CHEMISTRY_ATOM_TYPE_DATA_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atomic_orbital_types.h"
#include "bcl_chemistry_element_types.h"
#include "bcl_chemistry_hybrid_orbital_types.h"
#include "math/bcl_math_polynomial.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomTypeData
    //! @brief contains hybridization and bond geometry data, which is used in Atom
    //!
    //! @see @link example_chemistry_atom_type_data.cpp @endlink
    //! @author mueller, woetzen, mendenjl
    //! @date 08/23/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomTypeData :
      public util::ObjectInterface
    {
      friend class AtomTypes;
      friend class BondLengths;

    public:

    ///////////
    // Enums //
    ///////////

      //! enum properties for atom types
      // see Hinze, Jaffe: Electronegativity. I. Orbital Electronegativity of Neutral Atoms
      // each atom type has up to two different parameter sets for sigma- and pi-orbitals
      enum Properties
      {
        e_SigmaValenceStateIonizationPotential,    //!< sigma ionization potential
        e_SigmaValenceStateElectronAffinity,       //!< sigma electron affinity
        e_SigmaOrbitalElectronegativityMulliken,   //!< sigma orbital electronegativity, Mulliken scale
        e_SigmaOrbitalElectronegativityPauling,    //!< sigma orbital electronegativity, Pauli scale
        e_PiValenceStateIonizationPotential,       //!< pi ionization potential
        e_PiValenceStateElectronAffinity,          //!< pi electron affinity
        e_PiOrbitalElectronegativityMulliken,      //!< pi orbital electronegativity, Mulliken scale
        e_PiOrbitalElectronegativityPauling,       //!< pi orbital electronegativity, Pauli scale
        e_LonePairIonizationPotential,             //!< lone pair ionization potential
        e_LonePairElectronAffinity,                //!< lone pair electron affinity
        e_LonePairElectronegativityMulliken,       //!< lone pair electronegativity
        e_AdditiveAtomicPolarizability,            //!< additive atomic polarizability, see J. Am. Chem. Soc. 1990, 112, 8533-8542
        e_VdWaalsRadiusCSD,                        //!< Vdw radius for the atom type using data from the CSD; see notes below
        e_CovalentRadiusSingleBond,                //!< Covalent radius for single bond
        e_CovalentRadiusDoubleBond,                //!< Covalent radius for double bond
        e_CovalentRadiusTripleBond,                //!< Covalent radius for triple bond
        e_CovalentRadiusAromaticBond,              //!< Covalent radius for an aromatic bond
        e_CovalentRadiusTypical,                   //!< Nominal average covalent radius
        s_NumberOfProperties                       //!< Number of properties
      };

      //! @note e_VdWaalsRadiusCSD is only valid for atoms without hydrogens, and can only be used to
      //! @note detect bad geometries if both atoms have no H (otherwise H bonding may lead to close contacts, also
      //! @note could imply missing bonds), and which do not have opposite charges
      //! @note e_VdWaalsRadiusCSD was calculated using data from the cambridge structural database
      //! @note Vdw radii may be violated in certain bridged ring systems

      //! @brief element type property as string
      //! @param PROPERTY the property desired
      //! @return the property as string
      static const std::string &GetPropertyName( const Properties &PROPERTY);

      //! PropertyEnum simplifies the usage of the Properties enum of this class
      typedef util::WrapperEnum< Properties, &GetPropertyName, s_NumberOfProperties> PropertyEnum;

      //! how the atom type can contribute to a pi system
      enum PiContributionType
      {
        e_Zero,       //!< no pi electrons will be contributed; e.g. atom type with no non-singular bonds, no lone pairs, like B_TrTrTr
        e_One,        //!< exactly one pi electron will be contributed by this atom type (atom type with a double bond)
        e_Two,        //!< exactly two pi electrons will be contributed by this atom type (atom type with a triple bond or two 2x bonds)
        e_ZeroOrTwo,  //!< Either zero or two pi electrons will be contributed by this atom type (atom type with lone pairs but only singular bonds)
      };

    private:

    //////////
    // data //
    //////////

      ElementType m_ElementType;                 //!< element type

      HybridOrbitalType m_Hybridization;         //!< Unhybridized, SP, SP2, or SP3

      size_t m_NumberHybridOrbitalsSigmaBinding; //!< hybrid orbitals that are binding
      size_t m_NumberHybridOrbitalsNonbinding;   //!< hybrid orbitals that are non binding with their electrons
      size_t m_NumberElectronsInBonds;           //!< Number of e- bonds (calculated)
      size_t m_NumberBonds;                      //!< Number of bonds (calculated)
      size_t m_MaxEContributionToPiSystem;       //!< Number of e- in pi system (calculated)

      storage::Set< AtomicOrbitalTypesEnum> m_PiOrbitalsBinding; //!< pi orbitals involved in binding

      //! non hybridized non binding orbitals with their electrons
      storage::Set< AtomicOrbitalTypesEnum> m_AtomicOrbitalsNonbinding;

      double m_Properties[ int( s_NumberOfProperties)]; //!< real-valued properties
      double m_OrbitalENegPos; //!< orbital electronegativity associated with the charged state

       //! estimated stability of the atom type; only used to resolve clashes in atom types
      double m_StabilityMetric;

      //! charge of atom type = difference between electrons in bonds and non-bonds minus valence electrons
      short m_Charge;

      //! whether this atom is conjugated
      bool m_Conjugated;

      math::Polynomial m_SigmaChargeToEN; //!< sigma electronegativity = f(SigmaCharge)
      math::Polynomial m_PiChargeToEN;    //!< pi electronegativity = f(PiCharge)
      math::Polynomial m_SigmaChargeToLonePairEN; //!< lone pair electronegativity = f(SigmaCharge)

      // inflection points are needed for the above charge equations for the rare cases that the charge can end up on
      // the other side of the quadratic polynomial curve, in which case charges run off towards infinity
      // To avoid this problem, if the desired charge falls below the inflection point of the above quadratic or linear
      // functions, we use the value at the inflection point
      double m_SigmaChargeToENInflection;         //!< Inflection point for sigma charge to inflection
      double m_PiChargeToENInflection;            //!< Inflection point for pi charge to inflection
      double m_SigmaChargeToLonePairENInflection; //!< Inflection point for sigma charge to inflection

      bool m_IsGasteigerType; //!< Is the type a proper gasteiger type
      // before SmallMoleculeStandardizer is used on a molecule all atoms in it will not have a proper gasteiger type
      // after SmallMoleculeStandardizer is used, all atoms in it will have have a proper gasteiger atom type, provided
      // one exists

      PiContributionType m_PiContribution; //!< Contribution of this atom type to a pi-conjugated system

      //! bool for whether the atom type always forms nearly linear bonds
      bool m_FormsOnlyLinearBonds;

      //! Two letter code, used for atom interaction scores only at present
      std::string m_TwoLetterCode;

      //! Base atom type; only contains element type and charge
      util::SiPtr< const AtomType> m_BaseType;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AtomTypeData();

      //! @brief the usual constructor
      AtomTypeData
      (
        const ElementType &ELEMENT_TYPE,
        const HybridOrbitalType &HYBRIDIZATION,
        const size_t &HYBRID_ORBITALS_IN_SIGMA_BONDS,
        const size_t &HYBRID_ORBITALS_NONBINDING,
        const storage::Set< AtomicOrbitalTypesEnum> &PI_ORBITALS_IN_BONDS,
        const storage::Set< AtomicOrbitalTypesEnum> &ATOMIC_ORBITALS_NONBINDING,
        const double &SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL,
        const double &SIGMA_VALENCE_STATE_ELECTRON_AFFINITY,
        const double &PI_VALENCE_STATE_IONIZATION_POTENTIAL,
        const double &PI_VALENCE_STATE_ELECTRON_AFFINITY,
        const double &LONE_PAIR_IONIZATION_POTENTIAL,
        const double &LONE_PAIR_ELECTRON_AFFINITY,
        const double &ATOMIC_POLARIZABILITY
      );

      //! @brief constructor from just an element type and charge
      AtomTypeData
      (
        const ElementType &ELEMENT_TYPE,
        const short &CHARGE
      );

      //! @brief Clone function
      //! @return pointer to new AtomTypeData
      AtomTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! return ElementType
      const ElementType &GetElementType() const;

      //! @brief returns the hybridization of the atom type
      //! @return the type of hybrid orbital
      const HybridOrbitalType &GetHybridOrbitalType() const;

      //! return two letter code that indicates element name, and, for the primary organic elements BCNS, is followed by the
      //! number of bonds provided that number is at least 3, otherwise the letter X is used instead.
      const std::string &GetTwoLetterCode() const;

      //! @brief returns the number of hybridized orbitals
      //! @return the number of hybridized orbitals
      size_t GetNumberHybridOrbitals() const;

      //! @brief returns the number of lone pairs in hybrid orbitals
      //! @return the number of lone pairs in hybrid orbitals
      size_t GetNumberHybridLonePairs() const;

      //! @return Number of hybridized bonds
      size_t GetNumberHybridBonds() const;

      //! @return Number of bonds
      size_t GetNumberBonds() const;

      //! @return Number of Sigma orbitals that are not hybridized
      size_t GetNumberUnhybridizedSigmaOrbitals() const;

      //! @return Number of Sigma orbitals
      size_t GetNumberSigmaOrbitals() const;

      //! @return Number of electrons in p orbitals (whether hybridized or not)
      size_t GetNumberElectronsInPOrbitals() const;

      //! @return Number of pi-orbitals
      size_t GetNumberPiOrbitals() const;

      //! @return Charge
      const short &GetFormalCharge() const;

      //! @return valence electrons in sp orbitals
      size_t GetValenceElectronsSP() const;

      //! @brief GetNumberElectronsInBonds calculates the total number of electrons in pi-orbital and sigma bonds
      size_t GetNumberElectronsInBonds() const;

      //! @brief return Number of unhybridized lone pairs
      size_t GetNumberUnhybridizedLonePairs() const;

      //! @brief bool for whether the atom type always forms nearly linear bonds
      bool GetFormsOnlyLinearBonds() const;

      //! @brief atom type property as double
      //! @param PROPERTY the property desired
      //! @return the property as double
      double GetAtomTypeProperty( const AtomTypeData::Properties &PROPERTY) const;

      //! @brief given a partial sigma charge, determine the electronegativity
      //! @param CHARGE the partial sigma charge
      //! @return the resulting electronegativity
      double GetSigmaENFromCharge( const double &CHARGE) const;

      //! @brief given a partial pi charge, determine the electronegativity
      //! @param CHARGE the partial pi charge
      //! @return the resulting electronegativity
      double GetPiENFromCharge( const double &CHARGE) const;

      //! @brief given a partial sigma charge, determine the lone pair electronegativity
      //! @param CHARGE the partial sigma charge
      //! @return the resulting electronegativity
      double GetLonePairENFromSigmaCharge( const double &CHARGE) const;

      //! @brief get the function used to compute GetSigmaENFromCharge
      //! @return the function used to compute GetSigmaENFromCharge
      const math::Polynomial &GetSigmaENFromChargeFunction() const;

      //! @brief get the function used to compute GetPiENFromCharge
      //! @return the function used to compute GetPiENFromCharge
      const math::Polynomial &GetPiENFromChargeFunction() const;

      //! @brief get the function used to compute GetLonePairENFromCharge
      //! @return the function used to compute GetLonePairENFromCharge
      const math::Polynomial &GetLonePairENFromChargeFunction() const;

      //! @return the orbital electronegativity associated with the charged state
      double GetOrbitalENegPos() const;

      //! @brief determine if this atom type can participate in pi-bond conjugation
      //! @return true iff this atom type has any non-single bonds or lone pairs
      bool IsConjugated() const;

      //! @brief is this a well characterized gasteiger atom type
      //! @return true iff this atom type is this a well characterized gasteiger atom type
      bool IsGasteigerAtomType() const;

      //! @brief get the base atom type
      //! @return the atom type with only element type and charge information
      const AtomType &GetBaseAtomType() const;

      //! @brief Get the max number of electrons available for contribution to an aromatic ring
      //! @return the max electrons contributed by this atom type to a pi system
      size_t GetMaxEContributionToPiSystem() const;

      //! @brief Get the type of contribution this atom type can make to a pi system
      //! @return the type of contribution this atom type can make to a pi system
      PiContributionType GetPiElectronContributionType() const;

      //! @brief Get the stability metric.  Electronic stability is indicated by a larger number
      //! This is used to decide between atom types when no other means are possible
      double GetStabilityMetric() const;

    ////////////////
    // operations //
    ////////////////

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
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! Type difference specifies the difference between two atom types
      enum TypeDifference
      {
        e_None,                    //!< atom types are identical
        e_NumberBondingSOrbitals,  //!< different # of electrons in bonding s-orbitals
        e_NumberBondingPOrbitals,  //!< different # of electrons in bonding p-orbitals
        e_NumberLonePairOrbitals,  //!< different # of electrons in lone pairs
        e_Other,                   //!< multiple differences or other differences (e.g. element type),
        s_NumberTypeDifferences
      };

      //! @brief calculate the stability metric.  Electronic stability is indicated by a larger number
      //! This is used to decide between atom types when no other means are possible
      double CalculateStabilityMetric() const;

      //! @brief Estimate electronegativities of the anionic, cationic, or ground state given values from at least one other state
      //! @param ELECTRONEGATIVITIES electronegativities for the charge = -1,0, and +1 state
      //! @param TYPE_DIFFERENCE what orbital the electronegativity relates to
      void EstimateUndefinedElectronegativities
      (
        linal::Vector< double> &ELECTRONEGATIVITIES,
        const TypeDifference &TYPE_DIFFERENCE
      );

      //! @brief set orbital electronegativity for the charged species
      //! @param TYPE_DIFFERENCE whether the electronegativity is S/P Bonding/LonePair orbitals
      //! @param ELECTRONEGATIVITY electronegativity of the charged or neutral species, if available
      //! Used only in AtomTypes constructor
      void SetSimilarOrbitalElectronegativity
      (
        const TypeDifference &TYPE_DIFFERENCE,
        const double &ELECTRONEGATIVITY = util::GetUndefinedDouble()
      );

      //! @brief set orbital electronegativity from the neutral type
      //! @param TYPE_DIFFERENCE whether the electronegativity is S/P Bonding/LonePair orbitals
      //! @param NEUTRAL_TYPE neutral type of the atom, if available
      //! Used only in AtomTypes constructor
      void SetSimilarOrbitalElectronegativity
      (
        const TypeDifference &TYPE_DIFFERENCE,
        const AtomTypeData &NEUTRAL_TYPE
      );

      //! @brief get the average ionization potential ratio between cation and neutral atom type that differ by TYPE_DIFFERENCE
      //! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
      //! @return the ratio
      double GetAverageIPChangeCationToNeutral( const TypeDifference &TYPE_DIFFERENCE) const;

      //! @brief get the average ionization potential ratio between neutral and cation atom type that differ by TYPE_DIFFERENCE
      //! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
      //! @return the ratio
      double GetAverageIPChangeNeutralToAnion( const TypeDifference &TYPE_DIFFERENCE) const;

      //! @brief type difference as string
      //! @param TYPE_DIFFERENCE the type difference for which a string is desired
      //! @return the type difference as a string
      static const std::string &GetTypeDifferenceName( const TypeDifference &TYPE_DIFFERENCE);

      //! @brief determine the difference betweent his atom type data and another
      //! @param OTHER the atom type data to compare this atom type data to
      //! @return the corresponding TypeDifference
      TypeDifference DifferenceFrom( const AtomTypeData &OTHER);

      //! @brief get the electronegativity type corresponding to a TypeDifference
      //! @param TYPE_DIFFERENCE the type difference to get the corresponding electronegativity for
      //! @return the electronegativity type corresponding to TypeDifference
      double GetElectronegativity( const TypeDifference &TYPE_DIFFERENCE) const;

      //! @brief get the ionization potential type corresponding to a TypeDifference
      //! @param TYPE_DIFFERENCE the type difference to get the corresponding ionization potential for
      //! @return the ionization potential type corresponding to TypeDifference
      double GetIonizationPotential( const TypeDifference &TYPE_DIFFERENCE) const;

      //! @brief get the electron affinity type corresponding to a TypeDifference
      //! @param TYPE_DIFFERENCE the type difference to get the corresponding electron affinity for
      //! @return the electron affinity type corresponding to TypeDifference
      double GetElectronAffinity( const TypeDifference &TYPE_DIFFERENCE) const;

      //! @brief get the electronegativity from charge corresponding to a TypeDifference
      //! @param TYPE_DIFFERENCE the type difference to get the corresponding function for
      //! @return the electronegativity as a function of charge polynomial corresponding to TypeDifference
      math::Polynomial &GetElectronegativityFromChargeFunction( const TypeDifference &TYPE_DIFFERENCE);

      //! @brief set a particular data
      //! @param DATA the property to set
      //! @param VALUE the value to set the property to
      void SetProperty( const Properties &DATA, const double &VALUE);

      //! @brief set that this atom type forms only linear bonds
      void SetFormsOnlyLinearBonds();

      //! @brief GetAverageNeutralSigmaIVToEARatio helper function for AtomTypes::CalculateElectronegativityValues
      //! @return reference to a double, which returns the ratio of Average(e_SigmaValenceStateIonizationPotential) for neutral atoms vs. anions
      static double &GetAverageNeutralSigmaIPToAnionIPRatio();

      //! @brief GetAverageNeutralPiIPToAnionIPRatio helper function for AtomTypes::CalculateElectronegativityValues
      //! @return reference to a double, which returns the ratio of Average(e_PiValenceStateIonizationPotential) for neutral atoms vs. anions
      static double &GetAverageNeutralPiIPToAnionIPRatio();

      //! @brief GetAverageCationSigmaIPToNeutralIPRatio helper function for AtomTypes::CalculateElectronegativityValues
      //! @return reference to a double, which returns the ratio of Average(e_SigmaValenceStateIonizationPotential) for cations vs. neutral atoms
      static double &GetAverageCationSigmaIPToNeutralIPRatio();

      //! @brief GetAverageCationPiIPToNeutralIPRatio helper function for AtomTypes::CalculateElectronegativityValues
      //! @return reference to a double, which returns the ratio of Average(e_PiValenceStateIonizationPotential) for cations vs. neutral atoms
      static double &GetAverageCationPiIPToNeutralIPRatio();

    }; // class AtomTypeData

//    //! @brief Output operator for ATOM_TYPE_PROPERTY
//    //! @param OSTREAM stream to write output for
//    //! @param ATOM_TYPE_PROPERTY property to write
//    //! @return ostream which was written to
//    BCL_API
//    std::ostream &operator <<( std::ostream &OSTREAM, const AtomTypeData::Properties &ATOM_TYPE_PROPERTY);
//
//    //! @brief Input operator for ATOM_TYPE_PROPERTY
//    //! @param ISTREAM stream to read in from the stream
//    //! @param ATOM_TYPE_PROPERTY what property to read
//    //! @return istream which was read from
//    BCL_API
//    std::istream &operator >>( std::istream &ISTREAM, AtomTypeData::Properties &ATOM_TYPE_PROPERTY);

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_TYPE_DATA_H_

