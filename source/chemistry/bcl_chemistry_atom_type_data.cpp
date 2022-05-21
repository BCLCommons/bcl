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

#define BCL_ENUM_DECLARATION_ONLY

// include header of this class
#include "chemistry/bcl_chemistry_atom_type_data.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_types.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_message.h"

#undef BCL_ENUM_DECLARATION_ONLY
#include "util/bcl_util_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AtomTypeData::s_Instance
    (
      GetObjectInstances().AddInstance( new AtomTypeData())
    );

  ///////////
  // Enums //
  ///////////

    //! @brief element type property as string
    //! @param PROPERTY the property desired
    //! @return the property as string
    const std::string &AtomTypeData::GetPropertyName( const AtomTypeData::Properties &PROPERTY)
    {
      static const std::string s_properties[] =
      {
        "SigmaValenceStateIonizationPotential",
        "SigmaValenceStateElectronAffinity",
        "SigmaOrbitalElectronegativityMulliken",
        "SigmaOrbitalElectronegativityPauling",
        "PiValenceStateIonizationPotential",
        "PiValenceStateElectronAffinity",
        "PiOrbitalElectronegativityMulliken",
        "PiOrbitalElectronegativityPauling",
        "LonePairIonizationPotential",
        "LonePairElectronAffinity",
        "LonePairElectronegativity",
        "AdditiveAtomicPolarizability",
        "VdWaalsRadiusCSD",
        "CovalentRadiusSingleBond",
        "CovalentRadiusDoubleBond",
        "CovalentRadiusTripleBond",
        "CovalentRadiusAromaticBond",
        "CovalentRadiusAverage",
        GetStaticClassName< Properties>()
      };

      return s_properties[ PROPERTY];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomTypeData::AtomTypeData() :
      m_NumberHybridOrbitalsSigmaBinding( 0),
      m_NumberHybridOrbitalsNonbinding( 0),
      m_NumberElectronsInBonds( 0),
      m_NumberBonds( 0),
      m_MaxEContributionToPiSystem( 0),
      m_OrbitalENegPos( 0),
      m_StabilityMetric( 0),
      m_Charge( 0),
      m_Conjugated( false),
      m_SigmaChargeToENInflection( -std::numeric_limits< double>::infinity()),
      m_PiChargeToENInflection( -std::numeric_limits< double>::infinity()),
      m_SigmaChargeToLonePairENInflection( -std::numeric_limits< double>::infinity()),
      m_IsGasteigerType( false),
      m_PiContribution( e_Zero),
      m_FormsOnlyLinearBonds( false)
    {
      // make all the properties undefined
      for( int property_number( 0); property_number < s_NumberOfProperties; ++property_number)
      {
        m_Properties[ property_number] = util::GetUndefined< double>();
      }
    }

    //! @brief the usual constructor
    AtomTypeData::AtomTypeData
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
    ) :
      m_ElementType( ELEMENT_TYPE),
      m_Hybridization( HYBRIDIZATION),
      m_NumberHybridOrbitalsSigmaBinding( HYBRID_ORBITALS_IN_SIGMA_BONDS),
      m_NumberHybridOrbitalsNonbinding( HYBRID_ORBITALS_NONBINDING),
      m_PiOrbitalsBinding( PI_ORBITALS_IN_BONDS),
      m_AtomicOrbitalsNonbinding( ATOMIC_ORBITALS_NONBINDING),
      m_OrbitalENegPos( 0),
      m_StabilityMetric( 0),
      m_Charge
      (
        m_ElementType->GetElectronConfiguration().ValenceElectronsSP()
        - HYBRID_ORBITALS_IN_SIGMA_BONDS - m_PiOrbitalsBinding.GetSize()
        - 2 * ( HYBRID_ORBITALS_NONBINDING + ATOMIC_ORBITALS_NONBINDING.GetSize())
      ),
      m_SigmaChargeToENInflection( -std::numeric_limits< double>::infinity()),
      m_PiChargeToENInflection( -std::numeric_limits< double>::infinity()),
      m_SigmaChargeToLonePairENInflection( -std::numeric_limits< double>::infinity()),
      m_IsGasteigerType( true),
      m_FormsOnlyLinearBonds( false)
    {
      m_Properties[ e_SigmaValenceStateIonizationPotential] = SIGMA_VALENCE_STATE_IONIZATION_POTENTIAL;
      m_Properties[ e_SigmaValenceStateElectronAffinity] = SIGMA_VALENCE_STATE_ELECTRON_AFFINITY;
      m_Properties[ e_PiValenceStateIonizationPotential] = PI_VALENCE_STATE_IONIZATION_POTENTIAL;
      m_Properties[ e_PiValenceStateElectronAffinity] = PI_VALENCE_STATE_ELECTRON_AFFINITY;
      m_Properties[ e_LonePairIonizationPotential] = LONE_PAIR_IONIZATION_POTENTIAL;
      m_Properties[ e_LonePairElectronAffinity] = LONE_PAIR_ELECTRON_AFFINITY;
      m_Properties[ e_AdditiveAtomicPolarizability] = ATOMIC_POLARIZABILITY;
      // make all remaining properties undefined
      for( int property_number( e_VdWaalsRadiusCSD); property_number < s_NumberOfProperties; ++property_number)
      {
        m_Properties[ property_number] = util::GetUndefined< double>();
      }

      // calculate and store the mulliken electronegativities, (ElectronAffinity+IonizationPotential)/2
      m_Properties[ e_PiOrbitalElectronegativityMulliken]
        = 0.5 * ( m_Properties[ e_PiValenceStateElectronAffinity] + m_Properties[ e_PiValenceStateIonizationPotential]);
      m_Properties[ e_SigmaOrbitalElectronegativityMulliken]
        = 0.5 * (
                  m_Properties[ e_SigmaValenceStateElectronAffinity]
                  + m_Properties[ e_SigmaValenceStateIonizationPotential]
                );
      m_Properties[ e_LonePairElectronegativityMulliken]
        = 0.5 * (
                  3.0 * m_Properties[ e_LonePairElectronAffinity]
                  - m_Properties[ e_LonePairIonizationPotential]
                );

      // convert from the mulliken electronegativity scale to pauling, e.g. see Hinze, Jaffe, 1963
      m_Properties[ e_SigmaOrbitalElectronegativityPauling]
        = 0.336 * ( m_Properties[ e_SigmaOrbitalElectronegativityMulliken] - 0.615);
      m_Properties[ e_PiOrbitalElectronegativityPauling]
        = 0.336 * ( m_Properties[ e_PiOrbitalElectronegativityMulliken] - 0.615);

      m_NumberElectronsInBonds = m_NumberHybridOrbitalsSigmaBinding + m_PiOrbitalsBinding.GetSize();
      m_NumberBonds = GetNumberHybridBonds() == 0 ? m_NumberElectronsInBonds : GetNumberHybridBonds();
      if( GetNumberPiOrbitals() > GetNumberUnhybridizedSigmaOrbitals())
      {
        // case where there is at least one double or triple bond
        // then the # contributed is exactly = this expression since the electrons are directly on the bonds
        // in the ring
        m_MaxEContributionToPiSystem =
          std::min( size_t( 2), GetNumberPiOrbitals() - GetNumberUnhybridizedSigmaOrbitals());
      }
      else
      {
        // all singular bonds; contribution to pi system is 2 if there are any lone pairs, 0 otherwise
        m_MaxEContributionToPiSystem =
          2 * std::min( GetNumberUnhybridizedLonePairs() + GetNumberHybridLonePairs(), size_t( 1));
      }
      m_StabilityMetric = CalculateStabilityMetric();

      m_Conjugated =
        GetNumberPiOrbitals() > GetNumberUnhybridizedSigmaOrbitals()
         || !m_AtomicOrbitalsNonbinding.IsEmpty()
         || m_NumberHybridOrbitalsNonbinding
         || m_Hybridization == GetHybridOrbitalTypes().e_SP2;
      m_PiContribution = PiContributionType( m_MaxEContributionToPiSystem);
      if( m_MaxEContributionToPiSystem == size_t( 2) && m_NumberElectronsInBonds == m_NumberBonds)
      {
        // lone pair, no pi electrons in bonds, so the lone pair electrons are either all in or all out of the aromatic
        // system
        m_PiContribution = e_ZeroOrTwo;
      }
      m_TwoLetterCode = m_ElementType->GetChemicalSymbol();
      m_TwoLetterCode.resize( size_t( 2), 'X');
      if( m_TwoLetterCode[ 1] == 'X')
      {
        if( m_NumberBonds > size_t( 2))
        {
          if
          (
            m_ElementType == GetElementTypes().e_Boron
            || m_ElementType == GetElementTypes().e_Carbon
            || m_ElementType == GetElementTypes().e_Nitrogen
            || m_ElementType == GetElementTypes().e_Sulfur
          )
          {
            m_TwoLetterCode[ 1] = char( m_NumberBonds + '0');
          }
        }
        else if
        (
          m_NumberBonds == size_t( 2)
          && ( m_ElementType == GetElementTypes().e_Boron || m_ElementType == GetElementTypes().e_Carbon)
        )
        {
          m_TwoLetterCode[ 1] = '2';
        }
      }
    }

    //! @brief constructor from just an element type and charge
    AtomTypeData::AtomTypeData
    (
      const ElementType &ELEMENT_TYPE,
      const short &CHARGE
    ) :
      m_ElementType( ELEMENT_TYPE),
      m_NumberHybridOrbitalsSigmaBinding( util::GetUndefined< size_t>()),
      m_NumberHybridOrbitalsNonbinding( util::GetUndefined< size_t>()),
      m_NumberElectronsInBonds( 0),
      m_NumberBonds( 0),
      m_MaxEContributionToPiSystem( 0),
      m_OrbitalENegPos( 0),
      m_StabilityMetric( 0),
      m_Charge( CHARGE),
      m_Conjugated( false),
      m_SigmaChargeToENInflection( -std::numeric_limits< double>::infinity()),
      m_PiChargeToENInflection( -std::numeric_limits< double>::infinity()),
      m_SigmaChargeToLonePairENInflection( -std::numeric_limits< double>::infinity()),
      m_IsGasteigerType( false),
      m_PiContribution( e_Zero),
      m_FormsOnlyLinearBonds( false)
    {
      // make all the properties undefined
      for( int property_number( 0); property_number < s_NumberOfProperties; ++property_number)
      {
        m_Properties[ property_number] = util::GetUndefined< double>();
      }
    }

    //! @brief Clone function
    //! @return pointer to new AtomTypeData
    AtomTypeData *AtomTypeData::Clone() const
    {
      return new AtomTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AtomTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return ElementType
    const ElementType &AtomTypeData::GetElementType() const
    {
      return m_ElementType;
    }

    //! @brief returns the hybridization of the atom type
    //! @return the type of hybrid orbital
    const HybridOrbitalType &AtomTypeData::GetHybridOrbitalType() const
    {
      return m_Hybridization;
    }

    //! return two letter code that indicates element name, and, for the primary organic elements BCNS, is followed by the
    //! number of bonds provided that number is at least 3, otherwise the letter X is used instead.
    const std::string &AtomTypeData::GetTwoLetterCode() const
    {
      return m_TwoLetterCode;
    }

    //! @brief returns the number of hybridized orbitals
    //! @return the number of hybridized orbitals
    size_t AtomTypeData::GetNumberHybridOrbitals() const
    {
      return m_NumberHybridOrbitalsSigmaBinding + m_NumberHybridOrbitalsNonbinding;
    }

    //! @brief returns the number of lone pairs in hybrid orbitals
    //! @return the number of lone pairs in hybrid orbitals
    size_t AtomTypeData::GetNumberHybridLonePairs() const
    {
      return m_NumberHybridOrbitalsNonbinding;
    }

    //! return Number of bonds
    size_t AtomTypeData::GetNumberHybridBonds() const
    {
      return m_NumberHybridOrbitalsSigmaBinding == 0
             ? m_PiOrbitalsBinding.GetSize()
             : m_NumberHybridOrbitalsSigmaBinding;
    }

    //! return Number of Sigma orbitals that are not hybridized
    size_t AtomTypeData::GetNumberUnhybridizedSigmaOrbitals() const
    {
      return size_t( m_PiOrbitalsBinding.Contains( e_S));
    }

    //! return Number of Sigma orbitals
    size_t AtomTypeData::GetNumberSigmaOrbitals() const
    {
      // start with the number of pi orbitals binding
      return GetNumberUnhybridizedSigmaOrbitals() + m_NumberHybridOrbitalsSigmaBinding;
    }

    //! @brief return Number of unhybridized lone pairs
    size_t AtomTypeData::GetNumberUnhybridizedLonePairs() const
    {
      return m_AtomicOrbitalsNonbinding.GetSize();
    }

    //! return Number of occupied p orbitals
    size_t AtomTypeData::GetNumberElectronsInPOrbitals() const
    {
      return m_NumberElectronsInBonds - GetNumberUnhybridizedSigmaOrbitals();
    }

    //! return Number of pi-orbitals
    size_t AtomTypeData::GetNumberPiOrbitals() const
    {
      return m_PiOrbitalsBinding.GetSize();
    }

    //! return Charge
    const short &AtomTypeData::GetFormalCharge() const
    {
      return m_Charge;
    }

    //! @brief GetNumberElectronsInBonds calculates the total number of electrons in pi-orbital and sigma bonds
    size_t AtomTypeData::GetNumberElectronsInBonds() const
    {
      return m_NumberElectronsInBonds;
    }

    //! @brief bool for whether the atom type always forms nearly linear bonds
    bool AtomTypeData::GetFormsOnlyLinearBonds() const
    {
      return m_FormsOnlyLinearBonds;
    }

    //! @brief atom type property as double
    //! @param PROPERTY the property desired
    //! @return the property as double
    double AtomTypeData::GetAtomTypeProperty( const AtomTypeData::Properties &PROPERTY) const
    {
      return m_Properties[ PROPERTY];
    }

    //! @brief given a partial sigma charge, determine the electronegativity
    //! @param CHARGE the partial sigma charge
    //! @return the resulting electronegativity
    double AtomTypeData::GetSigmaENFromCharge( const double &CHARGE) const
    {
      return m_SigmaChargeToEN( std::max( CHARGE, m_SigmaChargeToENInflection));
    }

    //! @brief given a partial pi charge, determine the electronegativity
    //! @param CHARGE the partial pi charge
    //! @return the resulting electronegativity
    double AtomTypeData::GetPiENFromCharge( const double &CHARGE) const
    {
      return m_PiChargeToEN( std::max( CHARGE, m_PiChargeToENInflection));
    }

    //! @brief given a partial sigma charge, determine the lone pair electronegativity
    //! @param CHARGE the partial sigma charge
    //! @return the resulting electronegativity
    double AtomTypeData::GetLonePairENFromSigmaCharge( const double &CHARGE) const
    {
      return m_SigmaChargeToLonePairEN( std::max( CHARGE, m_SigmaChargeToLonePairENInflection));
    }

    //! @brief get the function used to compute GetSigmaENFromCharge
    //! @return the function used to compute GetSigmaENFromCharge
    const math::Polynomial &AtomTypeData::GetSigmaENFromChargeFunction() const
    {
      return m_SigmaChargeToEN;
    }

    //! @brief get the function used to compute GetPiENFromCharge
    //! @return the function used to compute GetPiENFromCharge
    const math::Polynomial &AtomTypeData::GetPiENFromChargeFunction() const
    {
      return m_PiChargeToEN;
    }

    //! @brief get the function used to compute GetLonePairENFromCharge
    //! @return the function used to compute GetLonePairENFromCharge
    const math::Polynomial &AtomTypeData::GetLonePairENFromChargeFunction() const
    {
      return m_SigmaChargeToLonePairEN;
    }

    //! @return the orbital electronegativity associated with the charged state
    double AtomTypeData::GetOrbitalENegPos() const
    {
      return m_OrbitalENegPos;
    }

    //! @return valence electrons in sp orbitals
    size_t AtomTypeData::GetValenceElectronsSP() const
    {
      return m_ElementType->GetElectronConfiguration().ValenceElectronsSP() - GetFormalCharge();
    }

    //! @return Number of bonds
    size_t AtomTypeData::GetNumberBonds() const
    {
      return m_NumberBonds;
    }

    //! @brief determine if this atom type can participate in pi-bond conjugation
    //! @return true iff this atom type has any e- in p-orbitals or lone pairs
    bool AtomTypeData::IsConjugated() const
    {
      return m_Conjugated;
    }

    //! @brief is this a well characterized gasteiger atom type
    //! @return true iff this atom type is this a well characterized gasteiger atom type
    bool AtomTypeData::IsGasteigerAtomType() const
    {
      return m_IsGasteigerType;
    }

    //! @brief get the base atom type
    //! @return the atom type with only element type and charge information
    const AtomType &AtomTypeData::GetBaseAtomType() const
    {
      return *m_BaseType;
    }

    //! @brief Get the max number of electrons available for contribution to an aromatic ring
    //! @return the max electrons contributed by this atom type to a pi system
    size_t AtomTypeData::GetMaxEContributionToPiSystem() const
    {
      return m_MaxEContributionToPiSystem;
    }

    //! @brief Get the type of contribution this atom type can make to a pi system
    //! @return the type of contribution this atom type can make to a pi system
    AtomTypeData::PiContributionType AtomTypeData::GetPiElectronContributionType() const
    {
      return m_PiContribution;
    }

    //! @brief Get the stability metric.  Electronic stability is indicated by a larger number
    //! This is used to decide between atom types when no other means are possible
    double AtomTypeData::GetStabilityMetric() const
    {
      return m_StabilityMetric;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AtomTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ElementType, ISTREAM);
      io::Serialize::Read( m_Hybridization, ISTREAM);
      io::Serialize::Read( m_NumberHybridOrbitalsSigmaBinding, ISTREAM);
      io::Serialize::Read( m_NumberHybridOrbitalsNonbinding, ISTREAM);
      io::Serialize::Read( m_PiOrbitalsBinding, ISTREAM);
      io::Serialize::Read( m_AtomicOrbitalsNonbinding, ISTREAM);
      io::Serialize::Read( m_Charge, ISTREAM);
      io::Serialize::Read( m_SigmaChargeToEN, ISTREAM);
      io::Serialize::Read( m_PiChargeToEN, ISTREAM);
      io::Serialize::Read( m_SigmaChargeToLonePairEN, ISTREAM);

      // Ensure that the number of properties is the same as when the file was written
      size_t properties_in_files;
      io::Serialize::Read( properties_in_files, ISTREAM);
      BCL_Assert
      (
        properties_in_files == size_t( s_NumberOfProperties),
        "Number of properties in files was incorrect"
      );

      for( size_t a = 0; a < size_t( s_NumberOfProperties); a++)
      {
        ISTREAM >> m_Properties[ a];
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AtomTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ElementType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Hybridization, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberHybridOrbitalsSigmaBinding, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberHybridOrbitalsNonbinding, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PiOrbitalsBinding, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomicOrbitalsNonbinding, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Charge, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SigmaChargeToEN, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PiChargeToEN, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SigmaChargeToLonePairEN, OSTREAM, INDENT) << '\n';

      // Write out the number of properties, if this changes, the read function will
      // fail via BCL_Assert
      io::Serialize::Write( size_t( s_NumberOfProperties), OSTREAM, INDENT) << '\n';

      for( size_t a = 0; a < size_t( s_NumberOfProperties); a++)
      {
        io::Serialize::Write( m_Properties[ a], OSTREAM, INDENT) << '\n';
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the stability metric.  Electronic stability is indicated by a larger number
    //! This is used to decide between atom types when no other means are possible
    double AtomTypeData::CalculateStabilityMetric() const
    {
      // Consider the markov chain model with this atom type
      // Allowed transitions are
      // 1. Gaining an electron in a sigma orbital
      //    This releases SigmaValenceStateElectronAffinity eV
      // 2. Removing an electron from a sigma orbital
      //    This requires SigmaValenceStateIonizationPotential eV
      // 3. Gaining an electron in a pi orbital
      //    This releases PiValenceStateElectronAffinity eV
      // 4. Removing an electron from a pi orbital
      //    This releases PiValenceStateIonizationPotential eV
      // Given this markov chain model, the effective time constant is
      // 1 / TimeConstant
      // = 1 / ( 1/SigmaValenceStateElectronAffinity + 1/SigmaValenceStateIonizationPotential
      //         + 1/PiValenceStateElectronAffinity + 1/PiValenceStateIonizationPotential)
      // The atom type with the longest time constant should be the most stable, provided that they both have the same
      // charge

      double sigma_ea( m_Properties[ AtomTypeData::e_SigmaValenceStateElectronAffinity]);
      double sigma_ip( m_Properties[ AtomTypeData::e_SigmaValenceStateIonizationPotential]);
      double pi_ea( m_Properties[ AtomTypeData::e_PiValenceStateElectronAffinity]);
      double pi_ip( m_Properties[ AtomTypeData::e_PiValenceStateIonizationPotential]);

      double reciprocal_time_constant( 0.0);
      if( sigma_ea > 0.0 && util::IsDefined( sigma_ea)) // add data on the energy released by adding 1 sigma-electron
      {
        reciprocal_time_constant += 1.0 / sigma_ea;
      }
      if( sigma_ip > 0.0 && util::IsDefined( sigma_ip)) // add data on the energy required to remove a sigma-electron
      {
        reciprocal_time_constant += 1.0 / sigma_ip;
      }
      if( pi_ea > 0.0 && util::IsDefined( pi_ea)) // add data on the energy released by adding 1 pi-electron
      {
        reciprocal_time_constant += 1.0 / pi_ea;
      }
      if( pi_ip > 0.0 && util::IsDefined( pi_ip)) // add data on the energy required to remove a pi-electron
      {
        reciprocal_time_constant += 1.0 / pi_ip;
      }

      // if the time constant is 0.0 or less, then the type is very unstable, so we should return 0 for the stability
      return ( reciprocal_time_constant <= 0.0 ? 0.0 : 1.0 / reciprocal_time_constant);
    }

    //! @brief GetAverageNeutralSigmaIVToEARatio helper function for AtomTypes::CalculateElectronegativityValues
    //! @return reference to a double, which returns the ratio of Average(e_SigmaValenceStateIonizationPotential) for neutral atoms vs. anions
    double &AtomTypeData::GetAverageNeutralSigmaIPToAnionIPRatio()
    {
      static double s_Data( util::GetUndefined< double>());
      return s_Data;
    }

    //! @brief GetAverageNeutralPiIPToAnionIPRatio helper function for AtomTypes::CalculateElectronegativityValues
    //! @return reference to a double, which returns the ratio of Average(e_PiValenceStateIonizationPotential) for neutral atoms vs. anions
    double &AtomTypeData::GetAverageNeutralPiIPToAnionIPRatio()
    {
      static double s_Data( util::GetUndefined< double>());
      return s_Data;
    }

    //! @brief GetAverageCationSigmaIPToNeutralIPRatio helper function for AtomTypes::CalculateElectronegativityValues
    //! @return reference to a double, which returns the ratio of Average(e_SigmaValenceStateIonizationPotential) for cations vs. neutral atoms
    double &AtomTypeData::GetAverageCationSigmaIPToNeutralIPRatio()
    {
      static double s_Data( util::GetUndefined< double>());
      return s_Data;
    }

    //! @brief GetAverageCationPiIPToNeutralIPRatio helper function for AtomTypes::CalculateElectronegativityValues
    //! @return reference to a double, which returns the ratio of Average(e_PiValenceStateIonizationPotential) for cations vs. neutral atoms
    double &AtomTypeData::GetAverageCationPiIPToNeutralIPRatio()
    {
      static double s_Data( util::GetUndefined< double>());
      return s_Data;
    }

    //! @brief Estimate electronegativities of the anionic, cationic, or ground state given values from at least one other state
    //! @param ELECTRONEGATIVITIES electronegativities for the charge = -1,0, and +1 state
    //! @param TYPE_DIFFERENCE what orbital the electronegativity relates to
    void AtomTypeData::EstimateUndefinedElectronegativities
    (
      linal::Vector< double> &ELECTRONEGATIVITIES,
      const AtomTypeData::TypeDifference &TYPE_DIFFERENCE
    )
    {
      BCL_Assert( ELECTRONEGATIVITIES.GetSize() == 3, "EstimateUndefinedElectronegativities needs a vector of size 3");

      double &anion( ELECTRONEGATIVITIES( 0));   //! reference for the electronegativity of the anionic atom type
      double &neutral( ELECTRONEGATIVITIES( 1)); //! reference for the electronegativity of the neutral atom type
      double &cation( ELECTRONEGATIVITIES( 2));  //! reference for the electronegativity of the cationic atom type

      // mulliken definition of electronegativity for bonding orbitals: (ionization_potenial + electron_affinity)/2
      //                                          for lone pair orbitals: (3*ionization_potenial - electron_affinity)/2
      // DIFFERENCE could be for either binding pairs or lone pairs, so store the coefficients of ionization potential
      // and electron affinity
      static const double ionization_potential_coeff[ s_NumberTypeDifferences] = { 0.0, 0.5, 0.5, 1.5, 0.0};
      static const double electron_affinity_coeff[ s_NumberTypeDifferences] = { 0.0, 0.5, 0.5, -0.5, 0.0};

      // select the appropriate coefficient from the array
      const double ip_coeff( ionization_potential_coeff[ TYPE_DIFFERENCE]);
      const double ea_coeff( electron_affinity_coeff[ TYPE_DIFFERENCE]);

      // for brevity, in the following code, abbreviations are used as follows
      // EA -> electron affinity
      // IP -> ionization potential
      // EN -> electronegativity
      if( m_Charge == 0 && util::IsDefined( neutral))
      {
        if( !util::IsDefined( cation))
        {
          // cation.EN = ( cation.IP + cation.EA ) / 2 = ( cation.IP + neutral.IP ) / 2
          // We estimate cation.IP = neutral.IP * GetAverageCationSigmaIPToNeutralIPRatio(), then
          // cation.EN = neutral.IP * ( GetAverageCationSigmaIPToNeutralIPRatio() + 1.0) / 2
          cation = GetIonizationPotential( TYPE_DIFFERENCE)
                     * ( ip_coeff * GetAverageIPChangeCationToNeutral( TYPE_DIFFERENCE) + ea_coeff);
        }
        if( !util::IsDefined( anion))
        {
          // anion.EN = ( anion.IP + anion.EA ) / 2 = ( neutral.EA + anion.EA ) / 2
          // anion.EA = 0 (see Gasteiger Marsili 1979), then
          // anion.EN = neutral.EA / 2.0
          anion = ip_coeff * GetElectronAffinity( TYPE_DIFFERENCE);
        }
      }
      else if( m_Charge == 1 && util::IsDefined( cation))
      {
        if( !util::IsDefined( neutral))
        {
          // neutral.EN = ( neutral.IP + neutral.EA ) / 2 = ( cation.EA + neutral.EA ) / 2
          // We estimate neutral.EA = cation.EA / GetAverageNeutralSigmaIPToAnionIPRatio(), then
          // neutral.EN = cation.EA * ( 1.0 + 1.0 / GetAverageNeutralSigmaIPToAnionIPRatio()) / 2
          neutral = GetElectronAffinity( TYPE_DIFFERENCE)
                      * ( ip_coeff + ea_coeff / GetAverageIPChangeNeutralToAnion( TYPE_DIFFERENCE));
        }
        if( !util::IsDefined( anion))
        {
          // anion.EN = ( anion.IP + anion.EA ) / 2 = ( neutral.EA + anion.EA ) / 2
          // anion.EA = 0 (see Gasteiger Marsili 1979)
          // We estimate neutral.EA = cation.EA / GetAverageNeutralSigmaIPToAnionIPRatio()^2
          // anion.EN = neutral.EA / 2.0
          anion = ip_coeff * GetElectronAffinity( TYPE_DIFFERENCE)
                    / math::Sqr( GetAverageIPChangeNeutralToAnion( TYPE_DIFFERENCE));
        }
      }
      else if( m_Charge == -1 && util::IsDefined( anion))
      {
        if( !util::IsDefined( neutral))
        {
          neutral = GetIonizationPotential( TYPE_DIFFERENCE)
                      * ( ip_coeff * GetAverageIPChangeNeutralToAnion( TYPE_DIFFERENCE) + ea_coeff);
        }
        if( !util::IsDefined( cation))
        {
          cation = GetIonizationPotential( TYPE_DIFFERENCE)
                      * GetAverageNeutralSigmaIPToAnionIPRatio()
                      * ( ip_coeff * GetAverageIPChangeCationToNeutral( TYPE_DIFFERENCE) + ea_coeff);
        }
      }
      else
      {
        BCL_MessageCrt
        (
          "Could not estimate an electronegativity because necessary values were undefined"
        );
      }
    }

    //! @brief set orbital electronegativity for the charged species
    //! @param TYPE_DIFFERENCE whether the electronegativity is S/P Bonding/LonePair orbitals
    //! @param ELECTRONEGATIVITY electronegativity of the charged or neutral species, if available
    //! Used only in AtomTypes constructor
    void AtomTypeData::SetSimilarOrbitalElectronegativity
    (
      const TypeDifference &TYPE_DIFFERENCE,
      const double &ELECTRONEGATIVITY
    )
    {
      // only handle -1,0,and +1 charges
      if( m_Charge || !util::IsDefined( GetElectronegativity( TYPE_DIFFERENCE)))
      {
        return;
      }

      linal::Vector< double> charges( 3), electronegativities( 3, util::GetUndefined< double>());

      charges( 0) = -1.0;
      charges( 1) = 0.0;
      charges( 2) = 1.0;

      electronegativities( 1) = GetElectronegativity( TYPE_DIFFERENCE);
      electronegativities( 2) = ELECTRONEGATIVITY;

      EstimateUndefinedElectronegativities( electronegativities, TYPE_DIFFERENCE);

      math::Polynomial &charge_func( GetElectronegativityFromChargeFunction( TYPE_DIFFERENCE));
      charge_func = math::Polynomial::MakeFromCoordinates( charges, electronegativities, 2);

      const double highest_order_coefficient( charge_func.GetCoefficients()( 2));
      if( TYPE_DIFFERENCE == e_NumberBondingSOrbitals)
      {
        m_OrbitalENegPos = electronegativities( 2);
      }
      double &inflection_point
      (
        TYPE_DIFFERENCE == e_NumberBondingSOrbitals
        ? m_SigmaChargeToENInflection
        : TYPE_DIFFERENCE == e_NumberBondingPOrbitals
          ? m_PiChargeToENInflection
          : m_SigmaChargeToLonePairENInflection
      );
      if( highest_order_coefficient < 0.0)
      {
        // charge function must never be concave
        // so estimation must have been wrong.  Use a linear function about the two points instead
        charge_func.SetCoefficients
        (
          storage::Vector< double>::Create( electronegativities( 1), electronegativities( 2) - electronegativities( 1))
        );
        inflection_point = -std::numeric_limits< double>::infinity();
      }
      else
      {
        // inflection point of a quadratic equation is -b / ( 2 * a)
        inflection_point = -charge_func.GetCoefficients()( 1) / ( 2.0 * charge_func.GetCoefficients()( 2));
      }
    }

    //! @brief set orbital electronegativity from the neutral type
    //! @param TYPE_DIFFERENCE whether the electronegativity is S/P Bonding/LonePair orbitals
    //! @param NEUTRAL_TYPE neutral type of the atom, if available
    //! Used only in AtomTypes constructor
    void AtomTypeData::SetSimilarOrbitalElectronegativity
    (
      const TypeDifference &TYPE_DIFFERENCE,
      const AtomTypeData &NEUTRAL_TYPE
    )
    {
      switch( TYPE_DIFFERENCE)
      {
        case e_NumberBondingPOrbitals:
          m_PiChargeToEN = NEUTRAL_TYPE.m_PiChargeToEN;
          m_PiChargeToENInflection = NEUTRAL_TYPE.m_PiChargeToENInflection;
          break;
        case e_NumberBondingSOrbitals:
          m_SigmaChargeToEN = NEUTRAL_TYPE.m_SigmaChargeToEN;
          m_OrbitalENegPos = m_SigmaChargeToEN( m_Charge + 1.0);
          m_SigmaChargeToENInflection = NEUTRAL_TYPE.m_SigmaChargeToENInflection;
          break;
        case e_NumberLonePairOrbitals:
          m_SigmaChargeToLonePairEN = NEUTRAL_TYPE.m_SigmaChargeToLonePairEN;
          m_SigmaChargeToLonePairENInflection = NEUTRAL_TYPE.m_SigmaChargeToLonePairENInflection;
          break;
        default:
          break;
      }
    }

    //! @brief type difference as string
    //! @param TYPE_DIFFERENCE the type difference for which a string is desired
    //! @return the type difference as a string
    const std::string &AtomTypeData::GetTypeDifferenceName( const AtomTypeData::TypeDifference &TYPE_DIFFERENCE)
    {
      static const std::string s_Properties[] =
      {
        "None",
        "NumberBondingSOrbitals",
        "NumberBondingPOrbitals",
        "NumberLonePairOrbitals",
        "Other",
        GetStaticClassName< TypeDifference>()
      };

      return s_Properties[ TYPE_DIFFERENCE];
    }

    //! @brief get the electronegativity type corresponding to TypeDifference
    //! @param TYPE_DIFFERENCE the type difference to get the corresponding electronegativity for
    //! @return the electronegativity type corresponding to TypeDifference
    double AtomTypeData::GetElectronegativity( const AtomTypeData::TypeDifference &TYPE_DIFFERENCE) const
    {
      switch( TYPE_DIFFERENCE)
      {
        case e_NumberBondingPOrbitals: return m_Properties[ e_PiOrbitalElectronegativityMulliken];
        case e_NumberBondingSOrbitals: return m_Properties[ e_SigmaOrbitalElectronegativityMulliken];
        case e_NumberLonePairOrbitals: return m_Properties[ e_LonePairElectronegativityMulliken];
        default:
          return util::GetUndefined< double>();
      }
      return util::GetUndefined< double>();
    }

    //! @brief set a particular data
    //! @param DATA the property to set
    //! @param VALUE the value to set the property to
    void AtomTypeData::SetProperty( const Properties &DATA, const double &VALUE)
    {
      m_Properties[ DATA] = VALUE;
      math::RunningAverage< double> ave;
      if( util::IsDefined( m_Properties[ e_CovalentRadiusSingleBond]))
      {
        ave += m_Properties[ e_CovalentRadiusSingleBond];
      }
      if( util::IsDefined( m_Properties[ e_CovalentRadiusDoubleBond]))
      {
        ave += m_Properties[ e_CovalentRadiusDoubleBond];
      }
      if( util::IsDefined( m_Properties[ e_CovalentRadiusTripleBond]))
      {
        ave += m_Properties[ e_CovalentRadiusTripleBond];
      }
      if( util::IsDefined( m_Properties[ e_CovalentRadiusAromaticBond]))
      {
        ave += m_Properties[ e_CovalentRadiusAromaticBond];
      }
      m_Properties[ e_CovalentRadiusTypical] = ave.GetAverage();
    }

    //! @brief get the average ionization potential ratio between cation and neutral atom type that differ by TYPE_DIFFERENCE
    //! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
    //! @return the ratio
    double AtomTypeData::GetAverageIPChangeCationToNeutral( const TypeDifference &TYPE_DIFFERENCE) const
    {
      switch( TYPE_DIFFERENCE)
      {
        case e_NumberBondingPOrbitals: return GetAverageCationPiIPToNeutralIPRatio();
        case e_NumberBondingSOrbitals: return GetAverageCationSigmaIPToNeutralIPRatio();
        case e_NumberLonePairOrbitals:
          return 0.5 * ( GetAverageCationPiIPToNeutralIPRatio() + GetAverageCationSigmaIPToNeutralIPRatio());
        default:
          return util::GetUndefined< double>();
      }
      return util::GetUndefined< double>();
    }

    //! @brief get the average ionization potential ratio between neutral and cation atom type that differ by TYPE_DIFFERENCE
    //! @param TYPE_DIFFERENCE the type difference to get the corresponding ratio for
    //! @return the ratio
    double AtomTypeData::GetAverageIPChangeNeutralToAnion( const TypeDifference &TYPE_DIFFERENCE) const
    {
      switch( TYPE_DIFFERENCE)
      {
        case e_NumberBondingPOrbitals: return GetAverageNeutralPiIPToAnionIPRatio();
        case e_NumberBondingSOrbitals: return GetAverageNeutralSigmaIPToAnionIPRatio();
        case e_NumberLonePairOrbitals:
          return 0.5 * ( GetAverageNeutralPiIPToAnionIPRatio() + GetAverageNeutralSigmaIPToAnionIPRatio());
        default:
          return util::GetUndefined< double>();
      }
      return util::GetUndefined< double>();
    }

    //! @brief get the ionization potential type corresponding to TypeDifference
    //! @param TYPE_DIFFERENCE the type difference to get the corresponding ionization potential for
    //! @return the ionization potential type corresponding to TypeDifference
    double AtomTypeData::GetIonizationPotential( const TypeDifference &TYPE_DIFFERENCE) const
    {
      switch( TYPE_DIFFERENCE)
      {
        case e_NumberBondingPOrbitals: return m_Properties[ e_PiValenceStateIonizationPotential];
        case e_NumberBondingSOrbitals: return m_Properties[ e_SigmaValenceStateIonizationPotential];
        case e_NumberLonePairOrbitals: return m_Properties[ e_LonePairIonizationPotential];
        default:
          return util::GetUndefined< double>();
      }
      return util::GetUndefined< double>();
    }

    //! @brief get the electron affinity type corresponding to TypeDifference
    //! @param TYPE_DIFFERENCE the type difference to get the corresponding electron affinity for
    //! @return the electron affinity type corresponding to TypeDifference
    double AtomTypeData::GetElectronAffinity( const TypeDifference &TYPE_DIFFERENCE) const
    {
      switch( TYPE_DIFFERENCE)
      {
        case e_NumberBondingPOrbitals: return m_Properties[ e_PiValenceStateElectronAffinity];
        case e_NumberBondingSOrbitals: return m_Properties[ e_SigmaValenceStateElectronAffinity];
        case e_NumberLonePairOrbitals: return m_Properties[ e_LonePairElectronAffinity];
        default:
          return util::GetUndefined< double>();
      }
      return util::GetUndefined< double>();
    }

    //! @brief get the electronegativity from charge corresponding to a TypeDifference
    //! @param TYPE_DIFFERENCE the type difference to get the corresponding function for
    //! @return the electronegativity as a function of charge polynomial corresponding to TypeDifference
    math::Polynomial &AtomTypeData::GetElectronegativityFromChargeFunction( const TypeDifference &TYPE_DIFFERENCE)
    {
      static math::Polynomial s_undefined;
      switch( TYPE_DIFFERENCE)
      {
        case e_NumberBondingPOrbitals: return m_PiChargeToEN;
        case e_NumberBondingSOrbitals: return m_SigmaChargeToEN;
        case e_NumberLonePairOrbitals: return m_SigmaChargeToLonePairEN;
        default:
          return s_undefined;
      }
      return s_undefined;
    }

    //! @brief determine the difference betweent his atom type data and another
    //! @param OTHER the atom type data to compare this atom type data to
    //! @return the corresponding TypeDifference
    AtomTypeData::TypeDifference AtomTypeData::DifferenceFrom( const AtomTypeData &OTHER)
    {
      if( &OTHER == this)
      {
        return e_None;
      }

      if( m_Charge == OTHER.m_Charge || m_ElementType != OTHER.m_ElementType || m_Hybridization != OTHER.m_Hybridization)
      {
        return e_Other;
      }

      const size_t number_sigma_orbitals( GetNumberSigmaOrbitals());
      const size_t number_sigma_orbitals_b( OTHER.GetNumberSigmaOrbitals());

      if
      (
        GetNumberHybridLonePairs() + GetNumberUnhybridizedLonePairs()
        != OTHER.GetNumberHybridLonePairs() + OTHER.GetNumberUnhybridizedLonePairs()
      )
      {
        const size_t number_unhybridized_lone_pairs( GetNumberUnhybridizedLonePairs());
        const size_t number_unhybridized_lone_pairs_b( OTHER.GetNumberUnhybridizedLonePairs());
        const size_t number_hybridized_lone_pairs( GetNumberHybridLonePairs());
        const size_t number_hybridized_lone_pairs_b( OTHER.GetNumberHybridLonePairs());

        const size_t lone_pairs( number_unhybridized_lone_pairs + number_hybridized_lone_pairs);
        const size_t lone_pairs_b( number_unhybridized_lone_pairs_b + number_hybridized_lone_pairs_b);

        // is the difference between the types == the difference in the # of lone pairs?
        if
        (
          lone_pairs + m_Charge != lone_pairs_b + OTHER.m_Charge
          || GetNumberHybridOrbitals() != OTHER.GetNumberHybridOrbitals()
        )
        {
          // the type is different in multiple ways
          return e_Other;
        }

        if
        (
          GetNumberBonds() - m_Charge != OTHER.GetNumberBonds() - OTHER.m_Charge
          && ( m_NumberElectronsInBonds - GetNumberBonds() - m_Charge) != ( OTHER.m_NumberElectronsInBonds - OTHER.GetNumberBonds() - OTHER.m_Charge)
        )
        {
          return e_Other;
        }

        // yes
        return e_NumberLonePairOrbitals;
      }

      if( number_sigma_orbitals + m_Charge == number_sigma_orbitals_b + OTHER.m_Charge)
      {
        return e_NumberBondingSOrbitals;
      }

      // if the number of sigma orbitals remained the same, then the electrons must have exclusively went into the p-orbitals
      if( number_sigma_orbitals == number_sigma_orbitals_b)
      {
        return e_NumberBondingPOrbitals;
      }

      // the type is different in multiple ways
      return e_Other;
    }

    //! @brief set that this atom type forms only linear bonds
    void AtomTypeData::SetFormsOnlyLinearBonds()
    {
      m_FormsOnlyLinearBonds = true;
    }

  } // namespace chemistry
} // namespace bcl

