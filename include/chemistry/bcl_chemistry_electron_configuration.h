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

#ifndef BCL_CHEMISTRY_ELECTRON_CONFIGURATION_H_
#define BCL_CHEMISTRY_ELECTRON_CONFIGURATION_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ElectronConfiguration
    //! @brief describes the electron configuration of atoms
    //! @details Describes the electron configuration of an atom on quantum chemistry level.
    //!
    //! @see @link example_chemistry_electron_configuration.cpp @endlink
    //! @author woetzen, karakam
    //! @date 10/28/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ElectronConfiguration :
      public util::ObjectInterface
    {
    public:

      enum PrincipalQuantumNumber
      {
        e_1,
        e_2,
        e_3,
        e_4,
        e_5,
        e_6,
        e_7,
        s_MaxPrincipleQuantumNumber
      };

      enum AngularMomentumQuantumNumber
      {
        e_S,
        e_P,
        e_D,
        e_F,
        s_MaxAngularMomentumQuantumNumber
      };

      //! @brief PrincipalQuantumNumber as string
      //! @param NUM the PrincipalQuantumNumber desired
      //! @return the PrincipalQuantumNumber as string
      static const std::string &GetDescriptor( const PrincipalQuantumNumber &NUM);

      //! @brief AngularMomentumQuantumNumber as string
      //! @param NUM the AngularMomentumQuantumNumber desired
      //! @return the AngularMomentumQuantumNumber as string
      static const std::string &GetDescriptor( const AngularMomentumQuantumNumber &NUM);

      //! PrincipalQuantumNumberEnum is used for I/O of PrincipalQuantumNumber
      typedef util::WrapperEnum< PrincipalQuantumNumber, &GetDescriptor, s_MaxPrincipleQuantumNumber>
                PrincipalQuantumNumberEnum;

      //! AngularMomentumQuantumNumberEnum is used for I/O of AngularMomentumQuantumNumbers
      typedef util::WrapperEnum< AngularMomentumQuantumNumber, &GetDescriptor, s_MaxAngularMomentumQuantumNumber>
                AngularMomentumQuantumNumberEnum;

    private:

    //////////
    // data //
    //////////

      size_t m_ValenceElectronsSP;     //!< number S and P valence electrons (main group 1,...,8)
      size_t m_ValenceElectronsSPD;    //!< number S, P, and D valence electrons (group 1,...,18)
      size_t m_Electrons[ s_MaxPrincipleQuantumNumber][ 4];      //!< ElectronConfiguration
      size_t m_ValenceQuantumNumber;   //!< last quantum number with a non-zero number of electrons

      static const size_t s_MaxElectronsInOrbital[ 7][ 4]; //!< maximum electron for each orbital

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ElectronConfiguration();

      //! construct from actual number of electrons
      ElectronConfiguration
      (
        const size_t VALENCE_ELECTRONS_SP,
        const size_t VALENCE_ELECTRONS_SPD,
        const size_t NUMBER_ELECTRONS_1S,
        const size_t NUMBER_ELECTRONS_1P,
        const size_t NUMBER_ELECTRONS_1D,
        const size_t NUMBER_ELECTRONS_1F,
        const size_t NUMBER_ELECTRONS_2S,
        const size_t NUMBER_ELECTRONS_2P,
        const size_t NUMBER_ELECTRONS_2D,
        const size_t NUMBER_ELECTRONS_2F,
        const size_t NUMBER_ELECTRONS_3S,
        const size_t NUMBER_ELECTRONS_3P,
        const size_t NUMBER_ELECTRONS_3D,
        const size_t NUMBER_ELECTRONS_3F,
        const size_t NUMBER_ELECTRONS_4S,
        const size_t NUMBER_ELECTRONS_4P,
        const size_t NUMBER_ELECTRONS_4D,
        const size_t NUMBER_ELECTRONS_4F,
        const size_t NUMBER_ELECTRONS_5S,
        const size_t NUMBER_ELECTRONS_5P,
        const size_t NUMBER_ELECTRONS_5D,
        const size_t NUMBER_ELECTRONS_5F,
        const size_t NUMBER_ELECTRONS_6S,
        const size_t NUMBER_ELECTRONS_6P,
        const size_t NUMBER_ELECTRONS_6D,
        const size_t NUMBER_ELECTRONS_6F,
        const size_t NUMBER_ELECTRONS_7S,
        const size_t NUMBER_ELECTRONS_7P,
        const size_t NUMBER_ELECTRONS_7D,
        const size_t NUMBER_ELECTRONS_7F
      );

      //! @brief virtual copy constructor
      ElectronConfiguration *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! return number ValenceElectrons in the sigma valence orbitals
      size_t ValenceElectronsS() const
      {
        return m_Electrons[ m_ValenceQuantumNumber][ e_S];
      }

      //! @brief return the number of valence electrons in SP orbitals
      size_t UnpairedValenceElectronsSP() const
      {
        return std::min( m_ValenceElectronsSP, MaxValenceElectronsSP() - m_ValenceElectronsSP);
      }

      //! @return number ValenceElectrons in the pi valence orbitals
      size_t ValenceElectronsP() const
      {
        return m_Electrons[ m_ValenceQuantumNumber][ e_P];
      }

      //! @return number ValenceElectronsSP
      size_t ValenceElectronsSP() const;

      //! @return the maximum number of electrons in SP orbitals for the noble gas in this period
      size_t MaxValenceElectronsSP() const;

      //! @return number ValenceElectronsSPD
      size_t ValenceElectronsSPD() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief number of electrons in that orbital
      //! @param PRINCIPAL_QUANTUM_NUMBER 1, 2, 3, 4, 5, 6. or 7
      //! @param ANGULAR_MOMENTUM_QUANTUM_NUMBER S, P, D, or F
      //! @return number of electrons in that particular orbital indicated by PRINCIPAL_QUANTUM_NUMBER and ANGULAR_MOMENTUM_QUANTUM_NUMBER
      size_t operator ()
      (
        const PrincipalQuantumNumber PRINCIPAL_QUANTUM_NUMBER,
        const AngularMomentumQuantumNumber ANGULAR_MOMENTUM_QUANTUM_NUMBER
      ) const
      {
        return m_Electrons[ PRINCIPAL_QUANTUM_NUMBER][ ANGULAR_MOMENTUM_QUANTUM_NUMBER];
      }

      //! @brief number of electrons in the valence orbital
      //! @param ANGULAR_MOMENTUM_QUANTUM_NUMBER S, P, D, or F
      //! @return number of electrons in the valence orbital indicated by ANGULAR_MOMENTUM_QUANTUM_NUMBER
      size_t operator ()
      (
        const AngularMomentumQuantumNumber ANGULAR_MOMENTUM_QUANTUM_NUMBER
      ) const
      {
        return m_Electrons[ m_ValenceQuantumNumber][ ANGULAR_MOMENTUM_QUANTUM_NUMBER];
      }

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

    }; // ElectronConfiguration
    
  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_ELECTRON_CONFIGURATION_H_

