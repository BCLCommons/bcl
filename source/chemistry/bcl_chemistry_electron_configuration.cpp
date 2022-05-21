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
#include "chemistry/bcl_chemistry_electron_configuration.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{

  namespace chemistry
  {

    //! @brief PrincipalQuantumNumber as string
    //! @param NUM the PrincipalQuantumNumber desired
    //! @return the PrincipalQuantumNumber as string
    const std::string &ElectronConfiguration::GetDescriptor( const PrincipalQuantumNumber &NUM)
    {
      static const std::string s_descriptors[] =
      {
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        GetStaticClassName< PrincipalQuantumNumber>()
      };

      return s_descriptors[ NUM];
    }

    //! @brief AngularMomentumQuantumNumber as string
    //! @param NUM the AngularMomentumQuantumNumber desired
    //! @return the AngularMomentumQuantumNumber as string
    const std::string &ElectronConfiguration::GetDescriptor( const AngularMomentumQuantumNumber &NUM)
    {
      static const std::string s_descriptors[] =
      {
        "s",
        "p",
        "d",
        "f",
        GetStaticClassName< AngularMomentumQuantumNumber>()
      };

      return s_descriptors[ NUM];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ElectronConfiguration::s_Instance
    (
      GetObjectInstances().AddInstance( new ElectronConfiguration())
    );

    const size_t ElectronConfiguration::s_MaxElectronsInOrbital[ 7][ 4] =
    {
      { 2, 0, 0, 0},
      { 2, 6, 0, 0},
      { 2, 6, 0, 0},
      { 2, 6, 10, 0},
      { 2, 6, 10, 0},
      { 2, 6, 10, 14},
      { 2, 6, 10, 14}
    };

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    ElectronConfiguration::ElectronConfiguration() :
      m_ValenceElectronsSP( util::GetUndefined< size_t>()),
      m_ValenceElectronsSPD( util::GetUndefined< size_t>()),
      m_ValenceQuantumNumber( util::GetUndefined< size_t>())
    {
      for
      (
        size_t shell_type( 0), num_shells( s_MaxPrincipleQuantumNumber);
        shell_type < num_shells;
        ++shell_type
      )
      {
        for
        (
          size_t orbital_type( 0), num_orbitals( s_MaxAngularMomentumQuantumNumber);
          orbital_type < num_orbitals;
          ++orbital_type
        )
        {
          m_Electrons[ shell_type][ orbital_type] = util::GetUndefined< size_t>();
        }
      }
    }

    //! construct from actual number of electrons
    ElectronConfiguration::ElectronConfiguration
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
    ) :
      m_ValenceElectronsSP( VALENCE_ELECTRONS_SP),
      m_ValenceElectronsSPD( VALENCE_ELECTRONS_SPD),
      m_ValenceQuantumNumber( 0)
    {
      m_Electrons[ e_1][ e_S] = NUMBER_ELECTRONS_1S;
      m_Electrons[ e_1][ e_P] = NUMBER_ELECTRONS_1P;
      m_Electrons[ e_1][ e_D] = NUMBER_ELECTRONS_1D;
      m_Electrons[ e_1][ e_F] = NUMBER_ELECTRONS_1F;
      m_Electrons[ e_2][ e_S] = NUMBER_ELECTRONS_2S;
      m_Electrons[ e_2][ e_P] = NUMBER_ELECTRONS_2P;
      m_Electrons[ e_2][ e_D] = NUMBER_ELECTRONS_2D;
      m_Electrons[ e_2][ e_F] = NUMBER_ELECTRONS_2F;
      m_Electrons[ e_3][ e_S] = NUMBER_ELECTRONS_3S;
      m_Electrons[ e_3][ e_P] = NUMBER_ELECTRONS_3P;
      m_Electrons[ e_3][ e_D] = NUMBER_ELECTRONS_3D;
      m_Electrons[ e_3][ e_F] = NUMBER_ELECTRONS_3F;
      m_Electrons[ e_4][ e_S] = NUMBER_ELECTRONS_4S;
      m_Electrons[ e_4][ e_P] = NUMBER_ELECTRONS_4P;
      m_Electrons[ e_4][ e_D] = NUMBER_ELECTRONS_4D;
      m_Electrons[ e_4][ e_F] = NUMBER_ELECTRONS_4F;
      m_Electrons[ e_5][ e_S] = NUMBER_ELECTRONS_5S;
      m_Electrons[ e_5][ e_P] = NUMBER_ELECTRONS_5P;
      m_Electrons[ e_5][ e_D] = NUMBER_ELECTRONS_5D;
      m_Electrons[ e_5][ e_F] = NUMBER_ELECTRONS_5F;
      m_Electrons[ e_6][ e_S] = NUMBER_ELECTRONS_6S;
      m_Electrons[ e_6][ e_P] = NUMBER_ELECTRONS_6P;
      m_Electrons[ e_6][ e_D] = NUMBER_ELECTRONS_6D;
      m_Electrons[ e_6][ e_F] = NUMBER_ELECTRONS_6F;
      m_Electrons[ e_7][ e_S] = NUMBER_ELECTRONS_7S;
      m_Electrons[ e_7][ e_P] = NUMBER_ELECTRONS_7P;
      m_Electrons[ e_7][ e_D] = NUMBER_ELECTRONS_7D;
      m_Electrons[ e_7][ e_F] = NUMBER_ELECTRONS_7F;

      for( size_t valence_number = 1, max_shells = size_t( e_7) + 1; valence_number < max_shells; valence_number++)
      {
        if( m_Electrons[ valence_number][ e_S] > 0)
        {
          m_ValenceQuantumNumber = ElectronConfiguration::PrincipalQuantumNumber( valence_number);
        }
      }
    }

    //! @brief virtual copy constructor
    ElectronConfiguration *ElectronConfiguration::Clone() const
    {
      return new ElectronConfiguration( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ElectronConfiguration::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @return number ValenceElectronsSP
    size_t ElectronConfiguration::ValenceElectronsSP() const
    {
      return m_ValenceElectronsSP;
    }

    //! @return number ValenceElectronsSP
    size_t ElectronConfiguration::ValenceElectronsSPD() const
    {
      return m_ValenceElectronsSPD;
    }

    //! @return the maximum number of electrons in SP orbitals for the noble gas in this period
    size_t ElectronConfiguration::MaxValenceElectronsSP() const
    {
      return
        s_MaxElectronsInOrbital[ m_ValenceQuantumNumber][ e_S] + s_MaxElectronsInOrbital[ m_ValenceQuantumNumber][ e_P];
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ElectronConfiguration::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ValenceElectronsSP, ISTREAM);
      io::Serialize::Read( m_ValenceElectronsSPD, ISTREAM);
      io::Serialize::Read( m_ValenceQuantumNumber, ISTREAM);

      // reset each state
      for
      (
        size_t shell_type( 0), num_shells( s_MaxPrincipleQuantumNumber);
        shell_type < num_shells;
        ++shell_type
      )
      {
        for
        (
          size_t orbital_type( 0), num_orbitals( s_MaxAngularMomentumQuantumNumber);
          orbital_type < num_orbitals;
          ++orbital_type
        )
        {
          m_Electrons[ shell_type][ orbital_type] = 0;
        }
      }

      // get how many occupied states there were
      size_t num_occupied_states( 0);
      io::Serialize::Read( num_occupied_states, ISTREAM);
      for( size_t curr_state( 0); curr_state < num_occupied_states; ++curr_state)
      {
        std::string state;
        ISTREAM >> state;

        PrincipalQuantumNumber principle_n( PrincipalQuantumNumberEnum( state.substr( 0, 1)));
        AngularMomentumQuantumNumber angular_n( AngularMomentumQuantumNumberEnum( state.substr( 1, 1)));
        size_t read_electrons( util::ConvertStringToNumericalValue< size_t>( state.substr( 2)));

        m_Electrons[ principle_n][ angular_n] = read_electrons;
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &ElectronConfiguration::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ValenceElectronsSP, OSTREAM, INDENT);
      io::Serialize::Write( m_ValenceElectronsSPD, OSTREAM, INDENT);
      io::Serialize::Write( m_ValenceQuantumNumber, OSTREAM, INDENT);

      // count the number of occupied states
      size_t num_occupied_states( 0);
      for
      (
        size_t shell_type( 0), num_shells( s_MaxPrincipleQuantumNumber);
        shell_type < num_shells;
        ++shell_type
      )
      {
        for
        (
          size_t orbital_type( 0), num_orbitals( s_MaxAngularMomentumQuantumNumber);
          orbital_type < num_orbitals;
          ++orbital_type
        )
        {
          if( m_Electrons[ shell_type][ orbital_type] > 0)
          {
            ++num_occupied_states;
          }
        }
      }

      // write how many occupied states there are, followed by each occupied state
      io::Serialize::Write( num_occupied_states, OSTREAM, INDENT) << ' ';
      for
      (
        size_t shell_type( 0), num_shells( s_MaxPrincipleQuantumNumber);
        shell_type < num_shells;
        ++shell_type
      )
      {
        for
        (
          size_t orbital_type( 0), num_orbitals( s_MaxAngularMomentumQuantumNumber);
          orbital_type < num_orbitals;
          ++orbital_type
        )
        {
          if( m_Electrons[ shell_type][ orbital_type] > 0)
          {
            OSTREAM << PrincipalQuantumNumberEnum( PrincipalQuantumNumber( shell_type))
                    << AngularMomentumQuantumNumberEnum( AngularMomentumQuantumNumber( orbital_type))
                    << m_Electrons[ shell_type][ orbital_type] << ' ';
          }
        }
      }

      OSTREAM << '\n';

      // return
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
