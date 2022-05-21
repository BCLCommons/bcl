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
#include "chemistry/bcl_chemistry_element_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  ///////////
  // Enums //
  ///////////

    //! @brief element type property as string
    //! @param PROPERTY the property desired
    //! @return the property as string
    const std::string &ElementTypeData::GetPropertyName( const ElementTypeData::Properties &PROPERTY)
    {
      static const std::string s_properties[] =
      {
        "Mass",
        "GyromagneticRatio",
        "CovalentRadius",
        "VDWaalsRadius",
        "MeltingPoint",
        "BoilingPoint",
        "ElectroNegativity",
        "IonizationPotential",
        "MainGroup",
        "HardVDWaalsRadius",
        "DaltonVdwRadius",
        "DaltonPvdw",
        "LJRadius",
        "LJEpsilon",
        GetStaticClassName< Properties>()
      };
      return s_properties[ PROPERTY];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ElementTypeData::s_Instance( GetObjectInstances().AddInstance( new ElementTypeData()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined element
    ElementTypeData::ElementTypeData() :
      m_AtomicNumber( util::GetUndefined< size_t>()),
      m_Period( util::GetUndefined< size_t>()),
      m_ChemicalSymbol( "X"),
      m_ChemicalName( "UNDEFINED_ELEMENT"),
      m_ElectronConfiguration(),
      m_StructureFactor()
    {
      // set all properties to undefined
      std::fill( m_Properties, m_Properties + size_t( s_NumberOfProperties), util::GetUndefined< double>());

      // set color to white
      m_PymolColorRGB[ 0] = 1.0f;
      m_PymolColorRGB[ 1] = 1.0f;
      m_PymolColorRGB[ 2] = 1.0f;
    }

    //! @brief construct element from all its data
    //! @param ATOMIC_NUMBER           - number in the PSE
    //! @param PERIOD                  - in which period is the element
    //! @param MAIN_GROUP              - in which of the main groups does the element belong (0 for transition metals)
    //! @param CHEMICAL_SYMBOL         - one or two letters as in international PSE, first letter capital
    //! @param CHEMICAL_NAME           - full international name (first letter capital)
    //! @param MASS                    - the atomic mass as a weighted avergage of all isotopes
    //! @param GYROMAGNETIC_RATIO      - gyromagnetic ratio
    //! @param COVALENT_RADIUS         - radius of atom with electrons
    //! @param VDW_RADIUS              - vdw radius
    //! @param MELTING_POINT           - melting point of the plain element in its natural form
    //! @param BOILING_POINT           - boiling point of the plain element in its natural form
    //! @param ELECTRO_NEGATIVITY      - electronegativitiy
    //! @param IONIZATION_POTENTIAL    - first ionization potential
    //! @param ELECTRON_CONFIGURATION  - the electron configuration
    //! @param PYMOL_COLOR_R           - red   component for pymol color
    //! @param PYMOL_COLOR_G           - green component for pymol color
    //! @param PYMOL_COLOR_B           - blue  component for pymol color
    //! @param STRUCTURE_FACTOR        - form factor calculated from Crommer Mann coefficients
    ElementTypeData::ElementTypeData
    (
      const size_t ATOMIC_NUMBER,
      const size_t PERIOD,
      const size_t MAIN_GROUP,
      const std::string &CHEMICAL_SYMBOL,
      const std::string &CHEMICAL_NAME,
      const double MASS,
      const double GYROMAGNETIC_RATIO,
      const double COVALENT_RADIUS,
      const double VDW_RADIUS,
      const double MELTING_POINT,
      const double BOILING_POINT,
      const double VDW_DALTON,
      const double PVDW_DALTON,
      const double ELECTRO_NEGATIVITY,
      const double IONIZATION_POTENTIAL,
      const ElectronConfiguration &ELECTRON_CONFIGURATION,
      const float PYMOL_COLOR_R,
      const float PYMOL_COLOR_G,
      const float PYMOL_COLOR_B,
      const ElementStructureFactor &STRUCTURE_FACTOR
    ) :
      m_AtomicNumber( ATOMIC_NUMBER),
      m_Period( PERIOD),
      m_MainGroup( MAIN_GROUP),
      m_ChemicalSymbol( CHEMICAL_SYMBOL),
      m_ChemicalName( CHEMICAL_NAME),
      m_ElectronConfiguration( ELECTRON_CONFIGURATION),
      m_StructureFactor( STRUCTURE_FACTOR)
    {
      m_Properties[ e_Mass] = MASS;
      m_Properties[ e_GyromagneticRatio]   = GYROMAGNETIC_RATIO;
      m_Properties[ e_CovalentRadius]      = COVALENT_RADIUS;
      m_Properties[ e_VDWaalsRadius]       = VDW_RADIUS;
      m_Properties[ e_MeltingPoint]        = MELTING_POINT;
      m_Properties[ e_BoilingPoint]        = BOILING_POINT;
      m_Properties[ e_ElectroNegativity]   = ELECTRO_NEGATIVITY;
      m_Properties[ e_IonizationPotential] = IONIZATION_POTENTIAL;
      m_Properties[ e_DaltonVdwRadius]     = VDW_DALTON;
      m_Properties[ e_DaltonPvdw]          = PVDW_DALTON;
      m_Properties[ e_MainGroup]           = MAIN_GROUP;
      m_PymolColorRGB[ 0]                  = PYMOL_COLOR_R;
      m_PymolColorRGB[ 1]                  = PYMOL_COLOR_G;
      m_PymolColorRGB[ 2]                  = PYMOL_COLOR_B;
      // Hard shell radius based on crystallographic open database, atoms at least 3 bonds away. Generally there is little
      // variation in main organic group elements which is what this will primarily be used for.
      switch( m_AtomicNumber)
      {
        case 1: // H
          m_Properties[ e_HardVDWaalsRadius] = 0.55;
          break;
        case 6: // C
          m_Properties[ e_HardVDWaalsRadius] = 1.1;
        case 7: // N
        case 8: // O
        case 9: // F
          m_Properties[ e_HardVDWaalsRadius] = 1.0;
          break;
        default:
          m_Properties[ e_HardVDWaalsRadius] = 1.2;
          break;
      }
      // LJ parameters, from rosetta
      switch( m_AtomicNumber)
      {
        case 1: // H
          m_Properties[ e_LJRadius] = 1.3;
          m_Properties[ e_LJEpsilon] = 0.03;
          break;
        case 6: // C
          m_Properties[ e_LJRadius] = 1.9;
          m_Properties[ e_LJEpsilon] = 0.09;
          break;
        case 7: // N
          m_Properties[ e_LJRadius] = 1.85;
          m_Properties[ e_LJEpsilon] = 0.2;
          break;
        case 8: // O
          m_Properties[ e_LJRadius] = 1.7;
          m_Properties[ e_LJEpsilon] = 0.12;
          break;
        case 9: // F
          m_Properties[ e_LJRadius] = 1.6;
          m_Properties[ e_LJEpsilon] = 0.11;
          break;
        case 16: // S
          m_Properties[ e_LJRadius] = 2.0;
          m_Properties[ e_LJEpsilon] = 0.4;
          break;
        case 17: // Cl
          m_Properties[ e_LJRadius] = 2.03;
          m_Properties[ e_LJEpsilon] = 0.15;
          break;
        default:
          m_Properties[ e_LJRadius] = 1.05 * VDW_DALTON + 0.081;
          m_Properties[ e_LJEpsilon] = std::max( 0.0533 + 0.331 * VDW_DALTON - 0.65 * PVDW_DALTON, 0.01);
          break;
      }
    }

    //! @brief virtual copy constructor
    ElementTypeData *ElementTypeData::Clone() const
    {
      return new ElementTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ElementTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ElementTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_AtomicNumber, ISTREAM);
      io::Serialize::Read( m_Period, ISTREAM);
      io::Serialize::Read( m_ChemicalSymbol, ISTREAM);
      io::Serialize::Read( m_ChemicalName, ISTREAM);
      io::Serialize::Read( m_ElectronConfiguration, ISTREAM);

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

      for( size_t a = 0; a < 3; a++)
      {
        io::Serialize::Read( m_PymolColorRGB[ a], ISTREAM);
      }

      io::Serialize::Read( m_StructureFactor, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &ElementTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_AtomicNumber, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Period, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChemicalSymbol, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChemicalName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ElectronConfiguration, OSTREAM, INDENT) << '\n';

      // Write out the number of properties, if this changes, the read function will
      // fail via BCL_Assert
      io::Serialize::Write( size_t( s_NumberOfProperties), OSTREAM, INDENT) << '\n';

      for( size_t a = 0; a < size_t( s_NumberOfProperties); a++)
      {
        io::Serialize::Write( m_Properties[ a], OSTREAM, INDENT) << '\n';
      }

      for( size_t a = 0; a < 3; a++)
      {
        io::Serialize::Write( m_PymolColorRGB[ a], OSTREAM, INDENT) << '\n';
      }

      io::Serialize::Write( m_StructureFactor, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl

