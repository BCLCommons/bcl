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

#ifndef BCL_CHEMISTRY_ELEMENT_TYPE_DATA_H_
#define BCL_CHEMISTRY_ELEMENT_TYPE_DATA_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_electron_configuration.h"
#include "bcl_chemistry_element_structure_factor.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ElementTypeData
    //! @brief stores element properties
    //! @details This is a low level helper class to store element properties
    //!
    //! @see @link example_chemistry_element_type_data.cpp @endlink
    //! @author meilerj, woetzen, mendenjl
    //! @date 08/31/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ElementTypeData :
      public util::ObjectInterface
    {

    public:

    ///////////
    // Enums //
    ///////////

      //! enum  properties for element types
      enum Properties
      {
        e_Mass,                //!< Mass
        e_GyromagneticRatio,   //!< GyromagneticRatio
        e_CovalentRadius,      //!< CovalentRadius
        e_VDWaalsRadius,       //!< VdWaalsRadius
        e_MeltingPoint,        //!< MeltingPoint
        e_BoilingPoint,        //!< BoilingPoint
        e_ElectroNegativity,   //!< ElectroNegativity
        e_IonizationPotential, //!< IonizationPotential
        e_MainGroup,           //!< Main group # (1 (alkali metals) - 8 (noble gases))
        e_HardVDWaalsRadius,   //!< Hard van der waals radius
        e_DaltonVdwRadius,     //!< A cartography of the van der Waals territories. Alvarez
        e_DaltonPvdw,          //!< Percent of molecules in the van-der waals range
        e_LJRadius,            //!< Leonard jones radius
        e_LJEpsilon,           //!< Epsilon parameter
        s_NumberOfProperties   //!< Number of properties
      };

      //! @brief element type property as string
      //! @param PROPERTY the property desired
      //! @return the property as string
      static const std::string &GetPropertyName( const Properties &PROPERTY);

      //! PropertyEnum simplifies the usage of the Properties enum of this class
      typedef util::WrapperEnum< Properties, &GetPropertyName, s_NumberOfProperties> PropertyEnum;

    private:

    //////////
    // data //
    //////////

      size_t                 m_AtomicNumber;                            //!< atomic number
      size_t                 m_Period;                                  //!< Period
      size_t                 m_MainGroup;                               //!< Group # in the main group (1-8)
      std::string            m_ChemicalSymbol;                          //!< ChemicalSymbol
      std::string            m_ChemicalName;                            //!< ChemicalName
      ElectronConfiguration  m_ElectronConfiguration;                   //!< electron configuration
      double                 m_Properties[ int( s_NumberOfProperties)]; //!< real-valued properties
      float                  m_PymolColorRGB[ 3];                       //!< rgb components for pymol color @see http://pymolwiki.org/index.php/Color_Values
      ElementStructureFactor m_StructureFactor;                         //!< Form Factor calculated from Crommer Mann coefficients

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined element
      ElementTypeData();

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
      ElementTypeData
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
      );

      //! @brief virtual copy constructor
      ElementTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief atomic number
      //! @return atomic number
      size_t GetAtomicNumber() const
      {
        return m_AtomicNumber;
      }

      //! @return Period
      size_t GetPeriod() const
      {
        return m_Period;
      }

      //! @return main Group #
      size_t GetMainGroup() const
      {
        return m_MainGroup;
      }

      //! @brief GetChemicalSymbol
      //! @return chemical symbol one or two letters as AtomName
      const std::string &GetChemicalSymbol() const
      {
        return m_ChemicalSymbol;
      }

      //! @brief GetChemicalName
      //! @return full chemical name
      const std::string &GetChemicalName() const
      {
        return m_ChemicalName;
      }

      //! @brief element type property as double
      //! @param PROPERTY the property desired
      //! @return the property as double
      double GetProperty( const ElementTypeData::Properties &PROPERTY) const
      {
        return m_Properties[ PROPERTY];
      }

      //! @brief electron configuration
      //! @return the ElectronConfiguration
      const ElectronConfiguration &GetElectronConfiguration() const
      {
        return m_ElectronConfiguration;
      }

      //! @brief rgb color for displaying the element
      //! @return a pointer to three float representing the rgb values
      const float *GetPymolColorRGB() const
      {
        return m_PymolColorRGB;
      }

      //! @brief FormFactor calculated from Crommer Mann coefficients
      //! @return FormFactor calculated from Crommer Mann coefficients
      const ElementStructureFactor &GetStructureFactor() const
      {
        return m_StructureFactor;
      }

      //! @brief tell whether this element type can participate in a conjugated system
      //! @return true if this element can participate in a common conjugated system
      //! Specifically tests if the element has 1-4 valence electrons in P orbitals
      bool IsConjugatable() const
      {
        const size_t n_valence_p( m_ElectronConfiguration.ValenceElectronsP());

        return n_valence_p > 0 && n_valence_p < 5;
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

    }; // class ElementTypeData

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_ELEMENT_TYPE_DATA_H_

