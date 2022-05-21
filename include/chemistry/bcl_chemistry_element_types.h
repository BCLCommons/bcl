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

#ifndef BCL_CHEMISTRY_ELEMENT_TYPES_H_
#define BCL_CHEMISTRY_ELEMENT_TYPES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_element_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ElementTypes
    //! @brief enumeration class for element types
    //! @details Allows access to the element type properties of a given atom
    //!
    //! Source: many sources give these values, including http://environmentalchemistry.com/yogi/periodic/
    //!         The original source for these values is unknown
    //! VDW radii source: Cambridge structural database: http://www.ccdc.cam.ac.uk/products/csd/radii/
    //!
    //! @see @link example_chemistry_element_types.cpp @endlink
    //! @author meilerj
    //! @date 08/31/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ElementTypes :
      public util::Enumerate< ElementTypeData, ElementTypes>
    {
      friend class util::Enumerate< ElementTypeData, ElementTypes>;

    public:

    //////////
    // data //
    //////////

      // declare all elements
      const ElementType e_Hydrogen;
      const ElementType e_Helium;
      // 2nd period
      const ElementType e_Lithium;
      const ElementType e_Beryllium;
      const ElementType e_Boron;
      const ElementType e_Carbon;
      const ElementType e_Nitrogen;
      const ElementType e_Oxygen;
      const ElementType e_Fluorine;
      const ElementType e_Neon;
      // 3rd period
      const ElementType e_Sodium;
      const ElementType e_Magnesium;
      const ElementType e_Aluminum;
      const ElementType e_Silicon;
      const ElementType e_Phosphorus;
      const ElementType e_Sulfur;
      const ElementType e_Chlorine;
      const ElementType e_Argon;
      // 4th period
      const ElementType e_Potassium;
      const ElementType e_Calcium;
      const ElementType e_Scandium;
      const ElementType e_Titanium;
      const ElementType e_Vanadium;
      const ElementType e_Chromium;
      const ElementType e_Manganese;
      const ElementType e_Iron;
      const ElementType e_Cobalt;
      const ElementType e_Nickel;
      const ElementType e_Copper;
      const ElementType e_Zinc;
      const ElementType e_Gallium;
      const ElementType e_Germanium;
      const ElementType e_Arsenic;
      const ElementType e_Selenium;
      const ElementType e_Bromine;
      const ElementType e_Krypton;
      // 5th period
      const ElementType e_Rubidium;
      const ElementType e_Strontium;
      const ElementType e_Yttrium;
      const ElementType e_Zirconium;
      const ElementType e_Niobium;
      const ElementType e_Molybdenum;
      const ElementType e_Technetium;
      const ElementType e_Ruthenium;
      const ElementType e_Rhodium;
      const ElementType e_Palladium;
      const ElementType e_Silver;
      const ElementType e_Cadmium;
      const ElementType e_Indium;
      const ElementType e_Tin;
      const ElementType e_Antimony;
      const ElementType e_Tellurium;
      const ElementType e_Iodine;
      const ElementType e_Xenon;
      // 6th period
      const ElementType e_Cesium;
      const ElementType e_Barium;
      const ElementType e_Lanthanum;
      const ElementType e_Cerium;
      const ElementType e_Praseodymium;
      const ElementType e_Neodymium;
      const ElementType e_Promethium;
      const ElementType e_Samarium;
      const ElementType e_Europium;
      const ElementType e_Gadolinium;
      const ElementType e_Terbium;
      const ElementType e_Dysprosium;
      const ElementType e_Holmium;
      const ElementType e_Erbium;
      const ElementType e_Thulium;
      const ElementType e_Ytterbium;
      const ElementType e_Lutetium;
      const ElementType e_Hafnium;
      const ElementType e_Tantalum;
      const ElementType e_Tungsten;
      const ElementType e_Rhenium;
      const ElementType e_Osmium;
      const ElementType e_Iridium;
      const ElementType e_Platinum;
      const ElementType e_Gold;
      const ElementType e_Mercury;
      const ElementType e_Thallium;
      const ElementType e_Lead;
      const ElementType e_Bismuth;
      const ElementType e_Polonium;
      const ElementType e_Astatine;
      const ElementType e_Radon;
      // 7th period
      const ElementType e_Francium;
      const ElementType e_Radium;
      const ElementType e_Actinium;
      const ElementType e_Thorium;
      const ElementType e_Protactinium;
      const ElementType e_Uranium;
      const ElementType e_Neptunium;
      const ElementType e_Plutonium;
      const ElementType e_Americium;
      const ElementType e_Curium;
      const ElementType e_Berkelium;
      const ElementType e_Californium;
      const ElementType e_Einsteinium;
      const ElementType e_Fermium;
      const ElementType e_Mendelevium;
      const ElementType e_Nobelium;
      const ElementType e_Lawrencium;
      const ElementType e_Rutherfordium;
      const ElementType e_Dubnium;
      const ElementType e_Seaborgium;
      const ElementType e_Bohrium;
      const ElementType e_Hassium;
      const ElementType e_Meitnerium;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct ElementTypes with all instances of the enums
      ElementTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! obtain ElementTypes from symbol
      //! @param SYMBOL the atomic symbol; may also contain isotopic information
      ElementType ElementTypeLookup( const std::string &SYMBOL) const;

    }; // class ElementTypes

    BCL_API
    const ElementTypes &GetElementTypes();

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< chemistry::ElementTypeData, chemistry::ElementTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_CHEMISTRY_ELEMENT_TYPES_H_
