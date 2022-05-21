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

#ifndef BCL_DESCRIPTOR_ATOM_PI_CHARGE_H_
#define BCL_DESCRIPTOR_ATOM_PI_CHARGE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomPiCharge
    //! @brief calculates the pi charge for every atom in a molecule
    //!
    //! @see @link example_descriptor_atom_pi_charge.cpp @endlink
    //! @author mendenjl, kothiwsk
    //! @date Feb 10, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomPiCharge :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      storage::Vector< float> m_Charges;    //!< Calculated charge values for the current molecule

      bool m_GetChargeOrElectronegativity;  //!< True if sigma charge is desired and false if electronegativity is desired

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param GET_CHARGE true if sigma charge is desired, false if sigma electronegativity is desired
      AtomPiCharge( const bool &GET_CHARGE = true);

      //! @brief virtual copy constructor
      AtomPiCharge *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

      //! @brief Recalculates the charges for the current molecule
      void RecalculateCharges();

      //! @return a unique string for the given connectivity information
      static std::string GetPiChargeTypeString
      (
        const chemistry::ElementType &ELEMENT_TYPE_1,
        const size_t      &BOND_ORDER_AND_AROMATIC,
        const chemistry::ElementType &ELEMENT_TYPE_2
      );

      //! @brief create the connectivity string to pi-charge parameter map
      static storage::Map< std::string, storage::VectorND< 2, float> > MakePiChargeTypeStringToPiChargeParameterMap();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class AtomPiCharge

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ATOM_PI_CHARGE_H_
