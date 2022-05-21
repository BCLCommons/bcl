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

#ifndef BCL_CHEMISTRY_HAS_PROPERTIES_H_
#define BCL_CHEMISTRY_HAS_PROPERTIES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from the bcl
#include "bcl_chemistry_has_properties_interface.h"
#include "bcl_chemistry_small_molecule_misc_properties.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HasProperties
    //! @brief Mix-in interface for objects that have properties
    //!
    //! @see @link example_chemistry_has_properties.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Jan 18, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_BaseClass>
    class HasProperties :
      public t_BaseClass
    {

    private:

    //////////
    // data //
    //////////

      SmallMoleculeMiscProperties m_StoredProperties; //!< Properties loaded or saved on this object

    public:

      //! @brief default constructor
      HasProperties() :
        m_StoredProperties()
      {
      }

      //! @brief constructor with stored properties
      HasProperties( const SmallMoleculeMiscProperties &STORED_PROPERTIES) :
        m_StoredProperties( STORED_PROPERTIES)
      {
      }

      //! @brief constructor with stored properties
      template< typename t_OtherBaseType>
      HasProperties( const HasPropertiesInterface< t_OtherBaseType> &PROPERTIES) :
        m_StoredProperties( PROPERTIES.GetStoredProperties())
      {
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief access to the misc properties
      //! @return ref to SmallMoleculeMiscProperties
      const SmallMoleculeMiscProperties &GetStoredProperties() const
      {
        return m_StoredProperties;
      }

      //! @brief access to the misc properties
      //! @return ref to SmallMoleculeMiscProperties
      SmallMoleculeMiscProperties &GetStoredPropertiesNonConst()
      {
        return m_StoredProperties;
      }

    }; // class HasProperties

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_HAS_PROPERTIES_H_
