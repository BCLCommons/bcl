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

#ifndef BCL_CHEMISTRY_HAS_PROPERTIES_INTERFACE_H_
#define BCL_CHEMISTRY_HAS_PROPERTIES_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from the bcl
#include "bcl_chemistry_small_molecule_misc_properties.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HasPropertiesInterface
    //! @brief Mix-in interface for objects that have stored properties
    //!
    //! @remarks example unnecessary
    //! @author kothiwsk, mendenjl
    //! @date Jan 18, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_BaseClass>
    class HasPropertiesInterface :
      virtual public t_BaseClass
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief access to the misc properties
      //! @return ref to SmallMoleculeMiscProperties
      virtual const SmallMoleculeMiscProperties &GetStoredProperties() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief store a property on the molecule such that it will be present when written out
      //! @param NAME name of property
      //! @param PROPERTY vector of properties
      template< typename t_VectorType>
      void StoreProperty( const std::string &NAME, const t_VectorType &PROPERTY)
      {
        GetStoredPropertiesNonConst().SetMDLProperty( NAME, PROPERTY);
      }

      //! @brief set a property
      //! @param NAME name of property
      //! @param PROPERTY property of the molecule
      void StoreProperty( const std::string &NAME, const std::string &PROPERTY)
      {
        GetStoredPropertiesNonConst().SetMDLProperty( NAME, PROPERTY);
      }

      //! @brief set all properties
      //! @param PROPERTIES something that implements this interface
      void StoreProperties( const HasPropertiesInterface< t_BaseClass> &BASE)
      {
        GetStoredPropertiesNonConst() = BASE.GetStoredProperties();
      }

      //! @brief remove any property with NAME
      //! @param NAME the name of the property to remove from storage
      void RemoveProperty( const std::string &NAME)
      {
        GetStoredPropertiesNonConst().RemoveProperty( NAME);
      }

      //! @brief remove any property with NAME
      //! @param NAME the name of the property to remove from storage
      void ResetStoredProperties()
      {
        GetStoredPropertiesNonConst() = SmallMoleculeMiscProperties();
      }

      //! @brief determine whether the given property exists in the storage
      //! @param NAME name of property
      //! @return true iff the property exists for this molecule either in the stored properties
      bool IsPropertyStored( const std::string &NAME) const
      {
        return GetStoredProperties().GetMDLProperties().Has( NAME);
      }

      //! @brief get an mdl (stored) property
      //! @param NAME name of property
      //! @return stored value of the property
      const std::string &GetMDLProperty( const std::string &NAME) const
      {
        return GetStoredProperties().GetMDLProperty( NAME);
      }

      //! @brief get an mdl (stored) property
      //! @param NAME name of property
      //! @return stored values for that property
      linal::Vector< float> GetMDLPropertyAsVector( const std::string &NAME) const
      {
        return GetStoredProperties().GetMDLPropertyAsVector( NAME);
      }

    private:

      //! @brief access to the misc properties
      //! @return ref to SmallMoleculeMiscProperties
      virtual SmallMoleculeMiscProperties &GetStoredPropertiesNonConst() = 0;

    }; // class HasPropertiesInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_HAS_PROPERTIES_INTERFACE_H_
