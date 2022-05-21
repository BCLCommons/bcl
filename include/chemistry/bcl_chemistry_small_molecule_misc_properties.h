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

#ifndef BCL_CHEMISTRY_SMALL_MOLECULE_MISC_PROPERTIES_H_
#define BCL_CHEMISTRY_SMALL_MOLECULE_MISC_PROPERTIES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SmallMoleculeMiscProperties
    //! @brief TODO: add an general comment to this class
    //!
    //! @see @link example_chemistry_small_molecule_misc_properties.cpp @endlink
    //! @author loweew, woetzen
    //! @date Feb 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmallMoleculeMiscProperties :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! map of properties, where key is the property identifier/name and the value is the actual unparsed string as
      //! written in an mdl file
      storage::Map< std::string, std::string> m_MDLProperties;

    public:

      //! typedef for the iterator into the map
      typedef storage::Map< std::string, std::string>::const_iterator const_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SmallMoleculeMiscProperties();

      //! @brief constructor from member
      SmallMoleculeMiscProperties( const storage::Map< std::string, std::string> &PROPERTIES_MAP);

      //! @brief merge properties from two small molecule misc properties, keeping values from PROPERTIES_A wherever
      //!        they differ
      //! @param PROPERTIES_A, PROPERTIES_B two small molecule misc properties to merge, preferring values from A
      SmallMoleculeMiscProperties
      (
        const SmallMoleculeMiscProperties &PROPERTIES_A,
        const SmallMoleculeMiscProperties &PROPERTIES_B
      );

      //! @brief Clone function
      //! @return pointer to new SmallMoleculeMiscProperties
      SmallMoleculeMiscProperties *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get all mdl properties
      //! @return properties map from names to values
      const storage::Map< std::string, std::string> &GetMDLProperties() const
      {
        return m_MDLProperties;
      }

      //! @brief merge properties from another misc properties, keeping values from *this wherever
      //!        they differ
      //! @param PROPERTIES small molecule misc properties to merge, preferring values from A
      void Merge( const SmallMoleculeMiscProperties &PROPERTIES);

      //! @brief set an mdl property, overwrite if it exists, insert if it does not exist yet
      //! @param NAME name of property
      //! @param VALUE the actual textual value of the property
      void SetMDLProperty( const std::string &NAME, const std::string &VALUE);

      //! @brief set a property
      //! @brief NAME name of property
      //! @param PROPERTY vector of properties
      template< typename t_VectorType>
      void SetMDLProperty( const std::string &NAME, const t_VectorType &PROPERTY)
      {
        // convert vector into string; use io::Serialize so that nan's are output as "nan" independent of machine
        std::stringstream value;
        for
        (
          auto itr( PROPERTY.Begin()), itr_end( PROPERTY.End());
          itr != itr_end;
          ++itr
        )
        {
          io::Serialize::Write( *itr, value) << ' ';
        }

        // set the actual string property
        this->SetMDLProperty( NAME, value.str());
      }

      //! @brief get a stored property
      //! @param NAME name of property
      //! @return vector of property for the molecule; empty vector if property is not known
      linal::Vector< float> GetMDLPropertyAsVector( const std::string &NAME) const;

      //! @brief access to the textual mdl property with given name
      //! @param NAME name of property
      //! @return string value for the property, empty if property is no available
      const std::string &GetMDLProperty( const std::string &NAME) const;

      //! @return a constant iterator to the beginning of the properties map
      const_iterator Begin() const;

      //! @return a constant iterator to the end of the properties map
      const_iterator End() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief remove any property with NAME
      //! @param NAME the name of the property to remove
      void RemoveProperty( const std::string &NAME);

      //! @brief remove all descriptors that appear to be conformationally dependent
      void RemoveConformationalDescriptors();

      //! @brief get the size of the properties
      size_t GetSize() const;

      //! @brief test whether there are any stored properties
      //! @return true if there are no stored properties
      bool IsEmpty() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SmallMoleculeMiscProperties

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SMALL_MOLECULE_MISC_PROPERTIES_H_

