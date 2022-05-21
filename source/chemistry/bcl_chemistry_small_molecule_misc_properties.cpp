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
#include "chemistry/bcl_chemistry_small_molecule_misc_properties.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeMiscProperties::s_Instance
    (
      GetObjectInstances().AddInstance( new SmallMoleculeMiscProperties())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SmallMoleculeMiscProperties::SmallMoleculeMiscProperties() :
      m_MDLProperties()
    {
    }

    //! @brief constructor from member
    SmallMoleculeMiscProperties::SmallMoleculeMiscProperties
    (
      const storage::Map< std::string, std::string> &PROPERTIES_MAP
    ) :
      m_MDLProperties( PROPERTIES_MAP)
    {
    }

    //! @brief merge properties from two small molecule misc properties, keeping values from PROPS_A where they differ
    //! @param PROPS_A, PROPS_B two small molecule misc properties to merge, preferring values from A
    //! @return a merged small molecule misc properties
    SmallMoleculeMiscProperties::SmallMoleculeMiscProperties
    (
      const SmallMoleculeMiscProperties &PROPS_A,
      const SmallMoleculeMiscProperties &PROPS_B
    ) :
      m_MDLProperties( PROPS_A.m_MDLProperties)
    {
      Merge( PROPS_B);
    }

    //! @brief Clone function
    //! @return pointer to new SmallMoleculeMiscProperties
    SmallMoleculeMiscProperties *SmallMoleculeMiscProperties::Clone() const
    {
      return new SmallMoleculeMiscProperties( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SmallMoleculeMiscProperties::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief merge properties from another misc properties, keeping values from *this wherever
    //!        they differ
    //! @param PROPERTIES small molecule misc properties to merge, preferring values from A
    void SmallMoleculeMiscProperties::Merge( const SmallMoleculeMiscProperties &PROPERTIES)
    {
      // insert elements of b into this; this will skip any properties already in a
      m_MDLProperties.InternalData().insert
      (
        PROPERTIES.m_MDLProperties.Begin(),
        PROPERTIES.m_MDLProperties.End()
      );
    }

    //! @brief set an mdl property, overwrite if it exists, insert if it does not exist yet
    //! @param NAME name of property
    //! @param VALUE the actual textual value of the property
    void SmallMoleculeMiscProperties::SetMDLProperty( const std::string &NAME, const std::string &VALUE)
    {
      m_MDLProperties[ NAME] = VALUE;
    }

    //! @brief access to the textual mdl property with given name
    //! @param NAME name of property
    //! @return string value for the property, empty if property is no available
    const std::string &SmallMoleculeMiscProperties::GetMDLProperty( const std::string &NAME) const
    {
      static const std::string s_empty;

      // find property with given name
      storage::Map< std::string, std::string>::const_iterator itr( m_MDLProperties.Find( NAME));

      // if no such property exist, return empty value
      if( itr == m_MDLProperties.End())
      {
        return s_empty;
      }

      // return the actual value
      return itr->second;
    }

    //! @brief get a stored property
    //! @param NAME name of property
    //! @return vector of property for the molecule; empty vector if property is not known
    linal::Vector< float> SmallMoleculeMiscProperties::GetMDLPropertyAsVector( const std::string &NAME) const
    {
      const storage::Vector< std::string> result
      (
        util::SplitString( util::TrimString( GetMDLProperty( NAME)), " \t\n\r,")
      );
      storage::Vector< float> converted;
      converted.AllocateMemory( result.GetSize());
      // iterate over all strings and convert them to float
      for
      (
        storage::Vector< std::string>::const_iterator itr_str( result.Begin()), itr_str_end( result.End());
        itr_str != itr_str_end;
        ++itr_str
      )
      {
        // read the input value; if float , could be a nan so we need to use io::Serialize
        std::istringstream input( *itr_str);
        float new_value;
        io::Serialize::Read( new_value, input);
        converted.PushBack( new_value);
      }

      // end
      return converted;
    }

    //! @return a constant iterator to the beginning of the properties map
    SmallMoleculeMiscProperties::const_iterator SmallMoleculeMiscProperties::Begin() const
    {
      return m_MDLProperties.Begin();
    }

    //! @return a constant iterator to the end of the properties map
    SmallMoleculeMiscProperties::const_iterator SmallMoleculeMiscProperties::End() const
    {
      return m_MDLProperties.End();
    }

    //! @brief remove any property with NAME
    //! @param NAME the name of the property to remove
    void SmallMoleculeMiscProperties::RemoveProperty( const std::string &NAME)
    {
      storage::Map< std::string, std::string>::iterator itr( m_MDLProperties.Find( NAME));
      if( itr != m_MDLProperties.End())
      {
        m_MDLProperties.RemoveElement( itr);
      }
    }

    namespace
    {
      //! @brief helper function to determine whether a misc property name appears to be conformationally dependent
      //! @param NAME name of the misc property
      //! @return true if the descriptor matches any of the internal strings or expressions
      bool IsKnownConformationalProperty( const std::string &NAME)
      {
        // Complete descriptors names that are usually dependent on conformation
        static const storage::Set< std::string> s_conformational_names
        (
          storage::Set< std::string>::Create
          (
            "Girth", "CovalentSurfaceArea", "CovalentVolume", "VdwSurfaceArea", "VdwVolume", "Position",
            "Atom_VDWSurfaceArea", "Atom_VDWVolume", "Atom_CovalentRadius", "Atom_CovalentVolume",
            "Rgyr", "RadiusGyration"
          )
        );
        // prefixes that usually imply a conformationally-dependent descriptor.
        // This list is comprehensive only for BCL and Adriana
        static const storage::Vector< std::string> s_conformational_prefixes
        (
          storage::Vector< std::string>::Create
          (
            "CoulombicForce",
            "3da", "3DA",
            "RDF", "Rdf",
            "MolecularAsymmetry",
            "Triangulator",
            "Planarity",
            "Atom_SurfaceArea",
            "Atom_Volume",
            "Prediction",
            "SurfA"
          )
        );
        if( s_conformational_names.Contains( NAME))
        {
          return true;
        }
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( s_conformational_prefixes.Begin()), itr_end( s_conformational_prefixes.End());
          itr != itr_end;
          ++itr
        )
        {
          if( util::StartsWith( NAME, *itr))
          {
            return true;
          }
        }
        return false;
      }

    }

    //! @brief remove all descriptors that appear to be conformationally dependent
    void SmallMoleculeMiscProperties::RemoveConformationalDescriptors()
    {
      for
      (
        storage::Map< std::string, std::string>::iterator itr( m_MDLProperties.Begin()), itr_end( m_MDLProperties.End());
        itr != itr_end;
      )
      {
        if( IsKnownConformationalProperty( itr->first))
        {
          m_MDLProperties.RemoveElement( itr++);
        }
        else
        {
          ++itr;
        }
      }
    }

    //! @brief get the size of the properties
    size_t SmallMoleculeMiscProperties::GetSize() const
    {
      return m_MDLProperties.GetSize();
    }

    //! @brief test whether there are any stored properties
    //! @return true if there are no stored properties
    bool SmallMoleculeMiscProperties::IsEmpty() const
    {
      return m_MDLProperties.IsEmpty();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SmallMoleculeMiscProperties::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SmallMoleculeMiscProperties::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
