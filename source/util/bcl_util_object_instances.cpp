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
#include "util/bcl_util_object_instances.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! @brief Singleton set of known objects, used to avoid costly list searches
    storage::Map< std::string, SiPtr< const ObjectInterface> > &GetKnownObjects()
    {
      static storage::Map< std::string, SiPtr< const ObjectInterface> > s_known_objects;

      return s_known_objects;
    }

    //! @brief vector of names of all known objects
    //! @return a vector of names of all known bcl objects
    storage::Vector< std::string> &GetKnownObjectsNamesNonConst()
    {
      static storage::Vector< std::string> s_objects;
      return s_objects;
    }

    //! @brief vector of names of all known objects
    //! @return a vector of names of all known bcl objects
    const storage::Vector< std::string> &ObjectInstances::GetKnownObjectNames() const
    {
      return GetKnownObjectsNamesNonConst();
    }

    //! @brief returns a pointer to object for the DESCRIPTOR / a NULL ptr if there is no such DESCRIPTOR in the map
    //! @param DESCRIPTOR the bcl descriptor as returned by GetStaticClassname
    //! @return ptr to ObjectInterface
    const SiPtr< const ObjectInterface> &ObjectInstances::GetPtrToObjectFromIdentifier
    (
      const std::string &DESCRIPTOR
    ) const
    {
      // find the corresponding enum
      storage::Map< std::string, SiPtr< const ObjectInterface> >::const_iterator
        itr( GetKnownObjects().Find( DESCRIPTOR));

      //check that entry with DESCRIPTOR was found
      BCL_Assert
      (
        itr != GetKnownObjects().End(),
        "entry with descriptor \"" + DESCRIPTOR + "\" was not found"
      );

      //if for some reason the map contains an empty pointer
      if( !itr->second.IsDefined())
      {
        static const SiPtr< const ObjectInterface> s_empty;
        BCL_MessageCrt( "Object instance for " + DESCRIPTOR + " was NULL!");
        return s_empty;
      }

      return itr->second;
    }

    //! @brief get the address used for writing in files when a shared pointer is not intended for combining with
    //!        other shared pointers upon reading
    const size_t &ObjectInstances::GetUncombinedPointerAddress()
    {
      static const size_t s_uncombined( 12345678);
      return s_uncombined;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief returns a pointer to object for the DESCRIPTOR / an empty ptr if there is no such DESCRIPTOR in the map
    //! @param DESCRIPTOR the bcl descriptor as returned by GetStaticClassname
    //! @return ptr to ObjectInterface (must be a base of the actual object, since a dynamic cast is used) to a copy of
    //!         the class instance on every bcl object
    ObjectInterface *ObjectInstances::GetNewPtrToObjectFromIdentifier( const std::string &DESCRIPTOR) const
    {
      return GetPtrToObjectFromIdentifier( DESCRIPTOR)->Clone();
    }

    //! @brief write out an object behind a shared pointer
    //! @param PTR the shared pointer to write out
    //! @param OSTREAM the stream to write the pointer out to
    //! @param INDENT the amount of indent to use
    void ObjectInstances::WriteShPtr
    (
      const ObjectInterface *PTR,
      std::ostream &OSTREAM,
      const size_t INDENT
    ) const
    {
      // write the address; completely handles null pointers
      WritePtrAddress( size_t( PTR), OSTREAM, INDENT);
      if( !PTR)
      {
        return;
      }

      // check that the descriptor map contains the classname to guarantee a safe reading
      if( !GetKnownObjects().Has( PTR->GetClassIdentifier()))
      {
        BCL_MessageCrt
        (
          "entry with name: \"" + PTR->GetClassIdentifier() +
          "\" does not exist => reading from ShPtr will be impossible"
        );
      }

      //end
      PTR->WriteObject( OSTREAM, INDENT);
    }

    //! @brief write out a pointer address
    //! @param ADDRESS addressed memory location
    //! @param STREAM the stream to write the pointer out to
    //! @param INDENT the amount of indent to use
    void ObjectInstances::WritePtrAddress( const size_t &ADDRESS, std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write address of pointer
      io::Serialize::Write( ADDRESS ? GetUncombinedPointerAddress() : 0, OSTREAM, INDENT) << '\n';
    }

  /////////////////
  // data access //
  /////////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief add pointer of class instance to enumeration of ObjectInstances
    //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
    //! @return INTERFACE
    SiPtr< const ObjectInterface> ObjectInstances::AddInstance( const ObjectInterface *INTERFACE)
    {
      if( INTERFACE == NULL)
      {
        return SiPtr< const ObjectInterface>();
      }
      return AddInstanceWithName( INTERFACE, INTERFACE->GetClassIdentifier());
    }

    //! @brief add pointer of class instance to enumeration of ObjectInstances
    //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
    //! @return SiPtr to the interface
    SiPtr< const ObjectInterface> ObjectInstances::AddInstanceWithName
    (
      const ObjectInterface *INTERFACE,
      const std::string &NAME
    )
    {
      if( GetKnownObjects().Insert( std::make_pair( NAME, SiPtr< const ObjectInterface>( INTERFACE))).second)
      {
        GetKnownObjectsNamesNonConst().PushBack( NAME);
      }
      else
      {
        BCL_MessageCrt( NAME + " was tried to insert twice into object instances, skipped the second time!");
      }

      return SiPtr< const ObjectInterface>( INTERFACE);
    }

    //! @brief add pointer of class instance to enumeration of ObjectInstances, if it is not already there
    //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
    //! @return true if the instance was added
    bool ObjectInstances::TryAddInstance( ObjectInterface *INTERFACE)
    {
      if( INTERFACE == NULL)
      {
        return false;
      }

      if( !GetKnownObjects().Has( INTERFACE->GetClassIdentifier()))
      {
        // add the object to object instances, if it is not already there
        AddInstance( INTERFACE);
        return true;
      }
      return false;
    }

  } // namespace util

  //! @brief get the only instance of ObjectInstances
  util::ObjectInstances &GetObjectInstances()
  {
    static util::ObjectInstances s_instances;
    return s_instances;
  }

} // namespace bcl
