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

#ifndef BCL_UTIL_OBJECT_INSTANCES_H_
#define BCL_UTIL_OBJECT_INSTANCES_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"
#include "bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectInstances
    //! @brief is an enum derived class, which enumerates one instances each of all bcl classes
    //! derived from ObjectInterface.
    //! @details each class that derives from ObjectInterface has this code:
    //!
    //! class TestClass
    //! {
    //! public:
    //!
    //! //////////
    //! // data //
    //! //////////
    //!
    //!   //! single instance of that class
    //!   static const util::SiPtr< const util::ObjectInterface> s_Instance;
    //! };
    //!
    //! which is then initialized in either in the header (for template classes) or cpp (for non template classes)
    //!
    //! //////////
    //! // data //
    //! //////////
    //!
    //!   // instantiate s_Instance
    //!   const util::SiPtr< const util::ObjectInterface> TestClass::s_Instance
    //!   (
    //!     GetObjectInstances().AddInstance( new TestClass())
    //!   );
    //!
    //! and at runtime, when those static members are initialized, they get added to the enums. For template classes
    //! it is usually necessary, if you need that instance i.e. for reading them from a stream, to make an explicit
    //! instantiation of that class:
    //!
    //! template class TestClass< double>;
    //!
    //! the name for the enum comes from the ClassIdentifier - so it is necessary to overwrite them properly.
    //!
    //! @see @link example_util_object_instances.cpp @endlink
    //! @author woetzen
    //! @date 01.05.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectInstances
    {

      template< typename t_DataType>
      friend class ShPtr;

    public:

    //////////
    // data //
    //////////

      //! @brief vector of names of all known objects
      //! @return a vector of names of all known bcl objects
      const storage::Vector< std::string> &GetKnownObjectNames() const;

      //! @brief returns a pointer to object for the DESCRIPTOR / a NULL ptr if there is no such DESCRIPTOR in the map
      //! @param DESCRIPTOR the bcl descriptor as returned by GetStaticClassname
      //! @return ptr to ObjectInterface
      const SiPtr< const ObjectInterface> &GetPtrToObjectFromIdentifier( const std::string &DESCRIPTOR) const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

    /////////////////
    // data access //
    /////////////////

      //! @brief get the address used for writing in files when a shared pointer is not intended for combining with
      //!        other shared pointers upon reading
      static const size_t &GetUncombinedPointerAddress();

      //! @brief returns a pointer to object for the DESCRIPTOR / a NULL ptr if there is no such DESCRIPTOR in the map
      //! @param DESCRIPTOR the bcl descriptor as returned by GetStaticClassname
      //! @return ptr to ObjectInterface
      ObjectInterface *GetNewPtrToObjectFromIdentifier( const std::string &DESCRIPTOR) const;

      //! @brief write out an object behind a shared pointer
      //! @param PTR the shared pointer to write out
      //! @param STREAM the stream to write the pointer out to
      //! @param INDENT the amount of indent to use
      void WriteShPtr( const ObjectInterface *PTR, std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief write out a pointer address
      //! @param ADDRESS addressed memory location
      //! @param STREAM the stream to write the pointer out to
      //! @param INDENT the amount of indent to use
      void WritePtrAddress( const size_t &ADDRESS, std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief add pointer of class instance to enumeration of ObjectInstances
      //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
      //! @return SiPtr to the interface
      SiPtr< const ObjectInterface> AddInstance( const ObjectInterface *INTERFACE);

      //! @brief add pointer of class instance to enumeration of ObjectInstances
      //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
      //! @return SiPtr to the interface
      SiPtr< const ObjectInterface> AddInstanceWithName( const ObjectInterface *INTERFACE, const std::string &NAME);

      //! @brief add pointer of class instance to enumeration of ObjectInstances, if it is not already there
      //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
      //! @return true if the instance was added
      bool TryAddInstance( ObjectInterface *INTERFACE);

    }; // class ObjectInstances

  } // namespace util

  // forward declaration of function GetObjectInstances
  BCL_API util::ObjectInstances &GetObjectInstances();

} // namespace bcl

#endif // BCL_UTIL_OBJECT_INSTANCES_H_
