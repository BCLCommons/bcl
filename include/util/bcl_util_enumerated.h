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

#ifndef BCL_UTIL_ENUMERATED_H_
#define BCL_UTIL_ENUMERATED_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_data_label.h"
#include "bcl_util_own_ptr.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Enumerated
    //! @brief contains functions used for selection and initialization of classes derived from a common base (via util::Implementation)
    //!
    //! @tparam t_Interface the interface that should be derived from
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Jan 5, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_Interface>
    class Enumerated
    {
    public:

    //////////////
    // typedefs //
    //////////////

      typedef typename storage::Map< std::string, OwnPtr< t_Interface> >::const_iterator const_iterator;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief add an implementation of an interface to the known instances
      //! @param INSTANCE a new implementation
      //! @param ALIAS how class will be called; if left empty, the value of the classes' GetAlias function is used instead
      //! @return true
      //! Derived classes should declare a static variable s_Instance
      //! and then add it to the corresponding Enumerated interface map, e.g.
      //! class MyClass : public Enumerated< MyInterface>
      //! {
      //!   static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< true once this class has been registered with the enumerated class
      //! };
      //!
      //! const util::SiPtr< const util::ObjectInterface> MyClass::s_Instance( Enumerated< MyInterface>::AddInstance( MyClass()));
      //! note that whatever instance is passed to the enumerated class becomes the default instance of that class
      static t_Interface *AddInstance( t_Interface *INSTANCE, const std::string &ALIAS = std::string())
      {
        const std::string &alias( ALIAS.empty() ? INSTANCE->GetAlias() : ALIAS);
        // try to insert the instance
        BCL_Assert
        (
          GetInstanceMap().Insert
          (
            std::pair< std::string, OwnPtr< t_Interface> >
            (
              alias,
              OwnPtr< t_Interface>( INSTANCE, false)
            )
          ).second,
          "Already had an object with alias " + alias + " in the instance map for "
          + GetStaticClassName< t_Interface>()
        );

        // validate GetSerializer
        BCL_Assert( INSTANCE->GetCompleteSerializer().IsFinalized(), "Instance named " + alias + " has invalid serializer");

        // Try to add the object to object interface
        GetObjectInstances().TryAddInstance( INSTANCE);

        return INSTANCE;
      }

      //! @brief get the default implementation, if one is known
      //! @return the default implementation, or undefined ptr otherwise
      static const_iterator GetDefaultImplementation()
      {
        return GetInstanceMap().Find( "");
      }

      //! @brief test whether a particular implementation is already in the map
      //! @param ALIAS alias for the implementation
      static bool HaveImplementationWithAlias( const std::string &ALIAS)
      {
        return GetInstanceMap().Find( ALIAS) != GetInstanceMap().End();
      }

      //! @brief get an iterator to the beginning of the instance map
      static const_iterator Find( const std::string &ALIAS)
      {
        return GetInstanceMap().Find( ALIAS);
      }

      //! @brief get an iterator to the beginning of the instance map
      static const_iterator Begin()
      {
        return GetInstanceMap().Begin();
      }

      //! @brief get an iterator to the end of the instance map
      static const_iterator End()
      {
        return GetInstanceMap().End();
      }

      //! @brief get the number of known instances
      static size_t GetSize()
      {
        return GetInstanceMap().GetSize();
      }

    protected:

    /////////////
    // friends //
    /////////////

      friend class Implementation< t_Interface>; //!< necessary so that this class can call GetInstanceMap

      //! @brief get the map from alias names to implementations
      //! @return the map from alias names to implementations
      static storage::Map< std::string, OwnPtr< t_Interface> > &GetInstanceMap()
      {
        static storage::Map< std::string, OwnPtr< t_Interface> > s_instance_map;
        return s_instance_map;
      }

    }; // class Enumerated

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_ENUMERATED_H_
