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

#ifndef BCL_UTIL_ENUMS_INSTANCES_H_
#define BCL_UTIL_ENUMS_INSTANCES_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <map>

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnumsInstances
    //! @brief the class collects all instances of Enumerate derived classes.
    //! @details it keeps a map of all Enumerate derived class instances and provides reading and writing mechanism for
    //!          changing the contents of exsting enums, erasing them or extending them by additional enumerated objects
    //!
    //! @see @link example_util_enums_instances.cpp @endlink
    //! @author woetzen
    //! @date Feb 24, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EnumsInstances
    {

    /////////////
    // friends //
    /////////////

      template< typename t_DataTypeNew, typename t_DerivedNew> friend class Enumerate;

    public:

    //////////
    // data //
    //////////

      //! @brief GetFlagEnumsFiles gives commandline flag to pass custom enums definitions
      //! @return ShPtr to a FlagInterface which is used to pass filenames for custom enum definitions
      static ShPtr< command::FlagInterface> &GetFlagEnumsFiles();

    private:

    //////////
    // data //
    //////////

      std::map< std::string, ObjectInterface *> m_EnumsInstances;

      typedef std::map< std::string, ObjectInterface*>::const_iterator const_iterator;
      typedef std::map< std::string, ObjectInterface*>::iterator iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      EnumsInstances();

    public:

      //! @brief return the number of instances
      size_t GetSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add enums from files that are given in the commandline
      static void AddEnumsFromCommandline();

      //! @brief the one and only instance of this class
      static EnumsInstances &GetEnumsInstances();

    private:

      //! @brief insert an instance of an Enums derived class
      void Insert( ObjectInterface &ENUMS_OBJECT, const std::string &ENUMS_NAME);

    }; // class EnumsInstances

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_ENUMS_INSTANCES_H_
