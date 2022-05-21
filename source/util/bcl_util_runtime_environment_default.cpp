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
#include "util/bcl_util_runtime_environment_default.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RuntimeEnvironmentDefault::RuntimeEnvironmentDefault()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @return the class name if this function is overwritten
    const std::string &RuntimeEnvironmentDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the process number or slot number
    //! @return the process or slot, this app is running in
    int RuntimeEnvironmentDefault::GetProcessNumber() const
    {
      return GetUndefined< int>();
    }

  ////////////////
  // operations //
  ////////////////

    //! @copydoc RuntimeEnvironmentInterface::ResolveFileName()
    std::string RuntimeEnvironmentDefault::ResolveFileName( const std::string &FILE_NAME) const
    {
      return FILE_NAME;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RuntimeEnvironmentDefault::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &RuntimeEnvironmentDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Initializes the environment
    //! @return A boolean representing whether the initialization was a success or failure.
    bool RuntimeEnvironmentDefault::Initialize() const
    {
      return ( true);
    }

    //! @brief Finalizes the environment
    //! @return A boolean representing whether finalize was a success or failure.
    //! @param STATUS indicates if an error occured
    bool RuntimeEnvironmentDefault::Finalize( const int STATUS) const
    {
      return STATUS == 0;
    }

  } // namespace util
} // namespace bcl
