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

#ifndef BCL_BIOL_ROTAMER_LIBRARY_H_
#define BCL_BIOL_ROTAMER_LIBRARY_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_rotamer.h"
#include "storage/bcl_storage_hash_map.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotamerLibrary
    //! @brief Stores a rotamer library for amino acid side chains
    //!
    //! @see @link example_biol_rotamer_library.cpp @endlink
    //! @author fischea
    //! @date Aug 08, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API RotamerLibrary :
      public util::SerializableInterface
    {

    ///////////////
    // tytpedefs //
    ///////////////

    public:

      typedef storage::Pair< double, Rotamer> RotamerProbability;

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! file name of the rotamer library
      std::string m_LibraryFileName;

      storage::HashMap< std::string, storage::Vector< RotamerProbability> > m_Rotamers;

      //! map of rotamer library and serialization keys
      static storage::HashMap< std::string, util::ShPtr< RotamerLibrary> > s_RotamerLibraries;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

      //! @brief default constructor
      RotamerLibrary();

      //! @brief construct from members
      //! @param LIBRARY_FILE_NAME path to the rotamer library
      RotamerLibrary( const std::string &LIBRARY_FILE_NAME);

    public:

      //! @brief returns a loop library from the given file name
      //! @param LIBRARY_FILE_NAME path to the rotamer library
      static util::ShPtr< RotamerLibrary> CreateRotamerLibrary( const std::string &LIBRARY_FILE_NAME);

      //! @brief clone function
      //! @return pointer to a new RotamerLibrary
      RotamerLibrary *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief read a loop template library from the given input stream
      //! @param ISTREAM input stream from which to read the loop template library
      //! @param ANGLE_BIN_WIDTH bin width for chi angles
      void ReadLibrary( std::istream &ISTREAM, double ANGLE_BIN_WIDTH = 10.0);

    }; // class RotamerLibrary

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_ROTAMER_LIBRARY_H_
