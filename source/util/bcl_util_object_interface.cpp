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
#include "util/bcl_util_object_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_assert.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    //! @brief virtual destructor
    ObjectInterface::~ObjectInterface()
    {
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read the identifier
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectInterface::ReadIdentifier( std::istream &ISTREAM) const
    {
      //extract the class name form the stream
      std::string identifier( ExtractIdentifier( ISTREAM));

       // assert that this is a stream for this object
      BCL_Assert
      (
        GetClassIdentifier() == identifier,
        "wrong stream; expected \"" + GetClassIdentifier() + "\" but encountered \"" + identifier + "\""
      );

      //end
      return ISTREAM;
    }

    //! @brief write the identifier
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &ObjectInterface::WriteIdentifier( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write class name as identifier of the string
      OSTREAM << std::string( INDENT * 2, ' ') << GetClassIdentifier() << '\n';

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectInterface::ReadObject( std::istream &ISTREAM)
    {
      //read and check stream identifier
      ReadIdentifier( ISTREAM);

      //call the Read function of the class
      Read( ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &ObjectInterface::WriteObject( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write identifier
      WriteIdentifier( OSTREAM, INDENT);

      // call the write function of the class
      Write( OSTREAM, INDENT + 1);

      // end
      return OSTREAM;
    }

    // extracts the identifier string for the given ISTREAM
    std::string ObjectInterface::ExtractIdentifier( std::istream &ISTREAM)
    {
      //instantiate identifier string that will be extracted
      std::string identifier;

      // repeat getline until eof is reached or the extracted stream is not empty
      while( !ISTREAM.eof() && ISTREAM.good())
      {
        ISTREAM >> std::ws;
        //get identifier string
        ISTREAM >> identifier;

        if( !identifier.empty())
        {
          ISTREAM >> std::ws;
          break;
        }
      }

      // return identifier
      return identifier;
    }

    //! @brief extract namespace name without the "bcl::" part from a given object identifier
    //! @param IDENTIFIER object identifier such as "bcl::util::ClassName"
    //! @return namespace name without the "bcl::" part such as "util" or "assemble"
    const std::string ObjectInterface::ExtractNamespaceName( const std::string &IDENTIFIER)
    {
      static const std::string s_name( bcl::GetNamespaceIdentifier() + "::");
      BCL_Assert
      (
        IDENTIFIER.size() >= s_name.size() + 1
        && std::equal( IDENTIFIER.begin(), IDENTIFIER.begin() + s_name.size(), s_name.begin()),
        "The provided identifier does not have at least one namespace name after BCL : " + IDENTIFIER
      );

      const size_t namespace_end_pos( IDENTIFIER.find( "::", s_name.size() + 1));

      // return the second one after "bcl"
      return IDENTIFIER.substr( s_name.size(), namespace_end_pos - s_name.size());
    }

    //! @brief operator >> to read ReadWrite derived Object from std::istream
    //! @param ISTREAM inpustream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the inputstream where the object was read from
    std::istream &operator >>
    (
      std::istream &ISTREAM,
      ObjectInterface &OBJECT
    )
    {
      return OBJECT.ReadObject( ISTREAM);
    }

    //! @brief operator >> to read const ReadWrite derived Object from std::istream which will assert, since this is not possible
    //! @param ISTREAM inpustream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the inputstream where the object was read from
    std::istream &operator >>
    (
      std::istream &ISTREAM,
      const ObjectInterface &OBJECT
    )
    {
      BCL_Exit( "do not call this function with a const OBJECT", -1);
      return ISTREAM;
    }

    //! @brief write ReadWrite derived Object to std::ostream
    //! @param OSTREAM outputstream the object is written to
    //! @param OBJECT object derived from ReadWrite
    //! @return outputstream the object was written to
    std::ostream &operator <<
    (
      std::ostream &OSTREAM,
      const ObjectInterface &OBJECT
    )
    {
      return OBJECT.WriteObject( OSTREAM, 0);
    }

  } // namespace util
} // namespace bcl

