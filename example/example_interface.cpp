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
#include "example_interface.h"

// includes from bcl - sorted alphabetically
#include "example.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{

//////////
// data //
//////////

  //! paths associated with each enum
  const std::string ExampleInterface::s_ExampleInputPaths[ s_NumberExampleInputTypes] =
  {
    "biology",
    "chemistry",
    "cluster",
    "density",
    "descriptor",
    "fold",
    "graphics",
    "math",
    "mc",
    "model",
    "opti",
    "quality",
    "random",
    "restraint",
    "scorestat"
  };

  //! static string to hold the input folder
  const std::string ExampleInterface::s_ExampleInputFolderName = "input";

  //! static string to hold the output folder
  const std::string ExampleInterface::s_ExampleOutputFolderName = "output";

  //! static string to hold the BCL Object folder
  const std::string ExampleInterface::s_ExampleBclObjectFolderName = "bcl_objects";

  //! static string to hold the source code output folder name
  const std::string ExampleInterface::s_ExampleSourceCodeFolderName = "source_code";

  //! static string to hold default extension for bcl object files
  const std::string ExampleInterface::s_ExampleBclObjectExtension = ".bcl";

//////////////////////
// input and output //
//////////////////////

  //! @brief read from std::istream
  //! @param ISTREAM input stream
  //! @return istream which was read from
  std::istream &ExampleInterface::Read( std::istream &ISTREAM)
  {
    return ISTREAM;
  }

  //! @brief write to std::ostream
  //! @param OSTREAM output stream to write to
  //! @return output stream which was written to
  std::ostream &ExampleInterface::Write( std::ostream &OSTREAM, const size_t INDENT) const
  {
    return OSTREAM;
  }

//////////////////////
// helper functions //
//////////////////////

  //! @brief function for adding path to filename
  std::string ExampleInterface::AddExamplePathToFilename( const std::string &FILENAME, const std::string &EXTENSION)
  {
    if( EXTENSION.empty())
    {
      return GetExamples().GetExamplePath() + FILENAME;
    }
    return GetExamples().GetExamplePath() + FILENAME + "." + EXTENSION;
  }

  //! @brief returns the input path to given filename for the given input type
  //! @param EXAMPLE_INPUT_TYPE type of example input
  //! @param FILENAME input filename
  //! @return the input path to given filename for the given input type
  std::string ExampleInterface::AddExampleInputPathToFilename
  (
    const ExampleInterface::ExampleInputTypes &EXAMPLE_INPUT_TYPE,
    const std::string &FILENAME
  )
  {
    return
      GetExamples().GetExamplePath() +
      s_ExampleInputFolderName + PATH_SEPARATOR +
      s_ExampleInputPaths[ EXAMPLE_INPUT_TYPE] + PATH_SEPARATOR +
      FILENAME;
  }

  //! @brief adds the example output path to filename for given sample object
  //! @param OBJECT object for which the example this function is being called was written for
  //! @param FILENAME name of the file to output
  std::string ExampleInterface::AddExampleOutputPathToFilename
  (
    const util::ObjectInterface &OBJECT,
    const std::string &FILENAME
  )
  {
    return
      GetExamples().GetExamplePath() +
      s_ExampleOutputFolderName + PATH_SEPARATOR +
      util::ObjectInterface::ExtractNamespaceName( OBJECT.GetClassIdentifier()) + PATH_SEPARATOR +
      FILENAME;
  }

  //! @brief adds the example output path to filename for given namespace
  //! @param NAMESPACE_IDENTIFIER namespace identifier of the object to which this example file belongs to
  //! @param FILENAME name of the file to output
  std::string ExampleInterface::AddExampleOutputPathToFilename
  (
    const std::string &NAMESPACE_IDENTIFIER,
    const std::string &FILENAME
  )
  {
    return
      GetExamples().GetExamplePath() +
      s_ExampleOutputFolderName + PATH_SEPARATOR +
      util::ObjectInterface::ExtractNamespaceName( NAMESPACE_IDENTIFIER) + PATH_SEPARATOR + FILENAME;
  }

  //! @brief returns the bcl object output path for the given object
  //! @param OBJECT object to be outputted
  //! @param EXTENSION extension of file ".bcl" by default
  //! @return the bcl object output path for the given namespace and the object
  std::string ExampleInterface::GetExampleOutputPathForBCLObject
  (
    const util::ObjectInterface &OBJECT,
    const std::string &EXTENSION
  )
  {
    return
      GetExamples().GetExamplePath() +
      s_ExampleBclObjectFolderName + PATH_SEPARATOR +
      util::ObjectInterface::ExtractNamespaceName( OBJECT.GetClassIdentifier()) + PATH_SEPARATOR +
      io::File::ConvertClassIdentifierToFilename( OBJECT.GetClassIdentifier()) + EXTENSION;
  }

  //! @brief reads the given bcl object from associated example path
  //! @param OBJECT object to be read
  //! @param EXTENSION extension of file ".bcl" by default
  //! @return true on success
  bool ExampleInterface::ReadBCLObject( util::ObjectInterface &OBJECT, const std::string &EXTENSION) const
  {
    // initialize stream and read the object file
    io::IFStream read;

    // first, try to open the file
    if
    (
      !ExampleClass::ExampleMustOpenInputFile
      (
        GetClassIdentifier(),
        __LINE__,
        GetExampleOutputPathForBCLObject( OBJECT, EXTENSION),
        read,
        std::ios::in
      )
    )
    {
      return false;
    }

    // read in the object
    read >> OBJECT;
    io::File::CloseClearFStream( read);
    return true;
  }

  //! @brief reads the given bcl object from associated example path
  //! @param OBJECT object to be read
  //! @param PATH path where the object to be read lives
  //! @return true on success
  bool ExampleInterface::ReadBCLObjectfromPath
  (
    util::ObjectInterface &OBJECT,
    const std::string &PATH
  ) const
  {
    // initialize stream and read the object file
    io::IFStream read;

    // first, try to open the file
    if
    (
      !ExampleClass::ExampleMustOpenInputFile
      (
        GetClassIdentifier(),
        __LINE__,
        PATH,
        read,
        std::ios::in
      )
    )
    {
      return false;
    }

    // read in the object
    read >> OBJECT;
    io::File::CloseClearFStream( read);
    return true;
  }

  //! @brief writes the given bcl object to associated example path
  //! @param OBJECT object to be written
  //! @param EXTENSION extension of file ".bcl" by default
  //! @return true if the object could be written
  bool ExampleInterface::WriteBCLObject( const util::ObjectInterface &OBJECT, const std::string &EXTENSION) const
  {
    // initialize ofstream
    io::OFStream write;

    // try to open the file first
    if
    (
      !ExampleClass::ExampleMustOpenOutputFile
      (
        GetClassIdentifier(),
        __LINE__,
        GetExampleOutputPathForBCLObject( OBJECT, EXTENSION),
        write,
        std::ios::out
      )
    )
    {
      return false;
    }
    write << OBJECT;
    io::File::CloseClearFStream( write);
    return true;
  }

  //! @brief tests whether Write functions for two objects return different strings
  //! @param OBJECT An object to be tested
  //! @param DIFF_OBJECT Any other object
  //! @return false if OBJECT and DIFF_OBJECT are identical
  bool ExampleInterface::TestBCLObjectOutputDiffers
  (
    const util::ObjectInterface &OBJECT,
    const util::ObjectInterface &DIFF_OBJECT
  )
  {
    std::ostringstream output_object;
    std::ostringstream output_diff_object;
    output_object << OBJECT;
    output_diff_object << DIFF_OBJECT;
    return ( output_object.str() != output_diff_object.str());
  }

  //! @brief Checks whether the output from one object can be read into another object
  //!        such that both objects return the same string when written out
  //! @param OBJECT Any object
  //! @param OBJECT_STORAGE an object (ideally a bna a blank object
  //! @return false if OBJECT was not written to or read into OBJECT_STORAGE symmetrically
  //! @note there is no way for this function to know whether the outputs are equivalent in some sense;
  //! @note this function is useful when there is a one-to-one mapping between output and object
  //! @note it cannot help you if the outputs do not have to be identical (e.g. writing addresses of objects)
  bool ExampleInterface::TestBCLObjectIOForSymmetry
  (
    const util::ObjectInterface &OBJECT,
    const util::ObjectInterface &OBJECT_STORAGE
  )
  {
    // write OBJECT to a stream
    std::stringstream io_object;
    io_object << OBJECT;

    // clone OBJECT_STORAGE to keep it const
    util::ShPtr< util::ObjectInterface> object_storage_clone( OBJECT_STORAGE.Clone());

    // read whatever OBJECT wrote back into the cloned object
    io_object >> *object_storage_clone;

    // write OBJECT_STORAGE to a different stream
    std::ostringstream object_storage_output;
    object_storage_output << *object_storage_clone;

    const bool result( io_object.str() == object_storage_output.str());

    if( result == false)
    {
      BCL_MessageCrt
      (
        "asymmetric output; object wrote out: " + io_object.str()
        + " but clone wrote out " + object_storage_output.str()
      );
    }
    return result;
  }

} // namespace bcl
