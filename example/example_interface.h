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

#ifndef EXAMPLE_INTERFACE_H_
#define EXAMPLE_INTERFACE_H_

// include forward header of this class

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //! @brief Interface for all examples
  //! every example class is derived from that interface. The Run function is like a main, that performs
  //! different functionality test on the class the example is written for.
  //! It provides an AddExamplePathToFilename function, that - depending on the environment - adds a path
  //! where all example files live and new files should be written to.
  //!
  //! Each example is then added to the ExampleClass - and stored with a ShPtr to an ExampleInterface.
  //! The execution of all examples is handled by the ExampleClass
  //!
  //! @author woetzen, karakam,
  //!
  //! @date August 2007
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleInterface :
    public util::ObjectInterface
  {

  public:

  //////////
  // data //
  //////////

    //! enumerator that defines the different example input paths
    enum ExampleInputTypes
    {
      e_Biology,
      e_Chemistry,
      e_Cluster,
      e_Density,
      e_Descriptor,
      e_Fold,
      e_Graphics,
      e_Math,
      e_Mc,
      e_Model,
      e_Opti,
      e_Quality,
      e_Random,
      e_Restraint,
      e_Scorestat,
      s_NumberExampleInputTypes
    };

    //! paths associated with each enum
    static const std::string s_ExampleInputPaths[];

    //! static string to hold the input folder
    static const std::string s_ExampleInputFolderName;

    //! static string to hold the output folder
    static const std::string s_ExampleOutputFolderName;

    //! static string to hold the BCL Object folder
    static const std::string s_ExampleBclObjectFolderName;

    //! static string to hold the source code output folder name
    static const std::string s_ExampleSourceCodeFolderName;

    //! static string to hold default extension for bcl object files
    static const std::string s_ExampleBclObjectExtension;

  ////////////////
  // operations //
  ////////////////

    //! @brief virtual run routine
    //! this is performing the execution of the example and needs to be overwritten by each derived example
    virtual int Run() const = 0;

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    virtual std::istream &Read( std::istream &ISTREAM);

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the input path to given filename for the given input type
    //! @param EXAMPLE_INPUT_TYPE type of example input
    //! @param FILENAME input filename
    //! @return the input path to given filename for the given input type
    static std::string AddExampleInputPathToFilename
    (
      const ExampleInputTypes &EXAMPLE_INPUT_TYPE,
      const std::string &FILENAME
    );

    //! @brief adds the example output path to filename for given sample object
    //! @param OBJECT object for which the example this function is being called was written for
    //! @param FILENAME name of the file to be outputted
    static std::string AddExampleOutputPathToFilename
    (
      const util::ObjectInterface &OBJECT,
      const std::string &FILENAME
    );

    //! @brief adds the example output path to filename for given namespace identifier
    //! @param NAMESPACE_IDENTIFIER namespace identifier of the object to which this example file belongs to
    //! @param FILENAME name of the file to be outputted
    static std::string AddExampleOutputPathToFilename
    (
      const std::string &NAMESPACE_IDENTIFIER,
      const std::string &FILENAME
    );

    //! @brief returns the bcl object output path for the given object
    //! @param OBJECT object to be outputted
    //! @param EXTENSION extension of file ".bcl" by default
    //! @return the bcl object output path for the given namespace and the object
    static std::string GetExampleOutputPathForBCLObject
    (
      const util::ObjectInterface &OBJECT,
      const std::string &EXTENSION = s_ExampleBclObjectExtension
    );

    //! @brief reads the given bcl object from associated example path
    //! @param OBJECT object to be read
    //! @param EXTENSION extension of file ".bcl" by default
    //! @return true on success
    bool ReadBCLObject( util::ObjectInterface &OBJECT, const std::string &EXTENSION = s_ExampleBclObjectExtension) const;

    //! @brief reads the given bcl object from associated example path
    //! @param OBJECT object to be read
    //! @param PATH path where the object to be read lives
    //! @return true on success
    bool ReadBCLObjectfromPath
    (
      util::ObjectInterface &OBJECT,
      const std::string &PATH = s_ExampleBclObjectFolderName
    ) const;

    //! @brief writes the given bcl object to associated example path
    //! @param EXTENSION extension of file ".bcl" by default
    //! @param OBJECT object to be written
    //! @return true on success
    bool WriteBCLObject( const util::ObjectInterface &OBJECT, const std::string &EXTENSION = s_ExampleBclObjectExtension) const;

    //! @brief tests whether Write functions for two objects return different strings
    //! @param OBJECT An object to be tested
    //! @param DIFF_OBJECT Any other object
    //! @return false if OBJECT and DIFF_OBJECT are identical
    static bool TestBCLObjectOutputDiffers
    (
      const util::ObjectInterface &OBJECT,
      const util::ObjectInterface &DIFF_OBJECT
    );

    //! @brief Checks whether the output from one object can be read into another object
    //!        such that both objects return the same string when written out
    //! @param OBJECT Any object
    //! @param OBJECT_STORAGE an object (ideally a blank object)
    //! @return false if OBJECT was not written to or read into OBJECT_STORAGE symmetrically
    //! @note there is no way for this function to know whether the outputs are equivalent in some sense;
    //! @note this function is useful when there is a one-to-one mapping between output and object
    //! @note it cannot help you if the outputs do not have to be identical (e.g. writing addresses of objects)
    static bool TestBCLObjectIOForSymmetry
    (
      const util::ObjectInterface &OBJECT,
      const util::ObjectInterface &OBJECT_STORAGE
    );

    //! @brief function for adding path to filename
    //! @param FILENAME used filename
    //! @param EXTENSION used extension, default is empty
    //! @return filename with path and extension
    static std::string AddExamplePathToFilename( const std::string &FILENAME, const std::string &EXTENSION = "");

  }; // class ExampleInterface

} // namespace bcl

#endif //EXAMPLE_INTERFACE_H_
