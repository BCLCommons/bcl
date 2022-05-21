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

#ifndef BCL_MODEL_H_
#define BCL_MODEL_H_

// include the namespace forward header
#include "bcl_model.fwd.hh"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_model.h
  //! @brief namespace for machine learning classes and data structures
  //!
  //! @see @link example_model.cpp @endlink
  //! @author loweew, mendenjl
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace model
  {
    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Model
    //! @brief keeps a AddModelPath function out of namespace scope
    //!
    //! @see @link example_model.cpp @endlink
    //! @author weinerbe
    //! @date May 14, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API Model
    {
    private:

      //! @brief default constructor
      Model();

    public:

      //! @brief given a FILENAME, the model path is prepended to the filename
      //! @param FILENAME the filename to a model that is used in one of the scores
      //! @return string with model/FILENAME
      static std::string AddModelPath( const std::string &FILENAME);

      //! @brief get the model path flag
      static util::ShPtr< command::FlagInterface> &GetModelPathFlag();
    }; // class Model

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_H_
