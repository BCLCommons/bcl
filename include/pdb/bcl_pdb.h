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

#ifndef BCL_PDB_H_
#define BCL_PDB_H_

// include the namespace forward header
#include "bcl_pdb.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_pdb.h
  //! @brief namespace for classes reading/writing and representing PDB files that contain protein structures
  //! @details This namespace provides a handler class for reading pdbs, a factory for creating a protein model or chain
  //! from a handler, as well as components used in this process to represents residues and pdb line formats
  //!
  //! @see @link example_pdb.cpp @endlink
  //! @author woetzen
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace pdb
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

    //! @brief file extension for pdb file
    //! @return the file extension of a pdb file
    BCL_API
    const std::string &GetDefaultFileExtension();

  } // namespace pdb
} // namespace bcl

#endif //BCL_PDB_H_
