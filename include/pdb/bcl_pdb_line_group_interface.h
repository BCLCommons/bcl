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

#ifndef BCL_PDB_LINE_GROUP_INTERFACE_H_
#define BCL_PDB_LINE_GROUP_INTERFACE_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LineGroupInterface
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Feb 20, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LineGroupInterface :
      public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief linetypes within group
      //! @return set of line types
      virtual const storage::Set< LineType> &GetTypesOfLines() const = 0;

      //! @brief access to lines of given type
      //! @param LINE_TYPE the desire line type
      //! @return lines of given type
      virtual util::ShPtrList< Line> GetLines( const LineType &LINE_TYPE) const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief locate lines of given criterium
      //! @param CRITERIUM unary predicate with operator return true for line, that should be considered
      //! @return lines that are considered by criterium
      virtual util::ShPtrList< Line> CollectLines( const util::FunctionInterface< Line, bool> &CRITERIUM) const = 0;

      //! @brief pushback a new line into that group
      //! @param LINE ShPtr to the line
      //! @return true, if it fits into that group (line type is eligible)
      virtual bool PushBack( const util::ShPtr< Line> &LINE) = 0;

      //! @brief reset the line group
      virtual void Reset() = 0;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @return outputstream which was written to
      virtual std::ostream &WriteLines( std::ostream &OSTREAM) const = 0;

    }; // class LineGroupInterface

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_LINE_GROUP_INTERFACE_H_ 
