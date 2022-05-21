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

#ifndef BCL_RESTRAINT_HANDLER_BODY_H_
#define BCL_RESTRAINT_HANDLER_BODY_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerBody
    //! @brief HandlerBody is for creating restraint::Bodies from a file.
    //!
    //! @see @link example_restraint_handler_body.cpp @endlink
    //! @author alexanns
    //! @date 07/06/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API HandlerBody :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! ShPtr to a function interface which determines how the occupancy of a restraint body is determined
      util::ShPtr
      <
        util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool>
      > m_DetermineOccupancy;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      HandlerBody();

      //! @brief constructor from a ShPtr to a FunctionInterface defining how occupancy is determined
      //! @param OCCUPANCY ShPtr to a FunctionInterface defining how occupancy is determined
      HandlerBody
      (
        const util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> > &OCCUPANCY
      );

      //! @brief copy constructor
      HandlerBody *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief CreateRestraintsBody is the function which creates the Restraints from an istream
      //! @param ISTREAM is the istream from which the restraints will be created
      //! @return returns a ShPtrVector of RestraintInterfaces
      util::ShPtrVector< Body>
      CreateRestraintsBody( std::istream &ISTREAM) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read restraint from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write restraint to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief CreateRestraintsBody is the function which creates the Restraints from an istream
      //! @param ISTREAM is the istream from which the restraints will be created
      //! @return returns a ShPtrVector of RestraintInterfaces
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > GenerateDensityRods( std::istream &ISTREAM) const;

    }; //class HandlerBody

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_HANDLER_BODY_H_
