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

#ifndef BCL_RESTRAINT_BODY_H_
#define BCL_RESTRAINT_BODY_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Body
    //! @brief This is a class for holding bodies for use as restraints.
    //!
    //! @see @link example_restraint_body.cpp @endlink
    //! @author alexanns, karakam
    //! @date 06/12/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Body :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! is the a ShPtr to the list of bodies that is the restraint
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > m_Bodies;

      //! is the method to determine whether a Body is occupied by a SSE
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
      Body();

      //! @brief construct from bodies
      //! @param BODIES the ShPtrVector of bodies which will serve as restraints
      //! @param DETERMINE_OCCUPANCY is a BinaryFunctionInterface which is used to determine if a restraint body is occupied by a SSE
      Body
      (
        const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &BODIES,
        const util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
          &DETERMINE_OCCUPANCY
      );

      //! @brief copy constructor
      Body *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetBody returns a const reference to "m_Bodies"
      //! @return returns a const reference to "m_Bodies"
      const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &GetBody() const;

      //! @brief SetBody "m_Bodies" to a new ShPtr of bodies
      //! @param BODIES is the ShPtrVector of bodies which "m_Bodies" will be changed to
      void SetBody( const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &BODIES);

    ///////////////
    // operation //
    ///////////////

      //! @brief GetUnoccupied determines which of the bodies in "m_Bodies" is not occupied
      //!        The argument could be the SSEs of a protein model, and GetUnoccupied will determine which restraint
      //!        bodies are not occupied by the bodies of the protein model and return this list of restraint bodies
      //!        which are not occupied by bodies of the protein model
      //! @param SSES a ShPtrVector of SSE which could occupy "m_Bodies"
      //! @return returns a ShPtrVector of bodies which are not occupied by any of "BODIES"
      util::ShPtrVector< assemble::SSEGeometryInterface> GetUnoccupied( const util::SiPtrVector< const assemble::SSE> &SSES) const;

      //! @brief get body that is occupied by given SSE
      //! @param SSE the sse to be considered to identify the body that is occupied by it
      //! @return ShPtr to occuoied body - will be undefined if there is non
      util::ShPtr< assemble::SSEGeometryInterface> GetOccupiedBody( const assemble::SSE &SSE) const;

      //! @brief GenerateAssignment creates the assignment of "m_Bodies" with other assembel::SSEs
      //! @param SSES SiPtrVector of SSE which will be assigned with "m_Bodies"
      //! @return returns an Assignment which assigns "m_Bodies" with the appropriate members of "SSES"
      const SSEAssignment
      GenerateAssignment( const util::SiPtrVector< const assemble::SSE> &SSES) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read distance from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write distance to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class restraint::Body

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_BODY_H_
