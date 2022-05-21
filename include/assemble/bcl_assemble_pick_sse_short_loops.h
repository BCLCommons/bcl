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

#ifndef BCL_ASSEMBLE_PICK_SSE_SHORT_LOOPS_H_
#define BCL_ASSEMBLE_PICK_SSE_SHORT_LOOPS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_pick_criteria_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickSSEShortLoops
    //! @brief class is used for picking sses with short loops from a list of sses to sses in the protein model
    //! @details This picker class collects, from a given list of SSEs (that are from the SSEPool), SSEs with short loops
    //! (less than m_MaxShortLoopLength) to any SSE in the given ProteinModel, and then picks one randomly from
    //! these collected SSEs with short loops.
    //!
    //! @see @link example_assemble_pick_sse_short_loops.cpp @endlink
    //! @author karakam
    //! @date 03/27/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickSSEShortLoops :
      public find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! maximum loop length for defining short loops
      size_t m_MaxShortLoopLength;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @brief MAX_LOOP_LENGTH maximum number of residues between SSEs to be classified as short loop ( 5 by defult)
      PickSSEShortLoops( const size_t MAX_LOOP_LENGTH = 5);

      //! virtual copy constructor
      PickSSEShortLoops *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
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

      //! @brief Picks one of SSEs with the Short Loops to the given domain
      //! @param SSE_LIST is the SiPtrList which provides the pool of assemble::SSE to pick from
      //! @param SSE_DOMAIN is the domain to which SSEs from the SSE_LIST will be compared for short loops
      //! @return returns SiPtr to the assemble::SSE object which has a short loops to one of the SSEs in SSE_DOMAIN
      util::SiPtr< const SSE>
      Pick( const util::SiPtrList< const SSE> &SSE_LIST, const DomainInterface &SSE_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PickSSEShortLoops

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_PICK_SSE_SHORT_LOOPS_H_
