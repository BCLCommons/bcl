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

#ifndef BCL_ASSEMBLE_PICK_SSE_RANDOM_H_
#define BCL_ASSEMBLE_PICK_SSE_RANDOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "find/bcl_find_pick_interface.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickSSERandom
    //! @brief class for picking a random SSE
    //! @details This class picks a single SSE randomly from a given SiPtrList or SiPtrVector of SSEs.
    //!
    //! @todo since this class is only a special case of @link bcl_assemble_pick_sses_random.cpp @endlink it should be
    //! replaced by the class for the general case
    //!
    //! @see @link example_assemble_pick_sse_random.cpp @endlink
    //! @author karakam
    //! @date 03/27/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickSSERandom :
      virtual public find::PickInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE> >
    {

    private:

    //////////
    // data //
    //////////

      //! set of SSTypes to be picked from, if not be used it is empty
      storage::Set< biol::SSType> m_SSTypes;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PickSSERandom();

      //! @brief constructor from a set of SSTypes
      //! @param SS_TYPES Set of SSTypes
      PickSSERandom( const storage::Set< biol::SSType> &SS_TYPES);

      //! virtual copy constructor
      PickSSERandom *Clone() const;

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

      //! Picks the SSE object
      //! @param SSE_LIST is the SiPtrList which provides the pool of SSE to pick from
      //! @return returns SiPtr to the SSE object
      util::SiPtr< const SSE> Pick( const util::SiPtrList< const SSE> &SSE_LIST) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PickSSERandom

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_PICK_SSE_RANDOM_H_
