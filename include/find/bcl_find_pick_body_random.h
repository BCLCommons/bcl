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

#ifndef BCL_FIND_PICK_BODY_RANDOM_H_
#define BCL_FIND_PICK_BODY_RANDOM_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_find_pick_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickBodyRandom
    //! @brief gets a random assemble::SSEGeometryInterface out of many assemble::SSEGeometryInterface.
    //!
    //! @see @link example_find_pick_body_random.cpp @endlink
    //! @author alexanns
    //! @date 06/30/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickBodyRandom :
      public PickInterface< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface> >
    {

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PickBodyRandom();

      //! virtual copy constructor
      PickBodyRandom *Clone() const;

      //! virtual destructor
      ~PickBodyRandom();

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

      //! Pick the util::ShPtr< assemble::SSEGeometryInterface> object from restraint::Body based on an SSE
      //! @param BODIES is the object which will provide the pool of util::ShPtr< assemble::SSEGeometryInterface> to pick from
      //! @return return a util::ShPtr< assemble::SSEGeometryInterface> which can accomodate "SSE_CRITERIA"
      util::ShPtr< assemble::SSEGeometryInterface>
      Pick( const util::ShPtrVector< assemble::SSEGeometryInterface> &BODIES) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PickBodyRandom

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_PICK_BODY_RANDOM_H_
