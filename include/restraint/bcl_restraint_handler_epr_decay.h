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

#ifndef BCL_RESTRAINT_HANDLER_EPR_DECAY_H_
#define BCL_RESTRAINT_HANDLER_EPR_DECAY_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_epr_decay.h"
#include "bcl_restraint_handler_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerEPRDecay
    //! @brief reads in EPR decay measurements from text files.
    //!
    //! @see @link example_restraint_handler_epr_decay.cpp @endlink
    //! @author fischea
    //! @date Nov 11, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API HandlerEPRDecay :
      public HandlerBase< storage::Vector< EPRDecay> >
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from defaults
      HandlerEPRDecay();

      //! @brief virtual copy constructor
      //! @return pointer to a new HandlerEPRDecay
      HandlerEPRDecay *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reads atom distance restraints from an input stream
      //! @brief input stream to read the restraints from
      //! @return the read in restraints
      storage::Vector< EPRDecay> ReadRestraints( std::istream &ISTREAM) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; //class HandlerEPRDecay

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_HANDLER_EPR_DECAY_H_
