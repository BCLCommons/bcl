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

#ifndef BCL_RESTRAINT_HANDLER_BASE_H_
#define BCL_RESTRAINT_HANDLER_BASE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_handler_interface.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerBase
    //! @brief base class for objects that read restraints of the specified template type
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Oct 31, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_RestraintType>
    class HandlerBase :
      public HandlerInterface
    {

    public:

      //! @brief constructor, accepts default file extension (passed to base class)
      //! @param DEFAULT_EXT the default file extension for this class instance
      HandlerBase( const std::string &DEFAULT_EXT = "") :
        HandlerInterface( DEFAULT_EXT)
      {
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief Clone function
      virtual HandlerBase *Clone() const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief reads restraints formatted for this restraint type from an istream
      //! @brief input stream to read the restraints from
      //! @return the read in restraints
      virtual t_RestraintType ReadRestraints( std::istream &ISTREAM) const = 0;

      //! @brief Reads the restraints from -restrain_file_prefix + file_extension
      //! @return the read in restraints
      virtual t_RestraintType ReadRestraintsFromFile() const
      {
        io::IFStream read;
        io::File::MustOpenIFStream( read, HandlerInterface::GetFilename());
        return ReadRestraints( read);
      }

    }; // class HandlerBase

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_HANDLER_BASE_H_
