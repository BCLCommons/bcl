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

#ifndef BCL_RESTRAINT_HANDLER_INTERFACE_H_
#define BCL_RESTRAINT_HANDLER_INTERFACE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "score/bcl_score.fwd.hh"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_protocol_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerInterface
    //! @brief interface class for objects that read restraints
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Oct 31, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HandlerInterface :
      public util::SerializableInterface
    {

    private:

      //! the extension used to identify files containing data for this restraint
      std::string m_Extension;
      std::string m_DefaultExtension;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief Constructor, sets the default extension
      HandlerInterface( const std::string &DEFAULT_EXT = "") :
        m_DefaultExtension( DEFAULT_EXT)
      {
      }

      //! @brief Clone function
      virtual HandlerInterface *Clone() const = 0;

      //! @brief gives the filename postfix that indicates the file contains data for this restraint type
      //! @return string - the filename postfix that indicates the file contains data for this restraint type
      const std::string &GetFilenameExtension() const
      {
        return m_Extension.empty() ? m_DefaultExtension : m_Extension;
      }

      //! @brief get the restraints filename
      //! @return filename for this restraints; uses -restraint_prefix flag and GetFilenameExtension()
      std::string GetFilename() const;

      //! @brief test for the existence of the restraints file
      bool Exists() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class HandlerInterface

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_HANDLER_INTERFACE_H_
