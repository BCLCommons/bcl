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

#ifndef BCL_OPTI_PRINT_INTERFACE_H_
#define BCL_OPTI_PRINT_INTERFACE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_tracker.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrintInterface
    //! @brief PrintInterface is a class, for writing information concerning the approximation process
    //! @details provides the common interface for any printer that can be used within optimization framework.
    //!
    //! @tparam t_ArgumentType argument to the optimization
    //! @tparam t_ResultType the type of the objective to optimize on
    //!
    //! @remarks example unnecessary
    //! @author fischea, mendenjl
    //! @date Dec 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class PrintInterface :
      public virtual util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new PrintInterface< t_ArgumentType, t_ResultType>
      virtual PrintInterface< t_ArgumentType, t_ResultType> *Clone() const = 0;

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetAlias() const
      {
        // this function should be removed once the remaining printeres are moved into the serializable
        // framework
        static const std::string s_alias( "ToBeRemoved");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      virtual void Print( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const = 0;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "This class should not be in the enum because it failed to override GetSerializer");
        return serializer;
      }

    }; // template class PrintInterface< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PRINT_INTERFACE_H_
