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

#ifndef BCL_FIND_PICK_INTERFACE_H_
#define BCL_FIND_PICK_INTERFACE_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickInterface
    //! @brief Class is used for picking a single item a collection of substituents (e.x. all CA atoms) from an
    //! argument (e.x. protein model).
    //!
    //! @tparam t_ReturnType is the type of substituent which will be chosen
    //! @tparam t_ArgumentType is the type of argument from which t_ReturnType will be chosen
    //!
    //! @remarks example unnecessary
    //! @author alexanns, karakam
    //! @date 03/27/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType>
    class PickInterface :
      virtual public util::SerializableInterface
    {

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief clone function
      //! @return pointer to a new PickInterface
      virtual PickInterface< t_ReturnType, t_ArgumentType> *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! Pick the t_ReturnType object from t_ArgumentType out of many t_ReturnType objects
      //! @param ARGUMENT entity that contains t_ReturnTypes
      //! @return returns SiPtr to the chosen t_ReturnType
      virtual t_ReturnType Pick( const t_ArgumentType &ARGUMENT) const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "This class needs an implementation of GetSerializer()");

        return serializer;
      }

    }; // class PickInterface

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_PICK_INTERFACE_H_
