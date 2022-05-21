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

#ifndef BCL_FIND_PICK_CRITERIA_INTERFACE_H_
#define BCL_FIND_PICK_CRITERIA_INTERFACE_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickCriteriaInterface
    //! @brief Class is used for picking a single item a collection of substituents (e.x. all CA atoms) from an
    //! argument (e.x. protein model).
    //!
    //! @tparam t_ReturnType is the type of substituent which will be chosen
    //! @tparam t_ArgumentType is the type of argument from which t_ReturnType will be chosen
    //! @tparam t_Criteria is the criteria to be used when choosing t_ReturnType
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date 03/27/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType, typename t_CriteriaType>
    class PickCriteriaInterface :
      public util::SerializableInterface
    {

    public:

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

      //! @brief virtual copy constructor
      virtual PickCriteriaInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      virtual const std::string &GetAlias() const
      {
        static const std::string s_alias( "PickCriteriaInterface");
        return s_alias;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Picks a single item from a collection.");
        // serializer.AddInitializer
        //   (
        //    "temperature",
        //    "temperature is kept constant at this value",
        //    io::Serialization::GetAgent( &m_Temperature)
        //    );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! Pick the t_ReturnType object from t_ArgumentType out of many t_ReturnType objects
      //! @param CRITERIA contains the picking criteria
      //! @param ARGUMENT entity that contains t_ReturnType
      //! @return returns SiPtr to the chosen t_ReturnType
      virtual t_ReturnType
      Pick( const t_ArgumentType &ARGUMENT, const t_CriteriaType &CRITERIA) const = 0;

    }; // template class PickCriteriaInterface

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_PICK_CRITERIA_INTERFACE_H_
