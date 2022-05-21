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

#ifndef BCL_FOLD_PLACEMENT_INTERFACE_H_
#define BCL_FOLD_PLACEMENT_INTERFACE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementInterface
    //! @brief template class for placing t_ObjectType with respect to t_ArgumentType
    //! @details This interface provides operators that calculates the TransformationMatrix3D for placing a given t_ObjectType
    //! with respect to the t_ArgumentType
    //!
    //! @tparam t_ObjectType type of the object that will be placed
    //! @tparam t_ArgumentType type of the argument in which t_ObjectType will be placed
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Apr 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ObjectType, typename t_ArgumentType>
    class PlacementInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual PlacementInterface< t_ObjectType, t_ArgumentType> *Clone() const = 0;

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
        serializer.SetClassDescription( "Places an object in a different object.");

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns placement for the t_ObjectType in t_ArgumentType
      //! @param OBJECT SiPtr to entity that contains a t_ObjectType
      //! @param ARGUMENT entity that contains a t_ObjectType
      //! @return the transformationmatrix3d to place OBJECT in ARGUMENT
      virtual storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const t_ObjectType &OBJECT,
        const t_ArgumentType &ARGUMENT
      ) const = 0;

      //! @brief returns placement for the t_ObjectType without taking t_ArgumentType into consideration
      //! @param OBJECT SiPtr to entity that contains a t_ObjectType
      //! @return the transformationmatrix3d to place OBJECT
      virtual storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const t_ObjectType &OBJECT
      ) const
      {
        // return
        return storage::Pair< math::TransformationMatrix3D, bool>( math::TransformationMatrix3D(), true);
      }

    }; // template class PlacementInterface

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_PLACEMENT_INTERFACE_H_
