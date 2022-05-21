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

#ifndef BCL_CHEMISTRY_ATOM_WITH_POSITION_INTERFACE_H_
#define BCL_CHEMISTRY_ATOM_WITH_POSITION_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_type_data.h"
#include "coord/bcl_coord_orientation_interface.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomWithPositionInterface
    //! @brief Basic methods implemented in atom conformation
    //! @details includes the 3D coordinate of atom
    //!
    //! @see @link example_chemistry_atom_with_position_interface.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 05, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomWithPositionInterface :
      public coord::OrientationInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AtomWithConfigurationWithPosition
      virtual AtomWithPositionInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns a reference of atom type data attribute
      //! @return a reference of atom type data attribute
      virtual const AtomType &GetAtomType() const = 0;

      //! @brief return a reference to the position vector
      //! @return a reference to the position vector
      virtual const linal::Vector3D &GetPosition() const = 0;

      //! @brief return a reference to the element type of the atom
      //! @return the element type of the atom
      const ElementType &GetElementType() const
      {
        return GetAtomType()->GetElementType();
      }

      //! @brief get the charge
      //! @return the charge of the atom
      const short &GetCharge() const
      {
        return GetAtomType()->GetFormalCharge();
      }

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const
      {
        return GetPosition();
      }

    }; // class AtomWithPositionInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_WITH_POSITION_INTERFACE_H_ 
