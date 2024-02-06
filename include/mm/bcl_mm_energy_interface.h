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

#ifndef BCL_MM_ENERGY_INTERFACE_H_
#define BCL_MM_ENERGY_INTERFACE_H_

// include the namespace header
#include "bcl_mm.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "util/bcl_util_function_interface_serializable.h"

namespace bcl
{
  namespace mm
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnergyInterface
    //! @brief This class is an interface class for classes that compute energies with molecular mechanics force fields
    //!
    //! @see @link example_mm_energy_interface.cpp @endlink
    //! @author brownbp1
    //! @date Oct 22, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EnergyInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      virtual EnergyInterface *Clone() const = 0;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      virtual const std::string &GetAlias() const = 0;

      //! @brief splits the molecule according to GetComponentVertices
      //! @param MOLECULE the molecule for which the energy will be computed
      //! @return the energy of MOLECULE
      virtual double CalculateEnergy( const chemistry::FragmentComplete &MOLECULE) const = 0;

    };

  } // namespace mm
} // namespace bcl

#endif // BCL_MM_ENERGY_INTERFACE_H_
