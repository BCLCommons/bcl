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

#ifndef BCL_MM_RDKIT_ENERGY_H_
#define BCL_MM_RDKIT_ENERGY_H_

// include the namespace header
#include "bcl_mm.h"
#include "mm/bcl_mm_rdkit_force_fields.h"

// include other forward headers - sorted alphabetically
#include "bcl_mm_energy_interface.h"
#include "io/bcl_io_serialization.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "util/bcl_util_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mm
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RDKitEnergy
    //! @brief This class computes molecular mechanics energies using the UFF, MMFF94, or MMFF94s force fields
    //! as implemented in RDKit. Working energy unit is kcal/mol. Any energy conversion is done after the calculation.
    //!
    //! @see @link example_mm_rdkit_energy.cpp @endlink
    //! @author brownbp1
    //! @date Feb 04, 2024
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RDKitEnergy :
      public EnergyInterface
    {
    //////////
    // data //
    //////////

    public:

    protected:

      //! The MMFF variant that we are using
      RdkitForceFieldsEnum m_ForceFieldEnum;
      std::string m_ForceFieldString;

      //! The threshold to be used in adding non-bonded terms to the force field.
      //! Any non-bonded contact whose current distance is greater than nonBondedThresh * the minimum value for that contact
      //! will not be included.
      double m_NonbondedThreshold;

      //!  If true, nonbonded terms will not be added between fragments
      bool m_IgnoreInterFragmentInteractions;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      RDKitEnergy();

      //! full constructor
      RDKitEnergy
      (
        const RdkitForceFieldsEnum &VARIANT,
        const double NON_BONDED_THRESHOLD,
        const bool IGNORE_INTER_FRAG_INTERACTIONS
      );

      //! virtual copy constructor
      RDKitEnergy *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the force field variant
      RdkitForceFieldsEnum GetForceFieldEnum() const;

      //! @brief returns the force field variant as a string
      std::string GetForceFieldString() const;

      //! @brief returns the non-bonded threshold
      double GetNonbondedThreshold() const;

      //! @brief returns whether to ignore fragment interactions
      bool GetIgnoreInterFragmentInteractions() const;

    ///////////////////
    //   operations  //
    ///////////////////

      //! @brief sets the force field variant
      void SetForceFieldFromEnum( const RdkitForceFieldsEnum &VARIANT);

      //! @brief sets the force field variant from a string
      void SetForceFieldFromString( const std::string &VARIANT);

      //! @brief sets the force field string
      void SetForceFieldString( const std::string &VARIANT);

      //! @brief sets the non-bonded threshold
      void SetNonbondedThreshold( const double THRESHOLD);

      //! @brief sets whether to ignore fragment interactions
      void SetIgnoreInterFragmentInteractions( const bool IGNORE_INTER_FRAGMENT_INTERACTIONS);

      //! @brief Computes the force field potential energy of a molecule
      //! @param MOLECULE the molecule for which the energy will be computed
      //! @return the energy of MOLECULE
      double CalculateEnergy( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief Computes the force field potential energy of a molecule
      //! @param MOLECULE the molecule for which the energy will be computed
      //! @param VARIANT whether to use UFF, MMFF94, or MMFF94s
      //! @param NON_BONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
      //! @param IGNORE_INTER_FRAG_INTERACTIONS If true, nonbonded terms will not be added between fragments
      //! @return the energy of MOLECULE
      static double CalculateEnergy
      (
        const chemistry::FragmentComplete &MOLECULE,
        const std::string &VARIANT,
        const double NON_BONDED_THRESHOLD,
        const bool IGNORE_INTER_FRAG_INTERACTIONS
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    };

  } // namespace mm
} // namespace bcl

#endif // BCL_MM_RDKIT_ENERGY_H_
