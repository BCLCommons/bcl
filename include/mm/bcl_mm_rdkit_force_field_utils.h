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

#ifndef BCL_MM_RDKIT_FORCE_FIELD_UTILS_H_
#define BCL_MM_RDKIT_FORCE_FIELD_UTILS_H_

// include the namespace header
#include "bcl_mm.h"
#include "mm/bcl_mm_rdkit_force_fields.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"
#include "GraphMol/RWMol.h"

namespace bcl
{
  namespace mm
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RdkitForceFieldUtils
    //! @brief This class contains functions to manage calculations involving molecular mechanics force fields
    //! available through the RDKit library.
    //!
    //! @see @link example_mm_rdkit_force_field_utils.cpp @endlink
    //! @author brownbp1
    //! @date Feb 04, 2024
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RdkitForceFieldUtils :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief delete the default constructor
      RdkitForceFieldUtils() = delete;

      //! @brief virtual copy constructor
      RdkitForceFieldUtils *Clone() const;

      //! @brief returns class name
      //! the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief return the enum corresponding to the force field name
      //! @param NAME the string corresponding to the force field enum
      //! @return the force field enum
      static RdkitForceFieldsEnum GetRdkitForceFieldsEnum( const std::string &NAME);

      //! @brief return an MMFF force field
      //! @param RDMOL the rdkit mol to be parameterized
      //! @param MMFF_VARIANT specify MMFF94 or MMFF94s
      //! @param NONBONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
      //! Any non-bonded contact whose current distance is greater than nonBondedThresh * the minimum value for that contact
      //! will not be included.
      //! @param IGNORE_INTER_FRAGMENT_INTERACTIONS if true, nonbonded terms will not be added between fragments
      //! @param INITIALIZE if true, initialize the force field for use
      //! @return a force field
      static ::ForceFields::ForceField* ConstructForceField
      (
        ::RDKit::RWMol &RDMOL,
        const std::string &MMFF_VARIANT = "UFF",
        const double NONBONDED_THRESHOLD = 100.0,
        const bool IGNORE_INTER_FRAGMENT_INTERACTIONS = true,
        const bool INITIALIZE = true
      );

      //! @brief add positional constraints to force field for geometry optimization
      //! @param FORCE_FIELD the force field that is modified with the new restraint term
      //! @param ATOM_INDICES indices that are restrained during minimization
      //! @param MAX_UNRESTRAINED DISPLACEMENT coordinate displacement above which restraint force is applied
      //! @param RESTRAINT_FORCE restraint force
      static void AddPositionalRestraintsUFF
      (
        ::ForceFields::ForceField *FORCE_FIELD, // raw pointer unconventional for BCL outside of Clone(), but this is what RDKit requires
        const storage::Vector< size_t> &ATOM_INDICES,
        const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
        const storage::Vector< double> &RESTRAINT_FORCE
      );

      //! @brief add positional constraints to force field for geometry optimization
      //! @param FORCE_FIELD the force field that is modified with the new restraint term
      //! @param ATOM_INDICES indices that are restrained during minimization
      //! @param MAX_UNRESTRAINED DISPLACEMENT coordinate displacement above which restraint force is applied
      //! @param RESTRAINT_FORCE restraint force
      static void AddPositionalRestraintsMMFF
      (
        ::ForceFields::ForceField *FORCE_FIELD, // raw pointer unconventional for BCL outside of Clone(), but this is what RDKit requires
        const storage::Vector< size_t> &ATOM_INDICES,
        const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
        const storage::Vector< double> &RESTRAINT_FORCE
      );

    private:

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class RdkitForceFieldUtils

  } // namespace mm
} // namespace bcl

#endif // BCL_MM_RDKIT_FORCE_FIELD_UTILS_H_
