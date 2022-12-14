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

#ifndef BCL_MM_RDKIT_ENERGY_MINIMIZE_MMFF94_H_
#define BCL_MM_RDKIT_ENERGY_MINIMIZE_MMFF94_H_

// include the namespace header
#include "bcl_mm.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "mm/bcl_mm_rdkit_energy_mmff94.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "util/bcl_util_function_interface_serializable.h"

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"
#include "GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h"
#include "GraphMol/ForceFieldHelpers/MMFF/Builder.h"

namespace bcl
{
  namespace mm
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RdkitEnergyMinimizeMmff94
    //! @brief This class performs energy minimization (geometry optimization) using the MMFF94 or MMFF94s force fields
    //! implemented in RDKit. Working energy unit is kcal/mol. Any energy conversion is done after the calculation.
    //!
    //! @see @link example_mm_rdkit_energy_minimize_mmff94.cpp @endlink
    //! @author brownbp1
    //! @date Oct 23, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RdkitEnergyMinimizeMmff94 :
      public mm::RDKitEnergyMMFF94
    {

      //////////
      // data //
      //////////

    private:

      //! @brief Maximum number of iterations to perform for geometry optimization
      size_t m_MaxIterations;

      //! @brief The convergence criterion for forces
      double m_ForceTolerance;

      //! @brief The convergence criterion for energies
      double m_EnergyTolerance;

      //! @brief The atoms subject to coordinate restraints
      storage::Vector< size_t> m_PositionalRestraintAtoms;
      std::string m_PositionalRestraintAtomsString;

      //! @brief The maximum allowed displacement per-atom
      storage::Vector< double> m_MaxUnrestrainedDisplacement;

      //! @brief The per-atom restraint force
      storage::Vector< double> m_RestraintForce;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RdkitEnergyMinimizeMmff94();

      //! @brief full constructor with restrained atoms string
      RdkitEnergyMinimizeMmff94
      (
        const size_t MAX_ITERATIONS,
        const double FORCE_TOLERANCE,
        const double ENERGY_TOLERANCE,
        const std::string &POSITION_RESTRAINED_ATOMS_STRING,
        const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
        const storage::Vector< double> &RESTRAINT_FORCE
      );

      //! @brief full constructor with directly specified restrained atom indices
      RdkitEnergyMinimizeMmff94
      (
        const size_t MAX_ITERATIONS,
        const double FORCE_TOLERANCE,
        const double ENERGY_TOLERANCE,
        const storage::Vector< size_t> &POSITION_RESTRAINED_ATOMS,
        const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
        const storage::Vector< double> &RESTRAINT_FORCE
      );

      //! virtual copy constructor
      RdkitEnergyMinimizeMmff94 *Clone() const;

      /////////////////
      // data access //
      /////////////////

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the maximum number of iterations to perform for geometry optimization
      //! @returns the number of iterations
      size_t GetMaxIterations() const;

      //! @brief returns the force tolerance for the geometry optimization
      //! @returns the force tolerance
      double GetForceTolerance() const;

      //! @brief returns the energy tolerance for the geometry optimization
      //! @returns the energy tolerance
      double GetEnergyTolerance() const;

      //! @brief returns the atoms (by index) to which positional restraints are applied
      storage::Vector< size_t> GetPositionalRestraintAtoms() const;

      //! @brief returns the atoms (by index) to which positional restraints are applied
      storage::Vector< size_t> GetPositionalRestraintAtomsFromString( std::string &ATOMS) const;

      //! @brief returns the atoms (by string) to which positional restraints are applied
      std::string GetPositionalRestraintAtomsString() const;

      //! @brief returns the maximum displacement each atom can experience before restraint force activates
      storage::Vector< double> GetMaxUnrestrainedDisplacement() const;

      //! @brief returns the restraint force applied to each atom
      storage::Vector< double> GetRestraintForce() const;

      ////////////////
      // operations //
      ////////////////

      //! @brief set the maximum number of iterations to perform for geometry optimization
      void SetMaxIterations( const size_t &MAX_ITERATIONS);

      //! @brief set the force tolerance for the geometry optimization
      void SetForceTolerance( const double FORCE_TOLERANCE);

      //! @brief set the energy tolerance for the geometry optimization
      void SetEnergyTolerance( const double ENERGY_TOLERANCE);

      //! @brief sets the atoms to which positional restraints are applied
      void SetPositionalRestraintAtoms( storage::Vector< size_t> &ATOMS);

      //! @brief sets the atoms to which positional restraints are applied
      void SetPositionalRestraintAtomsFromString( std::string &ATOMS);

      //! @brief sets the atoms to which positional restraints are applied
      void SetPositionalRestraintAtomsString( std::string &ATOMS);

      //! @brief sets the maximum displacement each atom can experience before restraint force activates
      void SetMaxUnrestrainedDisplacement( storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT);

      //! @brief sets the restraint force applied to each atom
      void SetRestraintForce( storage::Vector< double> &RESTRAINT_FORCE);

    private:

      //! @brief add positional constraints to force field for geometry optimization
      //! @param FORCE_FIELD the force field that is modified with the new restraint term
      //! @param ATOM_INDICES indices that are restrained during minimization
      //! @param MAX_UNRESTRAINED DISPLACEMENT coordinate displacement above which restraint force is applied
      //! @param RESTRAINT_FORCE restraint force
      void AddPositionalRestraints
      (
        ::ForceFields::ForceField* FORCE_FIELD, // raw pointer unconventional for BCL outside of Clone(), but this is what RDKit requires
        const storage::Vector< size_t> &ATOM_INDICES,
        const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
        const storage::Vector< double> &RESTRAINT_FORCE
      ) const;

    public:

      //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
      //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
      //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
      //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
      //! is the final energy of the optimized geometry
      storage::Pair< int, double> OptimizeGeometry( chemistry::FragmentComplete &MOLECULE) const;

      //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
      //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
      //! @returns a triplet where the first value is the new minimized molecule, the second value indicates
      //! if the minimization was a success (0), failed to converge within the maximum number of iterations (1),
      //! or had missing parameters (-1), and where the second value is the final energy of the optimized geometry
      storage::Triplet< chemistry::FragmentComplete, int, double> OptimizeGeometry( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
      //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
      //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
      //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
      //! is the final energy of the optimized geometry
      static storage::Pair< int, double> OptimizeGeometry
      (
        chemistry::FragmentComplete &MOLECULE,
        const std::string &MMFF_VARIANT,
        const double NON_BONDED_THRESHOLD = 10.0,
        const bool IGNORE_INTER_FRAG_INTERACTIONS = true,
        const size_t MAX_ITERATIONS = 1000,
        const double FORCE_TOLERANCE = 1e-4,
        const double ENERGY_TOLERANCE = 1e-4
      );

      //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
      //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
      //! @returns a triplet where the first value is the new minimized molecule, the second value indicates
      //! if the minimization was a success (0), failed to converge within the maximum number of iterations (1),
      //! or had missing parameters (-1), and where the second value is the final energy of the optimized geometry
      static storage::Triplet< chemistry::FragmentComplete, int, double> OptimizeGeometry
      (
        const chemistry::FragmentComplete &MOLECULE,
        const std::string &MMFF_VARIANT,
        const double NON_BONDED_THRESHOLD = 10.0,
        const bool IGNORE_INTER_FRAG_INTERACTIONS = true,
        const size_t MAX_ITERATIONS = 1000,
        const double FORCE_TOLERANCE = 1e-4,
        const double ENERGY_TOLERANCE = 1e-4
      );

      //////////////////////
      // helper functions //
      //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL; override from SerializableInterface
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    };

  } // namespace mm
} // namespace bcl

#endif // BCL_MM_RDKIT_ENERGY_MINIMIZE_MMFF94_H_
