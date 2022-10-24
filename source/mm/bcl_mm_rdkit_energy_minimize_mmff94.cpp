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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "mm/bcl_mm_rdkit_energy_minimize_mmff94.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"
#include "GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h"
#include "GraphMol/ForceFieldHelpers/MMFF/Builder.h"

namespace bcl
{
  namespace mm
  {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    //! virtual copy constructor
    RdkitEnergyMinimizeMmff94 *RdkitEnergyMinimizeMmff94::Clone() const
    {
      return new RdkitEnergyMinimizeMmff94( *this);
    }

    /////////////////
    // data access //
    /////////////////

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RdkitEnergyMinimizeMmff94::GetAlias() const
    {
      static const std::string s_name( m_MMFFVariant == e_MMFF94 ? "EnergyMinimization_MMFF94" : "EnergyMinimization_MMFF94s");
      return s_name;
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RdkitEnergyMinimizeMmff94::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the maximum number of iterations to perform for geometry optimization
    //! @returns the number of iterations
    size_t RdkitEnergyMinimizeMmff94::GetMaxIterations() const
    {
      return m_MaxIterations;
    }

    //! @brief returns the force tolerance for the geometry optimization
    //! @returns the force tolerance
    double RdkitEnergyMinimizeMmff94::GetForceTolerance() const
    {
      return m_ForceTolerance;
    }

    //! @brief returns the energy tolerance for the geometry optimization
    //! @returns the energy tolerance
    double RdkitEnergyMinimizeMmff94::GetEnergyTolerance() const
    {
      return m_EnergyTolerance;
    }
    ////////////////
    // operations //
    ////////////////

    //! @brief set the maximum number of iterations to perform for geometry optimization
    void RdkitEnergyMinimizeMmff94::SetMaxIterations( const size_t &MAX_ITERATIONS)
    {
      m_MaxIterations = MAX_ITERATIONS;
    }

    //! @brief set the force tolerance for the geometry optimization
    void RdkitEnergyMinimizeMmff94::SetForceTolerance( const double FORCE_TOLERANCE)
    {
      m_ForceTolerance = FORCE_TOLERANCE;
    }

    //! @brief set the energy tolerance for the geometry optimization
    void RdkitEnergyMinimizeMmff94::SetEnergyTolerance( const double ENERGY_TOLERANCE)
    {
      m_EnergyTolerance = ENERGY_TOLERANCE;
    }

    //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
    //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
    //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
    //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
    //! is the final energy of the optimized geometry
    storage::Pair< int, double> RdkitEnergyMinimizeMmff94::OptimizeGeometry( chemistry::FragmentComplete &MOLECULE) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, GetMMFFVariantAsString());
      if( !mmff_mol_properties.isValid())
      {
        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
        return storage::Pair< int, double>();
      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *rdkit_mol, m_NonbondedThreshold, -1, m_IgnoreInterFragmentInteractions);
      ff->initialize();

      // minimization
      int min_result( ff->minimize( m_MaxIterations, m_ForceTolerance, m_EnergyTolerance));
      double energy( ff->calcEnergy());

      // modify input molecule
      MOLECULE = *chemistry::RdkitMolUtils::RDKitROMolToFragmentComplete( *rdkit_mol);

      // return calculation result
      return storage::Pair< int, double>( min_result, energy);
    }

    //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
    //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
    //! @returns a triplet where the first value is the new minimized molecule, the second value indicates
    //! if the minimization was a success (0), failed to converge within the maximum number of iterations (1),
    //! or had missing parameters (-1), and where the second value is the final energy of the optimized geometry
    storage::Triplet< chemistry::FragmentComplete, int, double> RdkitEnergyMinimizeMmff94::OptimizeGeometry
    (
      const chemistry::FragmentComplete &MOLECULE
    ) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, GetMMFFVariantAsString());
      if( !mmff_mol_properties.isValid())
      {
        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
        return storage::Triplet< chemistry::FragmentComplete, int, double>();
      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *rdkit_mol, m_NonbondedThreshold, -1, m_IgnoreInterFragmentInteractions);
      ff->initialize();

      // minimization
      int min_result( ff->minimize( m_MaxIterations, m_ForceTolerance, m_EnergyTolerance));
      double energy( ff->calcEnergy());

      // generate a new molecule from the optimized rdkit molecule
      chemistry::FragmentComplete optimized_molecule( *chemistry::RdkitMolUtils::RDKitROMolToFragmentComplete( *rdkit_mol));

      // return calculation result
      return storage::Triplet< chemistry::FragmentComplete, int, double>( optimized_molecule, min_result, energy);
    }

    //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
    //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
    //! @param MMFF_VARIANT whether to use MMFF94 or MMFF94s
    //! @param NON_BONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
    //! @param IGNORE_INTER_FRAG_INTERACTIONS If true, nonbonded terms will not be added between fragments
    //! @param MAX_ITERATIONS maximum number of iterations in which to achieve convergence
    //! @param FORCE_TOLERANCE convergence criterion for force
    //! @param ENERGY_TOLERANCE convergence criterion for energy
    //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
    //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
    //! is the final energy of the optimized geometry
    storage::Pair< int, double> RdkitEnergyMinimizeMmff94::OptimizeGeometry
    (
      chemistry::FragmentComplete &MOLECULE,
      const std::string &MMFF_VARIANT,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE
    )
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, MMFF_VARIANT);
      if( !mmff_mol_properties.isValid())
      {
        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
        return storage::Pair< int, double>();
      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *rdkit_mol, NON_BONDED_THRESHOLD, -1, IGNORE_INTER_FRAG_INTERACTIONS);
      ff->initialize();

      // minimization
      int min_result( ff->minimize( MAX_ITERATIONS, FORCE_TOLERANCE, ENERGY_TOLERANCE));
      double energy( ff->calcEnergy());

      // modify input molecule
      MOLECULE = *chemistry::RdkitMolUtils::RDKitROMolToFragmentComplete( *rdkit_mol);

      // return calculation result
      return storage::Pair< int, double>( min_result, energy);
    }

    //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
    //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
    //! @returns a triplet where the first value is the new minimized molecule, the second value indicates
    //! if the minimization was a success (0), failed to converge within the maximum number of iterations (1),
    //! or had missing parameters (-1), and where the second value is the final energy of the optimized geometry
    storage::Triplet< chemistry::FragmentComplete, int, double> RdkitEnergyMinimizeMmff94::OptimizeGeometry
    (
      const chemistry::FragmentComplete &MOLECULE,
      const std::string &MMFF_VARIANT,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE
    )
    {

      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, MMFF_VARIANT);
      if( !mmff_mol_properties.isValid())
      {
        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
        return storage::Triplet< chemistry::FragmentComplete, int, double>();
      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *rdkit_mol, NON_BONDED_THRESHOLD, -1, IGNORE_INTER_FRAG_INTERACTIONS);
      ff->initialize();

      // minimization
      int min_result( ff->minimize( MAX_ITERATIONS, FORCE_TOLERANCE, ENERGY_TOLERANCE));
      double energy( ff->calcEnergy());

      // generate a new molecule from the optimized rdkit molecule
      chemistry::FragmentComplete optimized_molecule( *chemistry::RdkitMolUtils::RDKitROMolToFragmentComplete( *rdkit_mol));

      // return calculation result
      return storage::Triplet< chemistry::FragmentComplete, int, double>( optimized_molecule, min_result, energy);
    }

  } // namespace mm
} // namespace bcl
