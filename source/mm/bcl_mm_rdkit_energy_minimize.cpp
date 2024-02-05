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
#include "mm/bcl_mm_rdkit_energy_minimize.h"
#include "mm/bcl_mm_rdkit_force_field_utils.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"

namespace bcl
{
  namespace mm
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RdkitEnergyMinimize::RdkitEnergyMinimize() :
      RDKitEnergy(),
      m_MaxIterations( 1000),
      m_ForceTolerance( 1.0e-4),
      m_EnergyTolerance( 1.0e-4),
      m_PositionalRestraintAtoms( storage::Vector< size_t>()),
      m_PositionalRestraintAtomsString( ""),
      m_MaxUnrestrainedDisplacement( storage::Vector< double>()),
      m_RestraintForce( storage::Vector< double>())
    {
    }

    //! @brief full constructor with restrained atoms string
    RdkitEnergyMinimize::RdkitEnergyMinimize
    (
      const RDKitEnergy &ENERGY,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE,
      const std::string &POSITION_RESTRAINED_ATOMS_STRING,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    ) :
      RDKitEnergy( ENERGY),
      m_MaxIterations( MAX_ITERATIONS),
      m_ForceTolerance( FORCE_TOLERANCE),
      m_EnergyTolerance( ENERGY_TOLERANCE),
      m_PositionalRestraintAtoms( storage::Vector< size_t>()),
      m_PositionalRestraintAtomsString( POSITION_RESTRAINED_ATOMS_STRING),
      m_MaxUnrestrainedDisplacement( MAX_UNRESTRAINED_DISPLACEMENT),
      m_RestraintForce( RESTRAINT_FORCE)
    {
      // set the positionally restrained atom indices in RISH
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief full constructor with directly specified restrained atom indices
    RdkitEnergyMinimize::RdkitEnergyMinimize
    (
      const RDKitEnergy &ENERGY,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE,
      const storage::Vector< size_t> &POSITION_RESTRAINED_ATOMS,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    ) :
      RDKitEnergy( ENERGY),
      m_MaxIterations( MAX_ITERATIONS),
      m_ForceTolerance( FORCE_TOLERANCE),
      m_EnergyTolerance( ENERGY_TOLERANCE),
      m_PositionalRestraintAtoms( POSITION_RESTRAINED_ATOMS),
      m_PositionalRestraintAtomsString( ""),
      m_MaxUnrestrainedDisplacement( MAX_UNRESTRAINED_DISPLACEMENT),
      m_RestraintForce( RESTRAINT_FORCE)
    {
    }

    //! virtual copy constructor
    RdkitEnergyMinimize *RdkitEnergyMinimize::Clone() const
    {
      return new RdkitEnergyMinimize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RdkitEnergyMinimize::GetAlias() const
    {
      static const std::string s_name( "RDKitEnergyMinimizer"); // TODO: change name based on force field
      return s_name;
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RdkitEnergyMinimize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the maximum number of iterations to perform for geometry optimization
    //! @returns the number of iterations
    size_t RdkitEnergyMinimize::GetMaxIterations() const
    {
      return m_MaxIterations;
    }

    //! @brief returns the force tolerance for the geometry optimization
    //! @returns the force tolerance
    double RdkitEnergyMinimize::GetForceTolerance() const
    {
      return m_ForceTolerance;
    }

    //! @brief returns the energy tolerance for the geometry optimization
    //! @returns the energy tolerance
    double RdkitEnergyMinimize::GetEnergyTolerance() const
    {
      return m_EnergyTolerance;
    }

    //! @brief returns the atoms (by index) to which positional restraints are applied
    storage::Vector< size_t> RdkitEnergyMinimize::GetPositionalRestraintAtoms() const
    {
      return m_PositionalRestraintAtoms;
    }

    //! @brief returns the atoms (by string) to which positional restraints are applied
    std::string RdkitEnergyMinimize::GetPositionalRestraintAtomsString() const
    {
      return m_PositionalRestraintAtomsString;
    }

    //! @brief returns the maximum displacement each atom can experience before restraint force activates
    storage::Vector< double> RdkitEnergyMinimize::GetMaxUnrestrainedDisplacement() const
    {
      return m_MaxUnrestrainedDisplacement;
    }

    //! @brief returns the restraint force applied to each atom
    storage::Vector< double> RdkitEnergyMinimize::GetRestraintForce() const
    {
      return m_RestraintForce;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set the maximum number of iterations to perform for geometry optimization
    void RdkitEnergyMinimize::SetMaxIterations( const size_t &MAX_ITERATIONS)
    {
      m_MaxIterations = MAX_ITERATIONS;
    }

    //! @brief set the force tolerance for the geometry optimization
    void RdkitEnergyMinimize::SetForceTolerance( const double FORCE_TOLERANCE)
    {
      m_ForceTolerance = FORCE_TOLERANCE;
    }

    //! @brief set the energy tolerance for the geometry optimization
    void RdkitEnergyMinimize::SetEnergyTolerance( const double ENERGY_TOLERANCE)
    {
      m_EnergyTolerance = ENERGY_TOLERANCE;
    }

    //! @brief sets the atoms to which positional restraints are applied
    void RdkitEnergyMinimize::SetPositionalRestraintAtoms( storage::Vector< size_t> &ATOMS)
    {
      m_PositionalRestraintAtoms = ATOMS;
    }

    //! @brief sets the atoms to which positional restraints are applied
    void RdkitEnergyMinimize::SetPositionalRestraintAtomsFromString( std::string &ATOMS)
    {
      if( m_PositionalRestraintAtomsString.size())
      {
        m_PositionalRestraintAtoms.Reset();
        m_PositionalRestraintAtoms = util::SplitStringToNumerical< size_t>( m_PositionalRestraintAtomsString);
      }
      else
      {
        BCL_MessageStd
        (
          "[WARNING] RdkitEnergyMinimize::SetPositionalRestraintAtomsFromString "
          "string is empty - the restrained atoms set will not be modified!"
        );
      }
    }

    //! @brief sets the atoms to which positional restraints are applied
    void RdkitEnergyMinimize::SetPositionalRestraintAtomsString( std::string &ATOMS)
    {
      m_PositionalRestraintAtomsString = ATOMS;
    }

    //! @brief sets the maximum displacement each atom can experience before restraint force activates
    void RdkitEnergyMinimize::SetMaxUnrestrainedDisplacement( storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT)
    {
      m_MaxUnrestrainedDisplacement = MAX_UNRESTRAINED_DISPLACEMENT;
    }

    //! @brief sets the restraint force applied to each atom
    void RdkitEnergyMinimize::SetRestraintForce( storage::Vector< double> &RESTRAINT_FORCE)
    {
      m_RestraintForce = RESTRAINT_FORCE;
    }

    //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
    //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
    //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
    //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
    //! is the final energy of the optimized geometry
    storage::Pair< int, double> RdkitEnergyMinimize::OptimizeGeometry
    (
      chemistry::FragmentComplete &MOLECULE
    ) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // build force field
      ::ForceFields::ForceField *ff
      (
        RdkitForceFieldUtils::ConstructForceField
        (
          *rdkit_mol,
          GetForceFieldString(),
          m_NonbondedThreshold,
          m_IgnoreInterFragmentInteractions,
          true
        )
      );

      // add constraints
      if( m_PositionalRestraintAtoms.GetSize())
      {
        if( GetForceFieldEnum() == e_MMFF94 || GetForceFieldEnum() == e_MMFF94s)
        {
          RdkitForceFieldUtils::AddPositionalRestraintsMMFF( ff, m_PositionalRestraintAtoms, m_MaxUnrestrainedDisplacement, m_RestraintForce);
        }
        else // assume only other option is UFF
        {
          RdkitForceFieldUtils::AddPositionalRestraintsUFF( ff, m_PositionalRestraintAtoms, m_MaxUnrestrainedDisplacement, m_RestraintForce);
        }
      }

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
    storage::Triplet< chemistry::FragmentComplete, int, double> RdkitEnergyMinimize::OptimizeGeometry
    (
      const chemistry::FragmentComplete &MOLECULE
    ) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // build force field
      ::ForceFields::ForceField *ff
      (
        RdkitForceFieldUtils::ConstructForceField
        (
          *rdkit_mol,
          GetForceFieldString(),
          m_NonbondedThreshold,
          m_IgnoreInterFragmentInteractions,
          true
        )
      );

      // add constraints
      if( m_PositionalRestraintAtoms.GetSize())
      {
        if( GetForceFieldEnum() == e_MMFF94 || GetForceFieldEnum() == e_MMFF94s)
        {
          RdkitForceFieldUtils::AddPositionalRestraintsMMFF( ff, m_PositionalRestraintAtoms, m_MaxUnrestrainedDisplacement, m_RestraintForce);
        }
        else // assume only other option is UFF
        {
          RdkitForceFieldUtils::AddPositionalRestraintsUFF( ff, m_PositionalRestraintAtoms, m_MaxUnrestrainedDisplacement, m_RestraintForce);
        }
      }

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
    //! @param FF whether to use MMFF94 or MMFF94s
    //! @param NON_BONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
    //! @param IGNORE_INTER_FRAG_INTERACTIONS If true, nonbonded terms will not be added between fragments
    //! @param MAX_ITERATIONS maximum number of iterations in which to achieve convergence
    //! @param FORCE_TOLERANCE convergence criterion for force
    //! @param ENERGY_TOLERANCE convergence criterion for energy
    //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
    //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
    //! is the final energy of the optimized geometry
    storage::Pair< int, double> RdkitEnergyMinimize::OptimizeGeometry
    (
      chemistry::FragmentComplete &MOLECULE,
      const std::string &FF,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE,
      const storage::Vector< size_t> &POSITIONAL_RESTRAINT_ATOMS,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    )
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // build force field
      ::ForceFields::ForceField *ff
      (
        RdkitForceFieldUtils::ConstructForceField
        (
          *rdkit_mol,
          FF,
          NON_BONDED_THRESHOLD,
          IGNORE_INTER_FRAG_INTERACTIONS,
          true
        )
      );

      // add constraints
      if( POSITIONAL_RESTRAINT_ATOMS.GetSize())
      {
        if( RdkitForceFieldUtils::GetRdkitForceFieldsEnum( FF) == e_MMFF94 || RdkitForceFieldUtils::GetRdkitForceFieldsEnum( FF) == e_MMFF94s)
        {
          RdkitForceFieldUtils::AddPositionalRestraintsMMFF( ff, POSITIONAL_RESTRAINT_ATOMS, MAX_UNRESTRAINED_DISPLACEMENT, RESTRAINT_FORCE);
        }
        else // assume only other option is UFF
        {
          RdkitForceFieldUtils::AddPositionalRestraintsUFF( ff, POSITIONAL_RESTRAINT_ATOMS, MAX_UNRESTRAINED_DISPLACEMENT, RESTRAINT_FORCE);
        }
      }

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
    storage::Triplet< chemistry::FragmentComplete, int, double> RdkitEnergyMinimize::OptimizeGeometry
    (
      const chemistry::FragmentComplete &MOLECULE,
      const std::string &FF,
      const double NON_BONDED_THRESHOLD,
      const bool IGNORE_INTER_FRAG_INTERACTIONS,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE,
      const storage::Vector< size_t> &POSITIONAL_RESTRAINT_ATOMS,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    )
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // build force field
      ::ForceFields::ForceField *ff
      (
        RdkitForceFieldUtils::ConstructForceField
        (
          *rdkit_mol,
          FF,
          NON_BONDED_THRESHOLD,
          IGNORE_INTER_FRAG_INTERACTIONS,
          true
        )
      );

      // add constraints
      if( POSITIONAL_RESTRAINT_ATOMS.GetSize())
      {
        if( RdkitForceFieldUtils::GetRdkitForceFieldsEnum( FF) == e_MMFF94 || RdkitForceFieldUtils::GetRdkitForceFieldsEnum( FF) == e_MMFF94s)
        {
          RdkitForceFieldUtils::AddPositionalRestraintsMMFF( ff, POSITIONAL_RESTRAINT_ATOMS, MAX_UNRESTRAINED_DISPLACEMENT, RESTRAINT_FORCE);
        }
        else // assume only other option is UFF
        {
          RdkitForceFieldUtils::AddPositionalRestraintsUFF( ff, POSITIONAL_RESTRAINT_ATOMS, MAX_UNRESTRAINED_DISPLACEMENT, RESTRAINT_FORCE);
        }
      }

      // minimization
      int min_result( ff->minimize( MAX_ITERATIONS, FORCE_TOLERANCE, ENERGY_TOLERANCE));
      double energy( ff->calcEnergy());

      // generate a new molecule from the optimized rdkit molecule
      chemistry::FragmentComplete optimized_molecule( *chemistry::RdkitMolUtils::RDKitROMolToFragmentComplete( *rdkit_mol));

      // return calculation result
      return storage::Triplet< chemistry::FragmentComplete, int, double>( optimized_molecule, min_result, energy);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RdkitEnergyMinimize::GetSerializer() const
    {
      io::Serializer parameters( RDKitEnergy::GetSerializer());
      parameters.SetClassDescription( "Optimizes the geometry of a molecule using the MMFF94(s) force field with or without restraints.");
      parameters.AddInitializer
      (
        "max_iterations",
        "The maximum number of iterations to perform during geometry optimization, above which the optimization will terminate "
        "even if convergence is not yet reached.",
        io::Serialization::GetAgent( &m_MaxIterations),
        "1000"
      );
      parameters.AddInitializer
      (
        "force_tolerance",
        "The convergence criterion for forces",
        io::Serialization::GetAgent( &m_ForceTolerance),
        "1.0e-4"
      );
      parameters.AddInitializer
      (
        "energy_tolerance",
        "The convergence criterion for energies",
        io::Serialization::GetAgent( &m_EnergyTolerance),
        "1.0e-4"
      );
      parameters.AddInitializer
      (
        "coordinate_restraint_atoms",
        "Atoms that will be restrained to their current positions in Cartesian coordinate space",
        io::Serialization::GetAgent( &m_PositionalRestraintAtomsString),
        ""
      );

      return parameters;
    }

      //! @brief Set the members of this property from the given LABEL; override from SerializableInterface
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool RdkitEnergyMinimize::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        // parse string encoding atoms to be restrained in coordinate space
        SetPositionalRestraintAtomsFromString( m_PositionalRestraintAtomsString);
        return true;
      }

  } // namespace mm
} // namespace bcl
