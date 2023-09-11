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
#include "mm/bcl_mm_rdkit_energy_minimize_uff.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include "ForceField/UFF/PositionConstraint.h"

namespace bcl
{
  namespace mm
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RdkitEnergyMinimizeUff::RdkitEnergyMinimizeUff() :
      RDKitEnergyUFF(),
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
    RdkitEnergyMinimizeUff::RdkitEnergyMinimizeUff
    (
      const RDKitEnergyUFF &ENERGY,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE,
      const std::string &POSITION_RESTRAINED_ATOMS_STRING,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    ) :
      RDKitEnergyUFF( ENERGY),
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
    RdkitEnergyMinimizeUff::RdkitEnergyMinimizeUff
    (
      const RDKitEnergyUFF &ENERGY,
      const size_t MAX_ITERATIONS,
      const double FORCE_TOLERANCE,
      const double ENERGY_TOLERANCE,
      const storage::Vector< size_t> &POSITION_RESTRAINED_ATOMS,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    ) :
      RDKitEnergyUFF( ENERGY),
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
    RdkitEnergyMinimizeUff *RdkitEnergyMinimizeUff::Clone() const
    {
      return new RdkitEnergyMinimizeUff( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RdkitEnergyMinimizeUff::GetAlias() const
    {
      static const std::string s_name( "EnergyMinimization_UFF");
      return s_name;
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RdkitEnergyMinimizeUff::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the maximum number of iterations to perform for geometry optimization
    //! @returns the number of iterations
    size_t RdkitEnergyMinimizeUff::GetMaxIterations() const
    {
      return m_MaxIterations;
    }

    //! @brief returns the force tolerance for the geometry optimization
    //! @returns the force tolerance
    double RdkitEnergyMinimizeUff::GetForceTolerance() const
    {
      return m_ForceTolerance;
    }

    //! @brief returns the energy tolerance for the geometry optimization
    //! @returns the energy tolerance
    double RdkitEnergyMinimizeUff::GetEnergyTolerance() const
    {
      return m_EnergyTolerance;
    }

    //! @brief returns the atoms (by index) to which positional restraints are applied
    storage::Vector< size_t> RdkitEnergyMinimizeUff::GetPositionalRestraintAtoms() const
    {
      return m_PositionalRestraintAtoms;
    }

    //! @brief returns the atoms (by string) to which positional restraints are applied
    std::string RdkitEnergyMinimizeUff::GetPositionalRestraintAtomsString() const
    {
      return m_PositionalRestraintAtomsString;
    }

    //! @brief returns the maximum displacement each atom can experience before restraint force activates
    storage::Vector< double> RdkitEnergyMinimizeUff::GetMaxUnrestrainedDisplacement() const
    {
      return m_MaxUnrestrainedDisplacement;
    }

    //! @brief returns the restraint force applied to each atom
    storage::Vector< double> RdkitEnergyMinimizeUff::GetRestraintForce() const
    {
      return m_RestraintForce;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set the maximum number of iterations to perform for geometry optimization
    void RdkitEnergyMinimizeUff::SetMaxIterations( const size_t &MAX_ITERATIONS)
    {
      m_MaxIterations = MAX_ITERATIONS;
    }

    //! @brief set the force tolerance for the geometry optimization
    void RdkitEnergyMinimizeUff::SetForceTolerance( const double FORCE_TOLERANCE)
    {
      m_ForceTolerance = FORCE_TOLERANCE;
    }

    //! @brief set the energy tolerance for the geometry optimization
    void RdkitEnergyMinimizeUff::SetEnergyTolerance( const double ENERGY_TOLERANCE)
    {
      m_EnergyTolerance = ENERGY_TOLERANCE;
    }

    //! @brief sets the atoms to which positional restraints are applied
    void RdkitEnergyMinimizeUff::SetPositionalRestraintAtoms( storage::Vector< size_t> &ATOMS)
    {
      m_PositionalRestraintAtoms = ATOMS;
    }

    //! @brief sets the atoms to which positional restraints are applied
    void RdkitEnergyMinimizeUff::SetPositionalRestraintAtomsFromString( std::string &ATOMS)
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
          "[WARNING] RdkitEnergyMinimizeUff::SetPositionalRestraintAtomsFromString "
          "string is empty - the restrained atoms set will not be modified!"
        );
      }
    }

    //! @brief sets the atoms to which positional restraints are applied
    void RdkitEnergyMinimizeUff::SetPositionalRestraintAtomsString( std::string &ATOMS)
    {
      m_PositionalRestraintAtomsString = ATOMS;
    }

    //! @brief sets the maximum displacement each atom can experience before restraint force activates
    void RdkitEnergyMinimizeUff::SetMaxUnrestrainedDisplacement( storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT)
    {
      m_MaxUnrestrainedDisplacement = MAX_UNRESTRAINED_DISPLACEMENT;
    }

    //! @brief sets the restraint force applied to each atom
    void RdkitEnergyMinimizeUff::SetRestraintForce( storage::Vector< double> &RESTRAINT_FORCE)
    {
      m_RestraintForce = RESTRAINT_FORCE;
    }

    //! @brief add positional constraints to force field for geometry optimization
    //! @param FORCE_FIELD the force field that is modified with the new restraint term
    //! @param ATOM_INDICES indices that are restrained during minimization
    //! @param MAX_UNRESTRAINED DISPLACEMENT coordinate displacement above which restraint force is applied
    //! @param RESTRAINT_FORCE restraint force
    void RdkitEnergyMinimizeUff::AddPositionalRestraints
    (
      ::ForceFields::ForceField *FORCE_FIELD, // raw pointer unconventional for BCL outside of Clone(), but this is what RDKit requires
      const storage::Vector< size_t> &ATOM_INDICES,
      const storage::Vector< double> &MAX_UNRESTRAINED_DISPLACEMENT,
      const storage::Vector< double> &RESTRAINT_FORCE
    ) const
    {
      // sanity check on vector sizes
      if( ATOM_INDICES.GetSize() == MAX_UNRESTRAINED_DISPLACEMENT.GetSize() && ATOM_INDICES.GetSize() == RESTRAINT_FORCE.GetSize())
      {
        // loop over atom indices and add restraint forces to our force field
        for( size_t i( 0), sz( ATOM_INDICES.GetSize()); i < sz; ++i)
        {
          ::ForceFields::UFF::PositionConstraintContrib *coord_cst;
          coord_cst = new ::ForceFields::UFF::PositionConstraintContrib
              (
                FORCE_FIELD, ATOM_INDICES( i),
                MAX_UNRESTRAINED_DISPLACEMENT( i),
                RESTRAINT_FORCE( i)
              );
          FORCE_FIELD->contribs().push_back( ForceFields::ContribPtr( coord_cst));
        }
      }
      else if( ATOM_INDICES.GetSize() == MAX_UNRESTRAINED_DISPLACEMENT.GetSize() && RESTRAINT_FORCE.GetSize() == size_t( 1))
      {
        for( size_t i( 0), sz( ATOM_INDICES.GetSize()); i < sz; ++i)
        {
          ::ForceFields::UFF::PositionConstraintContrib *coord_cst;
          coord_cst = new ::ForceFields::UFF::PositionConstraintContrib
              (
                FORCE_FIELD, ATOM_INDICES( i),
                MAX_UNRESTRAINED_DISPLACEMENT( i),
                RESTRAINT_FORCE( 0)
              );
          FORCE_FIELD->contribs().push_back( ForceFields::ContribPtr( coord_cst));
        }
      }
      // do not kill, but inform user that positional restraints are not added
      else
      {
        BCL_MessageStd
        (
          "[WARNING] RdkitEnergyMinimizeUff::AddPositionalRestraints "
          "The number of atoms does not match the number of max displacements and/or the number of provided restraint forces; "
          "alternatively, if the number of restraint forces to be added is one, then the number of atoms simply does not match the number of "
          "max displacements. NO POSITIONAL RESTRAINT ADDED!"
          " Number of atoms to be restrained: " + util::Format()( ATOM_INDICES.GetSize()) +
          " Number of max unrestrained distances: " + util::Format()( MAX_UNRESTRAINED_DISPLACEMENT.GetSize()) +
          " Number of restraint forces specified: " + util::Format()( RESTRAINT_FORCE.GetSize())
        )
      }
    }

    //! @brief optimizes the geometry of a molecule based on a molecular mechanics force field
    //! @param MOLECULE the molecule to be optimized; passed as non-const reference to be modified directly
    //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
    //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
    //! is the final energy of the optimized geometry
    storage::Pair< int, double> RdkitEnergyMinimizeUff::OptimizeGeometry( chemistry::FragmentComplete &MOLECULE) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
//      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, GetMMFFVariantAsString());
//      if( !mmff_mol_properties.isValid())
//      {
//        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
//        return storage::Pair< int, double>( -1, util::GetUndefinedDouble());
//      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( *rdkit_mol, m_NonbondedThreshold, -1, m_IgnoreInterFragmentInteractions);
      ff->initialize();

      // add constraints
      if( m_PositionalRestraintAtoms.GetSize())
      {
        AddPositionalRestraints( ff, m_PositionalRestraintAtoms, m_MaxUnrestrainedDisplacement, m_RestraintForce);
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
    storage::Triplet< chemistry::FragmentComplete, int, double> RdkitEnergyMinimizeUff::OptimizeGeometry
    (
      const chemistry::FragmentComplete &MOLECULE
    ) const
    {
      // convert to rdkit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol;
      rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOLECULE);

      // check validity
//      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, GetMMFFVariantAsString());
//      if( !mmff_mol_properties.isValid())
//      {
//        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
//        return storage::Triplet< chemistry::FragmentComplete, int, double>( MOLECULE, -1, util::GetUndefinedDouble());
//      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( *rdkit_mol, m_NonbondedThreshold, -1, m_IgnoreInterFragmentInteractions);
      ff->initialize();

      // add constraints
      if( m_PositionalRestraintAtoms.GetSize())
      {
        AddPositionalRestraints( ff, m_PositionalRestraintAtoms, m_MaxUnrestrainedDisplacement, m_RestraintForce);
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
    //! @param MMFF_VARIANT whether to use MMFF94 or MMFF94s
    //! @param NON_BONDED_THRESHOLD the threshold to be used in adding non-bonded terms to the force field.
    //! @param IGNORE_INTER_FRAG_INTERACTIONS If true, nonbonded terms will not be added between fragments
    //! @param MAX_ITERATIONS maximum number of iterations in which to achieve convergence
    //! @param FORCE_TOLERANCE convergence criterion for force
    //! @param ENERGY_TOLERANCE convergence criterion for energy
    //! @returns a pair where the first value indicates if the minimization was a success (0), failed to converge
    //! within the maximum number of iterations (1), or had missing parameters (-1), and where the second value
    //! is the final energy of the optimized geometry
    storage::Pair< int, double> RdkitEnergyMinimizeUff::OptimizeGeometry
    (
      chemistry::FragmentComplete &MOLECULE,
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
//      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, MMFF_VARIANT);
//      if( !mmff_mol_properties.isValid())
//      {
//        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
//        return storage::Pair< int, double>( -1, util::GetUndefinedDouble());
//      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( *rdkit_mol, NON_BONDED_THRESHOLD, -1, IGNORE_INTER_FRAG_INTERACTIONS);
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
    storage::Triplet< chemistry::FragmentComplete, int, double> RdkitEnergyMinimizeUff::OptimizeGeometry
    (
      const chemistry::FragmentComplete &MOLECULE,
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
//      ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, MMFF_VARIANT);
//      if( !mmff_mol_properties.isValid())
//      {
//        BCL_MessageStd( "Invalid MMFF molecule properties. Returning null.");
//        return storage::Triplet< chemistry::FragmentComplete, int, double>( MOLECULE, -1, util::GetUndefinedDouble());
//      }

      // generate an initialized force field ready for use
      ::ForceFields::ForceField *ff = ::RDKit::UFF::constructForceField( *rdkit_mol, NON_BONDED_THRESHOLD, -1, IGNORE_INTER_FRAG_INTERACTIONS);
      ff->initialize();

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
    io::Serializer RdkitEnergyMinimizeUff::GetSerializer() const
    {
      io::Serializer parameters( RDKitEnergyUFF::GetSerializer());
      parameters.SetClassDescription( "Optimizes the geometry of a molecule using the UFF force field with or without restraints.");
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
      bool RdkitEnergyMinimizeUff::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        // parse string encoding atoms to be restrained in coordinate space
        SetPositionalRestraintAtomsFromString( m_PositionalRestraintAtomsString);
        return true;
      }

  } // namespace mm
} // namespace bcl
