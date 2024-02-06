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
#include "bcl_app_molecule_minimize.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_directory_entry.h"
#include "math/bcl_math_limits.h"
#include "mm/bcl_mm_rdkit_energy_minimize.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    const ApplicationType MoleculeMinimize::MoleculeMinimize_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeMinimize(), GetAppGroups().e_Molecule)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    MoleculeMinimize::MoleculeMinimize() :
      m_OutputFlag
      (
        new command::FlagStatic
        (
          "output",
          "output file for the minimized molecules",
          command::Parameter( "output", "name for output file containing geometry optimized molecules")
        )
      ),
      m_ForceFieldFlag
      (
        new command::FlagStatic
        (
          "force_field",
          "force field to use for energy minimization",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "UFF", "MMFF94", "MMFF94s") ),
            "UFF"
          )
        )
      ),
      m_NonbondedThresholdFlag
      (
        new command::FlagStatic
        (
          "nonbonded_threshold",
          "The threshold to be used in adding non-bonded terms to the force field; "
          "any non-bonded contact whose current distance is greater than nonBondedThresh * the minimum value "
          "for that contact will not be included",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckRanged< double>(),
            "10.0"
          )
        )
      ),
      m_IgnoreInterFragmentInteractionsFlag
      (
        new command::FlagStatic
        (
          "ignore_inter_fragment_interactions",
          "if enabled, nonbonded terms will not be added between fragments"
        )
      ),
      m_MaxIterationsFlag
      (
        new command::FlagStatic
        (
          "max_iterations",
          "the maximum number of iterations to perform during geometry optimization; "
          "if 0 then compute the energy without additional optimization",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            "1000"
          )
        )
      ),
      m_ForceToleranceFlag
      (
        new command::FlagStatic
        (
          "force_tolerance",
          "the convergence criterion for forces (kcal/(mol·A))",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckRanged< double>( 0.0, 1.0),
            "1.0e-4"
          )
        )
      ),
      m_EnergyToleranceFlag
      (
        new command::FlagStatic
        (
          "energy_tolerance",
          "the convergence criterion for energies (kcal/mol)",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckRanged< double>( 0.0, 1.0),
            "1.0e-4"
          )
        )
      ),
      m_PositionalRestraintsFlag
      (
        new command::FlagStatic
        (
          "position_restraints",
          "apply atomic coordinate space restraints; by default, enabling this flag will add restraints to all heavy atoms",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "minimum",
                "the lowest allowable similarity that a molecule is allowed to have with a scaffold in order to be a candidate for conformer generation",
                command::ParameterCheckRanged< float>( 0.0, 1.0),
                "0.0"
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "maximum",
                "the highest allowable similarity that a molecule is allowed to have with a scaffold in order to be a candidate for conformer generation",
                command::ParameterCheckRanged< float>(0.0, 1.0),
                "1.0"
              )
            )
          )
        )
      ),
      m_PositionalRestraintsMDLFlag
      (
        new command::FlagStatic
        (
          "position_restraints_mdl",
          "the MDL property name corresponding to the atom indices that are to be restrained; "
          "if no MDL property is specified but 'position_restraints' is enabled, then by default "
          "all heavy atoms will be restrained to their starting coordinates",
          command::Parameter
          (
            "",
            "",
            "SampleByParts"
          )
        )
      ),
      m_PositionalRestraintsMDLComplementFlag
      (
        new command::FlagStatic
        (
          "position_restraints_mdl_complement",
          "if enabled, use the complement atom indices of those specified by the MDL property specified by 'position_restraints_mdl'"
        )
      ),
      m_MaxUnrestrainedDisplacementMDLFlag
      (
        new command::FlagStatic
        (
          "max_unrestrained_displacement_mdl",
          "the MDL property that specifies the maximum allowed unrestrained displacement per-atom; "
          "if enabled, must be provided for all atoms in the molecule or else the default value "
          "will be applied to all atoms",
          command::Parameter
          (
            "",
            "",
            ""
          )
        )
      ),
      m_MaxUnrestrainedDisplacementDefaultFlag
      (
        new command::FlagStatic
        (
          "max_unrestrained_displacement",
          "scalar value specifying the maximum allowed unrestrained displacement in Angstroms for each atom during minimization; "
          "functions as the default value if 'max_unrestraiend_displacement_mdl' is not provided or is invalid",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckRanged< double>( 0.0, std::numeric_limits< double>::max()),
            "0.0"
          )
        )
      ),
      m_RestraintForceMDLFlag
      (
        new command::FlagStatic
        (
          "restraint_force_mdl",
          "the MDL property that specifies the per-atom restraint force "
          "if enabled, must be provided for all atoms in the molecule or else the default value "
          "will be applied to all atoms",
          command::Parameter
          (
            "",
            "",
            ""
          )
        )
      ),
      m_RestraintForceDefaultFlag
      (
        new command::FlagStatic
        (
          "restaint_force",
          "scalar value specifying the restraint force in kcal/(mol·A) for each atom during minimization; "
          "functions as the default value if 'restraint_force_mdl' is not provided or is invalid",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckRanged< double>( 0.0, std::numeric_limits< double>::max()),
            "10.0"
          )
        )
      )
    {
    }

    //! copy constructor, only copy the flags
    MoleculeMinimize::MoleculeMinimize( const MoleculeMinimize &PARENT) :
        m_OutputFlag( PARENT.m_OutputFlag),
        m_ForceFieldFlag( PARENT.m_ForceFieldFlag),
        m_NonbondedThresholdFlag( PARENT.m_NonbondedThresholdFlag),
        m_IgnoreInterFragmentInteractionsFlag( PARENT.m_IgnoreInterFragmentInteractionsFlag),
        m_MaxIterationsFlag( PARENT.m_MaxIterationsFlag),
        m_ForceToleranceFlag( PARENT.m_ForceToleranceFlag),
        m_EnergyToleranceFlag( PARENT.m_EnergyToleranceFlag),
        m_PositionalRestraintsFlag( PARENT.m_PositionalRestraintsFlag),
        m_PositionalRestraintsMDLFlag( PARENT.m_PositionalRestraintsMDLFlag),
        m_PositionalRestraintsMDLComplementFlag( PARENT.m_PositionalRestraintsMDLComplementFlag),
        m_MaxUnrestrainedDisplacementMDLFlag( PARENT.m_MaxUnrestrainedDisplacementMDLFlag),
        m_MaxUnrestrainedDisplacementDefaultFlag( PARENT.m_MaxUnrestrainedDisplacementDefaultFlag),
        m_RestraintForceMDLFlag( PARENT.m_RestraintForceMDLFlag),
        m_RestraintForceDefaultFlag( PARENT.m_RestraintForceDefaultFlag)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeMinimize
    MoleculeMinimize *MoleculeMinimize::Clone() const
    {
      return new MoleculeMinimize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeMinimize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeMinimize::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      sp_cmd->AddFlag( m_OutputFlag);
      sp_cmd->AddFlag( m_ForceFieldFlag);
      sp_cmd->AddFlag( m_NonbondedThresholdFlag);
      sp_cmd->AddFlag( m_IgnoreInterFragmentInteractionsFlag);
      sp_cmd->AddFlag( m_MaxIterationsFlag);
      sp_cmd->AddFlag( m_ForceToleranceFlag);
      sp_cmd->AddFlag( m_EnergyToleranceFlag);
      sp_cmd->AddFlag( m_PositionalRestraintsFlag);
      sp_cmd->AddFlag( m_PositionalRestraintsMDLFlag);
      sp_cmd->AddFlag( m_PositionalRestraintsMDLComplementFlag);
      sp_cmd->AddFlag( m_MaxUnrestrainedDisplacementMDLFlag);
      sp_cmd->AddFlag( m_MaxUnrestrainedDisplacementDefaultFlag);
      sp_cmd->AddFlag( m_RestraintForceMDLFlag);
      sp_cmd->AddFlag( m_RestraintForceDefaultFlag);

      sp_cmd->AddFlag( sdf::GetAddHydrogensFlag());
      sp_cmd->AddFlag( sdf::GetNeutralizeChargesFlag());
      sp_cmd->AddFlag( sdf::GetExplicitAromaticityFlag());
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd, storage::Set< command::FlagTypeEnum>( command::e_AppGeneric));
      return sp_cmd;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeMinimize::GetDescription() const
    {
      return "Performs molecular mechanics energy-based geometry optimization on molecules";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeMinimize::GetReadMe() const
    {
      static std::string s_read_me =
        "MoleculeMinimize performs energy minimization of molecular geometries.";
      return s_read_me;
    }

  ////////////////
  //    main    //
  ////////////////

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeMinimize::Main() const
    {
      // initialization
      InitializeOutputFiles();

      size_t feed_index( 0);
      for( chemistry::FragmentFeed feed; feed.NotAtEnd(); ++feed, ++feed_index)
      {
        chemistry::FragmentComplete mol( *feed);

        // add any restraints
        const storage::Vector< size_t> restraint_atoms( GetPositionalRestraintAtoms( mol));
        const storage::Vector< double> max_unrestrained_displacement( GetMaxUnrestrainedDisplacement( mol, restraint_atoms));
        const storage::Vector< double> restraint_force( GetRestraintForce( mol, restraint_atoms));

        // minimize
        if( m_MaxIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>() )
        {
          storage::Pair< int, double> minimization_result
          (
            mm::RdkitEnergyMinimize::OptimizeGeometry
            (
              mol,
              m_ForceFieldFlag->GetFirstParameter()->GetValue(),
              m_NonbondedThresholdFlag->GetFirstParameter()->GetNumericalValue< double>(),
              m_IgnoreInterFragmentInteractionsFlag->GetFlag(),
              m_MaxIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
              m_ForceToleranceFlag->GetFirstParameter()->GetNumericalValue< double>(),
              m_EnergyToleranceFlag->GetFirstParameter()->GetNumericalValue< double>(),
              restraint_atoms,
              max_unrestrained_displacement,
              restraint_force
            )
          );

          // handle properties relevant to minimization
          mol.StoreProperties( *feed);
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_Mode", util::Format()( "Minimize"));
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_FF", m_ForceFieldFlag->GetFirstParameter()->GetValue());
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_Success", minimization_result.First() == int( 0) ? util::Format()( "1") : util::Format()( "0"));
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_FailedConvergence", minimization_result.First() == int( 1) ? util::Format()( "1") : util::Format()( "0"));
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_MissingParams", minimization_result.First() == int( -1) ? util::Format()( "1") : util::Format()( "0"));
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_Energy", linal::Vector< float>( 1, minimization_result.Second() ));
        }
        // single point energy
        else
        {
          double energy_result
          (
            mm::RDKitEnergy::CalculateEnergy
            (
              mol,
              m_ForceFieldFlag->GetFirstParameter()->GetValue(),
              m_NonbondedThresholdFlag->GetFirstParameter()->GetNumericalValue< double>(),
              m_IgnoreInterFragmentInteractionsFlag->GetFlag()
            )
          );

          // handle properties relevant to single point energy calculation
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_Mode", util::Format()( "Energy"));
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_FF", m_ForceFieldFlag->GetFirstParameter()->GetValue());
          mol.GetStoredPropertiesNonConst().SetMDLProperty("MoleculeMinimize_Energy", linal::Vector< float>( 1, energy_result ));
        }
        mol.WriteMDL( m_Output);
      }

      // close and end
      io::File::CloseClearFStream( m_Output);
      return 0;
    }

  //////////////////////
  // helper functions //
  /////////////////////

    //! @brief prepare output streams
    void MoleculeMinimize::InitializeOutputFiles() const
    {
      // clear output file if it exists
      io::DirectoryEntry entry( m_OutputFlag->GetFirstParameter()->GetValue());
      if( entry.DoesExist())
      {
        entry.Remove();
      }

      if( m_OutputFlag->GetFlag())
      {
        BCL_MessageStd("Energy minimized molecules will be written to "
          + util::Format()( m_OutputFlag->GetFirstParameter()->GetValue()));
        io::File::MustOpenOFStream( m_Output, m_OutputFlag->GetFirstParameter()->GetValue());
      }
    }

    //! @brief identify atoms to which positional restraints need be applied
    //! @param MOLECULE the molecule being geometry optimized
    //! @return the atom indices in MOLECULE requiring restraints
    const storage::Vector< size_t> MoleculeMinimize::GetPositionalRestraintAtoms
    (
      const chemistry::FragmentComplete &MOLECULE
    ) const
    {
      // no restraints
      if( !m_PositionalRestraintsFlag->GetFlag())
      {
        return storage::Vector< size_t>();
      }

      else if( m_PositionalRestraintsMDLFlag->GetFlag() && !MOLECULE.GetStoredProperties().GetMDLProperty( m_PositionalRestraintsMDLFlag->GetFirstParameter()->GetValue() ).empty())
      {
        std::string restraint_atoms_str( MOLECULE.GetStoredProperties().GetMDLProperty( m_PositionalRestraintsMDLFlag->GetFirstParameter()->GetValue() ));
        if( m_PositionalRestraintsMDLComplementFlag->GetFlag())
        {
          storage::Vector< size_t> restraint_atoms( util::SplitStringToNumerical< size_t>( restraint_atoms_str) );
          storage::Set< size_t> restraint_atoms_set( restraint_atoms.Begin(), restraint_atoms.End());
          restraint_atoms.Reset();
          for( size_t mol_atom( 0), mol_sz( MOLECULE.GetSize()); mol_atom <  mol_sz; ++mol_atom)
          {
            if( restraint_atoms_set.Find(mol_atom) == restraint_atoms_set.End())
            {
              restraint_atoms.PushBack(mol_atom);
            }
          }
          return restraint_atoms;
        }
        return storage::Vector< size_t>( util::SplitStringToNumerical< size_t>( restraint_atoms_str) );
      }

      // default all heavy atom restraints
      else
      {
        storage::Vector< size_t> restraint_atoms;
        for( size_t atom_i( 0), sz( MOLECULE.GetSize()); atom_i < sz; ++atom_i)
        {
          const chemistry::AtomVector< chemistry::AtomComplete> &atoms( MOLECULE.GetAtomVector());
          if( atoms( atom_i).GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
          {
            restraint_atoms.PushBack( atom_i);
          }
        }
        return restraint_atoms;
      }
      return storage::Vector< size_t>();
    }

    //! @brief get the displacement allowed by each atom before the restraint force is applied
    //! @param MOLECULE the molecule being geometry optimized
    //! @return the per-atom allowed unrestrained displacement
    const storage::Vector< double> MoleculeMinimize::GetMaxUnrestrainedDisplacement
    (
      const chemistry::FragmentComplete &MOLECULE,
      const storage::Vector< size_t> &RESTRAINT_ATOMS
    ) const
    {
      // no restraints
      if( !m_PositionalRestraintsFlag->GetFlag())
      {
        return storage::Vector< double>();
      }
      // no valid MDL option for max unrestrained displacement
      else if
      (
          !m_MaxUnrestrainedDisplacementMDLFlag->GetFlag() ||
          MOLECULE.GetStoredProperties().GetMDLProperty( m_MaxUnrestrainedDisplacementMDLFlag->GetFirstParameter()->GetValue() ).empty() ||
          MOLECULE.GetStoredProperties().GetMDLPropertyAsVector( m_MaxUnrestrainedDisplacementMDLFlag->GetFirstParameter()->GetValue() ).GetSize() != RESTRAINT_ATOMS.GetSize()
      )
      {
        // use the default value
        return storage::Vector< double>( RESTRAINT_ATOMS.GetSize(), m_MaxUnrestrainedDisplacementDefaultFlag->GetFirstParameter()->GetNumericalValue< double>() );
      }
      // restraints enabled, valid MDL
      else
      {
        const std::string max_unrestrained_displacement_str
        (
          MOLECULE.GetStoredProperties().GetMDLProperty( m_MaxUnrestrainedDisplacementMDLFlag->GetFirstParameter()->GetValue() )
        );
        return storage::Vector< double>( util::SplitStringToNumerical< double>( max_unrestrained_displacement_str) );
      }
      return storage::Vector< double>();
    }

    //! @brief get the restraint force felt by each atom
    //! @param MOLECULE the molecule being geometry optimized
    //! @return the per-atom restraint force
    const storage::Vector< double> MoleculeMinimize::GetRestraintForce
    (
      const chemistry::FragmentComplete &MOLECULE,
      const storage::Vector< size_t> &RESTRAINT_ATOMS
    ) const
    {
      // no restraints
      if( !m_PositionalRestraintsFlag->GetFlag())
      {
        return storage::Vector< double>();
      }
      // no valid MDL option for restraint force
      else if
      (
          !m_RestraintForceMDLFlag->GetFlag() ||
          MOLECULE.GetStoredProperties().GetMDLProperty( m_RestraintForceMDLFlag->GetFirstParameter()->GetValue() ).empty() ||
          MOLECULE.GetStoredProperties().GetMDLPropertyAsVector( m_RestraintForceMDLFlag->GetFirstParameter()->GetValue() ).GetSize() != RESTRAINT_ATOMS.GetSize()
      )
      {
        // use the default value
        return storage::Vector< double>( RESTRAINT_ATOMS.GetSize(), m_RestraintForceDefaultFlag->GetFirstParameter()->GetNumericalValue< double>() );
      }
      // restraints enabled, valid MDL
      else
      {
        const std::string restraint_force_str
        (
          MOLECULE.GetStoredProperties().GetMDLProperty( m_RestraintForceMDLFlag->GetFirstParameter()->GetValue() )
        );
        return storage::Vector< double>( util::SplitStringToNumerical< double>( restraint_force_str) );
      }
      return storage::Vector< double>();
    }

  } // namespace app
} // namespace bcl

