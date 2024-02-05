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

#ifndef BCL_APP_MOLECULE_MINIMIZE_H_
#define BCL_APP_MOLECULE_MINIMIZE_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"

// includes from rdkit
//#include "GraphMol/RWMol.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeMinimize
    //! @brief Application for processing performing energy minimization of molecular geometries.
    //! @details This application makes use of the RDKit external library. It converts a BCL molecule into
    //! an RDKit molecule, utilizes the RDkit ForceFields framework to perform the energy minimization, and
    //! then converts the resultant molecule back to a BCL molecule before output.
    //!
    //! @see @link example_app_molecule_minimize.cpp @endlink
    //! @author brownbp1
    //! @date Aug 27, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeMinimize :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! Output filename base
      util::ShPtr< command::FlagInterface> m_OutputFlag;

      //! Force field to use
      util::ShPtr< command::FlagInterface> m_ForceFieldFlag;

      //! The threshold to be used in adding non-bonded terms to the force field.
      //! Any non-bonded contact whose current distance is greater than nonBondedThresh * the minimum value for that contact
      //! will not be included.
      util::ShPtr< command::FlagInterface> m_NonbondedThresholdFlag;

      //!  If true, nonbonded terms will not be added between fragments
      util::ShPtr< command::FlagInterface> m_IgnoreInterFragmentInteractionsFlag;

      //! @brief Maximum number of iterations to perform for geometry optimization
      util::ShPtr< command::FlagInterface> m_MaxIterationsFlag;

      //! @brief The convergence criterion for forces
      util::ShPtr< command::FlagInterface> m_ForceToleranceFlag;

      //! @brief The convergence criterion for energies
      util::ShPtr< command::FlagInterface> m_EnergyToleranceFlag;

      //! Position restraint flag
      util::ShPtr< command::FlagInterface> m_PositionalRestraintsFlag;

      //! MDL property that specifies which atoms are to be restrained
      util::ShPtr< command::FlagInterface> m_PositionalRestraintsMDLFlag;

      //! @brief MDL property that specifies the maximum allowed displacement per-atom (per-atom vector)
      util::ShPtr< command::FlagInterface> m_MaxUnrestrainedDisplacementMDLFlag;

      //! @brief The maximum allowed displacement per-atom (scalar applied to all atoms)
      util::ShPtr< command::FlagInterface> m_MaxUnrestrainedDisplacementDefaultFlag;

      //! @brief MDL property that specifies the per-atom restraint force (per-atom vector)
      util::ShPtr< command::FlagInterface> m_RestraintForceMDLFlag;

      //! @brief The per-atom restraint force (scalar applied to all atoms)
      util::ShPtr< command::FlagInterface> m_RestraintForceDefaultFlag;

      //! Output stream
      mutable io::OFStream m_Output;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeMinimize();

      //! copy constructor, only copy the flags
      MoleculeMinimize( const MoleculeMinimize &PARENT);

    public:

      // instantiate enumerator for MoleculeMinimize class
      static const ApplicationType MoleculeMinimize_Instance;

      //! @brief Clone function
      //! @return pointer to new MoleculeMinimize
      MoleculeMinimize *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

    ////////////////
    //    main    //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // helper functions //
    /////////////////////

      //! @brief prepare output streams
      void InitializeOutputFiles() const;

      //! @brief identify atoms to which positional restraints need be applied
      //! @param MOLECULE the molecule being geometry optimized
      //! @return the atom indices in MOLECULE requiring restraints
      const storage::Vector< size_t> GetPositionalRestraintAtoms( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief get the displacement allowed by each atom before the restraint force is applied
      //! @param MOLECULE the molecule being geometry optimized
      //! @return the per-atom allowed unrestrained displacement
      const storage::Vector< double> GetMaxUnrestrainedDisplacement( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief get the restraint force felt by each atom
      //! @param MOLECULE the molecule being geometry optimized
      //! @return the per-atom restraint force
      const storage::Vector< double> GetRestraintForce( const chemistry::FragmentComplete &MOLECULE) const;


    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    }; // MoleculeMinimize

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_MINIMIZE_H_
