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
#include "io/bcl_io_directory_entry.h"
#include "math/bcl_math_limits.h"

// external includes - sorted alphabetically
#include "GraphMol/ForceFieldHelpers/MMFF/MMFF.h"
#include "GraphMol/ForceFieldHelpers/FFConvenience.h"
#include "ForceField/ForceField.h"
#include "ForceField/MMFF/PositionConstraint.h"
#include "ForceField/MMFF/DistanceConstraint.h"
#include "ForceField/MMFF/AngleConstraint.h"
#include "ForceField/MMFF/TorsionConstraint.h"

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
      m_OutputFilenameBase
      (
        new command::FlagStatic
        (
          "output",
          "base name for output",
          command::Parameter( "output", "base name for output")
        )
      )
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

    ////////////////////
    // common options //
    ////////////////////

      // add flags for input
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // Output filename base
      sp_cmd->AddFlag( m_OutputFilenameBase);

      // molecule reading prefs
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags
      (
        *sp_cmd,
        storage::Set< command::FlagTypeEnum>( command::e_Pthread)
      );

      // return assembled Command object
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

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeMinimize::Main() const
    {
      // clear output file if it exists
      const std::string output_filename( m_OutputFilenameBase->GetFirstParameter()->GetValue() + util::Format()( ".sdf"));
      io::DirectoryEntry entry( output_filename);
      if( entry.DoesExist())
      {
        entry.Remove();
      }

      // I/O
      io::OFStream output;
      io::File::MustOpenOFStream( output, output_filename, std::ios::app);
      size_t feed_index( 0);
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol, rdkit_mol_min, rdkit_mol_min_ctrl, rdkit_mol_min_cst = nullptr;
      std::shared_ptr< chemistry::FragmentComplete> mol, min_mol_ctrl, min_mol, min_mol_cst = nullptr;
      for( chemistry::FragmentFeed feed; feed.NotAtEnd(); ++feed, ++feed_index)
      {
        // Convert the BCL fragment into an RDKit molecule
        rdkit_mol = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( *feed);
        rdkit_mol_min_ctrl = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( *feed);
        rdkit_mol_min = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( *feed);
        rdkit_mol_min_cst = chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( *feed);
        if( rdkit_mol == nullptr)
        {
          continue;
        }

        // minimize
        Minimize( rdkit_mol_min);
        MinimizeWithConstraints( rdkit_mol_min_cst);

        // required for force field construction
        std::string const mmff_variant("MMFF94s");
        size_t max_iters( 1000);
        ::RDKit::MMFF::MMFFMolProperties mmffMolProperties( *rdkit_mol_min_ctrl, mmff_variant);

  //      // output from minimizer
  //      // first: -1 if parameters were missing, 0 if the optimization converged, 1 if more iterations are required.
  //      // second: the energy
  //      std::pair<int, double> res = std::make_pair(-1, 0.0);
  //      if ( mmffMolProperties.isValid())
  //      {
  //        // construct force field
  //        ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( MOL, 10, -1, true);
  ////        ff->initialize();
  //
  //        // run minimizationMOL
  //        res = ::RDKit::ForceFieldsHelper::OptimizeMolecule(*ff, max_iters);
  //        delete ff;
  //      }
        ::RDKit::MMFF::MMFFOptimizeMolecule( *rdkit_mol_min_ctrl, max_iters, mmff_variant, 10, -1, true);


        // write complete
        mol = chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete(*rdkit_mol);
        min_mol_ctrl = chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete(*rdkit_mol_min_ctrl);
        min_mol = chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete(*rdkit_mol_min);
        min_mol_cst = chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete(*rdkit_mol_min_cst);
//        if( min_mol == nullptr)
//        {
//          BCL_MessageStd("Molecule could not be minimized! Returning input...");
//          min_mol = std::make_shared< chemistry::FragmentComplete>( *feed);
//        }

        mol->WriteMDL( output);
        min_mol_ctrl->WriteMDL( output);
        min_mol->WriteMDL( output);
        min_mol_cst->WriteMDL( output);
      }
      io::File::CloseClearFStream( output);
      return 0;
    }


    //! @brief minimize without any constraints
    void MoleculeMinimize::Minimize( std::shared_ptr< ::RDKit::RWMol> &MOL) const
    {
      // required for force field construction
      std::string const mmff_variant("MMFF94s");
      size_t max_iters( 1000);
      ::RDKit::MMFF::MMFFMolProperties mmffMolProperties( *MOL, mmff_variant);

      // output from minimizer
      // first: -1 if parameters were missing, 0 if the optimization converged, 1 if more iterations are required.
      // second: the energy
      std::pair<int, double> res = std::make_pair(-1, 0.0);
      if ( mmffMolProperties.isValid())
      {
        // construct force field
        ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *MOL, 10, -1, true);
        ff->initialize();

        // run minimizationMOL
        res = ::RDKit::ForceFieldsHelper::OptimizeMolecule(*ff, max_iters);
        ff->minimize( max_iters);
        double energy( ff->calcEnergy());
        BCL_MessageStd( "Minimization final energy: " + util::Format()( energy));
        delete ff;
      }
    }

    //! @brief minimize with constraints
    void MoleculeMinimize::MinimizeWithConstraints( std::shared_ptr< ::RDKit::RWMol> &MOL) const
    {
      // required for force field construction
      std::string const mmff_variant("MMFF94s");
      size_t max_iters( 1000);
      ::RDKit::MMFF::MMFFMolProperties mmffMolProperties( *MOL, mmff_variant);

      // output from minimizer
      // first: -1 if parameters were missing, 0 if the optimization converged, 1 if more iterations are required.
      // second: the energy
      std::pair<int, double> res = std::make_pair(-1, 0.0);
      if ( mmffMolProperties.isValid())
      {
        // construct force field
        ::ForceFields::ForceField *ff = ::RDKit::MMFF::constructForceField( *MOL, 10, -1, true);
        ff->initialize();

        // add constraint force field terms
        for( size_t i( 0); i < 11; ++i)
        {
          ::ForceFields::MMFF::PositionConstraintContrib *coord_cst;
          coord_cst = new ::ForceFields::MMFF::PositionConstraintContrib( ff, i, 0.0, 1.0e5);
          ff->contribs().push_back( ForceFields::ContribPtr( coord_cst));
        }

        // run minimization
//        res = ::RDKit::ForceFieldsHelper::OptimizeMolecule(*ff, max_iters);
        ff->minimize( max_iters);
        double energy( ff->calcEnergy());
        BCL_MessageStd( "Minimization with constraints final energy: " + util::Format()( energy));
        delete ff;
      }
    }

  } // namespace app
} // namespace bcl

