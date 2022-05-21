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
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_statistics.h"
#include "molecule/bcl_app_molecule_reorder.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    MoleculeReorder::MoleculeReorder() :
        m_InputFilenamesFlag
        (
          new command::FlagDynamic
          (
            "input_filenames",
            "filenames for sdf input files",
            command::Parameter
            (
              "input_filename",
              "an sdf input file",
              command::ParameterCheckFileExistence()
            ),
            1
          )
        ),
        m_RandomizeFlag
        (
          new command::FlagStatic
          (
            "randomize",
            "randomly reorder the molecules in the ensemble"
          )
        ),
        m_SortByPropertyFlag
        (
          new command::FlagStatic
          (
            "sort",
            "sort by any property",
            command::Parameter
            (
              "property",
              "small molecule property",
              command::ParameterCheckSerializable( descriptor::CheminfoProperty()),
              ""
            )
          )
        ),
        m_ReverseFlag
        (
          new command::FlagStatic
          (
            "reverse",
            "reverse the molecules (after sorting, if a sort property was given)"
          )
        ),
        m_MaxMoleculesPerSDFFlag
        (
          new command::FlagStatic
          (
            "max_molecules",
            "split incoming molecules into files of this size",
            command::Parameter
            (
              "filename",
              "sdf filename for where to write out molecules",
              command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
              util::Format()( std::numeric_limits< size_t>::max())
            )
          )
        ),
        m_Canonicalize
        (
          new command::FlagStatic
          (
            "canonicalize",
            "Set this flag to sort atoms in the molecule by descending cahn-ingold-prelog priority"
            "Useful for applications that require that the atoms in the molecule are aligned"
          )
        ),
        m_Atoms
        (
          new command::FlagDynamic
          (
            "atom_order",
            "Reorder atoms in the given molecules into the order specified (0-indexed"
            "e.g. if the order is 3 2, then the output molecules will contain the 4th atom (first), then the 3rd atom "
            "(second), with all other atoms omitted",
            command::Parameter
            (
              "order",
              "atom to use for this position",
              command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max())
            )
          )
        ),
        m_OutputFilenameFlag
        (
          new command::FlagStatic
          (
            "output",
            "sdf filename for where to write out molecules",
            command::Parameter
            (
              "filename",
              "sdf filename for where to write out molecules"
            )
          )
        ),
        m_OutputMaxFlag
        (
          new command::FlagDynamic
          (
            "output_max",
            "the maximum number of molecules to write out",
            command::Parameter
            (
              "max",
              "the maximum number of molecules to write to a file",
              command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
              util::Format()( std::numeric_limits< size_t>::max())
            )
          )
        ),
        m_OutputOrderFlag
        (
          new command::FlagDynamic
          (
            "output_order",
            "fine control over ordering of molecules using an explicit index list",
            command::Parameter
            (
              "ordering",
              "a 0-indexed list of how molecules should be reordered.",
              command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
              util::Format()( std::numeric_limits< size_t>::max())
            )
          )
        )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeReorder *MoleculeReorder::Clone() const
    {
      return new MoleculeReorder( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeReorder::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeReorder::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be edited
      sp_cmd->AddFlag( m_InputFilenamesFlag);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // randomly permute the molecules in the ensemble
      sp_cmd->AddFlag( m_RandomizeFlag);

      // property to sort by
      sp_cmd->AddFlag( m_SortByPropertyFlag);

      // flag for reversing the molecules in the ensemble, optionally after reversing
      sp_cmd->AddFlag( m_ReverseFlag);

      // The maximum number of molecules to store per sdf
      sp_cmd->AddFlag( m_MaxMoleculesPerSDFFlag);

      // whether to canonicalize atom ordering within molecules
      sp_cmd->AddFlag( m_Canonicalize);

      // explicit reordering of atoms within molecules
      sp_cmd->AddFlag( m_Atoms);

      //! File to write molecules into
      sp_cmd->AddFlag( m_OutputFilenameFlag);

      //! The maximum number of molecules to write out
      sp_cmd->AddFlag( m_OutputMaxFlag);

      //! explicit reordering of compounds using a list of indices
      sp_cmd->AddFlag( m_OutputOrderFlag);

    ///////////////////
    // default flags //
    ///////////////////

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief check that all the parameter choices were valid
    //! @return true if all the parameter choices were valid
    bool MoleculeReorder::CheckParametersAreAcceptable() const
    {
      // check that a reordering flag was given
      if( ( ( m_SortByPropertyFlag->GetFlag() && m_OutputOrderFlag->GetFlag()) || m_ReverseFlag->GetFlag()) && m_RandomizeFlag->GetFlag())
      {
        BCL_MessageCrt
        (
          "Choose only one reordering method for molecules: either sort by a property/reverse, randomize, or give an explicit ordering"
        );
        return false;
      }
      else if( m_Canonicalize->GetFlag() && m_Atoms->GetFlag())
      {
        BCL_MessageCrt
        (
          "Choose only one reordering method for atoms: either canonicalize or give an explicit atom order"
        );
        return false;
      }
      else if
      (
        !m_SortByPropertyFlag->GetFlag()
        && !m_OutputOrderFlag->GetFlag()
        && !m_ReverseFlag->GetFlag()
        && !m_RandomizeFlag->GetFlag()
        && !m_MaxMoleculesPerSDFFlag->GetFlag()
        && !m_Canonicalize->GetFlag()
        && !m_Atoms->GetFlag()
      )
      {
        // check that some method of reordering was requested
        BCL_MessageCrt
        (
          "Choose a reordering method: sort by a property and/or reverse, randomize, give an explicit ordering, or split up the ensemble"
        );
        return false;
      }

      return true;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeReorder::Main() const
    {
      if( !CheckParametersAreAcceptable())
      {
        return 0;
      }

      chemistry::FragmentEnsemble ensemble;

      // load input files; eventually, have options to load from db
      if( m_InputFilenamesFlag->GetFlag())
      {
        io::IFStream input;

        util::Stopwatch load_timer( false);

        const storage::Vector< std::string> files( m_InputFilenamesFlag->GetStringList());

        // load in the ensembles, one by one
        for
        (
          storage::Vector< std::string>::const_iterator itr_files( files.Begin()), itr_files_end( files.End());
          itr_files != itr_files_end;
          ++itr_files
        )
        {
          io::File::MustOpenIFStream( input, *itr_files);
          ensemble.ReadMoreFromMdl( input, sdf::GetCommandLineHydrogensPref());
        }
        // close the file stream
        io::File::CloseClearFStream( input);
        BCL_MessageStd
        (
          "Loaded ensembles with " + util::Format()( ensemble.GetSize())
          + " molecules total in " + load_timer.GetProcessDuration().GetTimeAsHourMinuteSecond()
        );
      }

      if( m_RandomizeFlag->GetFlag())
      {
        util::Stopwatch timer( false);
        ensemble.Shuffle();
        BCL_MessageStd
        (
          "Randomized the molecules in " + timer.GetProcessDuration().GetTimeAsHourMinuteSecond()
        );
      }
      else if( m_SortByPropertyFlag->GetFlag())
      {
        SortByNumericProperty
        (
          ensemble,
          descriptor::CheminfoProperty( m_SortByPropertyFlag->GetFirstParameter()->GetValue())
        );
      }
      else if( m_OutputOrderFlag->GetFlag())
      {
        storage::Vector< std::string> opts( m_OutputOrderFlag->GetStringList());
        storage::Vector< size_t> order( opts.GetSize());
        for( size_t i( 0), last( opts.GetSize()); i < last; ++i)
        {
          BCL_Assert
          ( 
            util::TryConvertFromString( order( i), opts( i), util::GetLogger()), 
            "Could not parse \"" + opts( i) + "\" as a number"
          );
        }
        storage::Vector< chemistry::FragmentComplete> mol_vector( ensemble.Begin(), ensemble.End());
        mol_vector.Reorder( order);
        ensemble = chemistry::FragmentEnsemble();
        for( size_t m( 0), end( mol_vector.GetSize()); m < end; ++m)
        {
          ensemble.PushBack( mol_vector( m));
        }
        BCL_MessageStd( "Reordered " + util::Format()( mol_vector.GetSize()) + " molecules");
      }

      if( m_ReverseFlag->GetFlag())
      {
        ensemble.GetMolecules().Reverse();
      }

      if( m_Atoms->GetFlag())
      {
        ReorderAtoms( ensemble);
      }
      else if( m_Canonicalize->GetFlag())
      {
        CanonicalizeAtoms( ensemble);
      }

      if( m_OutputFilenameFlag->GetFlag())
      {
        WriteEnsemble( ensemble);
      }

      // end
      return 0;
    }

    //! @brief reorder the atoms in the ensemble
    //! @param ENSEMBLE ensemble in which to reorder the atoms of the molecule
    void MoleculeReorder::ReorderAtoms( chemistry::FragmentEnsemble &ENSEMBLE) const
    {
      const storage::Vector< size_t> new_order( m_Atoms->GetNumericalList< size_t>());
      const size_t min_atoms( 1 + math::Statistics::MaximumValue( new_order.Begin(), new_order.End()));
      size_t mol_number( 0);
      for
      (
        chemistry::FragmentEnsemble::iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End());
        itr != itr_end;
        ++itr, ++mol_number
      )
      {
        if( itr->GetNumberAtoms() < min_atoms)
        {
          BCL_MessageCrt
          (
            "Could not reorder atoms on molecule #" + util::Format()( mol_number) + " because it had < atoms "
            + std::string( "than the max value given to -atoms: ") + util::Format()( min_atoms - 1)
          );
        }
        else
        {
          chemistry::AtomVector< chemistry::AtomComplete> atoms( itr->GetAtomInfo(), itr->GetBondInfo());
          atoms.Reorder( new_order);
          *itr = chemistry::FragmentComplete( atoms, itr->GetName(), itr->GetStoredProperties().GetMDLProperties());
        }
      }
    }

    //! @brief reorder the atoms in the ensemble
    //! @param ENSEMBLE ensemble in which to reorder the atoms of the molecule
    void MoleculeReorder::CanonicalizeAtoms( chemistry::FragmentEnsemble &ENSEMBLE) const
    {
      for
      (
        chemistry::FragmentEnsemble::iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->Canonicalize();
      }
    }

    //! @brief write out an ensemble, taking into account m_MaxMoleculesPerSDFFlag
    //! @param ENSEMBLE the ensemble to write
    void MoleculeReorder::WriteEnsemble( const chemistry::FragmentEnsemble &ENSEMBLE) const
    {
      util::Stopwatch timer( false);
      const std::string filename( m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      const size_t extension_position( filename.rfind( '.'));

      std::string filename_base
      (
        extension_position != std::string::npos
        ? filename.substr( 0, extension_position)
          : filename
      );

      std::string extension
      (
        extension_position != std::string::npos
        ? filename.substr( extension_position)
          : ""
      );

      const size_t max_molecules_total
      ( 
        m_OutputMaxFlag->GetFlag() ? 
          m_OutputMaxFlag->GetFirstParameter()->GetNumericalValue< size_t>()
          : std::numeric_limits< size_t>::max()
      );

      const size_t max_molecules_per_file
      (
        m_MaxMoleculesPerSDFFlag->GetFlag()
        ? m_MaxMoleculesPerSDFFlag->GetFirstParameter()->GetNumericalValue< size_t>()
          : std::numeric_limits< size_t>::max()
      );

      io::OFStream output;
      if( ENSEMBLE.GetSize() < max_molecules_per_file)
      {
        io::File::MustOpenOFStream( output, filename);
        size_t mol_no( 0);
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr_mol( ENSEMBLE.Begin()), itr_mol_end( ENSEMBLE.End());
          itr_mol != itr_mol_end && mol_no < max_molecules_total;
          ++itr_mol, ++mol_no
        )
        {
          itr_mol->WriteMDL( output);
        }
        io::File::CloseClearFStream( output);
        BCL_MessageStd
        (
          "Wrote ensemble containing " + util::Format()( mol_no)
          + " molecules to " + filename
          + " in " + timer.GetProcessDuration().GetTimeAsHourMinuteSecond()
        );
      }
      else
      {
        size_t n_to_write( std::min< size_t>( ENSEMBLE.GetSize(), max_molecules_total));
        const size_t number_files( ( n_to_write - 1) / max_molecules_per_file + 1);
        size_t required_width( 0);
        for
        (
          size_t number_files_dec( number_files - 1);
          number_files_dec > 0;
          number_files_dec /= 10, required_width++
        )
        {
        }

        util::Format file_id_formatter;
        file_id_formatter.Fill( '0').R().W( required_width);

        chemistry::FragmentEnsemble::const_iterator itr( ENSEMBLE.Begin());
        for
        (
          size_t molecule_id( 0), size( n_to_write);
          molecule_id < size;
        )
        {
          io::File::MustOpenOFStream
          (
            output,
            filename_base + file_id_formatter( molecule_id / max_molecules_per_file) + extension
          );
          for
          (
            const size_t max_molecule_id_this_file( std::min( size, molecule_id + max_molecules_per_file));
            molecule_id < max_molecule_id_this_file;
            ++itr, ++molecule_id
          )
          {
            itr->WriteMDL( output);
          }
          io::File::CloseClearFStream( output);
        }
        BCL_MessageStd
        (
          "Wrote ensemble containing " + util::Format()( n_to_write)
          + " molecules to file sequence beginning with "
          + filename_base + file_id_formatter( 0) + extension
          + " and ending with "
          + filename_base + file_id_formatter( number_files) + extension
          + " in " + timer.GetProcessDuration().GetTimeAsHourMinuteSecond()
        );
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeReorder::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoleculeReorder::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief function that sorts the ensemble by the numeric value of a particular property
    //! @param ENSEMBLE ensemble to sort
    //! @param PROPERTY property to sort by
    void MoleculeReorder::SortByNumericProperty
    (
      chemistry::FragmentEnsemble &ENSEMBLE,
      descriptor::CheminfoProperty PROPERTY
    ) const
    {
      util::Stopwatch timer( false);

      // old ensemble will hold the ensemble in its original order
      chemistry::FragmentEnsemble old_ensemble;
      swap( ENSEMBLE.GetMolecules().InternalData(), old_ensemble.GetMolecules().InternalData());
      // now ENSEMBLE is empty

      // create a map from property values to sorted ensembles with those values
      storage::Map< linal::Vector< float>, util::SiPtrList< chemistry::FragmentComplete> > sorted_ensembles;
      for
      (
        storage::List< chemistry::FragmentComplete>::iterator
        itr( old_ensemble.Begin()), itr_end( old_ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        sorted_ensembles[ PROPERTY->SumOverObject( *itr)].PushBack( util::SiPtr< chemistry::FragmentComplete>( &*itr));
      }

      // write the molecules into the ensemble in sorted order
      for
      (
        storage::Map< linal::Vector< float>, util::SiPtrList< chemistry::FragmentComplete> >::iterator
        itr( sorted_ensembles.Begin()), itr_end( sorted_ensembles.End());
        itr != itr_end;
        ++itr
      )
      {
        for
        (
          util::SiPtrList< chemistry::FragmentComplete>::iterator
          itr_ensemble( itr->second.Begin()), itr_ensemble_end( itr->second.End());
          itr_ensemble != itr_ensemble_end;
          ++itr_ensemble
        )
        {
          ENSEMBLE.PushBack( **itr_ensemble);
        }
      }

      BCL_MessageStd
      (
        "Sorted " + util::Format()( ENSEMBLE.GetSize())
        + " molecules by "
        + PROPERTY.GetString()
        + " in "
        + timer.GetProcessDuration().GetTimeAsHourMinuteSecond()
      );
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeReorder::GetReadMe() const
    {
      static std::string s_read_me =
        "MoleculeReorder reorders small molecules given in sdf files\n"
        " Molecule order can be randomized, reversed, and/or sorted by a property\n"
        " Additionally, atoms within molecules can be reordered canonically, or to an explicit order\n";

      return s_read_me;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeReorder::GetDescription() const
    {
      return "Reorder atoms in molecules, or reorder molecules in SD files by property";
    }

    const ApplicationType MoleculeReorder::MoleculeReorder_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeReorder(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
