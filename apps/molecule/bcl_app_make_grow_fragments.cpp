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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_groups.h"
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MakeGrowFragments
    //! @brief Application for removing hydrogens only from specific atoms
    //!
    //! @author geanesar
    //! @date Sept 25 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class MakeGrowFragments :
      public Interface
    {

    public:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      util::ShPtr< command::FlagInterface> m_NumberOpenValences;

      util::ShPtr< command::FlagInterface> m_FunctionalizationPointsFlag;

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      MakeGrowFragments();

      MakeGrowFragments *Clone() const
      {
        return new MakeGrowFragments( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_InputFilenamesFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_NumberOpenValences);
        sp_cmd->AddFlag( m_FunctionalizationPointsFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

    private:

      //! @brief enumerate molecules with open valences at each hydrogen atom position
      //! @param NUMBER_OPEN_VALENCES the number of open valences per molecule
      //! @param ENSEMBLE molecules corresponding to a different open valence (or open valence combo)
      util::ShPtr< chemistry::FragmentEnsemble> OpenMultipleValences
      (
        size_t NUMBER_OPEN_VALENCES,
        util::ShPtr< chemistry::FragmentEnsemble> ENSEMBLE
      ) const
      {
        size_t max_to_remove( NUMBER_OPEN_VALENCES);
        for( size_t iteration( 0); iteration < max_to_remove; ++iteration)
        {
          util::ShPtr< chemistry::FragmentEnsemble> new_ensemble( new chemistry::FragmentEnsemble());
          for
          (
              chemistry::FragmentEnsemble::const_iterator itr_mol( ENSEMBLE->Begin()), itr_mol_end( ENSEMBLE->End());
              itr_mol != itr_mol_end;
              ++itr_mol
          )
          {
            const chemistry::AtomVector< chemistry::AtomComplete> &atom_vector( itr_mol->GetAtomVector());
            size_t n_atoms( atom_vector.GetSize());
            storage::Vector< size_t> candidate_atoms;
            candidate_atoms.AllocateMemory( n_atoms);
            for( size_t atom_no( 0); atom_no < n_atoms; ++atom_no)
            {
              if( atom_vector( atom_no).GetElementType()->GetAtomicNumber() == 1)
              {
                continue;
              }
              if( atom_vector( atom_no).GetNumberCovalentlyBoundHydrogens())
              {
                candidate_atoms.PushBack( atom_no);
              }
            }

            // Go through all of the atoms that have hydrogens and create new molecules with one of those hydrogens removed
            for( size_t cand_no( 0); cand_no < candidate_atoms.GetSize(); ++cand_no)
            {
              size_t candidate_atom( candidate_atoms( cand_no));
              const storage::Vector< chemistry::BondConformational> &bonds( atom_vector( candidate_atom).GetBonds());
              size_t remove_atom( n_atoms);
              for( size_t bond_no( 0); bond_no < bonds.GetSize(); ++bond_no)
              {
                const chemistry::AtomConformationalInterface &target_atom( bonds( bond_no).GetTargetAtom());
                if( target_atom.GetElementType()->GetAtomicNumber() == 1)
                {
                  remove_atom = itr_mol->GetAtomIndex( target_atom);
                  break;
                }
              }
              if( remove_atom == n_atoms)
              {
                continue;
              }
              chemistry::AtomVector< chemistry::AtomComplete> new_vector( atom_vector);
              storage::Vector< size_t> keep_atoms;
              keep_atoms.AllocateMemory( n_atoms - 1);
              for( size_t i( 0); i < n_atoms; ++i)
              {
                if( i == remove_atom)
                {
                  continue;
                }
                keep_atoms.PushBack( i);
              }
              new_vector.Reorder( keep_atoms);
              new_ensemble->PushBack( chemistry::FragmentComplete( new_vector, ""));
            }
          }
          ENSEMBLE = new_ensemble;
        }
        return ENSEMBLE;
      }

      //! @brief functionalize a molecule at specific sites by removing specific hydrogen atoms
      //! @param ENSEMBLE molecules to have select hydrogen atoms removed
      util::ShPtr< chemistry::FragmentEnsemble> Functionalize
      (
        util::ShPtr< chemistry::FragmentEnsemble> ENSEMBLE
      ) const
      {
        storage::Vector< size_t> fxnl_pts; // functionalized points

        // convert the functionalization points to numeric values
        storage::Vector< size_t> fxnl_pts_str( m_FunctionalizationPointsFlag->GetNumericalList< size_t>());
        storage::Set< size_t> fxnl_pts_set;

        // For each molecule in input ensemble remove the hydrogen atoms individually
        util::ShPtr< chemistry::FragmentEnsemble> new_ensemble( new chemistry::FragmentEnsemble());
        for( chemistry::FragmentEnsemble::iterator ens_itr( ENSEMBLE->Begin()); ens_itr != ENSEMBLE->End(); ++ens_itr)
        {
          for( size_t i( 0), l( fxnl_pts_str.GetSize()); i < l; ++i)
          {
            // tell the user if they managed to pass a non-numerical character
            size_t point( fxnl_pts_str( i));
            if( point < ens_itr->GetSize())
            {
              fxnl_pts_set.Insert( point);
            }
          }
          fxnl_pts = storage::Vector< size_t>( fxnl_pts_set.Begin(), fxnl_pts_set.End());

          // You must have specified at least one functional atom
          BCL_Assert
          (
            !fxnl_pts.IsEmpty(),
            "No valid functionalization point(s) specified! Do better or use the 'open valences' option."
          );
          storage::Vector< size_t> indexing( storage::CreateIndexVector( ens_itr->GetSize()));
          chemistry::AtomVector< chemistry::AtomComplete> current_mol( ens_itr->GetAtomVector());
          for( storage::Vector< size_t>::iterator pts( fxnl_pts.Begin()); pts != fxnl_pts.End(); ++pts)
          {
            if( current_mol( *pts).GetNumberCovalentlyBoundHydrogens())
            {
              storage::Vector< chemistry::BondConformational> bonds( current_mol( *pts).GetBonds());
              for( auto bonds_itr( bonds.Begin()); bonds_itr != bonds.End(); ++bonds_itr)
              {
                const chemistry::AtomConformationalInterface &target_h( bonds_itr->GetTargetAtom());
                if( target_h.GetElementType()->GetAtomicNumber() == size_t( 1))
                {
                  // find where the hydrogen atom is in the atom indexer and remove it
                  size_t removal_index = current_mol.GetAtomIndex( target_h);
                  size_t h_index( indexing.Find( removal_index, 0));
                  indexing.RemoveElements( h_index, 1);
                }
              }
            }
            else
            {
              continue;
            }
          }
          current_mol.Reorder( indexing);
          new_ensemble->PushBack( chemistry::FragmentComplete( current_mol, ens_itr->GetName()));
        }
        return new_ensemble;
      }

    public:

      int Main() const
      {

        // Read in the molecules
        util::ShPtr< chemistry::FragmentEnsemble> ensemble( new chemistry::FragmentEnsemble());
        io::IFStream input;
        storage::Vector< std::string> filenames( m_InputFilenamesFlag->GetStringList());
        for( size_t fileno( 0); fileno < filenames.GetSize(); ++fileno)
        {
          io::File::MustOpenIFStream( input, filenames( fileno));
          ensemble->ReadMoreFromMdl( input, sdf::e_Saturate);
          io::File::CloseClearFStream( input);
        }

        // Produce molecules with multiple open valences
        if( m_FunctionalizationPointsFlag->GetFlag())
        {
          ensemble = Functionalize( ensemble);
        }
        else
        {
          ensemble = OpenMultipleValences( m_NumberOpenValences->GetFirstParameter()->GetNumericalValue< size_t>(), ensemble);
        }

        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
        ensemble->WriteMDL( output);
        io::File::CloseClearFStream( output);

        return 0;
      } // Main()

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

    private:

      static const ApplicationType MakeGrowFragments_Instance;

    }; // MakeGrowFragments

      //! @brief standard constructor
    MakeGrowFragments::MakeGrowFragments() :
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output", "flag selecting the output file name",
          command::Parameter
          (
            "filename", "filename for output sdf"
          )
        )
      ),
      m_InputFilenamesFlag
      (
        new command::FlagDynamic
        (
          "input_filenames",
          "filenames for input sdf",
          command::Parameter
          (
            "filenames of derived structure sdfs",
            "name of files containing derived structures",
            command::ParameterCheckFileExistence()
          ),
          1,
          21
        )
      ),
      m_NumberOpenValences
      (
        new command::FlagDynamic
        (
          "open_sites",
          "the number of open sites to leave on molecules",
          command::Parameter
          (
            "open valences",
            "the number of open valences that should be on each molecule",
            command::ParameterCheckRanged< size_t>( 1, 3),
            "1"
          )
        )
      ),
      m_FunctionalizationPointsFlag
      (
        new command::FlagDynamic
        (
          "functionalize",
          "which atoms (zero-indexed) of the input molecule should be functionalized with reactive groups",
          command::Parameter
          (
            "atoms",
            "the atom indices to modify during preprocessing",
            command::ParameterCheckRanged< size_t>( 0, 999)
          )
        )
      )
    {
    }

    const ApplicationType MakeGrowFragments::MakeGrowFragments_Instance
    (
      GetAppGroups().AddAppToGroup( new MakeGrowFragments(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
