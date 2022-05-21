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

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_site.h"
#include "storage/bcl_storage_table.hpp"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ExtractSites
    //! @brief this application extracts SITES form PDB files
    //! @details Sites define binding or catalytic locations within a protein structure. These are defined over by
    //!          contacting or involved residues that are scattered over the site. This application, seeks to extract
    //!          the topological fragment, which contains the majority of these residues.
    //!
    //! @author woetzen
    //! @date Feb 3, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ExtractSites :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! command line parameters
      util::ShPtr< command::FlagStatic> m_PdbListFlag;
      util::ShPtr< command::Parameter> m_PdbPathParam;
      util::ShPtr< command::Parameter> m_PdbListFileParam;
      util::ShPtr< command::Parameter> m_PdbPrefixParam;
      util::ShPtr< command::Parameter> m_PdbSuffixParam;
      util::ShPtr< command::Parameter> m_PdbPathHierarchyParam;

      //! table with sites
      mutable storage::Table< std::string> m_SitesTable;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ExtractSites();

    public:

      //! @brief Clone function
      //! @return pointer to new ExtractSites
      ExtractSites *Clone() const
      {
        return new ExtractSites( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType ExtractSites_Instance;

    }; // class ExtractSites

    // default constructor
    ExtractSites::ExtractSites() :
      m_PdbListFlag
      (
        new command::FlagStatic
        (
          "pdblist",
          "\tfrom a given path and a list of pdb codes, sites are derived"
        )
      ),
      m_PdbPathParam
      (
        new command::Parameter
        (
          "pdb_path",
          "path of pdb files, where for each code, there is a corresponding pdb file present"
        )
      ),
      m_PdbListFileParam
      (
        new command::Parameter( "list", "a list file of pdb codes", command::ParameterCheckFileExistence())
      ),
      m_PdbPrefixParam
      (
        new command::Parameter( "prefix", "prefix for each pdb file")
      ),
      m_PdbSuffixParam
      (
        new command::Parameter( "suffix", "suffix for each pdb file")
      ),
      m_PdbPathHierarchyParam
      (
        new command::Parameter
        (
          "hierarchy",
          "boolean to indicate whether a pdb hierarchy is used in input paths so for pdbtag 1abc.pdb it looks at {path}/ab/1abcA.pdb",
          command::ParameterCheckRanged< size_t>( 0, 1), "0"
        )
      )
    {
      // attach flags to parameters
      m_PdbListFlag->PushBack( m_PdbPathParam);
      m_PdbListFlag->PushBack( m_PdbListFileParam);
      m_PdbListFlag->PushBack( m_PdbPrefixParam);
      m_PdbListFlag->PushBack( m_PdbSuffixParam);
      m_PdbListFlag->PushBack( m_PdbPathHierarchyParam);
    }

    util::ShPtr< command::Command> ExtractSites::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      //flag and parameter for listfile
      sp_cmd->AddFlag( m_PdbListFlag);

      // pdb factory flags
//      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
//      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());
//      pdb::Factory::GetFlagAAClass()->GetParameterList()( 0)->SetDefaultParameter( biol::GetAAClasses().e_AABackBone);
//      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());
//      sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());
//      sp_cmd->AddFlag( pdb::Factory::GetFlagBiomolecule());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // end
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ExtractSites::Main() const
    {
      // setup table of sites
      const std::string header_strings[] =
      {
        "pdb_id",
        "site_name",
        "ligand_name",
        "ligand_fullname",
        "min_res_chain_id",
        "min_res_id",
        "max_res_chain_id",
        "max_res_id",
        "nr_residues",
        "motif_length",
        "nr_helices",
        "nr_strands",
        "nr_coils",
        "nr_sses",
        "length"
      };

      m_SitesTable = storage::Table< std::string>( storage::TableHeader( storage::Vector< std::string>( 15, header_strings)));

//      // create factory
//      pdb::Factory pdb_factory;

      // list of pdb files
      const std::string pdb_path( m_PdbPathParam->GetValue() + PATH_SEPARATOR);
      const std::string pdb_prefix( m_PdbPrefixParam->GetValue());
      const std::string pdb_suffix( m_PdbSuffixParam->GetValue());
      const bool pdb_hierarchy( m_PdbPathHierarchyParam->GetNumericalValue< bool>());

      io::IFStream read_pdb_list;
      io::File::MustOpenIFStream( read_pdb_list, m_PdbListFileParam->GetValue());

      const util::StringReplacement ligand_fullname_rectifier( util::StringReplacement::e_Any, " ", "_");

      while( read_pdb_list.good())
      {
        std::string pdb_code;
        read_pdb_list >> pdb_code;

        // sometimes it was not read to the end of line
        if( util::TrimString( pdb_code).empty())
        {
          continue;
        }

        std::transform( pdb_code.begin(), pdb_code.end(), pdb_code.begin(), tolower);

        // filename
        const std::string pdb_filename( pdb_path + ( pdb_hierarchy ? ( pdb_code.substr( 1, 2) + PATH_SEPARATOR) : "") + pdb_prefix + pdb_code + pdb_suffix);
        BCL_MessageCrt( "processing file: " + pdb_filename);

        io::IFStream read_pdb;
        if( !io::File::TryOpenIFStream( read_pdb, pdb_filename))
        {
          BCL_MessageCrt( "unable to open: " + pdb_filename);
          continue;
        }

        pdb::Handler handler( read_pdb);

        // does handler have at least one chain
        if( handler.GetProteinChains().empty())
        {
          BCL_MessageCrt( "no amino acid chain found: " + pdb_code);
          continue;
        }

        // get the sites
        util::ShPtrList< pdb::Site> sites( handler.GetSites());

        // iterate through sites
        for( util::ShPtrList< pdb::Site>::const_iterator site_itr( sites.Begin()), site_itr_end( sites.End()); site_itr != site_itr_end; ++site_itr)
        {
          // site
          const pdb::Site &current_site( **site_itr);

          // skipping sites without ligand
          if( !current_site.GetLigand().IsDefined())
          {
            continue;
          }

          // ligand
          const pdb::Ligand &current_ligand( *current_site.GetLigand());

          storage::Row< std::string> &row( m_SitesTable.InsertRow( pdb_code, true));
          row[ "pdb_id"] = pdb_code;
          row[ "site_name"] = current_site.GetName();
          row[ "ligand_name"] = current_ligand.GetResidueName();
          std::string new_fullname( current_ligand.GetFullname());
          ligand_fullname_rectifier.ReplaceEachIn( new_fullname);
          row[ "ligand_fullname"] = new_fullname;
          // residues
          const storage::List< pdb::ResidueSimple> &residues( current_site.GetChainResidues());
          if( !residues.IsEmpty())
          {
            const char min_chain_id( residues.FirstElement().GetChainID());
            const char max_chain_id( residues.LastElement().GetChainID());
            const int min_seq_id( residues.FirstElement().GetPDBID());
            const int max_seq_id( residues.LastElement().GetPDBID());
            row[ "min_res_chain_id"] = std::string( 1, min_chain_id);
            row[ "min_res_id"] = util::Format()( min_seq_id);
            row[ "max_res_chain_id"] = std::string( 1, max_chain_id);
            row[ "max_res_id"] = util::Format()( max_seq_id);
            row[ "nr_residues"] = util::Format()( residues.GetSize());
            int diff( 999);
            if( min_chain_id == max_chain_id)
            {
              diff = max_seq_id - min_seq_id;
            }
            row[ "motif_length"] = util::Format()( diff);
            row[ "nr_helices"] = "0";
            row[ "nr_strands"] = "0";
            row[ "nr_coils"]   = "0";
            row[ "nr_sses"]    = "0";
            row[ "length"]     = "0";
          }

          // locate the sses
          const storage::List< assemble::LocatorAA> aa_locators( current_site.AALocators());

          assemble::Domain site_domain;
          storage::Set< char> chain_ids;
          // create a model from the pdb
          const assemble::ProteinModel model( pdb::Factory().ProteinModelFromPDB( handler));

          // iterate through all sses in the model, and check if the sse contains one of the aas
          const util::SiPtrVector< const assemble::SSE> model_sses( model.GetSSEs());

          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( model_sses.Begin()), sse_itr_end( model_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // iterate through all aa locators
            for
            (
              storage::List< assemble::LocatorAA>::const_iterator
                loc_itr( aa_locators.Begin()), loc_itr_end( aa_locators.End());
              loc_itr != loc_itr_end;
              ++loc_itr
            )
            {
              // TODO: check for the first residue, how many residues are to the left; for the last how many residues are to the right;
              // if they are within SSEs (helix, strand) and there are not enough residues space to the border, take the next helix or strand, as attachment points for the fragment
              // if the amino acid can be found in this sse
              if( loc_itr->Locate( **sse_itr).IsDefined())
              {
                // find this sse in the model and insert it into the domain
                const util::ShPtr< assemble::SSE> sp_sse( model.FindSSE( **sse_itr));
                if( sp_sse.IsDefined())
                {
                  site_domain.Insert( sp_sse);
                  chain_ids.Insert( sp_sse->GetChainID());
                }

                // move on to next sse
                break;
              }
            } // locator
          } // sses

          // domain for that site
          // check that is only a single chain
          if( chain_ids.GetSize() > 1)
          {
            BCL_MessageVrb( "skipping site with residues in multiple chains: " + current_site.GetName());
            continue;
          }

          // also add the sse before the first and after the last SSE in that fragment, and add the missing SSEs in between
          if( site_domain.IsEmpty())
          {
            continue;
          }

          const util::ShPtr< assemble::SSE> &first_sse_domain( *site_domain.GetData().Begin());
          const util::ShPtr< assemble::SSE> &last_sse_domain( *site_domain.GetData().ReverseBegin());

          // iterate over model sses
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_itr( model_sses.Begin()), sse_itr_end( model_sses.End());

          util::SiPtr< const assemble::SSE> preceeding_sse( **sse_itr);

          // there is no preceeding sse, if the first sse in the model is the first sse in the site domain
          if( first_sse_domain != preceeding_sse)
          {
            while( sse_itr != sse_itr_end)
            {
              if( ( *sse_itr)->DoesPrecede( *first_sse_domain))
              {
                preceeding_sse = *sse_itr;
                break;
              }
              ++sse_itr;
            }
          }

          // fill the site with all sses in between
          while( sse_itr != sse_itr_end && *sse_itr != last_sse_domain)
          {
            // find this sse in the model and insert it into the domain
            const util::ShPtr< assemble::SSE> sp_sse( model.FindSSE( **sse_itr));
            if( sp_sse.IsDefined())
            {
              site_domain.Insert( sp_sse);
              chain_ids.Insert( sp_sse->GetChainID());
            }
            ++sse_itr;
          }
          // goto proceeding sse
          util::SiPtr< const assemble::SSE> proceeding_sse( last_sse_domain);
          if( sse_itr != sse_itr_end && ++sse_itr != sse_itr_end)
          {
            proceeding_sse = *sse_itr;
          }

          BCL_MessageStd
          (
            "preceeding_sse: " + preceeding_sse->GetIdentification() +
            "\nfirst_sse:      " + first_sse_domain->GetIdentification() +
            "\nlast_sse:       " + last_sse_domain->GetIdentification() +
            "\nproceeding_sse: " + proceeding_sse->GetIdentification()
          );

          if( preceeding_sse->GetType()->IsStructured())
          {
            // insert into domain
            // find this sse in the model and insert it into the domain
            const util::ShPtr< assemble::SSE> sp_sse( model.FindSSE( *preceeding_sse));
            if( sp_sse.IsDefined())
            {
              site_domain.Insert( sp_sse);
            }
          }
          if( proceeding_sse->GetType()->IsStructured())
          {
            // insert into domain
            // find this sse in the model and insert it into the domain
            const util::ShPtr< assemble::SSE> sp_sse( model.FindSSE( *proceeding_sse));
            if( sp_sse.IsDefined())
            {
              site_domain.Insert( sp_sse);
            }
          }

          row[ "nr_helices"] = util::Format()( site_domain.GetNumberSSE( biol::GetSSTypes().HELIX));
          row[ "nr_strands"] = util::Format()( site_domain.GetNumberSSE( biol::GetSSTypes().STRAND));
          row[ "nr_coils"]   = util::Format()( site_domain.GetNumberSSE( biol::GetSSTypes().COIL));
          row[ "nr_sses"]    = util::Format()( site_domain.GetNumberSSEs());
          row[ "length"]     = util::Format()( ( *site_domain.GetData().ReverseBegin())->GetLastAA()->GetSeqID() - ( *site_domain.GetData().Begin())->GetFirstAA()->GetSeqID());
        } // site
      } // pdbs
      io::File::CloseClearFStream( read_pdb_list);

      // print table
      m_SitesTable.WriteFormatted( util::GetLogger());

      // end
      return 0;
    }

    const ApplicationType ExtractSites::ExtractSites_Instance
    (
      GetAppGroups().AddAppToGroup( new ExtractSites(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
