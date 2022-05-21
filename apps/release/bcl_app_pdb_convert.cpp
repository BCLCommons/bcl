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
#include "assemble/bcl_assemble_aa_exposures.h"
#include "assemble/bcl_assemble_biomolecule.h"
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_printer_protein_model_multimer.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_moment_of_inertia.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_quality.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "biol/bcl_biol_dssp.h"
#include "biol/bcl_biol_membrane.h"
#include "biol/bcl_biol_protein_params.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_flag_static_and_dynamic.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_aa_sasa.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_sum_function.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_printer_membrane.h"
#include "score/bcl_score_sse_pair_packing.h"
#include "score/bcl_score_strand_pairing.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"
#include "sspred/bcl_sspred_mahssmi.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PDBConvert
    //! @brief Class converts pdb file into bcl-style pdb file
    //! @details Class reads a given pdb to and converts it to a bcl-style pdb with idealized SSEs,
    //! has options to print out the fasta and bcl style pdb output for the given pdb
    //!
    //! @see @link example_app_pdb_convert.cpp @endlink
    //! @author woetzen, alexanns, karakam, lib14
    //! @date 03/24/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PDBConvert :
      public InterfaceRelease
    {
    public:

    ///////////
    // types //
    ///////////

      //! @enum enumerate the different ways pdb files can be written
      enum PdbWriteFilesType
      {
        e_Split,
        e_Full,
        e_Both,
        s_NumberPdbWriteFilesType
      };

      //! @brief PdbWriteFilesType as string
      //! @param PDB_WRITE_FILES_TYPE the PdbWriteFilesType
      //! @return the string for the PdbWriteFilesType
      static const std::string &GetPdbWriteFilesTypeDescriptor( const PdbWriteFilesType &PDB_WRITE_FILES_TYPE)
      {
        static const std::string s_descriptors[] =
        {
          "Split",
          "Full",
          "Both",
          GetStaticClassName< PdbWriteFilesType>()
        };

        return s_descriptors[ PDB_WRITE_FILES_TYPE];
      }

      //! @typedef PdbWriteFilesTypeEnum is used for I/O of PdbWriteFilesType
      typedef util::WrapperEnum< PdbWriteFilesType, &GetPdbWriteFilesTypeDescriptor, s_NumberPdbWriteFilesType> PdbWriteFilesTypeEnum;

      //! @enum enumerate the different rosetta loop files
      enum RosettaLoopFileType
      {
        e_CCD,
        e_KIC,
        s_NumberRosettaLoopFileType
      };

      //! @brief RosettaLoopFileType as string
      //! @param ROSETTA_LOOP_FILE_TYPE the RosettaLoopFileType
      //! @return the string for the RosettaLoopFileType
      static const std::string &GetRosettaLoopFileTypeDescriptor( const RosettaLoopFileType &ROSETTA_LOOP_FILE_TYPE)
      {
        static const std::string s_descriptors[] =
        {
          "CCD",
          "KIC",
          GetStaticClassName< RosettaLoopFileType>()
        };

        return s_descriptors[ ROSETTA_LOOP_FILE_TYPE];
      }

      //! @typedef RosettaLoopFileType is used for I/O of RosettaLoopFileType
      typedef util::WrapperEnum< RosettaLoopFileType, &GetRosettaLoopFileTypeDescriptor, s_NumberRosettaLoopFileType> RosettaLoopFileTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! input pdb file parameter
      util::ShPtr< command::ParameterInterface> m_PdbFilename;

      //! flag for writing fasta
      util::ShPtr< command::FlagInterface> m_WriteFastaFlag;

      //! flag for input fasta to write SEQRES to pdb
      util::ShPtr< command::FlagInterface> m_PdbFromFastaFlag;

      //! flag for individual chain
      util::ShPtr< command::FlagInterface> m_ChainsFlag;

      //! flag for output prefix
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

      //! flag for writing bclpdb file
      util::ShPtr< command::FlagInterface> m_BclPdbFlag;

      //! flag to output dssp
      util::ShPtr< command::FlagInterface> m_DSSPFlag;

      //! flag to indicate whether body information should be written to output pdb
      util::ShPtr< command::FlagInterface> m_BodyOutputFlag;

      //! flag to idealize pdb
      util::ShPtr< command::FlagInterface> m_IdealizeFlag;

      //! flag for passing a pdbtm xml file containing the membrane and transformation matrix
      util::ShPtr< command::FlagInterface> m_PdbTmXmlFlag;

      //! flag for aligning protein against the moment of inertia
      util::ShPtr< command::FlagStatic> m_MomentOfInertiaFlag;
      util::ShPtr< command::ParameterInterface> m_MomentOfInertiaPropertyParam;
      util::ShPtr< command::ParameterInterface> m_MomentOfInertiaExposureParam;

      //! flag for writing a loops file
      util::ShPtr< command::FlagInterface> m_RosettaLoopFileFlag;

      //! flag for changing an empty chain id to a given chain id
      util::ShPtr< command::FlagDynamic> m_RenameChainIDsFlag;

      //! flag for adding side chains to a given pdb
      util::ShPtr< command::FlagInterface> m_AddSideChainsFlag;

      //! flag for defining the pdb to superimpose onto
      util::ShPtr< command::FlagStaticAndDynamic> m_SuperimposePDBFlag;
      util::ShPtr< command::Parameter> m_SuperimposPDBFileParam;
      util::ShPtr< command::Parameter> m_SuperimposeQualityMeasureParam;

      //! flag for applying a transformation to the protein
      util::ShPtr< command::FlagInterface> m_TransformationFlag;

      //! flag for writing topology graphviz file
      util::ShPtr< command::FlagInterface> m_TopologyGraphFlag;

      //! flag for writing MAHSSMI analysis of secondary structure & membrane topology
      util::ShPtr< command::FlagInterface> m_MahssmiFlag;

      //! flag for writing CIPhiPsi analysis of secondary structure & membrane topology
      util::ShPtr< command::FlagInterface> m_CIPhiPsiFlag;

      //! flag for specifying the native multimer
      util::ShPtr< command::FlagInterface> m_MultimerNativePDBFileFlag;

      //! flag to write protein params for each chain
      util::ShPtr< command::FlagStatic> m_ProteinParamsFlag;

      //! flag to generate a single biomolecule from symmetrically related chains in the protein
      util::ShPtr< command::FlagStatic> m_SuperimposeToBiomoleculeFlag;
      util::ShPtr< command::Parameter> m_SuperimposeToBiomoleculeParam;
      util::ShPtr< command::Parameter> m_SuperimposeToBiomoleculeThresholdParam;

      //! flag to merge broken transmembrane helices into a single helix
      util::ShPtr< command::FlagStatic> m_MergeBrokenTMHelixFlag;
      util::ShPtr< command::Parameter> m_MergeBrokenTMHelixSlopeParam;
      util::ShPtr< command::Parameter> m_MergeBrokenTMHelixHalfMembraneWidthParam;
      util::ShPtr< command::Parameter> m_MergeBrokenTMHelixNewTMMinSSESizeParam;
      util::ShPtr< command::Parameter> m_MergeBrokenTMHelixRedefineNonTMParam;

      //! flag for splitting an ensemble
      util::ShPtr< command::FlagInterface> m_SplitEnsembleFlag;

      //! flag for writing as .sdf file (for oligomers with complete coordinates only)
      util::ShPtr< command::FlagInterface> m_WriteSDFFlag;

      //! flag for writing .sasa file that contains neighbor count and neighbor vector
      util::ShPtr< command::FlagStatic> m_WriteSASAFlag;
      util::ShPtr< command::Parameter> m_MinimalSequenceSeparation;

      static const char s_ChainIDPairSeparator;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      PDBConvert();

    public:

      //! @brief Clone function
      //! @return pointer to new PDBConvert
      PDBConvert *Clone() const
      {
        return new PDBConvert( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Converts PDBs to/from FASTAs and performs various deterministic manipulations of PDBs";
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief write a file giving the loop boundaries for a given protein model
      //! @param CHAIN is the chain whose loop boundaries will be written to a file
      void WriteLoopsFile( const assemble::Chain &CHAIN) const;

      //! @brief write protein parameters for the given sequence into a params file
      //! @param SEQUENCE the sequence to use to derive the protein parameters
      void WriteProteinParamsFile( const biol::AASequence &SEQUENCE) const;

      //! @brief write protein parameters for each sequence into a params file
      //! @param MODEL model with sequences to use to derive the protein parameters
      void WriteProteinParamsFiles( const assemble::ProteinModel &MODEL) const;

      //! @brief determines helical transmembrane sequence stretches and merges broken segments
      //! @param MODEL the model whose helical tm segments will be merged into continuous tm helices if they are broken
      //! @return protein model with continuous transmembrane helices
      assemble::ProteinModel MergeHelixBasedTMStretches( const assemble::ProteinModel &MODEL) const;

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

    private:

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief process the model according to the command line flags
      //! @param PROTEIN_MODEL protein model to be processed
      //! @param IDENTIFIER optional string identifier to be used with output files
      void ProcessModel( assemble::ProteinModel &PROTEIN_MODEL, const std::string &IDENTIFIER = "") const;

      static const ApplicationType PDBConvert_Instance;

    }; // class PDBConvert

    //! @brief the Main function
    //! @return error code - 0 for success
    int PDBConvert::Main() const
    {
      if( m_PdbFilename->GetWasSetInCommandLine())
      {
        // instantiate pdb
        const pdb::Factory factory;

        // if splitting ensemble
        if( m_SplitEnsembleFlag->GetFlag())
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_PdbFilename->GetValue());
          pdb::Handler pdb( read);
          io::File::CloseClearFStream( read);

          // get the ensemble
          assemble::ProteinEnsemble ensemble( factory.ProteinEnsembleFromPDB( pdb));

          // iterate over the ensemble
          size_t model_no( 0);
          for
          (
            assemble::ProteinEnsemble::iterator itr( ensemble.Begin()), itr_end( ensemble.End());
            itr != itr_end; ++itr, ++model_no
          )
          {
            ProcessModel( **itr, "_" + util::Format()( model_no));
          }
        }
        // build the model and process it
        else
        {
          // use ProteinModelFromPDBFilename to retrieve filename; this sets ProteinModelData::e_PDBFile
          // unlike initializing a pdb::Handler with an istream, as is done above
          assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( m_PdbFilename->GetValue()));
          ProcessModel( model);
        }
      } // if pdb was given

      // if pdb from fasta flag was given
      if( m_PdbFromFastaFlag->GetFlag() && !m_PdbFromFastaFlag->GetParameterList().IsEmpty())
      {
        pdb::Factory factory( biol::GetAAClasses().e_AA);
        assemble::ProteinModel model;
        // iterate over all given fasta sequences and insert them as chains into the model
        for
        (
          util::ShPtrVector< command::ParameterInterface>::const_iterator
            itr( m_PdbFromFastaFlag->GetParameterList().Begin()),
            itr_end( m_PdbFromFastaFlag->GetParameterList().End());
          itr != itr_end;
          ++itr
        )
        {
          const std::string &fasta_file_name( ( *itr)->GetValue());
          const std::string::size_type extension_pos( fasta_file_name.find( ".fasta"));
          if( extension_pos == std::string::npos || extension_pos == 0)
          {
            BCL_MessageStd
            (
              "given file does not have extension .fasta " + fasta_file_name + "! Skipping!"
            );
            continue;
          }

          // chain id for sequence from last char of fasta filename without extension
          char chain_id( fasta_file_name[ extension_pos - 1]);
          if( chain_id == '_')
          {
            chain_id = ' ';
          }

          if( model.GetChain( chain_id).IsDefined())
          {
            BCL_MessageCrt( "already read chain with that identifier " + util::Format()( chain_id) + "! Skipping!");
            continue;
          }
          io::IFStream read;
          io::File::MustOpenIFStream( read, fasta_file_name);
          model.Insert( factory.ChainFromFastaStream( chain_id, read));
          io::File::CloseClearFStream( read);
        }
        const std::string model_file_name( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "fasta.pdb");
        io::OFStream write;
        io::File::MustOpenOFStream( write, model_file_name);
        factory.WriteModelToPDB( model, write);
        io::File::CloseClearFStream( write);

        // write protein params if desired
        if( m_ProteinParamsFlag->GetFlag())
        {
          WriteProteinParamsFiles( model);
        }
      }

      return 0;
    } // end Main

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &PDBConvert::GetReadMe() const
    {
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::PDBConvert, terms of use, appropriate citation, installation "
        "procedures, BCL::PDBConvert execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::PDBConvert?\n"
        "This application is designed to perform tasks in converting pdbs. Fasta sequences can be translated in pdb "
        "SEQRES information. Pdb files can be read which will issue warnings, if their is any inconsistency in the "
        "SEQRES vs. ATOM section and if secondary structure definitions are off. The pdb con then be renumbered, split "
        "by chain, or fastas can be generated. The relevant bio molecule can be generated or with the help of a "
        "pdbtm-xml file, the proper membrane protein structure can be generated and translated into the membrane. "
        "Rosetta loop-files can be generated and the according missing residues can be insterted into the pdb with "
        "zero coordinates. Unnatural amino acids can be converted into natural amino acids."
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::PDBConvert.\n"
        "When using BCL::PDBConvert in a publication, please cite the following publication describing the application's "
        "development:\n"
        "\n"
        "Woetzel, N., Karakas, M., Alexander, N., Straitzbichler, R. and Meiler, J. (2011). bcl::PDBConvert - a pdb "
        "file processing and correction tool. In preparation.\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::PDBConvert.\n"
        "Two main usages are of particular importance:\n"
        "1. Use with a fasta file:\n"
        "if used with -pdb_from_fasta it is possible to print a fasta sequence into a pdb (-bcl_pdb) with SEQRES lines "
        "and possibly zero coordinates for all residues (ATOM lines) using -loop.\n"
        "2. Use with an input pdb:\n"
        "If used with -pdb the input pdb will be read and messages on possible errors in the pdb will be generated. It "
        "is possible to restrict reading of chains by invoking (-chains), secondary structure elements with certain "
        "size (-min_sse_size) and the atoms to be read (-aaclass). A fasta file may be written for each chain by "
        "invoking -fasta. An output pdb is only written if -bcl_pdb is specified and the flag enables the following "
        "functionality:\n"
        "The output pdb resid will be renumbered, so that the first residue in the SEQRES will be 1 (except when "
        "-write_pdb_res_ids is used), and atom serials will be renumbered starting with 1 for the first atom in the "
        "output pdb (except when used with -write_pdb_atom_ids).  Side chains can be added (or overwritten) using "
        "-side_chains. The biomolecule can be generated using REMARK 350 of the pdb with flag -biomolecule. "
        "Transmembrane proteins can be handled by invoking -pdbtm_xml with a .xml file from the pdbtm (it will "
        "transform the protein into the virtual membrane (z-axis as normal on membrane) and will eventually use the "
        "given bio molecule information from that .xml file. Chain ids can be renamed using -rename_chain_id. The "
        "output pdb can be idealized (straight helices and strands - they are not bent  afterwards) with -idealize. "
        "Undefined coordinates can be printed as 0 coordinates if output will be used by rosetta or other programs "
        "with -write_zero_coordinates. If unnatural aatypes are used (not one of the 20 standard AA types) one can "
        "convert them to their natural amino acid \"parent\" with -convert_to_natural_aa_type. When loops need to be "
        "build for undefined amino acids, a rosetta loop file can be written with -loop_file_rosetta.\n"
        "Additional functionality is available through other flags, that are not mentioned here.\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl_pdb_convert.exe -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        "For more general information about the product, type bcl_pdb_convert.exe -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::PDBConvert.\n"
        "BCL::PDBConvert is under ongoing further development. Other logical errors that one can encounter will be "
        "checked in future versions. Specifying a set of residues that should be included in the loop file will be an "
        "option. More unnatural amino acids will be added, as they are encountered. A DSSP handler will be added, so "
        "that the program does not rely on the secondary structure defined in the pdb file.\n"
        "\n"
        + DefaultSectionSeparator()
      );

      return s_readme;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> PDBConvert::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddParameter( m_PdbFilename);

      sp_cmd->AddFlag( m_WriteFastaFlag);
      sp_cmd->AddFlag( m_PdbFromFastaFlag);
      sp_cmd->AddFlag( m_BclPdbFlag);
      sp_cmd->AddFlag( m_DSSPFlag);
      sp_cmd->AddFlag( m_BodyOutputFlag);
      sp_cmd->AddFlag( m_ChainsFlag);
      sp_cmd->AddFlag( m_OutputPrefixFlag);
      sp_cmd->AddFlag( m_IdealizeFlag);

      // biomolecule multimer
      sp_cmd->AddFlag( m_PdbTmXmlFlag);

      // moment of inertia
      sp_cmd->AddFlag( m_MomentOfInertiaFlag);

      sp_cmd->AddFlag( m_RenameChainIDsFlag);

      sp_cmd->AddFlag( m_AddSideChainsFlag);
      sp_cmd->AddFlag( m_SuperimposePDBFlag);
      sp_cmd->AddFlag( m_TransformationFlag);
      sp_cmd->AddFlag( m_TopologyGraphFlag);
      sp_cmd->AddFlag( m_MahssmiFlag);
      sp_cmd->AddFlag( m_CIPhiPsiFlag);

      pdb::Factory::GetFlagAAClass()->GetParameterList()( 0)->SetDefaultParameter( biol::GetAAClasses().e_AAComplete.GetName());
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
      sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());
      sp_cmd->AddFlag( pdb::Factory::GetFlagWritePDBResID());
      sp_cmd->AddFlag( pdb::Factory::GetFlagWritePDBAtomID());
      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());
      sp_cmd->AddFlag( pdb::Factory::GetFlagPDBAtomName());
      sp_cmd->AddFlag( pdb::Factory::GetFlagBiomolecule());

      // loop relevant flags
      sp_cmd->AddFlag( m_RosettaLoopFileFlag);
      sp_cmd->AddFlag( pdb::Factory::GetFlagWriteZeroCoordinatesForUndefinedAminoAcids());
      sp_cmd->AddFlag( pdb::Factory::GetFlagWriteHydrogens());
      sp_cmd->AddFlag( pdb::Handler::GetFlagHelixClasses());

      sp_cmd->AddFlag( m_MultimerNativePDBFileFlag);
      sp_cmd->AddFlag( m_ProteinParamsFlag);
      sp_cmd->AddFlag( m_SuperimposeToBiomoleculeFlag);

      // merge broken transmembrane helices into one transmembrane helix
      sp_cmd->AddFlag( m_MergeBrokenTMHelixFlag);

      sp_cmd->AddFlag( m_SplitEnsembleFlag);

      // add flag for conversion into sdf
      sp_cmd->AddFlag( m_WriteSDFFlag);

      // add flag for writring SASA file
      sp_cmd->AddFlag( m_WriteSASAFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief WriteLoopsFile writes a file giving the loop boundaries for a given protein model
    //! @param CHAIN is the chain whose loop boundaries will be written to a file
    void PDBConvert::WriteLoopsFile( const assemble::Chain &CHAIN) const
    {
      const RosettaLoopFileTypeEnum file_type( m_RosettaLoopFileFlag->GetFirstParameter()->GetValue());

      const char chain_id( CHAIN.GetChainID());
      // output filename
      const std::string loops_filename
      (
        m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ( chain_id == ' ' ? '_' : chain_id) + ".loops"
      );
      const int first_seq_id_chain( CHAIN.GetSequence()->GetFirstAA()->GetSeqID());
      const int last_seq_id_chain( CHAIN.GetSequence()->GetLastAA()->GetSeqID());

      // create io::OFStream write
      io::OFStream write;

      // open "write" and bind it to "loops_filename"
      io::File::MustOpenOFStream( write, loops_filename);

      // create SiPtrVector "coil_sses" and initialize with the COIL SSEs of "CHAIN"
      util::SiPtrVector< const assemble::SSE> coil_sses( CHAIN.GetSSEs( biol::GetSSTypes().COIL));

      // initialize loop segments
      storage::List< storage::VectorND< 2, int> > loop_segments;
      int previous_end_seq_id( 0);

      // iterate through "coil_sses"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( coil_sses.Begin()), sse_itr_end( coil_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        static const int s_min_loop_length( 5);

        int start_seq_id( std::max( first_seq_id_chain, ( *sse_itr)->GetFirstAA()->GetSeqID() - 1));
        int end_seq_id( std::min( last_seq_id_chain, ( *sse_itr)->GetLastAA()->GetSeqID() + 1));

        // iterate to increase loop size if too small
        while( end_seq_id - start_seq_id < s_min_loop_length - 1)
        {
          start_seq_id = std::max( first_seq_id_chain, start_seq_id - 1);
          end_seq_id = std::min( last_seq_id_chain, end_seq_id + 1);
        }

        // if this is the first entry
        if( loop_segments.IsEmpty())
        {
          loop_segments.PushBack( storage::VectorND< 2, int>( start_seq_id, end_seq_id));
        }
        // has previous loops
        else
        {
          // check for no overlap
          if( start_seq_id > previous_end_seq_id + 1)
          {
            loop_segments.PushBack( storage::VectorND< 2, int>( start_seq_id, end_seq_id));
          }
          // overlap w/ previous entry exists
          else
          {
            // combine the two loop segments
            loop_segments.Last()->Second() = end_seq_id;
          }
        }
        previous_end_seq_id = end_seq_id;
      }

      const util::Format format_seq_id( util::Format().R().Fill( ' ').ForceW().W( 3));

      // iterate through loops
      for
      (
        storage::List< storage::VectorND< 2, int> >::const_iterator sse_itr( loop_segments.Begin()),
          sse_itr_end( loop_segments.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // skip if we do expect undefined coordinates for adjacent aas
        // output the seq id of the first and last amino acid of the sse currently denoted by the sse denoted by sse_itr
        if( file_type == e_CCD)
        {
          write << "LOOP " << format_seq_id( sse_itr->First()) << " " << format_seq_id( sse_itr->Second())
                << " 0 " << " 0.0 " << " 0" << '\n';
        }
        else if( file_type == e_KIC)
        {
          if( sse_itr->First() == first_seq_id_chain || sse_itr->Second() == last_seq_id_chain)
          {
            BCL_MessageCrt
            (
              "loop at the C or N terminus. Rosetta KIC cannot build coordinates for those - use loop file at your own "
              "risk or consider using CCD!"
            )
            continue;
          }
          write << "LOOP " << format_seq_id( sse_itr->First()) << " "
                << format_seq_id( sse_itr->Second()) << " 0 " << " 0.0 " << " 1" << '\n';
        }
      }
    }

    //! @brief write protein parameters for the given sequence into a params file
    //! @param SEQUENCE the sequence to use to derive the protein parameters
    void PDBConvert::WriteProteinParamsFile( const biol::AASequence &SEQUENCE) const
    {
      // open file for write
      const std::string out_file_name( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + SEQUENCE.GetChainID() + ".params");
      io::OFStream ostream;

      // try opening
      if( io::File::TryOpenOFStream( ostream, out_file_name))
      {
        // calculate the protein parameters for the given sequence and write the result table to the open file
        biol::ProteinParams protein_params;
        const storage::Table< double> params( protein_params( SEQUENCE));
        params.WriteFormatted( ostream);
        io::File::CloseClearFStream( ostream);
      }
    }

    //! @brief write protein parameters for each sequence into a params file
    //! @param MODEL model with sequences to use to derive the protein parameters
    void PDBConvert::WriteProteinParamsFiles( const assemble::ProteinModel &MODEL) const
    {
      // iterate through all chains
      for
      (
        assemble::ProteinModel::const_iterator itr( MODEL.GetChains().Begin()), itr_end( MODEL.GetChains().End());
        itr != itr_end;
        ++itr
      )
      {
        WriteProteinParamsFile( *( *itr)->GetSequence());
      }
    }

    //! @brief determines helical transmembrane sequence stretches and merges broken segments
    //! @param MODEL the model whose helical tm segments will be merged into continuous tm helices if they are broken
    //! @return protein model with continuous transmembrane helices
    assemble::ProteinModel PDBConvert::MergeHelixBasedTMStretches( const assemble::ProteinModel &MODEL) const
    {
      // to hold the tm segments for each chain and the sses that make up each tm segment
      storage::Map< char, storage::Vector< util::SiPtrVector< const assemble::SSE> > > chain_tmsegments_helices;

      // iterate over the chains of the model to determine the transmembrane segments
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( MODEL.GetChains().Begin()), chain_itr_end( MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // reference to current chain
        const assemble::Chain &current_chain( **chain_itr);

        // current chain id
        const char chain_id( current_chain.GetChainID());

        // add in the current chain id
        chain_tmsegments_helices[ chain_id] = storage::Vector< util::SiPtrVector< const assemble::SSE> >();

        // reference to the current vector of tm segments for this chain
        storage::Vector< util::SiPtrVector< const assemble::SSE> > &tmsegments_helices
        (
          chain_tmsegments_helices[ chain_id]
        );

        // add initial sse vector to hold the sses for the first tm segment
        tmsegments_helices.PushBack( util::SiPtrVector< const assemble::SSE>());

        // keep track of slope of previous sse
        double previous_slope( util::GetUndefinedDouble());

        // get all sses
        util::SiPtrVector< const assemble::SSE> helix_sses( current_chain.GetSSEs());

        // iterate over the sses of the current chain
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_itr( helix_sses.Begin()), sse_itr_end( helix_sses.End());
            sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          const assemble::SSE &current_sse( **sse_itr);

          // calculate slope of this sse
          const int seq_id_change( current_sse.GetLastAA()->GetSeqID() - current_sse.GetFirstAA()->GetSeqID());
          const double first_aa_zcoord( current_sse.GetFirstAA()->GetCA().GetCoordinates().Z());
          const double last_aa_zcoord( current_sse.GetLastAA()->GetCA().GetCoordinates().Z());
          const double ca_z_coord_change( last_aa_zcoord - first_aa_zcoord);

          // minimum size an sse needs to be in order to be able to cause a new tm segment - helps to keep very short
          // sses from giving bad slopes and causing the current tm segment to be interrupted
          const int new_tm_min_sse_size( m_MergeBrokenTMHelixNewTMMinSSESizeParam->GetNumericalValue< int>());

          // calculate current slope of this sse (z coordinate change / seq id chain) i.e. how quickly the sse is
          // traversing the membrane. amphipathic helices will have small slope
          const double slope( seq_id_change > new_tm_min_sse_size ? ca_z_coord_change / seq_id_change : previous_slope);

          // true if the slope has not yet been defined for this chain
          if( !util::IsDefined( slope))
          {
            // go to next sse
            continue;
          }

          BCL_MessageStd
          (
            "slope of sse " + ( *sse_itr)->GetIdentification() + " is "
            + util::Format()( slope)
          );

          // true if not sloped enough, don't want amphipathics
          if( math::Absolute( slope) < math::Absolute( m_MergeBrokenTMHelixSlopeParam->GetNumericalValue< double>()))
          {
            BCL_MessageStd
            (
              "sse not sloped enough : " + current_sse.GetIdentification() +
              " has slope " + util::Format()( slope)
            );

            // start a new tm segment since this one is now ended
            tmsegments_helices.PushBack( util::SiPtrVector< const assemble::SSE>());

            // go to next sse
            continue;
          }

          // membrane thickness
          static const size_t membrane_half_width
          (
            m_MergeBrokenTMHelixHalfMembraneWidthParam->GetNumericalValue< double>()
          );

          // true if beginning or ending CA coordinates are not within membrane
          if
          (
            math::Absolute( first_aa_zcoord) > membrane_half_width &&
            math::Absolute( last_aa_zcoord) > membrane_half_width
          )
          {
            BCL_MessageStd
            (
              "sse not within membrane : " + current_sse.GetIdentification() +
              " first aa zcoord " + util::Format()( first_aa_zcoord) + " last_aa_zcoord "
              + util::Format()( last_aa_zcoord)
            );

            // start a new tm segment since this one is now ended
            tmsegments_helices.PushBack( util::SiPtrVector< const assemble::SSE>());

            // go to next sse
            continue;
          }

          // true if there are no tm segments yet for the current chain
          if( tmsegments_helices.IsEmpty() || tmsegments_helices.LastElement().IsEmpty())
          {
            BCL_MessageStd
            (
              "tmsegments or new segment empty. Previous slope was "
              + util::Format()( previous_slope) + " and is now being set to " + util::Format()( slope)
            );

            // set previous slope
            previous_slope = slope;
          }

          BCL_Assert
          (
            util::IsDefined( previous_slope), "previous slope is not defined.tmsegments_helices\n"
          );

          // true if slopes have different signs, i.e. new transmembrane segment
          if( slope / previous_slope < 0.0)
          {
            BCL_MessageStd
            (
              "adding new tmsegment since slope of " + util::Format()( slope) +
              " has different slope that previous slope of " + util::Format()( previous_slope)
            );

            // add vector for new tm segment
            tmsegments_helices.PushBack( util::SiPtrVector< const assemble::SSE>());
          }

          BCL_MessageStd( "adding sse " + ( *sse_itr)->GetIdentification());

          // add this sse to the current tm segment for this chain
          tmsegments_helices.LastElement().PushBack( *sse_itr);

          BCL_MessageStd
          (
            "setting previous slope of " + util::Format()( previous_slope) + " to " + util::Format()( slope)
          );

          // set slope
          previous_slope = slope;
        }
      }

      // copy the protein model
      assemble::ProteinModel new_model( *util::ShPtr< assemble::ProteinModel>( MODEL.HardCopy()));

      // combine any multi-sse tm segments into a single sse tm segment and replace them in the protein model
      for
      (
        storage::Map< char, storage::Vector< util::SiPtrVector< const assemble::SSE> > >::const_iterator
          chain_itr( chain_tmsegments_helices.Begin()), chain_itr_end( chain_tmsegments_helices.End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        // reference to current tm segments in the current chain
        const storage::Vector< util::SiPtrVector< const assemble::SSE> > &tm_segments( chain_itr->second);

        // to hold the new tm sses inserted into the protein model
        storage::Set< util::SiPtr< const assemble::SSE> > tm_sses;

        // iterate over the tm segments
        for
        (
          storage::Vector< util::SiPtrVector< const assemble::SSE> >::const_iterator
            tm_segments_itr( tm_segments.Begin()), tm_segments_itr_end( tm_segments.End());
          tm_segments_itr != tm_segments_itr_end;
          ++tm_segments_itr
        )
        {
          // sses making up the current tm segment
          const util::SiPtrVector< const assemble::SSE> &sses( *tm_segments_itr);

          // to keep track if the first sse of this segment has been found or not, don't want to start with anything
          // except for helix
          bool first_sse_not_found( true);

          // iterators to point to the beginning sse and ending sse of this tm segment
          util::SiPtrVector< const assemble::SSE>::const_iterator first_helix_itr, last_helix_itr( sses.End());

          // iterate over the sses of the current tm segment to find the first and last helices as start and ends
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( sses.Begin()), sse_itr_end( sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // find first helix as start of tm region - true if sse is type helix
            if( first_sse_not_found && ( *sse_itr)->GetType() == biol::GetSSTypes().HELIX)
            {
              first_helix_itr = sse_itr;
              first_sse_not_found = false;
              BCL_MessageStd
              (
                "first helix of tm region is " + ( *first_helix_itr)->GetIdentification()
              );
            }
            // find last helix as end of tm region
            if( ( *sse_itr)->GetType() == biol::GetSSTypes().HELIX)
            {
              last_helix_itr = sse_itr;
            }
          } //< end iteration through helices of this tm segment

          // true if no helices in tm segment
          if( last_helix_itr == sses.End())
          {
            // go to next tm segment
            continue;
          }

          BCL_MessageStd
          (
            "last helix of tm region is " + ( *last_helix_itr)->GetIdentification()
          );

          // make sure the tm segment is actually spanning the membrane as defined by simply checking if the
          // changes sign between the starting residue of the first helix to the ending residue of the last helix.
          // change in sign of z-coordinate indicates it is crossing over the center of the membrane
          const double z_coord_first( ( *first_helix_itr)->GetFirstAA()->GetCA().GetCoordinates().Z());
          const double z_coord_last( ( *last_helix_itr)->GetLastAA()->GetCA().GetCoordinates().Z());

          // true if the tm segment does not cross the center of the membrane
          if( z_coord_first / z_coord_last > 0.0)
          {
            BCL_MessageStd
            (
              "possible tm segment starting with helix "
              + ( *first_helix_itr)->GetIdentification() + " and ending with helix "
              + ( *last_helix_itr)->GetIdentification() +
              " does not cross the center of the membrane, so it is ignored. Start residue CA z coord is "
              + util::Format()( z_coord_first) + " last residue CA z coord is " + util::Format()( z_coord_last)
            );

            // go to next tm segment
            continue;
          }

          // to hold new aa sequence of combined sses for this tm segment
          biol::AASequence new_sequence;

          // iterate over the sses to fill the new sequence
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( first_helix_itr), sse_itr_end( last_helix_itr + 1);
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            new_sequence.AppendSequence( **sse_itr);
          }

          // make new sse and replace it in the protein model
          util::ShPtr< assemble::SSE> new_sse( new assemble::SSE( new_sequence, biol::GetSSTypes().HELIX));
          new_model.ReplaceWithOverlapping( new_sse);

          // add new sse to the set of continuous tm sses
          tm_sses.Insert( new_sse);
        }

        // if only want tm segment SSEs, remove other sses in chain
        if( m_MergeBrokenTMHelixRedefineNonTMParam->GetNumericalValue< size_t>())
        {
          // get all the SSEs
          util::SiPtrVector< const assemble::SSE> chain_sses( new_model.GetChain( chain_itr->first)->GetSSEs());

          // iterate over the sses
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_itr( chain_sses.Begin()), sse_itr_end( chain_sses.End());
            sse_itr != sse_itr_end;
            ++sse_itr
          )
          {
            // true if the current sse is a tm segment
            if( tm_sses.Contains( *sse_itr))
            {
              // do nothing, go to next sse
              continue;
            }

            // replace the sse with the same sse but as a coil in the model
            util::ShPtr< assemble::SSE> new_sse( new assemble::SSE( **sse_itr, biol::GetSSTypes().COIL));
            new_model.ReplaceWithOverlapping( new_sse);
          }
        }
      }

      // return new model with continous tm helices
      return new_model;
    }

    //! @brief process the model according to the command line flags
    //! @param PROTEIN_MODEL protein model to be processed
    //! @param IDENTIFIER optional string identifier to be used with output files
    void PDBConvert::ProcessModel( assemble::ProteinModel &PROTEIN_MODEL, const std::string &IDENTIFIER) const
    {
      // instantiate the pdb factory
      pdb::Factory factory;

      // process model with dssp
      if( m_DSSPFlag->GetFlag())
      {
        biol::DSSP model_dssp;
        math::MutateResult< assemble::ProteinModel> new_model( model_dssp( PROTEIN_MODEL));
        if( new_model.GetArgument().IsDefined())
        {
          io::OFStream write;
          io::File::MustOpenOFStream( write, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".dssp");
          model_dssp.WriteToFile( write);
          io::File::CloseClearFStream( write);
          PROTEIN_MODEL = *new_model.GetArgument();
        }
      }

      // process model with mahssmi
      if( m_MahssmiFlag->GetFlag() || m_CIPhiPsiFlag->GetFlag())
      {
        if( m_MahssmiFlag->GetFlag())
        {
          sspred::MethodHandler::ReadAllPredictionsForProteinModel
          (
            PROTEIN_MODEL,
            m_OutputPrefixFlag->GetFirstParameter()->GetValue(),
            ""
          );
        }
        sspred::PDB::SetEnvironmentTypes( PROTEIN_MODEL, true);
        io::OFStream write;
        if( m_MahssmiFlag->GetFlag())
        {
          // compute and write out mahssmi analysis
          sspred::Mahssmi::Calculate( PROTEIN_MODEL, true);
          io::File::MustOpenOFStream( write, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".mahssmi");
          sspred::Mahssmi::WriteAnalysis( write, PROTEIN_MODEL);
          io::File::CloseClearFStream( write);
        }
        if( m_CIPhiPsiFlag->GetFlag())
        {
          // compute and write out context-independent phi/psi analysis
          sspred::CIPhiPsi().Calculate( PROTEIN_MODEL, true);
          io::File::MustOpenOFStream( write, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".ciphipsi");
          sspred::CIPhiPsi().WriteAnalysis( write, PROTEIN_MODEL);
          io::File::CloseClearFStream( write);
        }
      }

      // true if broken tm helices should be merged to form a continous tm helix
      if( m_MergeBrokenTMHelixFlag->GetFlag())
      {
        // update model with new helices
        PROTEIN_MODEL = MergeHelixBasedTMStretches( PROTEIN_MODEL);
      }

      const std::string model_chain_ids( PROTEIN_MODEL.GetChainIDs());
      io::IFStream read;

      // write protein params if desired
      if( m_ProteinParamsFlag->GetFlag())
      {
        WriteProteinParamsFiles( PROTEIN_MODEL);
      }

      // generate multimer
      if( pdb::Factory::GetFlagBiomolecule()->GetFlag())
      {
        // initialize multiplier
        util::ShPtr< assemble::ProteinModelMultiplier> sp_multiplier;

        if( m_PdbTmXmlFlag->GetFirstParameter()->GetWasSetInCommandLine())
        {
          // get membrane and matrix from xml file
          io::File::MustOpenIFStream( read, m_PdbTmXmlFlag->GetFirstParameter()->GetValue());
          // chain ids with matrices that are applied to those chains to generate parts of the multimer
          storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > bio_transformation_matrices
          (
            biol::Membrane::BioTransformationMatricesFromPDBTMXML( read, model_chain_ids)
          );
          io::File::CloseClearFStream( read);

          if( bio_transformation_matrices.IsEmpty())
          {
            BCL_MessageStd( "no biomatrices in xml given, using pdb REMARK 350");
          }
          else
          {
            BCL_MessageStd( "biomatrices in xml given");
            sp_multiplier = util::ShPtr< assemble::ProteinModelMultiplier>
            (
              new assemble::ProteinModelMultiplier( bio_transformation_matrices, PROTEIN_MODEL)
            );
          }
        }

        // either no pdbtmxml given, or no biomatrix given in xml
        if( !sp_multiplier.IsDefined())
        {
          // see if the multiplier was set by the pdb factory
          sp_multiplier = util::ShPtr< assemble::ProteinModelMultiplier>
          (
            PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Multiplier)
          );
        }

        // if there are matrices given somehow
        if( sp_multiplier.IsDefined())
        {
          // build the multimer
          PROTEIN_MODEL = sp_multiplier->operator ()( PROTEIN_MODEL);

          // write multiplier information
          BCL_MessageStd( "mapping of chain ids after applying biomolecule multiplier:");
          if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Standard))
          {
            sp_multiplier->ChainIDMapping().WriteFormatted( util::GetLogger());
          }
        }
        else
        {
          BCL_MessageCrt
          (
            "could not generate bio multimer, since no matrices were given"
          );
        }

        if( m_MultimerNativePDBFileFlag->GetFlag())
        {
          const std::string native_multimer_filename( m_MultimerNativePDBFileFlag->GetFirstParameter()->GetValue());
          pdb::Factory factory;
          const assemble::ProteinModel native_multimer
          (
            factory.ProteinModelFromPDBFilename( native_multimer_filename)
          );
          const quality::Measure quality
          (
            quality::GetMeasures().GetEnumFromName
            (
              m_SuperimposeQualityMeasureParam->GetValue()
            )
          );
          BCL_MessageStd( "determining best multimer chain assignment");
          PROTEIN_MODEL = assemble::PrinterProteinModelMultimer::CalculateBestMultimer
          (
            PROTEIN_MODEL,
            native_multimer,
            quality,
            sp_multiplier
          );
        }
      }

      // transform according to moment of inertia
      if( m_MomentOfInertiaFlag->GetFlag())
      {
        const assemble::ProteinModelMomentOfInertia moment_of_inertia_calculator
        (
          biol::AATypeData::PropertyTypeEnum( m_MomentOfInertiaPropertyParam->GetValue()),
          *assemble::AAExposure( m_MomentOfInertiaExposureParam->GetValue()),
          false
        );

        const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> transformation_moments
        (
          moment_of_inertia_calculator.TransformationAndMoments( PROTEIN_MODEL)
        );
        BCL_MessageStd( "transformation:\t" + util::Format()( transformation_moments.First()));
        BCL_MessageStd( "moments:\t" + util::Format()( transformation_moments.Second()));

        // transform the model
        PROTEIN_MODEL.Transform( transformation_moments.First());
      }
      // transform according to given pdbtm xml
      else if( m_PdbTmXmlFlag->GetFirstParameter()->GetWasSetInCommandLine())
      {
        // get membrane and matrix from xml file
        io::File::MustOpenIFStream( read, m_PdbTmXmlFlag->GetFirstParameter()->GetValue());
        const storage::Pair< biol::Membrane, math::TransformationMatrix3D>
          membrane_matrix( biol::Membrane::MembraneAndTransformationFromPDBTMXML( read, double( 10), double( 2.5)));
        io::File::CloseClearFStream( read);

        // transform protein according to given matrix
        PROTEIN_MODEL.Transform( membrane_matrix.Second());

        // insert the membrane
        util::ShPtr< assemble::ProteinModelData> sp_data( PROTEIN_MODEL.GetProteinModelData());
        sp_data->Insert( assemble::ProteinModelData::e_Membrane, util::CloneToShPtr( membrane_matrix.First()));
        PROTEIN_MODEL.SetProteinModelData( sp_data);

        // include membrane printer so membrane info is added to pdb output
        factory.AppendPrinter( util::ShPtr< pdb::PrinterMembrane>( new pdb::PrinterMembrane()));
      }

      if( m_IdealizeFlag->GetFlag())
      {
        PROTEIN_MODEL.SetToIdealConformation();
      }

      // check if chains should be superimposed and written as biomolecule entries
      if( m_SuperimposeToBiomoleculeFlag->GetFlag())
      {
        const assemble::Biomolecule biomolecule( biol::GetAtomTypes().GetBackBoneAtomTypes(), *quality::SuperimposeMeasure( m_SuperimposeToBiomoleculeParam->GetValue()), m_SuperimposeToBiomoleculeThresholdParam->GetNumericalValue< double>());
        math::MutateResult< assemble::ProteinModel> result( biomolecule( PROTEIN_MODEL));
        if( result.GetArgument().IsDefined())
        {
          PROTEIN_MODEL = *result.GetArgument();
        }
      }

      // change chain id of empty identifier to given chain id
      const storage::Vector< std::string> chain_id_pairs( m_RenameChainIDsFlag->GetObjectList< std::string>());
      for
      (
        storage::Vector< std::string>::const_iterator
          chain_pair_itr( chain_id_pairs.Begin()), chain_pair_itr_end( chain_id_pairs.End());
        chain_pair_itr != chain_pair_itr_end;
        ++chain_pair_itr
      )
      {
        std::string chain_pair( *chain_pair_itr);
        if( std::count( chain_pair.begin(), chain_pair.end(), s_ChainIDPairSeparator) != 1)
        {
          BCL_MessageCrt
          (
            "cannot use the given chain pair for renaming: " + chain_pair
          );
          continue;
        }
        if( chain_pair.size() > 3)
        {
          BCL_MessageCrt
          (
            "too many characters in chain pair for renaming: " + chain_pair
          );
          continue;
        }

        // replace separator with space
        chain_pair[ chain_pair.find( s_ChainIDPairSeparator)] = ' ';

        // original and new chain id
        const char original_chain_id( chain_pair[ 0]);
        const char new_chain_id( chain_pair[ chain_pair.size() - 1]);

        if( original_chain_id == new_chain_id)
        {
          continue;
        }

        // assert that chain with new id does not exist already
        BCL_Assert
        (
          !PROTEIN_MODEL.GetChain( new_chain_id).IsDefined(),
          "chain with new chain id already exists for pair: " + ( *chain_pair_itr)
        );

        // original chain
        util::ShPtr< assemble::Chain> original_chain( PROTEIN_MODEL.GetChain( original_chain_id));
        BCL_Assert
        (
          original_chain.IsDefined(),
          "chain with original chain id does not exist for pair: " + ( *chain_pair_itr)
        );

        // rename
        original_chain->SetChainID( new_chain_id);
      }

      // list of chain ids of interest
      std::string chain_ids;

      // if specific chains are provided
      if( m_ChainsFlag->GetFlag())
      {
        // iterate over provided chains
        for
        (
          util::ShPtrVector< command::ParameterInterface>::const_iterator param_itr( m_ChainsFlag->GetParameterList().Begin()),
            param_itr_end( m_ChainsFlag->GetParameterList().End());
          param_itr != param_itr_end;
          ++param_itr
        )
        {
          // push back the first letter
          chain_ids.push_back( ( *param_itr)->GetValue()[ 0]);
        }
      }
      // else if no specific chains were given, get all the chains from the pdb
      else
      {
        chain_ids = PROTEIN_MODEL.GetChainIDs();
      }

      // if fasta flag is given
      if( m_WriteFastaFlag->GetFlag())
      {
        // iterate over chain_ids
        for
        (
          std::string::const_iterator chain_id_itr( chain_ids.begin()), chain_id_itr_end( chain_ids.end());
          chain_id_itr != chain_id_itr_end;
          ++chain_id_itr
        )
        {
          // get the chain id, use underscore if equal to space
          const char chain_id_for_protein_model
          (
            ( *chain_id_itr)
          );

          // output filename
          const std::string fasta_filename
          (
            m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
            ( chain_id_for_protein_model == ' ' ? '_' : chain_id_for_protein_model) + ".fasta"
          );

          // open output stream
          io::OFStream write;
          io::File::MustOpenOFStream( write, fasta_filename);

          // output the fasta to the opened ofstream
          PROTEIN_MODEL.GetChain( chain_id_for_protein_model)->GetSequence()->WriteFasta( write);

          io::File::CloseClearFStream( write);
        }
      }

      // graphvis file
      if( m_TopologyGraphFlag->GetFlag())
      {
        const assemble::Topology protein_topology
        (
          assemble::CollectorTopologyCombined().CalculateTopology( PROTEIN_MODEL.GetSSEs())
        );

        io::OFStream write;
        // output filename
        const std::string graphviz_filename
        (
          m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".viz"
        );

        // scoring function used for edge coloring
        math::SumFunction< assemble::SSEGeometryPacking, double> score_sum;
        score_sum.NewOperand( score::SSEPairPacking());
        score_sum.NewOperand( score::StrandPairing());

        io::File::MustOpenOFStream( write, graphviz_filename);
        protein_topology.WriteGraphVizScript( write, score_sum);
        io::File::CloseClearFStream( write);

        BCL_MessageStd
        (
          "create topology graph with the following commandline:\n"
          "neato -Tpng " + graphviz_filename + " -o " + m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".png"
        );
      }

      // true if sdfs should be written
      if( m_WriteSDFFlag->GetFlag())
      {
        // iterate over chain_ids
        for
        (
          std::string::const_iterator chain_id_itr( chain_ids.begin()), chain_id_itr_end( chain_ids.end());
          chain_id_itr != chain_id_itr_end;
          ++chain_id_itr
        )
        {
          // get the chain id, use underscore if equal to space
          const char chain_id_for_protein_model
          (
            ( *chain_id_itr)
          );

          // create a molecule
          util::ShPtr< biol::AASequence> sequence( PROTEIN_MODEL.GetChain( chain_id_for_protein_model)->GetSequence());
          chemistry::AAFragmentComplete aa_fragment( sequence->GetMembers(), true);

          // check if molecule is too large
          if( aa_fragment.GetNumberAtoms() > 999)
          {
            BCL_MessageCrt
            (
              "Could not write chain " + std::string( size_t( 1), chain_id_for_protein_model)
              + " as an sdf because it contains more atoms (" + util::Format()( aa_fragment.GetNumberAtoms()) + ") than"
              " the MDL/SDF format support (999)"
            );
            continue;
          }

          // output filename
          const std::string sdf_filename
          (
            m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
            ( chain_id_for_protein_model == ' ' ? '_' : chain_id_for_protein_model) + ".sdf"
          );

          // open output stream
          io::OFStream write;
          io::File::MustOpenOFStream( write, sdf_filename);

          // output the sdf to the opened ofstream
          aa_fragment.WriteMDL( write);

          io::File::CloseClearFStream( write);
        }
      }

      // when SASA is requested
      if( m_WriteSASAFlag->GetFlag())
      {
        // iterate over all chains
        for
        (
          std::string::const_iterator chain_id_itr( chain_ids.begin()), chain_id_itr_end( chain_ids.end());
          chain_id_itr != chain_id_itr_end;
          ++chain_id_itr
        )
        {
          // get the chain id, use underscore if equal to space
          const char chain_id_for_protein_model
          (
            ( *chain_id_itr)
          );

          // TODO
          // check if each chain has coordinates
          size_t num_res_with_coordinates( 0);
          for
          (
              util::ShPtrVector< biol::AABase>::const_iterator
                aa_itr( PROTEIN_MODEL.GetChain( chain_id_for_protein_model)->GetSequence()->Begin()),
                aa_itr_end( PROTEIN_MODEL.GetChain( chain_id_for_protein_model)->GetSequence()->End());
              aa_itr != aa_itr_end;
              ++aa_itr
          )
          {
            if( ( *aa_itr)->HasDefinedCoordinates())
            {
              ++num_res_with_coordinates;
            }
          }

          if( num_res_with_coordinates == 0)
          {
            BCL_MessageStd
            (
              m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
                " chain " + chain_id_for_protein_model +
                " does not have defined coordinates."
            );
          }

          // output filename
          const std::string ncnv_filename
          (
            m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
            ( chain_id_for_protein_model == ' ' ? '_' : chain_id_for_protein_model) + ".ncnv"
          );

          // open output stream
          io::OFStream write;
          io::File::MustOpenOFStream( write, ncnv_filename);

          descriptor::AASasa::WriteNCNV( m_MinimalSequenceSeparation->GetNumericalValue< size_t>(), write, chain_id_for_protein_model, PROTEIN_MODEL);

          io::File::CloseClearFStream( write);
        }
      }

      // when written pdb is requested
      if( m_BclPdbFlag->GetFlag())
      {
        const PdbWriteFilesTypeEnum write_modus( m_BclPdbFlag->GetFirstParameter()->GetValue());

        // true if loops should be added to the protein or zero coordinates are requested
        if( m_RosettaLoopFileFlag->GetFlag() || pdb::Factory::GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()->GetFlag())
        {
          // add loops to "single_chain_model"
          PROTEIN_MODEL.AddLoops( true, false);

          // print warning if modus is full and loop file is requested
          if( m_RosettaLoopFileFlag->GetFlag() && write_modus == e_Full)
          {
            BCL_MessageCrt
            (
              "Rosetta loop files not available for multi-chain models, use \"Split\" or \"Both\" with -bcl_pdb flag"
            );
          }
        }

        // only write single chain pdbs if desired
        if( ( write_modus == e_Split || write_modus == e_Both))
        {
          // iterate over chain_ids
          for
          (
            std::string::const_iterator chain_id_itr( chain_ids.begin()), chain_id_itr_end( chain_ids.end());
            chain_id_itr != chain_id_itr_end;
            ++chain_id_itr
          )
          {
            // get the chain id, use underscore if equal to space
            const char chain_id( ( *chain_id_itr));

            // output filename
            const std::string pdb_filename
            (
              m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
                ( chain_id == ' ' ? '_' : chain_id) + IDENTIFIER + ".pdb"
            );

            // open output stream
            io::OFStream write;
            io::File::MustOpenOFStream( write, pdb_filename);

            // output the pdb to the opened ofstream
            assemble::ProteinModel single_chain_model( PROTEIN_MODEL.GetChain( chain_id));

            // true if loops file should be written
            if( m_RosettaLoopFileFlag->GetFlag())
            {
              BCL_Assert
              (
                !pdb::Factory::GetFlagWritePDBResID()->GetFlag(),
                "flag -" + pdb::Factory::GetFlagWritePDBResID()->GetName() + " was selected together with rosetta loop "
                "file flag. Written rosetta loop file will not make sense if the pdb is not renumbered!"
              );
              // write loop file
              WriteLoopsFile( *PROTEIN_MODEL.GetChain( chain_id));
            }

            factory.WriteModelToPDB( single_chain_model, write, m_BodyOutputFlag->GetFlag());

            io::File::CloseClearFStream( write);
          }
        }

        if( write_modus == e_Full || write_modus == e_Both)
        {
          // if AddSideChains flag was given
          if( m_AddSideChainsFlag->GetFlag())
          {
            // Protein Model should contain at least AABackBones (default flag) for the SideChainFactory to work properly
            BCL_Assert
            (
                biol::AAClass( pdb::Factory::GetFlagAAClass()->GetFirstParameter()->GetValue())
              >= biol::GetAAClasses().e_AABackBone,
              "Protein model is not instantiated with at least AABackbone"
            );

            //instantiate AASideChainFactory no hydrogens include backbone atoms for sidechain placement superimposition
            biol::AASideChainFactory side_chain_factory( false, true);

            // add side chains to model and set it to model
            PROTEIN_MODEL = *side_chain_factory.ProteinModelWithSideChains( PROTEIN_MODEL);
          }

          // if the superimposition flag is set
          if( m_SuperimposePDBFlag->GetFlag())
          {
            // read the native pdb in
            io::File::MustOpenIFStream( read, m_SuperimposePDBFlag->GetFirstParameter()->GetValue());
            pdb::Handler native_pdb( read);
            io::File::CloseClearFStream( read);

            // instantiate protein model
            const assemble::ProteinModel native_model( factory.ProteinModelFromPDB( native_pdb));

            // initialize the atom types
            // iterate over dynamic parameters
            const util::SiPtrVector< const command::ParameterInterface> dynamic_parameters_atom_types( m_SuperimposePDBFlag->GetDynamicParameterList());
            storage::Set< biol::AtomType> atom_types;
            for
            (
              util::SiPtrVector< const command::ParameterInterface>::const_iterator
                itr( dynamic_parameters_atom_types.Begin()), itr_end( dynamic_parameters_atom_types.End());
              itr != itr_end;
              ++itr
            )
            {
              atom_types.Insert( biol::AtomType( ( *itr)->GetValue()));
            }

            // get the enum from name
            const quality::SuperimposeMeasure measure
            (
              quality::GetSuperimposeMeasures().GetEnumFromName
              (
                m_SuperimposeQualityMeasureParam->GetValue()
              )
            );

            // transform the model
            assemble::Quality::SuperimposeModel( measure, PROTEIN_MODEL, native_model, atom_types);
          } // if superimpose and full or both

          // if transformation flag is given
          if( m_TransformationFlag->GetFlag())
          {
            // read in the transformation matrix
            io::IFStream read;
            io::File::MustOpenIFStream( read, m_TransformationFlag->GetFirstParameter()->GetValue());
            math::TransformationMatrix3D transform;
            read >> transform;
            io::File::CloseClearFStream( read);

            // apply the transformation
            PROTEIN_MODEL.Transform( transform);
          }

          // write complete bcl pdb
          // output filename
          const std::string pdb_filename( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + IDENTIFIER + "bcl.pdb");

          // open output stream
          io::OFStream write;
          io::File::MustOpenOFStream( write, pdb_filename);

          // output the pdb to the opened ofstream
          factory.WriteModelToPDB( PROTEIN_MODEL, write, m_BodyOutputFlag->GetFlag());

          io::File::CloseClearFStream( write);
        }
      } // if write full of both bcl pdb

      io::File::CloseClearFStream( read);
    }

    //! default constructor
    PDBConvert::PDBConvert() :
      m_PdbFilename
      (
        new command::Parameter
        (
          "pdb_filename",
          "filename for input pdb",
          ""
        )
      ),
      m_WriteFastaFlag
      (
        new command::FlagStatic
        (
          "fasta",
          "write fasta file"
        )
      ),
      m_PdbFromFastaFlag
      (
        new command::FlagDynamic
        (
          "pdb_from_fasta",
          "list of fasta files of which the seqres will be written to a file",
          command::Parameter
          (
            "fasta_file",
            "name of fasta file with chain extension like: 1ubiA.fasta",
            command::ParameterCheckExtension( ".fasta"),
            ""
          ),
          0, 100
        )
      ),
      m_ChainsFlag
      (
        new command::FlagDynamic
        (
          "chains",
          "list of chains to be used for output",
          command::Parameter
          (
            "chain_id",
            "one letter chain identifier"
          ),
          0,
          util::GetUndefined< size_t>()
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "prefix for output to be used if any output flag is selected",
          command::Parameter
          (
            "prefix",
            "output prefix",
            "out"
          )
        )
      ),
      m_BclPdbFlag
      (
        new command::FlagStatic
        (
          "bcl_pdb",
          "write bcl pdb file",
          command::Parameter
          (
            "write_modus",
            "select if the pdb should be split up (write individual files for each chain), if it should be remain combined as one pdb, or if both should be written",
            command::ParameterCheckSerializable( PdbWriteFilesTypeEnum()),
            GetPdbWriteFilesTypeDescriptor( e_Full)
          )
        )
      ),
      m_DSSPFlag
      (
        new command::FlagStatic
        (
          "dssp",
          "reassign secondary structure based on DSSP and write dssp file"
        )
      ),
      m_BodyOutputFlag
      (
        new command::FlagStatic
        (
          "body",
          "write body information to bcl pdb file"
        )
      ),
      m_IdealizeFlag
      (
        new command::FlagStatic
        (
          "idealize",
          "idealize bcl pdb file"
        )
      ),
      m_PdbTmXmlFlag
      (
        new command::FlagStatic
        (
          "pdbtm_xml",
          "transform according to the given transformation matrix given in file and if available, use BIOMATRIX instead of the one given in PDB to generate multimer",
          command::Parameter
          (
            "pdbtm_xml_file_name",
            "filename of pdbtm xml file containing membrane normal and transformationmatrix",
            ""
          )
        )
      ),
      m_MomentOfInertiaFlag
      (
        new command::FlagStatic
        (
          "align_moment_of_inertia",
          "align the protein model with its moment of interia -> largest align with z-axis, smallest with x-axis"
        )
      ),
      m_MomentOfInertiaPropertyParam
      (
        new command::Parameter
        (
          "property", "an amino acid property to use for the moment of inertia",
          command::ParameterCheckAllowed
          (
            storage::Vector< std::string>
            (
              biol::AATypeData::s_NumberPropertyTypes - 1,
              &biol::AATypeData::GetPropertyDescriptor( biol::AATypeData::e_NaturalPrevalence)
            )
          ),
          biol::AATypeData::GetPropertyDescriptor( biol::AATypeData::e_FreeEnergyCore)
        )
      ),
      m_MomentOfInertiaExposureParam
      (
        new command::Parameter
        (
          "amino acid exposure",
          "multiply the property with the exposure - more exposed residues will have higher contribution",
          command::ParameterCheckEnumerate< assemble::AAExposures>(),
          "none"
        )
      ),
      m_RosettaLoopFileFlag
      (
        new command::FlagStatic
        (
          "loop_file_rosetta",
          "write coil regions to rosetta-format loop-file; if generated pdb file is used as input to rosetta loop "
          "building, also use flag: -" + pdb::Factory::GetFlagWriteZeroCoordinatesForUndefinedAminoAcids()->GetName(),
          command::Parameter
          (
            "loop_file_type",
            "choice of rosetta loop file type depending on used protocol",
            command::ParameterCheckSerializable( RosettaLoopFileTypeEnum()),
            GetRosettaLoopFileTypeDescriptor( e_CCD)
          )
        )
      ),
      m_RenameChainIDsFlag
      (
        new command::FlagDynamic
        (
          "rename_chain_id",
          "change the first given chain id to be the second given chain id - unlimited number of pair possible\n"
          "  -rename_chain_id A,B  ,D C, will rename chain A to B, \' \' to D and C to \' \'",
          command::Parameter
          (
            "chain_id_pair",
            "original,new",
            std::string( "A") + s_ChainIDPairSeparator + "B"
          ),
          0, 26
        )
      ),
      m_AddSideChainsFlag
      (
        new command::FlagStatic
        (
          "side_chains",
          "add side chains to pdb file"
        )
      ),
      m_SuperimposePDBFlag
      (
        new command::FlagStaticAndDynamic
        (
          "superimpose",
          "superimpose onto a given pdb file, using the given quality measure and defined atom types",
          command::Parameter
          (
            "superimpose_atom_type",
            "backbone atom type to use for superimposition",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>( biol::GetAtomTypes().GetBackBoneAtomNames())
            )
          ), 0, biol::GetAtomTypes().GetBackBoneAtomNames().GetSize()
        )
      ),
      m_SuperimposPDBFileParam
      (
        new command::Parameter
        (
          "superimpose_pdb_filename",
          "parameter for the filename of the pdbfile onto which given pdbs will be superimposed",
          command::ParameterCheckExtension( ".pdb"),
          ""
        )
      ),
      m_SuperimposeQualityMeasureParam
      (
        new command::Parameter
        (
          "superimpose_measure",
          "\tquality measure to use for superimposition",
          command::ParameterCheckEnumerate< quality::SuperimposeMeasures>(),
          quality::GetSuperimposeMeasures().e_RMSD.GetName()
        )
      ),
      m_TransformationFlag
      (
        new command::FlagStatic
        (
          "transform",
          "apply given transformation",
          command::Parameter( "transformation_file", "file containing BCL transformation matrix", "")
        )
      ),
      m_TopologyGraphFlag
      (
        new command::FlagStatic
        (
          "topology",
          "write topology graph file for graphviz(neato)"
        )
      ),
      m_MahssmiFlag
      (
        new command::FlagStatic
        (
          "mahssmi",
          "Use the membrane-aware hybrid secondary structure and membrane topology identification algorithm and write "
          "out .mahssmi files for each chain.  This option requires that StrideDSSP and PALSSE files are available. "
          "and in the same directory as -output_prefix.  For membrane proteins, a pdbtm .xml file should also be "
          "in the same location"
        )
      ),
      m_CIPhiPsiFlag
      (
        new command::FlagStatic
        (
          "ciphipsi",
          "Use the context-insensitive, phi-psi based secondary structure and membrane topology identification algorithm"
          "and write out .ciphipsi files for each chain. For membrane proteins, a pdbtm .xml file should also be "
          "in the same location"
        )
      ),
      m_MultimerNativePDBFileFlag
      (
        new command::FlagStatic
        (
          "multimer_native",
          "multimer pdb filename. If a multimer is being generated, pass a native file here so that the generated "
          "multimer is a close as possible to the native according to the superimpose measure provided",
          command::Parameter( "pdb_filename", "pdb filename", "multimer_native.pdb")
        )
      ),
      m_ProteinParamsFlag
      (
        new command::FlagStatic
        (
          "protein_params",
          "write protein parameters like mol weight, extinction coefficient, AA counts and other metrices for each chain in a \".params\" file"
        )
      ),
      m_SuperimposeToBiomoleculeFlag
      (
        new command::FlagStatic
        (
          "superimpose_chain_biomolecule",
          "generate a biomolecule entry REMARK 350 by testing superimpositions of chains"
        )
      ),
      m_SuperimposeToBiomoleculeParam
      (
        new command::Parameter
        (
          "superimpose_measure",
          "\tquality measure to use for superimposition",
          command::ParameterCheckEnumerate< quality::SuperimposeMeasures>(),
          quality::GetSuperimposeMeasures().e_RMSD.GetName()
        )
      ),
      m_SuperimposeToBiomoleculeThresholdParam
      (
        new command::Parameter
        (
          "superimpose_thresold",
          "\tthreshold for superimposition measure to consider two chains as related",
          command::ParameterCheckRanged< double>( 0.0, 1000.0),
          ""
        )
      ),
      m_MergeBrokenTMHelixFlag
      (
        new command::FlagStatic
        (
          "merge_broken_tm_helices",
          "Merge transmembrane helices that are broken into a single transmembrane helix"
        )
      ),
      m_MergeBrokenTMHelixSlopeParam
      (
        new command::Parameter
        (
          "zcoord_change_per_residue",
          "the slope at which a helix must be traversing the membrane to be considered transmembrane"
          " (i.e. to avoid amphipathic helices)",
          "0.75"
        )
      ),
      m_MergeBrokenTMHelixHalfMembraneWidthParam
      (
        new command::Parameter
        (
          "membrane_half_width",
          "Maximum distance a helix can be from the center of the membrane and be considered transmembrane. Either the"
          " first or last residue CA z-coordinate has to be within this value (absolute values are used)",
          "25"
        )
      ),
      m_MergeBrokenTMHelixNewTMMinSSESizeParam
      (
        new command::Parameter
        (
          "new_tm_min_sse_size",
          "The minimum size ( = amount of sequence difference = (#resi - 1)) of an sse for it to be possible for it"
          " to start a new tm region. Useful for avoid one or two residues in from breaking a tm stretch.",
          command::ParameterCheckRanged< int>( 1, util::GetUndefined< int>()),
          "1"
        )
      ),
      m_MergeBrokenTMHelixRedefineNonTMParam
      (
        new command::Parameter
        (
          "redefine_non_tm",
          "If set to true (1), only the tm regions will be defined secondary structure elements after merging has been "
          "done. This gives the protein model with only transmembrane regions defined as SSEs. Must be 0 for false or "
          "1 for true",
          command::ParameterCheckRanged< size_t>( 0, 1),
          "0"
        )
      ),
      m_SplitEnsembleFlag
      (
        new command::FlagStatic
        (
          "split_ensemble",
          "split an ensemble (several models in one PDB file) into separate files"
        )
      ),
      m_WriteSDFFlag
      (
        new command::FlagStatic
        (
          "sdf",
          "Writes model(s) out as (MDL) sdf entries in addition to the usual PDB entries. "
          "Input files must contain coordinates for all non-H atoms if this option is used. "
          "The MDL format is limited to molecules with 999 or fewer atoms; larger models will not be written out"
        )
      ),
      m_WriteSASAFlag
      (
        new command::FlagStatic
        (
          "sasa",
          "Write neighbor counts as well as neighbor vectors into a single .ncnv file."
        )
      ),
      m_MinimalSequenceSeparation
      (
        new command::Parameter
        (
          "min_seq_separation",
          "minimal sequence separation between two residues for them to consider as spatial neighbors",
          command::ParameterCheckRanged< size_t>( 0, 50),
          "0"
        )
      )
    {
      m_SuperimposePDBFlag->PushBack( m_SuperimposPDBFileParam);
      m_SuperimposePDBFlag->PushBack( m_SuperimposeQualityMeasureParam);

      m_MomentOfInertiaFlag->PushBack( m_MomentOfInertiaPropertyParam);
      m_MomentOfInertiaFlag->PushBack( m_MomentOfInertiaExposureParam);

      m_MergeBrokenTMHelixFlag->PushBack( m_MergeBrokenTMHelixSlopeParam);
      m_MergeBrokenTMHelixFlag->PushBack( m_MergeBrokenTMHelixHalfMembraneWidthParam);
      m_MergeBrokenTMHelixFlag->PushBack( m_MergeBrokenTMHelixNewTMMinSSESizeParam);
      m_MergeBrokenTMHelixFlag->PushBack( m_MergeBrokenTMHelixRedefineNonTMParam);

      m_SuperimposeToBiomoleculeFlag->PushBack( m_SuperimposeToBiomoleculeParam);
      m_SuperimposeToBiomoleculeFlag->PushBack( m_SuperimposeToBiomoleculeThresholdParam);

      m_WriteSASAFlag->PushBack( m_MinimalSequenceSeparation);
    }

    const char PDBConvert::s_ChainIDPairSeparator = ',';

    const ApplicationType PDBConvert::PDBConvert_Instance
    (
      GetAppGroups().AddAppToGroup( new PDBConvert(), GetAppGroups().e_Protein)
    );

  } // namespace app
} // namespace bcl
