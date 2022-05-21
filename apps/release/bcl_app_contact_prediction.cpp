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
#include "bcl_app_contact_prediction.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "contact/bcl_contact_map.h"
#include "contact/bcl_contact_prediction_map.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "sspred/bcl_sspred_jufo.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

    //! @brief return the bcl::commons name
    //! @return string for the bcl::commons name of that application
    std::string ContactPrediction::GetBCLScopedName() const
    {
      static const std::string s_bcl_commons_name( "BCL::Contact");
      return s_bcl_commons_name;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ContactPrediction::GetDescription() const
    {
      return "Prediction of amino acid contacts for a given sequence";
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ContactPrediction::Main() const
    {
      //initialize write and read stream objects
      io::IFStream read;
      io::OFStream write;

      // if a pdb list is not provided then assume it is a fasta file
      if( !m_PDBListFlag->GetFlag())
      {
        // get the fasta filename
        const std::string fasta_filename( m_InputFilenameParam->GetValue());

        // create a sequence
        io::File::MustOpenIFStream( read, fasta_filename);
        util::ShPtr< biol::AASequence> sp_sequence
        (
          new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read))
        );
        io::File::CloseClearFStream( read);

        // read the blast profile
        storage::VectorND< 2, std::string> path_and_tag( pdb::Handler::ExtractPathAndPDBTag( fasta_filename));
        const std::string blast_filename( path_and_tag.First() + PATH_SEPARATOR + path_and_tag.Second() + ".ascii");

        io::File::MustOpenIFStream( read, blast_filename);
        biol::BlastProfileHandler::ReadProfileForAASequence( read, *sp_sequence);
        io::File::CloseClearFStream( read);

        // calculate jufo for the given sequence
        sspred::JUFO::Calculate( *sp_sequence);

        // create a chain from this sequence
        assemble::Chain chain( sp_sequence);

        // create prediction map
        contact::PredictionMap prediction_map( chain);

        // prepare the output file name
        const std::string prediction_map_filename
        (
          m_OutputPathFlag->GetFirstParameter()->GetValue() + PATH_SEPARATOR + path_and_tag.Second() + ".contact"
        );

        // initialize write stream with the filename
        io::File::MustOpenOFStream( write, prediction_map_filename);

        // output the map to file
        prediction_map.WritePredictionMap( write, m_ThresholdFlag->GetFirstParameter()->GetNumericalValue< double>());

        // clean the write stream
        io::File::CloseClearFStream( write);
      }

      // if pdb flag is set
      else
      {
        // read the list of pdbs from file
        BCL_MessageStd( "read list of pdbs from file");
        io::File::MustOpenIFStream( read, m_InputFilenameParam->GetValue());

        //initialize vector of string for file names
        storage::Vector< storage::VectorND< 2, std::string> > pdb_tags;

        std::string full_path;

        // iterate over names
        while( !read.eof())
        {
          // read path
          read >> full_path;

          // extract the path and pdb tag from the full path provided
          storage::VectorND< 2, std::string> this_tag( pdb::Handler::ExtractPathAndPDBTag( full_path));

          // output the tag
          BCL_MessageStd
          (
            this_tag.Second() + " at " + this_tag.First()
          );

          // insert into to the vector of tags
          pdb_tags.PushBack( this_tag);
        }

        io::File::CloseClearFStream( read);
        BCL_MessageStd( "pdb list read");

        // iterate over proteins, read them in and fills up the data vector
        for
        (
          storage::Vector< storage::VectorND< 2, std::string> >::const_iterator pdb_itr( pdb_tags.Begin()),
            pdb_itr_end( pdb_tags.End());
          pdb_itr != pdb_itr_end;
          ++pdb_itr
        )
        {
          // initializes pdb tag
          const std::string path( pdb_itr->First());
          const std::string pdb_tag( pdb_itr->Second());
          const std::string pdb_filename( path + PATH_SEPARATOR + pdb_tag + ".pdb");

          BCL_MessageStd( pdb_tag + " at " + path);

          // initialize the pdb handler and factory
          pdb::Factory pdb_factory( biol::GetAAClasses().e_AABackBone);

          io::File::MustOpenIFStream( read, pdb_filename);
          pdb::Handler pdb_handler( read);
          io::File::CloseClearFStream( read);

          // read in the sequences from pdb file
          storage::Map< biol::SSType, size_t> ssetype_min_size;
          ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
          ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
          assemble::ProteinModel protein_model( pdb_factory.ProteinModelFromPDB( pdb_handler, ssetype_min_size));

          // if no prediction flag is not set
          if( !m_NoPredictionFlag->GetFlag())
          {
            // if blast profile is being read, read the ascii
            biol::BlastProfileHandler::TryReadProfileForProteinModel( protein_model, path + PATH_SEPARATOR + pdb_tag);

            // calculate jufo for the given protein model
            sspred::JUFO::Calculate( protein_model);
          }

          // iterate over all the sequences
          for
          (
            util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( protein_model.GetChains().Begin()),
              chain_itr_end( protein_model.GetChains().End());
            chain_itr != chain_itr_end;
            ++chain_itr
          )
          {
            // output chain information
            BCL_MessageStd
            (
              "\t chain :" + util::Format()( ( *chain_itr)->GetChainID()) +
              " : " + util::Format()( ( *chain_itr)->GetSequence()->GetSize()) + "residues"
            );

            // if this chain lacks JUFO,PSIPRED of Blast profile skip it
            if
            (
              !m_NoPredictionFlag->GetFlag() &&
              (
                !( *chain_itr)->GetSequence()->GetFirstAA()->GetSSPrediction( sspred::GetMethods().e_JUFO).IsDefined() ||
                !( *chain_itr)->GetSequence()->GetFirstAA()->GetBlastProfilePtr().IsDefined()
              )
            )
            {
              BCL_MessageStd
              (
                "The following chain does not have JUFO or BlastProfile and therefore is being skipped: " +
                util::Format()( pdb_tag) + " " + util::Format()( ( *chain_itr)->GetChainID())
              );
              continue;
            }
            // else create and output contact map
            else
            {
              CreateContactMap( *chain_itr, pdb_tag);
            }

          } // end iterate over chains
        } // end iterate over pdbs
      }
      BCL_MessageStd( "contact prediction has been completed");

      //end
      return 0;

    } // end Main

    //! @brief creates the contact map and real contact map( if requested) and outputs them
    //! @brief for the provided chain and pdb_tag
    //! @param SP_CHAIN ShPtr to the chain for which the contact map is going to be generated
    //! @param PDB_TAG pdb tag of the protein of interest
    void ContactPrediction::CreateContactMap
    (
      const util::ShPtr< assemble::Chain> &SP_CHAIN,
      const std::string &PDB_TAG
    ) const
    {
      // output stream
      io::OFStream write;

      // if no prediction flag is not set
      if( !m_NoPredictionFlag->GetFlag())
      {
        // initialize the prediction map
        contact::PredictionMap prediction_map( *SP_CHAIN);

        // prepare the output file name
        const std::string prediction_map_filename
        (
          m_OutputPathFlag->GetFirstParameter()->GetValue() + PATH_SEPARATOR +
          PDB_TAG + SP_CHAIN->GetChainID() + ".contact"
        );

        // initialize write stream with the filename
        io::File::MustOpenOFStream( write, prediction_map_filename);

        // output the map to file
        prediction_map.WritePredictionMap( write, m_ThresholdFlag->GetFirstParameter()->GetNumericalValue< double>());

        // clean the write stream
        io::File::CloseClearFStream( write);
      }

      // if m_RealContactsFlag is set
      if( m_RealContactsFlag->GetFlag())
      {
        // create the real contacts map
        contact::Map contact_map( SP_CHAIN);

        // prepare the output file name
        const std::string contact_map_filename
        (
          m_OutputPathFlag->GetFirstParameter()->GetValue() + PATH_SEPARATOR +
          PDB_TAG + SP_CHAIN->GetChainID() + ".contacts"
        );

        // create the output stream with the filename
        // make sure the file has been opened
        io::File::MustOpenOFStream( write, contact_map_filename);

        // output the map to file
        contact_map.WriteMap( write, true);

        // clean the write stream
        io::File::CloseClearFStream( write);
      }
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ContactPrediction::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddParameter( m_InputFilenameParam);
      sp_cmd->AddFlag( m_PDBListFlag);
      sp_cmd->AddFlag( m_ThresholdFlag);
      sp_cmd->AddFlag( m_RealContactsFlag);
      sp_cmd->AddFlag( m_NoPredictionFlag);
      sp_cmd->AddFlag( m_OutputPathFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! default constructor
    ContactPrediction::ContactPrediction() :
      m_InputFilenameParam
      (
        new command::Parameter( "input_filename", "\name of the input file")
      ),
      m_PDBListFlag
      (
        new command::FlagStatic
        (
          "pdb_list",
          "\t provides the functionality to input multiple pdb tags"
        )
      ),
      m_ThresholdFlag
      (
        new command::FlagStatic
        (
          "threshold",
          "\tthe merged contact prediction threshold, residue pairs below this prediction do not get printed out",
          command::Parameter
          (
            "threshold",
            "\tthreshold value",
            command::ParameterCheckRanged< double>( 0.0, 1.0),
            "0.4"
          )
        )
      ),
      m_RealContactsFlag
      (
        new command::FlagStatic( "real_contacts", "\toutput real contacts")
      ),
      m_NoPredictionFlag
      (
        new command::FlagStatic( "no_prediction", "\tdo not read ssfiles and do not generate predictions")
      ),
      m_OutputPathFlag
      (
        new command::FlagStatic
        (
          "output_path",
          "\tdefine the output path",
          command::Parameter( "output_path", "\tany valid output path", ".")
        )
      )
    {
    }

    const ApplicationType ContactPrediction::ContactPrediction_Instance
    (
      GetAppGroups().AddAppToGroup( new ContactPrediction(), GetAppGroups().e_BioInfo)
    );

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ContactPrediction::GetReadMe() const
    {
      // construct static readme
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::Contact, terms of use, appropriate citation, installation "
        "procedures, BCL::Contact execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::Contact?\n"
        "BCL::Contact is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part\n"
        " a larger library of applications called BCL::Commons.  BCL::Contact predicts residue-residue pair contacts \n"
        "for a given protein primary sequence. It utilizes an Artificial Neural Network (ANN) as prediction method. "
        "Position-specific scoring matrices from PsiBlast, amino acid properties and secondary structure predictions "
        "from Bcl::Jufo are used as input. BCL::Contact relies solely on sequence information and has  individual ANNs "
        "specialized for helix-helix, helix-strand, strand-helix, strand-strand, and sheet-sheet contacts. Output is a "
        "three-state (helix, strand, undefined) secondary structure probability profile for each amino acid in the "
        "sequence of interest.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::Contact.\n"
        "When using BCL::Contact in a publication, please cite the following publications describing the application's "
        "development:\n"
        "\n"
        "Karakas M., Woetzen N., Meiler J. BCL::Contact - low confidence fold recognition hits boost protein contact "
        "prediction and de novo structure determination J. Comput. Biol., 17, (2), 153-168. ; 2010\n"
        "Link:  www.meilerlab.org/index.php/publications/showPublication/pub_id/81\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        "1) BCL::Contact:\n"
        + DefaultInstallationProcedure() +
        "\n"
        "2) PsiBlast:\n"
        "In order to run BCL::Contact, PsiBlast must be installed.  Free downloads can be found at "
        "http://www.ncbi.nlm.nih.gov/BLAST/download.shtml.  Please make sure to download the version for the correct "
        "platform and follow all instructions on the screen.\n"
        "\n"
        "The Basic Local Alignment Search Tool (Blast) finds regions of local similarity between sequences. The "
        "program compares nucleotide or protein sequences to sequence databases and calculates the statistical "
        "significance of matches. Blast can be used to infer functional and evolutionary relationships between "
        "sequences as well as help identify members of gene families.  For more information please visit "
        "http://blast.ncbi.nlm.nih.gov/.\n"
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::Contact.\n"
        "Running BCL::Contact consists of three main steps.\n"
        "\n"
        "1) Create the fasta sequence file for the protein to be studied:\n"
        "You will need the protein sequence in fasta format for both BCL::Contact and PsiBlast, and it should be "
        "stored in a <.fasta> file.  An example is given below. For more information about fasta formats, please visit "
        "http://www.ncbi.nlm.nih.gov/Blast/fasta.shtml.\n"
        "\n"
        "2) Create the PsiBlast position-specific scoring matrix:\n"
        "Run PsiBlast on the fasta sequence to produce a PsiBlast position-specific scoring matrix with extension "
        "\".ascii.\"  More information on running PsiBlast and adjusting various parameters can be found in the "
        "documentation which accompanied the download.  An example run of PsiBlast at the command line could look like "
        "the following, where $MyDataBase and $MyBlastProfile are the names of the database used and desired PsiBlast "
        "position-specific scoring matrix output file respectively:\n"
        "\n"
        "blastpgp -b 0 -j 3 -h 0.001 -d $MyDataBase -i MyFastaSequence.fasta -C MyCheckPoint.chk "
        "-Q $MyBlastProfile.ascii\n"
        "\n"
        "3) Run BCL::Contact:\n"
        "At a command prompt, navigate to the location of your BCL::Contact executable program. The syntax for running "
        "the application looks like the following:\n"
        "\n"
        "bcl.exe ContactPrediction <fasta_file> -output_path <output_path> -threshold <threshold>\n"
        "\n"
        "BCL::Contact needs a fasta sequence with the extension <.fasta> and a PsiBlast position-specific scoring "
        "matrix file with the extension <.ascii> to exist in the same directory.\n"
        "\n"
        "FLAGS:\n"
        "-output_path <output_path>   This flag defines the folder where the output files will be written to if not "
        "provided, the output files will be written to the folder to the current working directory.\n"
        "\n"
        "-threshold <threshold> This flag determines which residue pairs gets written to contact file. Only the "
        "residue pairs with a normalized prediction larger than the given threshold will be output. The default value "
        "for the threshold is 0.4. This value needs to between 0 and 1, the higher the value is, less residue pairs "
        "with higher confidence are output.\n"
        "\n"
        "INPUT AND OUTPUT.\n"
        "\n"
        "BCL::Contact requires two inputs, a fasta file and a corresponding PsiBlast position-specific scoring matrix. "
        "The fasta file uses one letter codes for protein sequence and looks like the following:\n"
        "\n"
        ">1UBI:A|PDBID|CHAIN|SEQUENCE\n"
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n"
        "\n"
        "The PsiBlast position-specific scoring matrix on the other hand looks like following:\n"
        "\n"
        "Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per "
        "position, and relative weight of gapless real matches to pseudocounts\n"
        "\n"
        "        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V\n"
        " 1 M   -2 -2 -3 -4 -2 -1 -3 -4 -3  1  2 -2  9 -1 -4 -3 -2 -2 -2  0   0   0   0   0   0   0   0   0   0   1   4   0  94   0   0   0   0   0   0   0  1.03 0.31\n"
        " 2 Q   -2  0 -1 -1 -4  7  1 -3  1 -4 -3  0  1 -4 -2 -1 -2 -3 -2 -3   0   0   0   0   0  92   0   0   4   0   0   0   4   0   0   0   0   0   0   0  1.07 0.31\n"
        "\n"
        "\n"
        "The output file is formatted as following:\n"
        "\n"
        "CHAINS A A\n"
        "   5 V    9 T 0.038 0.057 0.410 0.163 0.117 0.075\n"
        "   5 V   10 G 0.038 0.047 0.063 0.273 0.470 0.044\n"
        "   5 V   11 K 0.046 0.082 0.141 0.493 0.370 0.219\n"
        "   5 V   12 T 0.011 0.408 0.050 0.830 0.507 0.659\n"
        "   5 V   13 I 0.013 0.130 0.035 0.808 0.715 0.748\n"
        "\n"
        "The first line specifies between which chains the contacts are predicted"
        "The individual columns represent following:\n"
        "Column 1: Amino acid position of the first residue\n"
        "Column 2: One letter code for amino acid type of the first residue\n"
        "Column 3: Amino acid position of the second residue\n"
        "Column 4: One letter code for amino acid type of the second residue\n"
        "Column 5: Prediction of specified two residues to be in a helix-helix contact (between 0 and 1)\n"
        "Column 6: Prediction of specified two residues to be in a helix-sheet contact (between 0 and 1)\n"
        "Column 7: Prediction of specified two residues to be in a sheet-helix contact (between 0 and 1)\n"
        "Column 8: Prediction of specified two residues to be in a strand-strand contact (between 0 and 1)\n"
        "Column 9: Prediction of specified two residues to be in a sheet-sheet contact (between 0 and 1)\n"
        "Column 10: Normalized prediction of specified two residues to be in a contact(between 0 and 1)\n"
        "This is calculated by summing over previous 5 columns after normalizing each column with the likelihood of "
        "specified residues to be in secondary structures specified by the contact-type. For example, for the first "
        "line, the helix-sheet prediction of 0.0057 is normalized by probability of first residue (#5) to be in a "
        "helix and probability of the second residue (#9) to be in a strand. These probabilities are retrieved from "
        "three-state predictions given by BCL::Jufo secondary structure prediction method\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing\n"
        "bcl.exe ContactPrediction  -help\n"
        "\n"
        "For more general information about the product, type\n"
        "bcl.exe ContactPrediction -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::Contact.\n"
        "BCL::Contact is under ongoing further development. For current research please refer to www.meilerlab.org and "
        "navigate to research\n"
        + DefaultSectionSeparator()
      );

      // return
      return s_readme;
    }

  } // namespace app
} // namespace bcl
