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
#include "bcl_app_jufo.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_with_cache_storage_file.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "sspred/bcl_sspred_jufo.h"
#include "sspred/bcl_sspred_jufo9d.h"
#include "sspred/bcl_sspred_method_handler.h"

namespace bcl
{
  namespace app
  {
  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief default constructor
    Jufo::Jufo() :
      m_FastaFileParam( new command::Parameter( "fasta", "fasta file to use for Jufo predictions")),
      m_AsciiFileFlag
      (
        new command::FlagStatic
        (
          "ascii",
          "blast ascii6 (just .ascii for 3d) file to use for Jufo predictions, if flag is not given, use ascii file in fasta file path",
          command::Parameter( "ascii", "blast ascii file to use for Jufo predictions", "")
        )
      ),
      m_OutputFlag
      (
        new command::FlagStatic
        (
          "output",
          "optional output prefix - default will be the given pdb id",
          command::Parameter( "output_filename", "output prefix", "")
        )
      ),
      m_Jufo3dFlag
      (
        new command::FlagStatic( "3d", "use original Jufo algorithm (only 3-state predictions)")
      ),
      m_MultimerFlag
      (
        new command::FlagStatic( "multimer", "use this flag if the protein is a multimer")
      )
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Jufo::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string Jufo::GetDescription() const
    {
      return "predicts secondary structure for a given protein primary sequence";
    }

  ////////////////
  // operations //
  ////////////////

    util::ShPtr< command::Command> Jufo::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add member parameters and flags
      sp_cmd->AddParameter( m_FastaFileParam);
      sp_cmd->AddFlag( m_AsciiFileFlag);
      sp_cmd->AddFlag( m_OutputFlag);
      sp_cmd->AddFlag( m_Jufo3dFlag);
      sp_cmd->AddFlag( m_MultimerFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! Main
    int Jufo::Main() const
    {
      // read fasta
      const std::string fasta_filename( m_FastaFileParam->GetValue());
      assemble::ProteinModelWithCache model
      (
        assemble::ProteinWithCacheStorageFile::GenerateProteinModelFromFile
        (
          fasta_filename,
          false,
          biol::GetAAClasses().e_AA
        ),
        false
      );

      // determine the prefix (everything that preceeds ".fasta"
      const size_t pos_period( fasta_filename.rfind( "."));
      const std::string prefix( fasta_filename.substr( 0, pos_period));

      // read ascii, use parameter if given, if not replace ".fasta" with ".ascii" in fasta_filename and try that
      BCL_Assert
      (
        biol::BlastProfileHandler::TryReadProfileForProteinModel( model, prefix, m_Jufo3dFlag->GetFlag() ? ".ascii" : ""),
        "Required blast profile file could not be located!"
      );

      // initialize output
      io::OFStream write;
      std::string output_name
      (
        m_OutputFlag->GetFlag() ?
          m_OutputFlag->GetFirstParameter()->GetValue() :
          prefix
      );

      // if using old Jufo
      if( m_Jufo3dFlag->GetFlag())
      {
        // calculate Jufo
        sspred::JUFO::Calculate( model);

        const biol::AASequence &seq( *model.GetChains().FirstElement()->GetSequence());

        // write predictions
        io::File::MustOpenOFStream( write, output_name + ( *sspred::GetMethods().e_JUFO)->GetFileExtension());
        sspred::MethodHandler::WritePredictionsForAASequence( write, seq, sspred::GetMethods().e_JUFO);
        io::File::CloseClearFStream( write);
      }
      // use new 9D jufo
      else
      {
        // calculate the predictions
        sspred::JUFO9D::Calculate( model, m_MultimerFlag->GetFlag());

        const biol::AASequence &seq( *model.GetChains().FirstElement()->GetSequence());

        // output the 9-state predictions
        output_name += ( *sspred::GetMethods().e_JUFO9D)->GetFileExtension();
        io::File::MustOpenOFStream( write, output_name);
        sspred::MethodHandler::WritePredictionsForAASequence
        (
          write,
          seq,
          sspred::GetMethods().e_JUFO9D,
          sspred::MethodHandler::e_NineState
        );
        // clear stream
        io::File::CloseClearFStream( write);

        // output the 3-state SS predictions
        io::File::MustOpenOFStream( write, output_name + "_ss");
        sspred::MethodHandler::WritePredictionsForAASequence
        (
          write,
          seq,
          sspred::GetMethods().e_JUFO9D,
          sspred::MethodHandler::e_ThreeState
        );
        // clear stream
        io::File::CloseClearFStream( write);

        // write out 3-state TM predictions
        io::File::MustOpenOFStream( write, output_name + "_tm");
        sspred::JUFO9D::WriteThreeStateTMPredictions( write, seq);
        // clear stream
        io::File::CloseClearFStream( write);

        // write out 2-state predictions for helix
        io::File::MustOpenOFStream( write, output_name + "_tmh");
        sspred::JUFO9D::WriteTwoStateTMPredictions( write, seq, biol::GetSSTypes().HELIX);
        // clear stream
        io::File::CloseClearFStream( write);

        // write out 2-state predictions for strand
        io::File::MustOpenOFStream( write, output_name + "_tms");
        sspred::JUFO9D::WriteTwoStateTMPredictions( write, seq, biol::GetSSTypes().STRAND);
        // clear stream
        io::File::CloseClearFStream( write);
      }

      // end
      return 0;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &Jufo::GetReadMe() const
    {
      // create a static string to hold readme information
      static const std::string s_readme_text
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::Jufo, terms of use, appropriate citation, installation "
        "procedures, BCL::Jufo execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::Jufo?\n"
        "\n"
        "BCL::Jufo is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part of "
        "a larger library of applications called BCL::Commons.  BCL::Jufo predicts secondary structure for a given "
        "protein primary sequence. It utilizes an Artificial Neural Network as prediction method.  Position-specific "
        "scoring matrices from PsiBlast as well as amino acid properties are used as input.  Output is a three-state"
        "(helix, strand, undefined) secondary structure probability profile for each amino acid in the sequence of "
        "interest.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        "\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::Jufo.\n"
        "\n"
        "When using BCL::Jufo in a publication, please cite the following publications describing the application's "
        "development:\n"
        "J. K. Leman, R. Mueller, M. Karakas, N. Woetzel, and J. Meiler, Simultaneous prediction of protein secondary "
        "structure and trans-membrane spans., Proteins, Jan. 2013.\n"
        "Link:  http://meilerlab.org/index.php/publications/showPublication/pub_id/142\n"
        "\n"
        "Meiler J., Baker D. Coupled Prediction of Protein Secondary and Tertiary Structure. PNAS, 100, (21), "
        "12105-12110. ; 2003\n"
        "Link:  http://meilerlab.org/index.php/publications/showPublication/pub_id/27\n"
        "\n"
        "Meiler J., Mueller M., Zeidler A., Schmaeschke F. Generation and Evaluation of Dimension Reduced Amino Acid "
        "Parameter Representations by Artificial Neural Networks. J. Mol. Model., 7, (9), 360-369. ; 2001\n"
        "Link:  http://meilerlab.org/index.php/publications/showPublication/pub_id/37\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        "1) BCL::Jufo:\n"
        + DefaultInstallationProcedure() +
        "\n"
        "2) PsiBlast:\n"
        "In order to run BCL::Jufo, PsiBlast must be installed.  Free downloads can be found at "
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
        "VI. RUNNING BCL::Jufo.\n"
        "\n"
        "Running BCL::Jufo consists of three main steps.\n"
        "\n"
        "1) Create the fasta sequence file for the protein to be studied:\n"
        "You will need the protein sequence in fasta format for both BCL::Jufo and PsiBlast, and it should be stored "
        "in a <.fasta> file.  An example is given below.  For more information about fasta formats, please visit "
        "http://www.ncbi.nlm.nih.gov/Blast/fasta.shtml.\n"
        "\n"
        "2) Create the PsiBlast position-specific scoring matrix:\n"
        "Run PsiBlast on the fasta sequence to produce a PsiBlast position-specific scoring matrix with extension "
        "\".ascii.\"  More information on running PsiBlast and adjusting various parameters can be found in the "
        "documentation which accompanied the download.  An example run of PsiBlast at the command line could look like "
        "the following, where $MyDataBase and $MyBlastProfile are the names of the database used and desired PsiBlast "
        "position-specific scoring matrix output file respectively:\n"
        "\n"
        "blastpgp -b 0 -j 3 -h 0.001 -d $MyDataBase -i MyFastaSequence.fasta -C MyCheckPoint.chk -Q $MyBlastProfile.ascii\n"
        "\n"
        "3) Run BCL::Jufo:\n"
        "At a command prompt, navigate to the location of your BCL::Jufo executable program.  The syntax for running "
        "the application looks like the following:\n"
        "\n"
        "bcl_jufo.exe -pdbid <prefix> -output <output_filename>\n"
        "\n"
        "BCL::Jufo needs a fasta sequence with the extension <.fasta> and a PsiBlast position-specific scoring matrix "
        "file with the extension <.ascii> to exist in the same directory.\n"
        "\n"
        "FLAGS:\n"
        "\n"
        "-pdbid <tag>   The application finds fasta and Blast needed to run BCL::Jufo by adding <.fasta> and "
        "<.ascii> extensions to the  <tag> value provided. The <tag> value is usually a four letter pdb code such as "
        "\"1UBI\" but also can include path \"/home/user/sequences/1UBI\".\n"
        "\n"
        "-output  <output_filename>   By default the output file is created with name <tag>  + \".Jufo\". The user can "
        "specify a different output name by providing it at <output_filename> such as \"sequence.Jufo” \" or "
        "\"/home/user/sequences/1UBI.Jufo.\"\n"
        "\n"
        "INPUT AND OUTPUT.\n"
        "\n"
        "BCL::Jufo requires two inputs, a fasta file and a corresponding PsiBlast position-specific scoring matrix. "
        "The fasta file uses one letter codes for protein sequence and looks like the following:\n"
        "\n"
        ">1UBI:A|PDBID|CHAIN|SEQUENCE\n"
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n"
        "\n"
        "The PsiBlast position-specific scoring matrix on the other hand looks like following:\n"
        "\n"
        "Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per "
        "position, and relative weight of gap less real matches to pseudo counts\n"
        "\n"
        "           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A \n"
        " R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T \n"
        " W   Y   V\n"
        "    1 M   -4 -4 -5 -6 -4 -4 -5 -6 -5  4  2 -4  9 -3 -5 -4 -3 -4 -4 -1    0\n"
        "  0   0   0   0   0   0   0   0  24  13   0  63   0   0   0   0   0  \n"
        "0   0  1.72 1.27\n"
        "    2 Q   -4  3  1 -1 -6  5  1 -5 -1 -2 -3  3 -1 -2 -4 -1 -2  2 -2 -5    0\n"
        " 15   6   3   0  32   6   0   1   3   2  20   1   1   0   4   2   3  \n"
        "1   0  0.74 1.36\n"
        "\n"
        "The output file is formatted as following:\n"
        "\n"
        "   1 M C   0.603  0.006  0.391\n"
        "   2 Q E   0.214  0.015  0.771\n"
        "   3 I E   0.225  0.015  0.760\n"
        "   4 F E   0.063  0.003  0.934\n"
        "   5 V E   0.077  0.019  0.905\n"
        "   6 K E   0.263  0.059  0.678\n"
        "   7 T E   0.324  0.049  0.627\n"
        "   8 L C   0.586  0.031  0.383\n"
        "   9 T C   0.901  0.032  0.067\n"
        "  10 G C   0.898  0.042  0.060\n"
        "\n"
        "The individual columns represent following:\n"
        "Column 1: Amino acid position in sequence\n"
        "Column 2: One letter code for amino acid type\n"
        "Column 3: One letter code for the predicted secondary structure;\n"
        "C=Undefined/Loop/Coil, E=Strand, H=Helix\n"
        "Column 4: Undefined/Loop/Coil probability (between 0 and 1) for this amino acid\n"
        "Column 5: Helix probability (between 0 and 1) for this amino acid\n"
        "Column 6: Strand probability (between 0 and 1) for this amino acid\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing bcl_jufo.exe -help\n"
        "\n"
        "For more general information about the product, type bcl_jufo.exe -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        "\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::Jufo.\n"
        "\n"
        "BCL::Jufo is under ongoing further development.  For current research please refer to www.meilerlab.org and "
        "navigate to research.  The information can be found under Simultaneous Prediction of Protein Secondary "
        "Structure and Trans-Membrane Spanning Regions.\n"
        + DefaultSectionSeparator()
      );

      // return readme information
      return s_readme_text;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &Jufo::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::Jufo: Simultaneous Prediction of Protein Secondary Structure and Trans-Membrane Spans\n\n"
        "A first step towards protein tertiary structure prediction is the identification of secondary structure "
        "elements from the sequence. In addition, the identification of trans-membrane spans is required for membrane "
        "proteins.\n\n"
        "The aim of this project is to simultaneously predict secondary structure and trans-membrane spans with a "
        "single tool. The rationale for this approach is the hypothesis that both phenomena are interrelated: When a "
        "nascent polypeptide reaches the membrane interface the altered dielectric environment (described by the free"
        "energy) leads to the formation of hydrogen bonds and therefore  secondary structure. To date there is no "
        "single tool available that is suitable for the prediction of both trans-membrane &alpha;-helices and "
        "&szlig;-strands at the same time. \n"
        "!jufo_ss_tm.gif!\n"
        "Fig. 1:\n"
        "\n\nAn Artificial Neural Network (ANN) will be trained on non-redundant datasets of both soluble and membrane "
        "proteins (sequence similarity &lt;25%) to accomplish this task (Fig. 1). As input serve several amino acid "
        "properties, position-specific scoring matrices from PSIBLAST, and knowledge-based free energies for the "
        "secondary structure types helix, strand, or coil and the regions trans-membrane, interface, or solution. A "
        "matrix of these numbers over a window of 31 amino acids forms the input for the ANN. The output is a "
        "nine-state prediction for the central residue in the sequence window.\n"
        "!jufo_schematic.png!\n"
        "Fig. 2:\n"
        "\n\nScheme of BCL::Jufo. As amino acid properties the steric parameter, polarizability, volume, iso-electric "
        "point, and the solvent-accessible surface area are used as input. The free energies are derived from the "
        "membrane protein database.\n"
      );
      return s_web_text;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const ApplicationType Jufo::Jufo_Instance
    (
      GetAppGroups().AddAppToGroup( new Jufo(), GetAppGroups().e_BioInfo)
    );

  } // namespace app
} // namespace bcl
