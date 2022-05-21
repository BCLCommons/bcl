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
#include "bcl_app_pdb_compare.h"
#include "align/bcl_align_alignment_node_reference.h"
#include "sched/bcl_sched_thunk_job.h"

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace app
  {
  //////////
  // data //
  //////////

    const ApplicationType ProteinCompare::PDBCompare_Instance
    (
      GetAppGroups().AddAppToGroup( new ProteinCompare(), GetAppGroups().e_Protein)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // default constructor
    ProteinCompare::ProteinCompare() :
      m_ReferencePDB
      (
        new command::FlagDynamic
        (
          "reference_pdb", "pdb which other pdbs will be compared to",
           command::Parameter
          (
            "reference_pdb_file_name",
            "filename of pdb which other pdbs will be compared to"
          ), 0, 1
        )
      ),
      m_TestPDB
      (
        new command::FlagDynamic
        (
          "test_pdb", "pdb which will be compared to the reference pdb",
           command::Parameter
          (
            "test_pdb_file_name",
            "filename of pdb which will be compared to the reference pdb"
          ), 0, 1
        )
      ),
      m_PDBList
      (
        new command::FlagDynamic
        (
          "pdb_list",
          "If reference_pdb flag is set then pdbs in the list will be compared to the reference pdb. If reference_pdb flag is not set then qualities between all pdbs will be calculated",
           command::Parameter
          (
            "pdb_list",
            "filename of list of pdbs for quality calculations"
          ), 0, 1
        )
      ),
      m_OutFile
      (
        new command::FlagStatic
        (
          "prefix",
          "\tprefix for output files",
          command::Parameter
          (
            "string",
            "\tprefix for output files",
            ""
          )
        )
      ),
      m_QualityNormalization100
      (
        new command::FlagStatic
        (
          "norm100",
          "adds additional entry for each quality measure that is normalized to 100 residues"
        )
      ),
      m_NoLabels
      (
        new command::FlagStatic
        (
          "no_labels",
          "indicates that no labeling for the matrix is desired"
        )
      ),
      m_InputPDBPrefix
      (
        new command::FlagStatic
        (
          "input_pdb_prefix",
          "\tindicates that a string must be prepended to the PDB strings in order to access the actual PDB files",
          command::Parameter
          (
            "string",
            "\tprefix needed for accessing PDB files",
            ""
          )
        )
      ),
      m_InputPDBPostfix
      (
        new command::FlagStatic
        (
          "input_pdb_postfix",
          "\tindicates that a string must be appended to the PDB strings in order to access the actual PDB files",
          command::Parameter
          (
            "string",
            "\tprefix needed for accessing PDB files",
            ""
          )
        )
      ),
      m_OutputDirectory
      (
        new command::FlagStatic
        (
          "output_dir",
          "\tthe output directory into which the quality measure matrices will be deposited",
          command::Parameter
          (
            "string",
            "\t indicates the output directory into which the quality measure matrices will be deposited",
            ""
          )
        )
      ),
      m_Topology
      (
        new command::FlagStatic
        (
          "topology",
          "use topology score as metric for comparison"
        )
      ),
      m_ConsideredAtoms(),
      m_Alignment(),
      m_Qualities()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Quality
    ProteinCompare *ProteinCompare::Clone() const
    {
      return new ProteinCompare( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ProteinCompare::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the original name of the application, as it was used in license files
    //! @return string for the bcl::commons name of that application
    //! This is necessary so that, if release application names change, licenses will continue to work
    std::string ProteinCompare::GetLicensedName() const
    {
      return "PDBCompare";
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ProteinCompare::GetDescription() const
    {
      return "calculates the structural differences between protein models using methods such as RMSD and GDT";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProteinCompare::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::PDBCompare, terms of use, appropriate citation, installation "
        "procedures, BCL::PDBCompare execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::PDBCompare?\n"
        "BCL::PDBCompare is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::PDBCompare is a utility for calculating the "
        "structural differences between protein models. A variety of comparison algorithms are available such as RMSD.\n"
        "\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::PDBCompare.\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::PDBCompare.\n"
        "Running BCL::PDBCompare consists of two main steps.\n"
        "\n"
        "1) Get at least two protein structures in PDB format.\n"
        "\n"
        "2) Run BCL::PDBCompare to calculate any number of comparisons between the models\n"
        "\n"
        "At a command prompt, navigate to the location of your BCL::PDBCompare executable program.\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe PDBCompare -help\n"
        "\n"
        "For more general information about the product, type bcl.exe PDBCompare -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::PDBCompare.\n"
        "BCL::PDBCompare is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );

      return readme;
    }

    //! @brief return command line flag specifying specific residues that should be used for quality calculations
    //! @return command line flag for specifying specific residues to use during quality calculations
    util::ShPtr< command::FlagInterface> &ProteinCompare::GetFlagSpecifyQualityResidues()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "specify_residues",
          "\tonly specified residues will be used for quality calculations. File format for specifying residues is <chain id> <resi seq id> per line. For example, 'A' 32"
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_file
      (
        new command::Parameter
        (
          "filename", "\tname of the file containing the list of residues", "residues.ls"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_file);
      }

      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ProteinCompare::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      //input pdb names
      sp_cmd->AddFlag( m_ReferencePDB);
      sp_cmd->AddFlag( m_TestPDB);

      // input quality flag
      sp_cmd->AddFlag( quality::Measures::GetFlagQualityMeasures());

      // input atom flag
      sp_cmd->AddFlag( biol::AtomTypes::GetFlagAtomTypes());

      // input list of pdbs
      sp_cmd->AddFlag( m_PDBList);

      //input alignment
      sp_cmd->AddFlag( assemble::Quality::GetFlagAlignments());

      // input specific residues
      sp_cmd->AddFlag( GetFlagSpecifyQualityResidues());

      // input output file name
      sp_cmd->AddFlag( m_OutFile);

      // normalization flag
      sp_cmd->AddFlag( m_QualityNormalization100);

      // no labeling the matrix flag
      sp_cmd->AddFlag( m_NoLabels);

      // pdb string prefix flag
      sp_cmd->AddFlag( m_InputPDBPrefix);

      // pdb string postfix flag
      sp_cmd->AddFlag( m_InputPDBPostfix);

      // output directory flag
      sp_cmd->AddFlag( m_OutputDirectory);

      // flag to use topology score as metric
      sp_cmd->AddFlag( m_Topology);

      // aa class type
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // sse min size for creating protein models from pdb
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // use backbone conformation if no sse definitions are in the pdb
      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());

      // add flag that allows conversion of non-natural aa types to their natural equivalents (e.g. Sel-Met to Met)
      sp_cmd->AddFlag( pdb::Factory::GetFlagConvertToNaturalAAType());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ProteinCompare::Main() const
    {
      // make sure test_pdb and pdb_list flags are not set simultaneously
      BCL_Assert
      (
        !( m_TestPDB->GetFlag() && m_PDBList->GetFlag()), "do not set test_pdb and pdb_list flags simultaneously"
      );

      // make sure one of test_pdb or pdb_list flags is set
      BCL_Assert
      (
        m_TestPDB->GetFlag() || m_PDBList->GetFlag(), "either test_pdb or pdb_list flags must be set"
      );

      // make sure that if the reference_pdb flag is not set that the pdb_list flag is set
      if( !m_ReferencePDB->GetFlag())
      {
        BCL_Assert
        (
          m_PDBList->GetFlag(), "if the reference_pdb flag is not set then the pdb_list flag must be provided"
        );
      }

      // make sure that the alignment flag is only used with the reference flag
      BCL_Assert
      (
        !( assemble::Quality::GetFlagAlignments()->GetFlag() && !m_ReferencePDB->GetFlag()),
        "alignment flag must be used with reference_pdb flag"
      );

      // create Vector "string_vector_a" to hold first list of pdbs
      storage::Vector< std::string> string_vector_a;

      // create Vector "string_vector_b" to hold second list of pdbs
      storage::Vector< std::string> string_vector_b;

      // bool for determining if one only needs to calculate half of the matrix
      bool half( false);

      // true if pdb list and reference pdb is given; calculate quality between reference pdb and each pdb in list
      if( m_PDBList->GetFlag() && m_ReferencePDB->GetFlag())
      {
        // create Vector "string_vector_a" and fill with the pdbs in the pdb list
        string_vector_a = StringVectorFromPDBListParameter( m_PDBList);

        // create Vector "string_vector_b" and add the reference pdb
        string_vector_b.PushBack( m_ReferencePDB->GetFirstParameter()->GetValue());
      }

      // true if test pdb and reference pdb is give; calculates quality between these two pdbs
      else if( m_TestPDB->GetFlag() && m_ReferencePDB->GetFlag())
      {
        // create Vector "string_vector_a" and add the reference pdb
        string_vector_a.PushBack( m_TestPDB->GetFirstParameter()->GetValue());

        // create Vector "string_vector_b" and add the reference pdb
        string_vector_b.PushBack( m_ReferencePDB->GetFirstParameter()->GetValue());
      }

      // pdb list is given but no reference pb is given so calculate pairwise quality between all pdbs in pdb list
      else if( m_PDBList->GetFlag() && !m_ReferencePDB->GetFlag())
      {
        // create Vector "string_vector_a" and fill with the pdbs in the pdb list
        string_vector_a = StringVectorFromPDBListParameter( m_PDBList);

        // create Vector "string_vector_b" with string_vector_a because no reference pdb was given
        string_vector_b = string_vector_a;

        // only calculate the half of the matrix for reasons above
        half = true;
      }
      else
      {
        BCL_Exit( "unable to do quality calculations", -1);
      }

      // get atoms for quality calculation
      m_ConsideredAtoms = biol::AtomTypes::GetCommandLineAtoms();

      // if an alignment is used
      if( assemble::Quality::GetFlagAlignments()->GetFlag())
      {
        // get the alignment files
        const storage::Vector< io::DirectoryEntry> align_files
        (
          assemble::Quality::GetFlagAlignments()->GetObjectList< io::DirectoryEntry>()
        );

        m_Alignment = assemble::Quality::GetChainAlignments( align::HandlerPIR< biol::AABase>(), align_files);
      }

      // columns for each pdb in vector b
      const util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( string_vector_b));

      // names of rows in table
      const storage::Vector< std::string> &row_names( string_vector_a);

      // if topology flag is given, fill table with values from topology calculation
      if( m_Topology->GetFlag())
      {
        // read in first pdb from list
        std::string first_pdb( *( string_vector_a.Begin()));
        pdb::Factory factory;
        io::IFStream read;
        io::File::MustOpenIFStream
        (
          read,
          first_pdb
        );
        
        // get protein model from PDB
        pdb::Handler pdb( read, true);
        io::File::CloseClearFStream( read);
        const assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb));

        // initialize string that will represent table header
        storage::Vector< std::string> header_names;

        //initialize table
        const util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( row_names));
        storage::Table< double> topology_table( sp_table_header);

        // insert empty rows
        for
        (
            storage::Vector< std::string>::const_iterator row_itr( row_names.Begin()), row_itr_end( row_names.End());
            row_itr != row_itr_end;
            ++row_itr
        )
        {
          topology_table.InsertRow( *row_itr);
        }

        // computation of comparisons of protein models
        storage::Table< double> table( ProteinModelCompareTopology( string_vector_a, topology_table));

        //output matrix to file "topology.txt"
        std::string out_file_name
        (
            m_OutputDirectory->GetFirstParameter()->GetValue()
          + m_OutFile->GetFirstParameter()->GetValue()
          + "_topology.txt"
        );

        BCL_MessageCrt( "outputting to file " + out_file_name);

        // write table to file
        io::OFStream write;
        io::File::MustOpenOFStream( write, out_file_name);
        table.WriteFormatted( write);

        // clear the stream
        io::File::CloseClearFStream( write);
      }
      // if using quality metrics other than topology score
      else
      {
        // columns for each pdb in vector b
        const util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( string_vector_b));

        // initialize map of quality measures
        storage::Map< quality::Measure, storage::Table< double> > quality_measures;

        // initialize map of normalized quality measures
        storage::Map< quality::Measure, storage::Table< double> > norm100_quality_measures;

        // construct the quality measure set
        m_Qualities = quality::Measures::GetCommandLineQualityMeasures();

        // iterate over list of quality measures passed to build up the tables of quality measures
        for
        (
          storage::Set< quality::Measure>::const_iterator
          quality_itr( m_Qualities.Begin()), quality_itr_end( m_Qualities.End());
          quality_itr != quality_itr_end; ++quality_itr
        )
        {
          // insert a table for "quality_measure" into each quality map
          quality_measures[ *quality_itr] = storage::Table< double>( sp_table_header);
          norm100_quality_measures[ *quality_itr] = storage::Table< double>( sp_table_header);

          // iterate over row names
          storage::Table< double> &current_table( quality_measures[ *quality_itr]);
          storage::Table< double> &current_table_norm100( norm100_quality_measures[ *quality_itr]);

          const double optimal_value( ( **quality_itr)->OptimalValue());
          // iterate through row names to insert empty rows, with optimal measure as the diagonal elements
          size_t row_count( 0);
          for
          (
            storage::Vector< std::string>::const_iterator row_itr( row_names.Begin()), row_itr_end( row_names.End());
            row_itr != row_itr_end;
            ++row_itr, ++row_count
          )
          {
            storage::Vector< double> row_values( sp_table_header->GetSize(), 0.0);
            if( half)
            {
              row_values( row_count) = optimal_value;
            }
            current_table.InsertRow( *row_itr, row_values);
            current_table_norm100.InsertRow( *row_itr, row_values);
          }
        }

        // print the sizes of the two pdb lists
        BCL_MessageStd( "size of pdb list a " + util::Format()( string_vector_a.GetSize()));
        BCL_MessageStd( "size of pdb list b " + util::Format()( string_vector_b.GetSize()));

        // do all of the quality calculations
        ProteinModelCompareAll( quality_measures, norm100_quality_measures, string_vector_a, string_vector_b, half);

        // output all of the quality measures in table format
        MatrixOutput( quality_measures, norm100_quality_measures);
      }
      //successful end
      return 0;
    }

    //! @brief StringVectorFromPDBListParameter gives all the strings in a file as a Vector
    //! @param FLAG the FlagInterface which has the string indicating the name of the list file
    //! @return return a Vector which has all of the strings contained within the file denoted by "FLAG"
    storage::Vector< std::string> ProteinCompare::StringVectorFromPDBListParameter
    (
      const util::ShPtr< command::FlagInterface> &FLAG
    ) const
    {
      // initialize write and read stream object
      io::IFStream pdb_list;

      // open pdb list file
      // make sure the file opened
      io::File::MustOpenIFStream( pdb_list, FLAG->GetFirstParameter()->GetValue());

      // Vector filled with the strings in the pdb list
      const storage::Vector< std::string> pdb_vector
      (
        util::StringLineListFromIStream( pdb_list)
      );
      // close and clear stream
      io::File::CloseClearFStream( pdb_list);

      return pdb_vector;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ThreadManager
    //! @brief manages threads for computing quality measures
    //!
    //! @remarks example unncessary
    //! @author mendenjl
    //! @date May 22, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ProteinModelThreadedComparer
    {

    public:

      typedef storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > ProteinList;

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Worker
      //! @brief calculates quality measures on pairs of proteins; requests additional work from Manager-type parent class
      //!
      //! @remarks example unncessary
      //! @author mendenjl
      //! @date May 22, 2014
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct Worker
      {
        //! Iterator to the first list of proteins
        ProteinList::const_iterator m_IteratorProteinA;

        //! Iterator to the second list of proteins
        ProteinList::const_iterator m_IteratorProteinB;

        //! Simple pointer to the output vector of quality measures
        util::SiPtr< linal::Vector< double> > m_QualityMeasureStorage;

        //! Simple pointer to the output vector of quality measures normalized per 100 residues
        util::SiPtr< linal::Vector< double> > m_QualityMeasureNorm100Storage;

        //! Simple pointer to the quality measures set
        util::SiPtr< const storage::Vector< quality::Measure> > m_QualityMeasures;

        //! Simple pointer to the alignment object, if it was given
        util::SiPtr< const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > m_Alignment;

        //! Simple pointer to the atoms to consider
        util::SiPtr< const storage::Set< biol::AtomType> > m_AtomsToConsider;

        //! simple pointer to the manager, needed so that the worker can get updated
        util::OwnPtr< ProteinModelThreadedComparer> m_Manager;

        //! @brief start threads to find all fragments that match a given molecule
        void RunThread();
      };

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief function to update the Worker object
      //! @param WORKER the worker to update
      //! @return true if the iterators are not at the end
      bool UpdateWorker( Worker &WORKER);

    //////////
    // data //
    //////////

    private:

      //! Bool: true if iterating over a symmetric list of proteins
      bool m_Symmetric;

      //! Iterator to m_IteratorProteinA of the last updated thread
      ProteinList::const_iterator m_IteratorProteinA;

      //! Iterator to m_IteratorProteinB of the last updated thread
      ProteinList::const_iterator m_IteratorProteinB;

      //! Begin iterator for protein B
      ProteinList::const_iterator m_IteratorProteinBBegin;

      //! End iterator for protein B
      ProteinList::const_iterator m_IteratorProteinBEnd;

      //! Size_t of current position in the proteins A list
      size_t m_PositionA;

      //! Size_t of current position in the proteins B list
      size_t m_PositionB;

      //! Size of list A; cached because list.GetSize() is an O(N) operation
      size_t m_ListASize;

      //! Size of list B; cached because list.GetSize() is an O(N) operation
      size_t m_ListBSize;

      //! Vector of quality measures
      storage::Vector< quality::Measure> m_QualityMeasures;

      //! Bool: true to calculate norm-100 values as well
      bool m_DoNorm100;

      //! Output vector of values per protein in list a, per protein in list b, and per quality measure,
      //! both the raw and normalized by per-100 residues
      storage::Vector< storage::Vector< linal::Vector< double> > > m_QualityValues;
      storage::Vector< storage::Vector< linal::Vector< double> > > m_QualityValuesNorm100;

      //! filenames
      storage::Vector< std::string> m_FilenamesA;

      //! Mutex to control access to thread state
      sched::Mutex m_WorkerStateMutex;

      //! Total number of comparisons that will be performed; cached for performance
      size_t m_TotalComparisons;

      //! Comparison id
      size_t m_ComparisonNumber;

      //! Status update interval
      size_t m_StatusUpdateInterval;

    public:

      //! @brief constructor from members; runs all comparisons
      //! @param STRING_VECTOR_A first vector of pdb names
      //! @param STRING_VECTOR_B second vector of pdb names
      //! @param LIST_A, LIST_B the lists of proteins to calculate quality measures over
      //!        If &LIST_A == &LIST_B, only the lower triangle of the matrix will be calculated
      //! @param ALIGNMENT non-empty to use a defined alignment for aligning all proteins
      //! @param INCLUDE_NORM100 whether to also calculate norm100 values for each quality measure
      //! @param ATOMS_TO_CONSIDER atoms to consider in the quality
      //! @param N_THREADS number of threads to launch
      //! @param QUALITY_MEASURES map that contains table with double qualitymeasure for each pair of pdbs mapped to the quality measure
      //! @param NORM100_QUALITY_MEASURES map that containes normalized quality measure
      ProteinModelThreadedComparer
      (
        const storage::Vector< std::string> &STRING_VECTOR_A,
        const storage::Vector< std::string> &STRING_VECTOR_B,
        const storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > &LIST_A,
        const storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > &LIST_B,
        const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENT,
        const bool INCLUDE_NORM100,
        const storage::Set< biol::AtomType> &ATOMS_TO_CONSIDER,
        const size_t N_THREADS,
        storage::Map< quality::Measure, storage::Table< double> > &QUALITY_MEASURES,
        storage::Map< quality::Measure, storage::Table< double> > &NORM100_QUALITY_MEASURES
      );

      //! @brief dummy clone function, needed so this object can be put in an own ptr
      ProteinModelThreadedComparer *Clone() const
      {
        BCL_Exit( "ProteinModelThreadedComparer has simple pointers to interally-held objecs and cannot be cloned", -1);
        return NULL;
      }
    };

    //! @brief constructor from members; runs all comparisons
    //! @param STRING_VECTOR_A first vector of pdb names
    //! @param STRING_VECTOR_B second vector of pdb names
    //! @param LIST_A, LIST_B the lists of proteins to calculate quality measures over
    //!        If &LIST_A == &LIST_B, only the lower triangle of the matrix will be calculated
    //! @param ALIGNMENT non-empty to use a defined alignment for aligning all proteins
    //! @param INCLUDE_NORM100 whether to also calculate norm100 values for each quality measure
    //! @param ATOMS_TO_CONSIDER atoms to consider in the quality
    //! @param N_THREADS number of threads to launch
    //! @param QUALITY_MEASURES map that contains table with double qualitymeasure for each pair of pdbs mapped to the quality measure
    //! @param NORM100_QUALITY_MEASURES map that containes normalized quality measure
    ProteinModelThreadedComparer::ProteinModelThreadedComparer
    (
      const storage::Vector< std::string> &STRING_VECTOR_A,
      const storage::Vector< std::string> &STRING_VECTOR_B,
      const storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > &LIST_A,
      const storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > &LIST_B,
      const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > &ALIGNMENT,
      const bool INCLUDE_NORM100,
      const storage::Set< biol::AtomType> &ATOMS_TO_CONSIDER,
      const size_t N_THREADS,
      storage::Map< quality::Measure, storage::Table< double> > &QUALITY_MEASURES,
      storage::Map< quality::Measure, storage::Table< double> > &NORM100_QUALITY_MEASURES
    ) :
      m_Symmetric( &LIST_A == &LIST_B),
      m_IteratorProteinA( LIST_A.Begin()),
      m_IteratorProteinB( LIST_B.Begin()),
      m_IteratorProteinBBegin( LIST_B.Begin()),
      m_IteratorProteinBEnd( m_Symmetric ? LIST_B.Begin() : LIST_B.End()),
      m_PositionA( 0),
      m_PositionB( 0),
      m_ListASize( STRING_VECTOR_A.GetSize()),
      m_ListBSize( STRING_VECTOR_B.GetSize()),
      m_QualityMeasures( QUALITY_MEASURES.GetKeysAsVector()),
      m_DoNorm100( INCLUDE_NORM100),
      m_QualityValues
      (
        m_ListASize,
        storage::Vector< linal::Vector< double> >( m_ListBSize, linal::Vector< double>( m_QualityMeasures.GetSize()))
      ),
      m_QualityValuesNorm100( m_DoNorm100 ? m_QualityValues : storage::Vector< storage::Vector< linal::Vector< double> > >()),
      m_FilenamesA( STRING_VECTOR_A),
      m_TotalComparisons( m_Symmetric ? m_ListASize * ( m_ListASize - 1) / 2 : m_ListASize * m_ListBSize),
      m_ComparisonNumber( 0)
    {
      // nothing to do if either list is empty
      if( !m_ListASize || !m_ListBSize)
      {
        m_TotalComparisons = 0;
        return;
      }
      BCL_Assert
      (
        !m_Symmetric || STRING_VECTOR_A == STRING_VECTOR_B,
        "Symmetric calculation should have referred to the same files"
      );
      BCL_Assert
      (
        LIST_A.GetSize() == STRING_VECTOR_A.GetSize(),
        "First list had different size from corresponding string vector"
      );

      // get the number of processors; start as many threads as allowed, unless there are too few calculations to perform
      size_t n_threads
      (
        std::min
        (
          std::max( sched::GetNumberCPUs(), size_t( 1)),
          std::min( N_THREADS, m_TotalComparisons)
        )
      );

      // skip status update if the message level is too high
      if( !util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Standard))
      {
        m_StatusUpdateInterval = m_TotalComparisons + 1;
      }
      else
      {
        // update the status bar at least once every 1000 comparisons, or every 0.1%, whichever is smaller
        m_StatusUpdateInterval = std::min( size_t( 999), size_t( m_TotalComparisons / 1000)) + 1;
      }

      // create a vector of workers, each of which will perform single protein-pair alignments and comparisons until
      // all comparisons are complete
      std::vector< Worker> workers( n_threads);
      for( std::vector< Worker>::iterator itr( workers.begin()), itr_end( workers.end()); itr != itr_end; ++itr)
      {
        // get the worker
        Worker &worker( *itr);

        // setup all the constant members of the worker
        worker.m_QualityMeasures = m_QualityMeasures;
        worker.m_Alignment = ALIGNMENT;
        worker.m_AtomsToConsider = ATOMS_TO_CONSIDER;
        worker.m_Manager = util::OwnPtr< ProteinModelThreadedComparer>( this, false);

        // call the update function
        UpdateWorker( worker);
      }

      // timer to see how long all the quality calculations take
      util::Stopwatch timer( "quality calculations", util::Time( std::numeric_limits< size_t>::max(), 0), util::Message::e_Critical);

      // make the job vector big enough to hold all the jobs

      // create a vector to hold the jobs
      util::ShPtrVector< sched::JobInterface> jobs;
      jobs.AllocateMemory( n_threads);
      const size_t group_id( 0);
      for( size_t processor_number( 0); processor_number < n_threads; ++processor_number)
      {
        // create the worker
        Worker &worker_reference( workers[ processor_number]);

        // create the job
        jobs.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::ThunkJob< Worker, void>
            (
              group_id,
              worker_reference,
              &Worker::RunThread,
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );
        // submit it to the scheduler
        sched::GetScheduler().RunJob( jobs.LastElement());
      }

      // join all the jobs
      for( size_t processor_number( 0); processor_number < n_threads; ++processor_number)
      {
        sched::GetScheduler().Join( jobs( processor_number));
      }

      // lastly, copy all values into the passed-in maps of tables.  The values in each linal::Vector< double>,
      // which correspond to the different quality measures, are guaranteed to be in the same order as the maps
      // because the vector of quality metrics is created from the iterated-over keys of the maps
      size_t quality_measure_id( 0);
      for
      (
        storage::Map< quality::Measure, storage::Table< double> >::iterator
          itr_table_raw( QUALITY_MEASURES.Begin()), itr_table_raw_end( QUALITY_MEASURES.End()),
          itr_table_norm( NORM100_QUALITY_MEASURES.Begin());
        itr_table_raw != itr_table_raw_end;
        ++itr_table_raw, ++itr_table_norm, ++quality_measure_id
      )
      {
        size_t protein_a_id( 0);
        for
        (
          storage::Table< double>::iterator
            itr_row_raw( itr_table_raw->second.Begin()), itr_row_raw_end( itr_table_raw->second.End()),
            itr_row_norm( itr_table_norm->second.Begin());
          itr_row_raw != itr_row_raw_end;
          ++itr_row_raw, ++itr_row_norm, ++protein_a_id
        )
        {
          size_t protein_b_id( 0);
          for
          (
            storage::Row< double>::iterator
              itr_col_raw( itr_row_raw->Second().Begin()),
              itr_col_raw_end( m_Symmetric ? itr_row_raw->Second().Begin() + protein_a_id : itr_row_raw->Second().End()),
              itr_col_norm( itr_row_norm->Second().Begin());
            itr_col_raw != itr_col_raw_end;
            ++itr_col_raw, ++itr_col_norm, ++protein_b_id
          )
          {
            *itr_col_raw = m_QualityValues( protein_a_id)( protein_b_id)( quality_measure_id);
            if( m_DoNorm100)
            {
              *itr_col_norm = m_QualityValuesNorm100( protein_a_id)( protein_b_id)( quality_measure_id);
            }
          }
        }
      }
    }

    //! @brief function to update the Worker object
    //! @param WORKER the worker to update
    //! @return true if the iterators are not at the end
    bool ProteinModelThreadedComparer::UpdateWorker( Worker &WORKER)
    {
      m_WorkerStateMutex.Lock();

      // test whether A is at the end or B is at the end and the next A is the end

      // check for reaching the end of the fragment list
      if( m_PositionA == m_ListASize || ( m_IteratorProteinB == m_IteratorProteinBEnd && ++m_PositionA == m_ListASize))
      {
        // write status so that user knows that comparisons are finished
        util::GetLogger().LogStatus
        (
          "=std=bcl::app=> Total comparisons done " + util::Format()( m_ComparisonNumber)
          + " / " + util::Format()( m_TotalComparisons) + " 100.00%"
        );
        m_WorkerStateMutex.Unlock();
        return false;
      }
      else if( m_IteratorProteinB == m_IteratorProteinBEnd)
      {
        // move to next protein in list A, reset iterator on list B
        ++m_IteratorProteinA;

        // already incremented m_PositionA in previous _if_ check
        m_IteratorProteinB = m_IteratorProteinBBegin;
        m_PositionB = 0;

        // for symmetric comparisons, lists A and B are the same, so stop iterator B once it reaches iterator A
        // so that only the lower triangular matrix is calculated
        if( m_Symmetric)
        {
          m_IteratorProteinBEnd = m_IteratorProteinA;
        }
      }

      // update worker's iterators and storage
      WORKER.m_IteratorProteinA = m_IteratorProteinA;
      WORKER.m_IteratorProteinB = m_IteratorProteinB++;
      WORKER.m_QualityMeasureStorage = m_QualityValues( m_PositionA)( m_PositionB);
      if( !m_QualityValuesNorm100.IsEmpty())
      {
        WORKER.m_QualityMeasureNorm100Storage = m_QualityValuesNorm100( m_PositionA)( m_PositionB);
      }

      ++m_ComparisonNumber;
      ++m_PositionB;

      // write status, if message level is high enough every 1000 comparisons
      if( ( m_ComparisonNumber % m_StatusUpdateInterval) == 0)
      {
        util::GetLogger().LogStatus
        (
          "=std=bcl::app=> Quality involving " + m_FilenamesA( m_PositionA) + " File "
          + util::Format()( m_PositionA) + " / " + util::Format()( m_ListASize)
          + " Total comparisons done " + util::Format()( m_ComparisonNumber)
          + " / " + util::Format()( m_TotalComparisons) + " "
          + util::Format().FFP( 2)( float( m_ComparisonNumber) * 100.0 / float( m_TotalComparisons))
          + "%"
        );
      }

      // unlock the mutex and return
      m_WorkerStateMutex.Unlock();
      return true;
    }

    //! @brief start threads to find all fragments that match a given molecule
    void ProteinModelThreadedComparer::Worker::RunThread()
    {
      do
      {
        // initialize alignment map
        storage::List< align::AlignmentNodeReference< biol::AABase> > alignment;
        // initialize coordinates vector
        storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > common_atom_coordinates;

        // initialize counters
        size_t pruned_ctr( 0);

        // if an alignment was read in
        if( !m_Alignment->IsEmpty())
        {
          // iterate over the map
          for
          (
            storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
            map_itr( m_Alignment->Begin()), map_itr_end( m_Alignment->End());
            map_itr != map_itr_end; ++map_itr
          )
          {
            // get the corresponding alignments from the protein models
            storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator find_itr_a
            (
              m_IteratorProteinA->Find( map_itr->first)
            );

            // get the corresponding alignments from the protein models
            storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator find_itr_b
            (
              m_IteratorProteinB->Find( map_itr->first)
            );

            // assert that the alignments were found
            BCL_Assert
            (
              find_itr_a != m_IteratorProteinA->End() && find_itr_b != m_IteratorProteinB->End(),
              "Unable to find chain ID, " + util::Format()( map_itr->first)
            );

            // copy shared pointers due to interface requirements; ultimately due to shared-pointer const-incorrectness
            alignment.PushBack( align::AlignmentNodeReference< biol::AABase>( find_itr_a->second, find_itr_b->second));
            // perform the alignment
            biol::AlignByAAData().AlignPairWithNode( alignment.LastElement(), map_itr->second);
          }
        }
        // no alignment was read in
        else
        {
          // iterate over the model a map
          for
          (
            storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator
            seq_itr( m_IteratorProteinA->Begin()), seq_itr_end( m_IteratorProteinA->End());
            seq_itr != seq_itr_end; ++seq_itr
          )
          {
            // try to find the corresponding sequence from the sse map
            storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > >::const_iterator find_itr
            (
              m_IteratorProteinB->Find( seq_itr->first)
            );

            // assert that the alignments were found
            BCL_Assert
            (
              find_itr != m_IteratorProteinB->End(),
              "Unable to find chain ID, " + util::Format()( seq_itr->first)
            );

            alignment.PushBack( align::AlignmentNodeReference< biol::AABase>( seq_itr->second, find_itr->second));
            // perform the alignment
            biol::AlignByAAData().AlignPairWithNode( alignment.LastElement(), biol::AACompareData());
          }
        }

        for
        (
          storage::List< align::AlignmentNodeReference< biol::AABase> >::iterator
          seq_itr( alignment.Begin()), seq_itr_end( alignment.End());
          seq_itr != seq_itr_end; ++seq_itr
        )
        {
          // get the common coordinates
          const storage::Pair< storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> >, size_t>
          this_coordinates_pair_and_size
          (
            assemble::Quality::CoordinatesAndResidueCountFromAlignment( *seq_itr, *m_AtomsToConsider)
          );
          pruned_ctr += this_coordinates_pair_and_size.Second();

          // update the total coordinates
          common_atom_coordinates.First().Append( this_coordinates_pair_and_size.First().First());
          common_atom_coordinates.Second().Append( this_coordinates_pair_and_size.First().Second());
        }

        // iterate over the quality measures
        linal::Vector< double>::iterator quality_storage_itr( m_QualityMeasureStorage->Begin());
        linal::Vector< double>::iterator quality_norm100_storage_itr
        (
          m_QualityMeasureNorm100Storage.IsDefined()
          ? m_QualityMeasureNorm100Storage->Begin()
          : m_QualityMeasureStorage->Begin()
        );

        for
        (
          storage::Vector< quality::Measure>::const_iterator
          itr( m_QualityMeasures->Begin()), itr_end( m_QualityMeasures->End());
          itr != itr_end;
          ++itr, ++quality_storage_itr, ++quality_norm100_storage_itr
        )
        {
          // "quality_value" initialize with current quality measure
          const double quality_value
          (
            ( **itr)->CalculateMeasure( common_atom_coordinates.First(), common_atom_coordinates.Second())
          );

          BCL_MessageDbg( util::Format()( quality_value));
          // assign the quality value
          ( *quality_storage_itr) = quality_value;

          // if normalized quality measure is desired then calculate it
          if( m_QualityMeasureNorm100Storage.IsDefined())
          {
            const double norm_value( assemble::Quality::RMSD100( quality_value, pruned_ctr));
            BCL_MessageVrb
            (
              "Calculating norm100_" + itr->GetName() + "= " + util::Format()( norm_value)
              + " from quality_value=" + util::Format()( quality_value)
              + " and amino_acid_count=" + util::Format()( pruned_ctr)
            );

            // "normalized_quality_value" initialize with the normalized current quality measure
            ( *quality_norm100_storage_itr) = norm_value;
          }
        }
      } while( m_Manager->UpdateWorker( *this));
    }

    //! @brief ProteinModelCompareAll calculates all quality measures between all necessary PDBs
    //! @param QUALITY_MEASURES map that contains table with double qualitymeasure for each pair of pdbs mapped to the quality measure
    //! @param NORM100_QUALITY_MEASURES map that containes normalized quality measure
    //! @param STRING_VECTOR_A first vector of pdb names
    //! @param STRING_VECTOR_B second vector of pdb names
    //! @param HALF whether only half the matrix will be calculated, applicable when string vectors are  the same
    void ProteinCompare::ProteinModelCompareAll
    (
      storage::Map< quality::Measure, storage::Table< double> > &QUALITY_MEASURES,
      storage::Map< quality::Measure, storage::Table< double> > &NORM100_QUALITY_MEASURES,
      const storage::Vector< std::string> &STRING_VECTOR_A,
      const storage::Vector< std::string> &STRING_VECTOR_B,
      const bool HALF
    ) const
    {
      // assert that given string vectors are of same size
      BCL_Assert
      (
        !HALF || STRING_VECTOR_A.GetSize() == STRING_VECTOR_B.GetSize(),
        "if only half of the matrix is calculated, the size of the protein list have to match"
      );

      // create "pdb_list_a" which holds the first pdb list
      storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > pdb_list_a
      (
        GenerateAlignments( STRING_VECTOR_A)
      );

      storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > pdb_list_b;
      if( !HALF)
      {
        pdb_list_b = GenerateAlignments( STRING_VECTOR_B);
      }

      // create the comparer object and do the comparisons
      ProteinModelThreadedComparer manager
      (
        STRING_VECTOR_A,
        STRING_VECTOR_B,
        pdb_list_a,
        HALF ? pdb_list_a : pdb_list_b,
        m_Alignment,
        m_QualityNormalization100->GetFlag(),
        m_ConsideredAtoms,
        sched::GetNumberCPUs(),
        QUALITY_MEASURES,
        NORM100_QUALITY_MEASURES
      );
    }

    //! @brief ProteinModelCompareTopology calculates topology metric between all necessary PDBs
    //! @param STRING_VECTOR vector of pdb names
    //! @param TOPOLOGY_TABLE table that will be filled with topology scores
    storage::Table< double> ProteinCompare::ProteinModelCompareTopology
    (
      const storage::Vector< std::string> &STRING_VECTOR,
      const storage::Table< double> &TOPOLOGY_TABLE
    ) const
    {
      // initialize table
      storage::Table< double> topology_table( TOPOLOGY_TABLE);

      // contains all protein models read in from pdb files
      util::ShPtrVector< assemble::ProteinModel> model_list;

      // iterate over pdb names and compare pairwise
      for
      (
        storage::Vector< std::string>::const_iterator a_itr( STRING_VECTOR.Begin()), 
        a_itr_end( STRING_VECTOR.End());
        a_itr != a_itr_end; a_itr++
      )
      {
        storage::Row< double> &row( topology_table[ *a_itr]);

        // read in protein model for this row
        const std::string proteinA( *a_itr);
        pdb::Factory factory;
        io::IFStream read;
        io::File::MustOpenIFStream
        (
          read,
          proteinA
        );

        pdb::Handler pdbA( read, true);
        io::File::CloseClearFStream( read);
       
        // get model from pdb and add it to protein list
        assemble::ProteinModel modelA( factory.ProteinModelFromPDB( pdbA));
        model_list.PushBack( util::ShPtr< assemble::ProteinModel>( new assemble::ProteinModel( modelA))); 

        // iterate over current model list to perform pairwise comparisons
        util::ShPtrVector< assemble::ProteinModel>::const_iterator model_itr( model_list.Begin());

        // iterate over all pdbs with which to perform pairwise calculations
        for
        (
          storage::Vector< std::string>::const_iterator b_itr( STRING_VECTOR.Begin());
          b_itr != a_itr; b_itr++, model_itr++
        )
        {
          // second model for comparison
          const assemble::ProteinModel modelB( **model_itr);

          // calculate topology score between these two protein models
          score::ProteinModelFragmentTopology topology_score;
          const double quality_value
          (
            topology_score.operator()( modelA, modelB)
          );

          BCL_MessageCrt( "quality value " + util::Format()( quality_value));

          // fill in matrix with this value
          double &score( row[ *b_itr]);
          score = quality_value;
        }
      }

      // return lower triangular matrix of pairwise protein model comparisons
      return topology_table;
    }

    //! @brief creates an alignment for each protein model in the string vector
    //! @param PDB_VECTOR is the string vector which gives the pdbs to create protein models from
    //! @return  an alignment for each protein model in the string vector
    storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > >
    ProteinCompare::GenerateAlignments( const storage::Vector< std::string> &PDB_VECTOR) const
    {
      // timer to see how long it takes to create the protein model information
      util::Stopwatch timer
      (
        "create protein models and defined aas",
        util::Time( std::numeric_limits< size_t>::max(), 0),
        util::Message::e_Critical
      );

      // create "alignment_list" an alignment of amino acids in the protein model
      storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > > alignment_list;

      // iterate though the vector of pdb filenames to create the protein models and get their defined amino acids
      for
      (
        storage::Vector< std::string>::const_iterator
          pdb_name_itr( PDB_VECTOR.Begin()), pdb_name_itr_end( PDB_VECTOR.End());
        pdb_name_itr != pdb_name_itr_end;
        ++pdb_name_itr
      )
      {
        // create "factory" to create protein model with amino acids of type AABackBone
        pdb::Factory factory;

        // initialize write and read stream objects
        io::IFStream read;

        // open pdb file
        io::File::MustOpenIFStream
        (
          read,
            m_InputPDBPrefix->GetFirstParameter()->GetValue()
          + *pdb_name_itr
          + m_InputPDBPostfix->GetFirstParameter()->GetValue()
        );

        // create pdb handler "pdb" and pass "read" and true to indicate ignore clashes
        pdb::Handler pdb( read, true);

        // close and clear "read"
        io::File::CloseClearFStream( read);

        // read the protein model
        assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb));

        // create the alignment and push it back
        if( GetFlagSpecifyQualityResidues()->GetFlag())
        {
          static const assemble::CollectorAASpecified s_collector
          (
            GetFlagSpecifyQualityResidues()->GetFirstParameter()->GetValue()
          );
          alignment_list.PushBack( assemble::Quality::CreateAlignmentFromCollectorAASpecified( model, s_collector));
        }
        else //< use all residues in sses
        {
          alignment_list.PushBack( assemble::Quality::CreateAlignmentFromProteinModelSSEs( model));
        }
      }

      // return "alignment_list"
      return alignment_list;
    }

    //! @brief MatrixOutput outputs the quality measures in matrix format to files of quality name (ex. RMSD, ROC3)
    //! @param QUALITY_MEASURES map that contains table with double quality measure for each pair of pdbs mapped to the quality measure
    //! @param NORM100_QUALITY_MEASURES map that contains normalized quality measures
    //! @return void outputs a matrix for each quality measure to a separate file in directory application is executed
    void ProteinCompare::MatrixOutput
    (
      const storage::Map< quality::Measure, storage::Table< double> > &QUALITY_MEASURES,
      const storage::Map< quality::Measure, storage::Table< double> > &NORM100_QUALITY_MEASURES
    ) const
    {
      // timer to measure how long it takes to output the tables
      util::Stopwatch timer
      (
        "Table Outputting",
        util::Time( std::numeric_limits< size_t>::max(), 0),
        util::Message::e_Critical
      );

      // initialize write and read stream objects for normalized and regular output
      io::OFStream write, write_norm100;

      // iterate over quality measures; one table per quality measure
      for
      (
        storage::Map< quality::Measure, storage::Table< double> >::const_iterator
          quality_itr( QUALITY_MEASURES.Begin()), norm100_quality_itr( NORM100_QUALITY_MEASURES.Begin()),
          quality_itr_end( QUALITY_MEASURES.End()), norm100_quality_itr_end( NORM100_QUALITY_MEASURES.End());
        quality_itr != quality_itr_end && norm100_quality_itr != norm100_quality_itr_end;
        ++quality_itr, ++norm100_quality_itr
      )
      {
        // "out_file_name" for unnormalized values
        std::string out_file_name
        (
            m_OutputDirectory->GetFirstParameter()->GetValue()
          + m_OutFile->GetFirstParameter()->GetValue()
          + quality_itr->first.GetName() + ".txt"
        );

        BCL_MessageTop( "outputting to file " + out_file_name);

        // open output file
        io::File::MustOpenOFStream( write, out_file_name);

        // write the matrices to files with or without labels depending on the flags given
        if( m_NoLabels->GetFlag())
        {
          quality_itr->second.WriteFormattedWithoutNames( write);
        }
        else
        {
          quality_itr->second.WriteFormatted( write);
        }
        // clear the stream
        io::File::CloseClearFStream( write);

        // do the same for the normalized values if requested
        if( m_QualityNormalization100->GetFlag())
        {
          // "norm100_out_file_name" for normalized output
          std::string norm100_out_file_name
          (
              m_OutputDirectory->GetFirstParameter()->GetValue()
            + m_OutFile->GetFirstParameter()->GetValue()
            + quality_itr->first.GetName() + "_norm100.txt"
          );
          // open output file
          io::File::MustOpenOFStream( write_norm100, norm100_out_file_name);

          if( m_NoLabels->GetFlag())
          {
            norm100_quality_itr->second.WriteFormattedWithoutNames( write_norm100);
          }
          else
          {
            // write out the normalized values
            norm100_quality_itr->second.WriteFormatted( write_norm100);
          }
          // reset the file
          io::File::CloseClearFStream( write_norm100);
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinCompare::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ProteinCompare::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
