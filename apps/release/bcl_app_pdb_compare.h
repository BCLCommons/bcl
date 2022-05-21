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

#ifndef BCL_APP_PDB_COMPARE_H_
#define BCL_APP_PDB_COMPARE_H_

// include header of this class
#include "app/bcl_app_apps.h"

// align
#include "align/bcl_align_alignment_node.h"
#include "align/bcl_align_handler_pir.h"

// assemble
#include "assemble/bcl_assemble_collector_aa_specified.h"
#include "assemble/bcl_assemble_collector_common_aa.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"

// biol
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_align_by_aa_data.h"

// command
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_ranged.h"

// coord
#include "coord/bcl_coord.h"

// io
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"

// pdb
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_line.h"
#include "pdb/bcl_pdb_residue.h"

// quality
#include "quality/bcl_quality_rmsd.h"

// sched
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"

//score
#include "score/bcl_score_protein_model_fragment_topology.h"
#include "score/bcl_score_protein_model_topology.h"

// storage
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector_nd.h"

// util
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinCompare
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_app_pdb_compare.cpp @endlink
    //! @author alexanns
    //! @note threading by mendenjl May 22, 2014
    //! @date Jan 14, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinCompare :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! first pdb
      util::ShPtr< command::FlagInterface> m_ReferencePDB;

      //! second pdb
      util::ShPtr< command::FlagInterface> m_TestPDB;

      //! list of pdbs
      util::ShPtr< command::FlagInterface> m_PDBList;

      //! stream to write the results to
      util::ShPtr< command::FlagInterface> m_OutFile;

      //! normalization flag
      util::ShPtr< command::FlagInterface> m_QualityNormalization100;

      //! flag for indicating that no labeling for the output matrix is desired
      util::ShPtr< command::FlagInterface> m_NoLabels;

      //! flag for indicating a prefix which must be prepended to the PDB strings in order to access the actual files
      util::ShPtr< command::FlagInterface> m_InputPDBPrefix;

      //! flag for indicating a prefix which must be appended to the PDB strings in order to access the actual files
      util::ShPtr< command::FlagInterface> m_InputPDBPostfix;

      //! flag for indicating the directory into which the quality measure matrices should be deposited
      util::ShPtr< command::FlagInterface> m_OutputDirectory;

      //! flag for use of topology score as quality metric
      util::ShPtr< command::FlagInterface> m_Topology;

      //! atoms as given on m_AtomList to be used for rmsd calculation
      mutable storage::Set< biol::AtomType> m_ConsideredAtoms;

      //! alignment
      mutable storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > m_Alignment;

      mutable storage::Set< quality::Measure> m_Qualities;

    public:

      static const ApplicationType PDBCompare_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinCompare();

      //! @brief Clone function
      //! @return pointer to new Quality
      ProteinCompare *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the original name of the application, as it was used in license files
      //! @return string for the bcl::commons name of that application
      //! This is necessary so that, if release application names change, licenses will continue to work
      std::string GetLicensedName() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief return command line flag specifying specific residues that should be used for quality calculations
      //! @return command line flag for specifying specific residues to use during quality calculations
      static util::ShPtr< command::FlagInterface> &GetFlagSpecifyQualityResidues();

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief creates an alignment for each protein model in the string vector
      //! @param PDB_VECTOR is the string vector which gives the pdbs to create protein models from
      //! @return  an alignment for each protein model in the string vector
      storage::List< storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > >
      GenerateAlignments( const storage::Vector< std::string> &PDB_VECTOR) const;

      //! @brief StringVectorFromPDBListParameter gives all the strings in a file as a Vector
      //! @param FLAG the FlagInterface which has the string indicating the name of the list file
      //! @return return a Vector which has all of the strings contained within the file denoted by "FLAG"
      storage::Vector< std::string> StringVectorFromPDBListParameter
      (
        const util::ShPtr< command::FlagInterface> &FLAG
      ) const;

      //! @brief ProteinModelCompareAll calculates all quality measures between all necessary PDBs
      //! @param QUALITY_MEASURES map that contains table with double qualitymeasure for each pair of pdbs mapped to the quality measure
      //! @param NORM100_QUALITY_MEASURES map that containes normalized quality measure
      //! @param STRING_VECTOR_A first vector of pdb names
      //! @param STRING_VECTOR_B second vector of pdb names
      //! @param HALF whether only half the matrix will be calculated, applicable when string vectors are  the same
      void ProteinModelCompareAll
      (
        storage::Map< quality::Measure, storage::Table< double> > &QUALITY_MEASURES,
        storage::Map< quality::Measure, storage::Table< double> > &NORM100_QUALITY_MEASURES,
        const storage::Vector< std::string> &STRING_VECTOR_A,
        const storage::Vector< std::string> &STRING_VECTOR_B,
        const bool HALF
      ) const;

      //! @brief ProteinModelCompareTopology calculates topology metric between all necessary PDBs
      //! @param STRING_VECTOR vector of pdb names
      //! @param TOPOLOGY_TABLE table that will be filled with topology scores
      storage::Table< double> ProteinModelCompareTopology
      (
        const storage::Vector< std::string> &STRING_VECTOR,
        const storage::Table< double> &TOPOLOGY_TABLE
      ) const;

      //! @brief MatrixOutput outputs the quality measures in matrix format to files of quality name (ex. RMSD, ROC3)
      //! @param QUALITY_MEASURES map that contains table with double qualitymeasure for each pair of pdbs mapped to the quality measure
      //! @param NORM100_QUALITY_MEASURES map that contains nomralized quality measures
      //! @return void outputs a matrix for each quality measure to a separate file in directory application is executed
      void MatrixOutput
      (
        const storage::Map< quality::Measure, storage::Table< double> > &QUALITY_MEASURES,
        const storage::Map< quality::Measure, storage::Table< double> > &NORM100_QUALITY_MEASURE
      ) const;

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
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PDBCompare

  } // namespace app
} // namespace bcl

#endif // BCL_APP_PDB_COMPARE_H_
